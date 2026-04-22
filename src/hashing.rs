//simplified hashing
// we asume we just have DNA and use ACGT

// 16 mb needs to be divisible by 4 for dna4
const CHUNK_SIZE: usize = 16 * 1024 * 1024;
// assert!(CHUNK_SIZE.is_multiple_of(4));

const TO_DNA4: [[u8; 4]; 256] = dna4_table();
const TO_CHAR: [u8; 4] = [b'A',b'C',b'G',b'T'];

const fn dna4_table() -> [[u8; 4]; 256] {
    //in case i hit something else it will still be ACGT and skip the last macro step
    let mut look_up_table = [[0b0000_0000, 0b0000_0001, 0b0000_0010, 0b0000_0011]; 256];
    macro_rules! add4 {
        ($base:expr, $arr:expr) => {
            look_up_table[$base as usize] = $arr;
        };
    }
    let a: u8 = 0b0000_0000;
    let c: u8 = 0b0000_0001;
    let g: u8 = 0b0000_0010;
    let t: u8 = 0b0000_0011;

    // derived from https://en.wikipedia.org/wiki/Nucleic_acid_notation
    // Standard Bases
    add4!(b'A', [a,a,a,a]); add4!(b'C', [c,c,c,c]); add4!(b'G', [g,g,g,g]); add4!(b'T', [t,t,t,t]); add4!(b'U', [t,t,t,t]);
    add4!(b'a', [a,a,a,a]); add4!(b'c', [c,c,c,c]); add4!(b'g', [g,g,g,g]); add4!(b't', [t,t,t,t]); add4!(b'u', [t,t,t,t]);

    // Weak,Strong,Amino,Ketone,Purine,Pyrimide
    add4!(b'W', [a,t,a,t]); add4!(b'S', [c,g,c,g]); add4!(b'M', [a,c,a,c]); add4!(b'K', [g,t,g,t]); add4!(b'R', [a,g,a,g]); add4!(b'Y', [c,t,c,t]);
    add4!(b'w', [a,t,a,t]); add4!(b's', [c,g,c,g]); add4!(b'm', [a,c,a,c]); add4!(b'k', [g,t,g,t]); add4!(b'r', [a,g,a,g]); add4!(b'y', [c,t,c,t]);

    // not ACGT, assuming they appear in similar proportions
    // each letter has a bias but they correspond to their complement B<>V, D<>H
    add4!(b'B', [c,g,t,c]); add4!(b'D', [a,g,t,a]); add4!(b'H', [a,c,t,t]); add4!(b'V', [a,c,g,g]);
    add4!(b'b', [c,g,t,c]); add4!(b'd', [a,g,t,a]); add4!(b'h', [a,c,t,t]); add4!(b'v', [a,c,g,g]);

    // literally no information and already filled up
    add4!(b'N', [a,c,g,t]); add4!(b'n', [a,c,g,t]);
    look_up_table
}


//some helpers
fn get0(kmer: u8) -> u8 { (kmer & 0b1100_0000) >> 6 }
fn get1(kmer: u8) -> u8 { (kmer & 0b0011_0000) >> 4 }
fn get2(kmer: u8) -> u8 { (kmer & 0b0000_1100) >> 2 }
fn get3(kmer: u8) -> u8 { kmer & 0b0000_0011 }


pub fn hash_to_dna4(chunk: &mut [u8], seed: &mut u64) {
    let mut i: usize;
    for byte in chunk.iter_mut() {
        //xorshift64 for pseudo randomness in translating other Amino Acids
        //from https://en.wikipedia.org/wiki/Xorshift
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed >> 17;
        i = (*seed & 0b11) as usize;
        *byte = TO_DNA4[*byte as usize][i];
    }
}

// //to not always overcount the same ACGT
// pub fn give_pseudo_random(some_bits: &[u8]) -> u8 {
pub fn rolling(kmer: &mut Vec<u8>, next_bit: u8) {
    for i in 1..kmer.len() {
        kmer[i-1] = (kmer[i-1] << 2) | get0(kmer[i]);
    }
    if let Some(last_byte) = kmer.last_mut() {
        *last_byte =  (*last_byte | next_bit) << 6
    }
}

// safe version, for the end_part
// here we calculate the correct kmer_length, to get the last ones we otherwise miss
// there is probably a smarter and faster way, but i don'T know
pub fn rolling_kmer_end_hash(char_stream: &[u8], seed: &mut u64, k_order: usize) -> KmerTrie {
    let mut end_trie: KmerTrie = KmerTrie::new(k_order as u8, 0); //this is another order, from the rest

    // we fill the next byte 
    let fake_order: usize = k_order.next_multiple_of(4) + 1;
    let mut kmer: Vec<u8> = vec![0b0000_0000, fake_order as u8]; 
    let mut i: usize;
    let mut bit;

    for byte in char_stream.iter().take(k_order){
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed >> 17;
        i = (*seed & 0b11) as usize;
        bit = TO_DNA4[*byte as usize][i];
        rolling(&mut kmer, bit);
    }
    end_trie.insert_kmer_counting(&kmer[..]);

    for byte in char_stream.iter().skip(k_order){
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed >> 17;
        i = (*seed & 0b11) as usize;
        bit = TO_DNA4[*byte as usize][i];
        rolling(&mut kmer, bit);
        end_trie.insert_kmer_counting(&kmer[..]);
    }

    for _ in k_order..fake_order{
        rolling(&mut kmer, 0b0000_0000); //to count the remainding kmer
        end_trie.insert_kmer_counting(&kmer[..]);
    }
    end_trie.reduce_to_order(k_order as u8);
    end_trie
}

pub fn compress_to_dna4(char_stream: &[u8]) -> [u8; CHUNK_SIZE/4]  {
    let mut hashed_chunk = [0b0000_0000; CHUNK_SIZE / 4];
    let mut i: usize = 0;
    for base in char_stream.chunks_exact(4) {
        hashed_chunk[i] = base[0] | (base[1]<<2) | (base[2]<<4) | (base[3]<<6);
        i += 1;
    }
    hashed_chunk
}

pub fn decompress_to_dna4(hashed_chunk: &[u8], len: usize) -> Vec<u8> {
    let mut rehashed_chunk = vec![0; len.next_multiple_of(4)];
    let mut i: usize = 0;
    for &kmer in hashed_chunk.iter() {
        rehashed_chunk[i] = get0(kmer);
        rehashed_chunk[i+1] = get1(kmer);
        rehashed_chunk[i+2] = get2(kmer);
        rehashed_chunk[i+3] = get3(kmer);
        i += 4;
    }
    rehashed_chunk.truncate(len);
    rehashed_chunk
}


fn kmer0(compressed_chunk: &[u8], block_length: usize) -> Vec<u8> {
    let mut return_vec: Vec<u8> = Vec::from(compressed_chunk);
    let b0: u8 = 0b1100_0000;
    for i in (block_length-1..return_vec.len()).step_by(block_length) {
        return_vec[i] &= b0;
    }
    return_vec
}

//probably some simd magic would be nice here
fn kmer_shift(compressed_chunk: &mut [u8]) {
    for i in 1..compressed_chunk.len() {
        compressed_chunk[i-1] = (compressed_chunk[i-1] << 2) | get0(compressed_chunk[i]);
    }
}

fn kmer_next(compressed_chunk: &mut [u8], block_length: usize) -> Vec<u8> {
    kmer_shift(compressed_chunk);
    kmer0(compressed_chunk, block_length)
}

// fn rolling_kmer(compressed_chunk: &mut [u8], order_k: u8) -> Kmer_Trie {
//     return_trie = Kmer_Trie::new(order_k, 0u8);
//     block_length = (order_k as usize).next_multiple_of(4);
//     let kmer: Vec<u8> = Vec::with_capacity[block_length];
//     let b0: u8 = 0b0000_0011;
//     for i, chunk in compressed_chunk.iter().take(block_length).enumerate(){
//         kmer[i] = chunk;
//     }
//     if let Some(last_byte) = compressed_chunk.last_mut() {
//         *last =  (*last | next_bit) << 6
//     }
//     for 
// }

fn kmer_count(compressed_chunk: &mut [u8], order_k: u8) -> KmerTrie {
    let mut return_trie = KmerTrie::new(order_k, 0u8);
    let block_length = (order_k as usize).next_multiple_of(4);
    return_trie.insert_kmer_counting(&kmer0(compressed_chunk, block_length)[..]);
    return_trie.insert_kmer_counting(&kmer_next(compressed_chunk, block_length)[..]);
    return_trie.insert_kmer_counting(&kmer_next(compressed_chunk, block_length)[..]);
    return_trie.insert_kmer_counting(&kmer_next(compressed_chunk, block_length)[..]);
    return_trie
}

//fn alt_run_markov(args : &SimulateArgs) {
    // #[arg(short, long)]
    // pub input: String,

    // #[arg(short, long, default_value_t = String::from("seq.fasta"))]
    // pub output: String,

    // #[arg(short, long)]
    // pub lens: Vec<usize>,

    // #[arg(long, default_value_t = 3)]
    // pub order: usize,

    // #[arg(long, default_value_t = 42)]
    // pub seed: u64,

    // #[arg(short, long, default_value_t = false)]
    // pub verbose: bool,

    
}
// fn main() {
//     // Collect command-line arguments into a vector of strings
//     let args: Vec<String> = args().collect();

//     // Ensure the path to the file is provided as the second argument
//     // Extracts a reference to the path string.
//     let path: &String = &args.get(1).expect("File path not specified.");

//     // Retrieve the metadata for the open file to get the length in `usize` bytes
//     // If the metadata retrieval fails, program will panic
//     let file_length: usize = metadata(&path)
//         .expect("Unable to query file details")
//         .len()
//         .try_into()
//         .expect("Cannot convert file length");

//     // Size of each block to be read from the file (16MB)
//     const BLOCK_SIZE: usize = 16_777_216;

//     // Number of threads to be used for reading the file
//     const THREADS: usize = 4;

//     // Determine the portion of the file each thread will handle
//     let division: usize = (file_length + THREADS - 1) / THREADS;

//     // Use scoped threads to ensure all threads are joined before the main thread exits
//     thread::scope(|scope: &thread::Scope| {
//         // Create `THREADS` number of threads
//         for i in 0..THREADS {
//             scope.spawn(move || {
//                 // Open a file handle per thread
//                 // Ensuring each thread works on its own file descriptor
//                 let mut thread_file: File = File::open(&path).expect("Unable to open file");

//                 // Initialize a 16MB buffer to hold data read by this thread
//                 let mut contents: Vec<u8> = vec![0_u8; BLOCK_SIZE];

//                 // Counter for the number of bytes read in each operation
//                 // Init'd as 1, as zero would indicate an EOF
//                 let mut read_length: usize = 1;

//                 // Counter to keep track of the total nuber of bytes read by this thread
//                 let mut read_total: usize = 0;

//                 // Determing the offset in the file where this thread should start reading.
//                 // Each thred reads a different portion of the file
//                 let offset: u64 = (i * division) as u64;

//                 // Seek to the starting position in the file for this thread
//                 // If it couldn't start from there, the program panics
//                 thread_file
//                     .seek(SeekFrom::Start(offset))
//                     .expect("Couldn't seek to position in file");

//                 // Read data in iterations until the thread's portion is read or EOF is reached
//                 // i.e there is no more bytes to be read from the file
//                 while (read_total < division) && (read_length != 0) {
//                     // Adjust the contents buffer size if the remaining bytes to be read are less than the block size
//                     if read_total + BLOCK_SIZE > division {
//                         // Reduce the vector to ensure we don't read beyond the division
//                         // If this doesn't happen, data for the next thread will be read causing offset seek errors
//                         contents.truncate(division - read_total);
//                     }

//                     // Read data into the contents buffer
//                     // `read_length` is updated with the number of bytes read
//                     read_length = thread_file.read(&mut contents).expect("Couldn't read file");
//                     read_total += read_length;
//                 }

//                 // Verify the entire file was read correctly by visualizing total number of bytes read from the file
//                 // by each thread and the original file length
//                 println!("Thread {i}: Total bytes read: {read_total} bytes || Expected: {division} bytes");
//             });
//         }
//     });
