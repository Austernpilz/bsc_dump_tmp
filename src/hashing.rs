//simplified hashing
// we asume we just have DNA and use ACGT

// 16 mb
const CHUNK_SIZE: u64 = 16 * 1024 * 1024;

const TO_DNA4: [[u8; 4]; 256] = {
    //in case i hit something else it will still be ACGT and skip the last macro step
    let mut look_up_table = [[0,1,2,3]; 256];
    macro_rules! add4 {
        ($c:expr, $dna:expr) => {
            look_up_table[$c as usize] = $dna;
        }
    };
    // derived from https://en.wikipedia.org/wiki/Nucleic_acid_notation
    // Standard Bases
    add4!(b'A', [0,0,0,0]); add4!(b'C', [1,1,1,1]); add4!(b'G', [2,2,2,2]); add4!(b'T', [3,3,3,3]); add4!(b'U', [3,3,3,3]);
    add4!(b'a', [0,0,0,0]); add4!(b'c', [1,1,1,1]); add4!(b'g', [2,2,2,2]); add4!(b't', [3,3,3,3]); add4!(b'u', [3,3,3,3]);

    // Weak,Strong,Amino,Ketone,Purine,Pyrimide
    add4!(b'W', [0,3,0,3]); add4!(b'S', [1,2,1,2]); add4!(b'M', [0,1,0,1]); add4!(b'K', [2,3,2,3]); add4!(b'R', [0,2,0,2]); add4!(b'Y', [1,3,1,3]);
    add4!(b'w', [0,3,0,3]); add4!(b's', [1,2,1,2]); add4!(b'm', [0,1,0,1]); add4!(b'k', [2,3,2,3]); add4!(b'r', [0,2,0,2]); add4!(b'y', [1,3,1,3]);

    // not ACGT, assuming they appear in similar proportions
    // each letter has a bias but they correspond to their complement B<>V, D<>H
    add4!(b'B', [1,2,3,1]); add4!(b'D', [0,2,3,0]); add4!(b'H', [0,1,3,3]); add4!(b'V', [0,1,2,2]);
    add4!(b'b', [1,2,3,1]); add4!(b'd', [0,2,3,0]); add4!(b'h', [0,1,3,3]); add4!(b'v', [0,1,2,2]);

    // literally no information and already filled up
    //add4!(b'N', [0,1,2,3]); add4!(b'n', [0,1,2,3]);
    look_up_table
}


pub fn translate_to_dna4(chunk: &mut [u8], seed: &mut u64){
    for byte in chunk.iter_mut() {
        //xorshift64 for pseudo randomness
        //from https://en.wikipedia.org/wiki/Xorshift
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed >> 17;

        let i = (*seed & 3) as usize;
        *byte = TO_DNA4[*byte as usize][i];
    }
}

//hash to compress the size
pub fn hash_to_4mer(chunk: &[u8]) -> vec[u8] {

    let hashed_chunk = Vec::with_capacity((chunk.len() / 4) as usize);

    for &base in chunk.chunks(4); {
        kmer: u8 = 0;
        for &base in chunk.iter() {
            kmer = (kmer<<2) | base;
        }
        hashed_chunk.push(kmer)
    }
    hashed_chunk
}

// for the last kmer, to still be alignable with longer structures
pub fn split_at(kmer: u8, pos: u8) -> (u8, u8) {
    if pos > 7 { return 0, kmer; }
    let mask = (1 << pos) - 1
    let front: u8 = kmer & !mask;
    let back: u8 = kmer & mask;
    (front, back)
}

pub fn rehash_from_4mer(chunk: vec[u8]) -> vec[u8] {
    let rehashed_chunk = Vec::with_capacity(chunk.len()*4);
    for kmer in chunk.iter() {
        for i in 0..4{
            rehashed_chunk.push(get_base_at_x(kmer, 0));
        }
    }
    rehashed_chunk
};

pub fn get_base_at_x(kmer: u8, pos: u8) -> u8 {
    let shift = (3 - pos) * 2;
    (kmer >> shift) & 3
}

pub fn give_pos(size_of_alphabet: u8) -> u8 {
    match (size_of_alphabet) {
        0..=2 => return 1,
        3..=4 => return 2,
        5..=8 => return 3,
        9..=16 => return 4,
        17..=32 => return 5,
        33..=64 => return 6,
        65..=128 => return 7,
        _ => return 8,
    }
}

pub fn give_k(size_of_alphabet: u8) -> u8 {
    let i: u8 = 1;
    let bit_size = give_pos(size_alphabet);
    while (bit_size * i <= 8) {
        i += 1;
    }
    i -= 1
}

pub fn chunked_reading_t()