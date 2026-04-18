use std::env::args;
use std::fs::{metadata, File};
use std::io::SeekFrom;
use std::io::{Read, Seek};
use std::thread;

#[derive(Debug, Default)]
enum kmer_trie_leaf {
    counts: [usize],
    previous_kmer_len: u8,
    is_leaf: bool,
}

impl kmer_trie_leaf {
    fn new(len_b: u8, len_previous_kmer: u8) -> Self {
        Self {
            counts: [0usize; (1<<len_b) as usize],
            len_previous_kmer,
            is_leaf : true,
        }
    }
}

#[derive(Debug, Default)]
enum kmer_trie_node {
    Node(Box<kmer_trie>),
    Leaf(Box<kmer_trie_leaf>),
}

#[derive(Debug, Default)]
struct kmer_trie {
    order: u8,
    depth: u8,
    len_b: u8, //size_alphabet in bit (2 in our case)
    sub_k: u8, //how many hashed chars can fit in one u8
    is_leaf: bool,
    children: HashMap<u8, kmer_trie_node>,
}

impl kmer_trie {
    fn new(order: u8, depth: u8, size_alphabet: u8) -> Self {
        Self {
            order,
            depth,
            len_b: give_pos(size_alphabet),
            sub_k: give_k(size_alphabet), //(bit_size:kmer) 1:8 2:4 3:2 4:2 5:1 6:1 7:1 8:1 9geht nicht :D
            leaf: False,
            children: HashMap::new(),
        }
    }

    fn new_node(&self) -> kmer_trie_node {
        kmer_trie_node::Node(Box::new(kmer_trie::new(
            self.order,
            self.depth + self.sub_k,
            self.len_b,
        )))
    }

    fn new_leaf(&self, len_last_kmer: u8) -> Vec<usize> {
        kmer_trie_node::Leaf(Box::new(kmer_trie_leaf(
            self.len_b,
            len_last_kmer,
        )))
    }

    //the hashing assumes trailing zeroes
    fn insert_kmer(&mut self, kmer: &[u8])
    {
        if self.depth > self.order || kmer.is_empty() { return; }
        // base == alphabet size in bits, 
        let to_fill = self.order - self.depth - 1; //adding the 1 here makes the rest easier
        if to_fill > self.sub_k{ //more than one base is still in kmer[1]
            let node = self.children
                .entry(kmer[0])
                .or_insert(self.new_node());
            *node.insert_kmer(&kmer[1..]);
        } else {
            let border = 8 - self.to_fill * self.len_b
            let (front, back) = split_at(kmer[1], border)
            back >>= (border - self.len_b)
            let node = self.children
                .entry(front)
                .or_insert(self.new_leaf(self.len_b));
            *node.counts[back as usize] += 1
    }

    //count all appearences in this branch
    fn get_sum(&self, count: &usize) {
        for (_, node) in self.children.iter() {
            match node {
                kmer_trie_node::Node(trie) => trie.get_sum(count),
                kmer_trie_node::Leaf(leaf) => count + leaf.iter().sum(),
            }
        }
    }

    fn get_transitions(&self, kmer: &[u8], counts: &mut [usize], k: u8){
        if kmer.is_empty() || k == 0 { 
            let c: usize = 0;
            self.get_sum(c);
            for count in counts.iter() {
                counts += c.strict_div(counts.len())
            }
        }

        if kmer.len() == 1 && k < self.sub_k {
            let (last_kmer, _) = split_at(kmer[0], 8 - k * self.len_b)
            for (key, node) in self.children.iter() {
                match node {
                    kmer_trie_node::Node(trie) => trie.get_sum(count),
                    kmer_trie_node::Leaf(leaf) => count + leaf.iter().sum(),
                }
        }
            if let Some(node) = self.children.get(kmer[0]) {
                match node {
                    kmer_trie_node::Node(trie) => trie.get_sum(count),
                    kmer_trie_node::Leaf(leaf) => count + leaf.iter().sum(),
                    _ => return,
                }

        }
        if kmer.len() > 1 && k > self.depth + self.sub_k + self.len_b {
            if !self.children.contains_key(kmer[0]) { return; }
            if let Some(node) = self.children.get(kmer[0]) {
            match node {
                kmer_trie_node::Node(next_tree) => next_tree.get_transitions(kmer[1..], counts),
                _ => return,
            }
        }
        if kmer.len() == 1 {
            let rest = self.order - k 
        }&& k < 8 {
            let (front, _) = split_at(kmer[1], 8 - k)
        }
    }



    fn get_leafs(&self, kmer: u8, k, counts: &mut [usize]){

    }

    fn get_(&self, kmer: &[u8], k: u8, counts: &mut [usize]) {
        if kmer.len() == 0 { return self.get_count(); }
        if kmer.len() == 1 { 
            return self.get_leafs(kmer[0], k, counts: &mut [usize]); 
        }
        if let Some(node) = self.children.get(kmer[0]) {
            match node {
                kmer_trie_node::Node(next_tree) => next_tree.get_counts(kmer[1..], counts),
                kmer_trie_node::Leaf(leaf) => 
            }
        }
        if kmer.len() == 1 { 

        }
        if let Some(node) = self.children.get(&target_kmer) {
            match node {
                kmer_trie_node::Node(next_tree) => next_tree.extract_path(target_kmer, path),
                kmer_trie_node::Leaf(_) => return, // Found the end
            }
        }
    }


}


 // //need to specify the length of of the kmer asked for
    // fn get_count(&self, kmer&[u8], k: u8, count: &usize) {
    //     if kmer.is_empty() ||
    //     if k < self.depth + self.sub_k && !kmer.is_empty() { 

    //     }
    //     if k > self.depth ||  { return self.get_sum(count); }
    //     if kmer.len() == 0 { return self.get_sum(count); }
    //     if kmer.len() == 1 {
    //         let 
    //         if k == self.order {
    //             let first, last = split_at(kmer[0], );
    //             if !self.children.contains_key(first) { return; }
    //             let Some(leaf) = self.children.get(first);
    //             if !leaf.is_leaf { return; }
    //             count + *leaf[last as usize];
    //         } else {
    //             let pos = k - self.depth;
    //             let first, last = split_at(kmer[0], pos);
    //             let Some(node) = self.children.get(first);
    //                 match node {
    //                     kmer_trie_node::Leaf(leaf) => count + *leaf[last as usize],
    //                     _ => return,
    //                 }
    //         }
    //         let pos = k - self.depth;
    //         for (key, node) in self.children.iter() {
    //             let pos = k-self.depth;
    //             if 
    //             match node {
    //                 kmer_trie_node::Node(trie) => trie.get_sum(count),
    //                 kmer_trie_node::Leaf(leaf) => count + leaf.iter().sum(),
    //             }
    //     }
    //     if k > self.order { return; }
    // 
    //     }
    //     if let Some(node) = self.children.get(kmer[0]) {
    //         match node {
    //             kmer_trie_node::Node(next_tree) => next_tree.get_counts(kmer[1..], counts),
    //             kmer_trie_node::Leaf(leaf) => 
    //         }
    //     }
    //     for (_, node) in self.children.iter() {
    //         match node {
    //             kmer_trie_node::Node(trie) => count + trie.get_sum(),
    //             kmer_trie_node::Leaf(leaf) => count + leaf.iter().sum(),
    //         }
    //     }
    // }
fn main() {
    // Collect command-line arguments into a vector of strings
    let args: Vec<String> = args().collect();

    // Ensure the path to the file is provided as the second argument
    // Extracts a reference to the path string.
    let path: &String = &args.get(1).expect("File path not specified.");

    // Retrieve the metadata for the open file to get the length in `usize` bytes
    // If the metadata retrieval fails, program will panic
    let file_length: usize = metadata(&path)
        .expect("Unable to query file details")
        .len()
        .try_into()
        .expect("Cannot convert file length");

    // Size of each block to be read from the file (16MB)
    const BLOCK_SIZE: usize = 16_777_216;

    // Number of threads to be used for reading the file
    const THREADS: usize = 4;

    // Determine the portion of the file each thread will handle
    let division: usize = (file_length + THREADS - 1) / THREADS;

    // Use scoped threads to ensure all threads are joined before the main thread exits
    thread::scope(|scope: &thread::Scope| {
        // Create `THREADS` number of threads
        for i in 0..THREADS {
            scope.spawn(move || {
                // Open a file handle per thread
                // Ensuring each thread works on its own file descriptor
                let mut thread_file: File = File::open(&path).expect("Unable to open file");

                // Initialize a 16MB buffer to hold data read by this thread
                let mut contents: Vec<u8> = vec![0_u8; BLOCK_SIZE];

                // Counter for the number of bytes read in each operation
                // Init'd as 1, as zero would indicate an EOF
                let mut read_length: usize = 1;

                // Counter to keep track of the total nuber of bytes read by this thread
                let mut read_total: usize = 0;

                // Determing the offset in the file where this thread should start reading.
                // Each thred reads a different portion of the file
                let offset: u64 = (i * division) as u64;

                // Seek to the starting position in the file for this thread
                // If it couldn't start from there, the program panics
                thread_file
                    .seek(SeekFrom::Start(offset))
                    .expect("Couldn't seek to position in file");

                // Read data in iterations until the thread's portion is read or EOF is reached
                // i.e there is no more bytes to be read from the file
                while (read_total < division) && (read_length != 0) {
                    // Adjust the contents buffer size if the remaining bytes to be read are less than the block size
                    if read_total + BLOCK_SIZE > division {
                        // Reduce the vector to ensure we don't read beyond the division
                        // If this doesn't happen, data for the next thread will be read causing offset seek errors
                        contents.truncate(division - read_total);
                    }

                    // Read data into the contents buffer
                    // `read_length` is updated with the number of bytes read
                    read_length = thread_file.read(&mut contents).expect("Couldn't read file");
                    read_total += read_length;
                }

                // Verify the entire file was read correctly by visualizing total number of bytes read from the file
                // by each thread and the original file length
                println!("Thread {i}: Total bytes read: {read_total} bytes || Expected: {division} bytes");
            });
        }
    });
