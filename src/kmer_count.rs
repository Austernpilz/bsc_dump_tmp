use std::collections::HashMap;

#[derive(Debug, Default)]
struct KmerTrieLeaf {
    counts: [usize; 4],
}

#[derive(Default)]
enum KmerTrieNode {
    Trie(Box<KmerTrie>),
    Leaf(KmerTrieLeaf),
    #[default]
    Empty
}

#[derive(Default)]
struct KmerTrie {
    order: u8,
    depth: u8,
    children: HashMap<u8, KmerTrieNode>,
}

impl KmerTrieLeaf {
    fn new() -> Self {
        Self {counts: [1,1,1,1]}
    }

    fn add_node(&mut self, node: &KmerTrieNode) {
        match node {
            KmerTrieNode::Leaf(leaf) => self.add_leaf(leaf.counts),
            KmerTrieNode::Trie(trie) => self.add_trie(trie),
            KmerTrieNode::Empty => return,
        }
    }

    fn add_trie(&mut self, trie: &KmerTrie) {
        self.add_leaf(trie.give_first_base_transition());
    }

    fn add_leaf(&mut self, count: [usize; 4]) {
        self.counts[0] += count[0] - 1;
        self.counts[1] += count[1] - 1;
        self.counts[2] += count[2] - 1;
        self.counts[3] += count[3] - 1;
        // for (base, c) in self.counts.iter_mut().zip(count) {
        //     *base += c - 1;
        // } // if leaf extend to other alphabets
    }

    fn get_sum(&self) -> usize {
        self.counts[0] + self.counts[1] + self.counts[2] + self.counts[3] - 4 
        //self.counts.iter().sum::<usize>() - 4
    }
}

impl KmerTrieNode {
    fn get_sum(&self) -> usize {
        return match self {
            KmerTrieNode::Leaf(leaf) => leaf.get_sum(),
            KmerTrieNode::Trie(trie) => trie.get_sum(),
            KmerTrieNode::Empty => 0,
        }
    }

    fn give_first_base_transition(&self) -> [usize; 4] {
        return match self {
            KmerTrieNode::Leaf(leaf) => leaf.counts,
            KmerTrieNode::Trie(trie) => trie.give_first_base_transition(),
            KmerTrieNode::Empty => [1,1,1,1],
        }
    }

    fn give_transitions(&self, kmer: &[u8]) -> [usize; 4] {
        return match self {
            KmerTrieNode::Leaf(leaf) => leaf.counts,
            KmerTrieNode::Trie(trie) => trie.give_transitions(&kmer),
            KmerTrieNode::Empty => [1,1,1,1],
        }
    }

    fn merge(&mut self, trie: &KmerTrie) {
        match self {
            KmerTrieNode::Leaf(leaf) => leaf.add_trie(trie),
            KmerTrieNode::Trie(my_trie) => my_trie.merge(trie),
            KmerTrieNode::Empty => return, 
            //would be really bad if it lands here, should insert before :D
        }
    }


    fn is_leaf(&self) -> bool {
        match self {
            KmerTrieNode::Leaf(_) => true,
            _ => false,
        }
    }

    fn is_trie(&self) -> bool {
        match self {
            KmerTrieNode::Trie(_) => true,
            _ => false,
        }
    }
}

fn new_trie(order: u8, depth: u8) -> KmerTrieNode {
        KmerTrieNode::Trie(Box::new(KmerTrie::new(order, depth + 4,)))
    }

fn new_leaf() -> KmerTrieNode {
    KmerTrieNode::Leaf(KmerTrieLeaf::new())
}

impl KmerTrie {
    fn new(order: u8, depth: u8) -> Self {
        Self {order, depth, children: HashMap::new()}
    }

    fn insert_kmer_counting(&mut self, kmer: &[u8])
    {
        //counting_step shoulnd't reach kmer.len() == 1, because i use k.next_multiple_of(4)+1 so only order 0 gives me 1
        if kmer.len() < 2 { return; } //for safeguarding
        if kmer.len() == 2 {
            if let KmerTrieNode::Leaf(leaf) = self.children.entry(kmer[0]).or_insert_with(|| new_leaf()) {
                leaf.counts[get3(kmer[0]) as usize] += 1;
            }
        } else if let KmerTrieNode::Trie(trie) = self.children.entry(kmer[0]).or_insert_with( || new_trie(self.order, self.depth)) {
            trie.insert_kmer_counting(&kmer[1..]);
        }
    }

    fn count_kmer(&mut self, kmer: &[u8], block_length: usize){
        //chunks didn't work, i don't know
        for i in (block_length-1..kmer.len()).step_by(block_length) {
            self.insert_kmer_counting(&kmer[i-block_length..i]);
        }
    }

    //count all appearences in this branch
    fn get_sum(&self) -> usize {
        // let mut count: usize = 0;
        // for (_, node_mer) in self.children.iter() {
        //     count += node_mer.get_sum();
        // //     match node_mer {
        // //         KmerTrieNode::Trie(trie) => count += trie.get_sum(),
        // //         KmerTrieNode::Leaf(leaf) => count += leaf.get_sum(), //leaf.get_sum() calcs -4 to reverse the over counting
        // //         KmerTrieNode::Empty => continue,
        // //     };
        // }
        // count
        self.children.values().fold(0, |acc, node| acc + node.get_sum())
    }

    fn give_first_base_transition(&self) -> [usize; 4] {
        //count every first base on this branch
        let mut transition = [1,1,1,1];
        for (key_mer, node_mer) in self.children.iter() {
            transition[get0(*key_mer)as usize] += node_mer.get_sum();
            // match (node_mer) {
            //     KmerTrieNode::Trie(trie) => {
            //         transition[get0(*key_mer)as usize] += trie.get_sum();
            //     },
            //     KmerTrieNode::Leaf(leaf) => {
            //         transition[get0(*key_mer)as usize] += leaf.get_sum();
            //     },
            // }
        }
        transition
        // self.children
        //     .iter()
        //     .fold(
        //         [1,1,1,1], |acc, (key_mer, node_mer)| acc[] + node.get_sum())
    }

    // fn get_matches_at_pos(&self, 4kmer: u8, pos: u8) -> (bool , [u8; 4]) {
    //     //pos means bit poisition, counting from right to left (76543210)
    //     let a: u8 = 0b0000_0000; // is here for reference, but used as the check
    //     let c: u8 = !(0b0000_0001 << shift);
    //     let g: u8 = !(0b0000_0010 << shift);
    //     let t: u8 = !(0b0000_0011 << shift);

    //     kmer &= t; //this is forcing this part to be A
    //     //make sure we get 4 disdinct solutions otherwise rust panics or we need to be unsafe
    //     let potential_matches = self.children.get_disjoint_mut([kmer, kmer | c, kmer | g, kmer | c | g]);

    // }

    fn find_last_matching_pos(&self, kmer: u8) -> u8 {
        //obviously i din't found the the full match :D 
        let mask3: u8 = 0b1111_1100; //3base match
        let mask2: u8 = 0b1111_0000; //2base match
        let mask1: u8 = 0b1100_0000; //1base match
        let mut pos2: bool = false;
        let mut pos1: bool = false;
        //would be faster with an array, but design choices in the beginning :(
        // beceause we can search in a trie, there must be something 
        for key_mer in self.children.keys() {
            if (key_mer & mask3) == (kmer & mask3) {
                return 2;
            }
            pos2 |= (key_mer & mask2) == (kmer & mask2); //look for pos 3 at the same time
            pos1 |= (key_mer & mask1) == (kmer & mask1);
        }
        if pos2 { return 4; }
        if pos1 { return 6; }
        8
    }

    fn give_last_transitions(&self, kmer: u8) -> [usize; 4] {
        //I didn't found a leaf so now i gotta find the others that fit
        let pos = self.find_last_matching_pos(kmer);
        if pos == 8 { return self.give_first_base_transition(); }

        let mut transitions: [usize; 4] = [1,1,1,1];
        let mask: u8 = 0b1111_1111 << pos;
        //here as well, the hashmap doesn't help
        for (key_mer, node_mer) in self.children.iter() {
            if (key_mer & mask) == (kmer & mask) {
                transitions[get3(key_mer >> pos) as usize] += node_mer.get_sum();
            }
        }
        transitions
    }

    fn give_transitions(&self, kmer: &[u8]) -> [usize; 4] {
        if kmer.len() == 0 {  return self.give_first_base_transition() }
        //we assume that the kmer was asked for in the right block size, with trailing zeroes
        if let Some(node) = self.children.get(&kmer[0]) {
            if kmer.len() == 1 {
                return node.give_first_base_transition();
            } else {
                return node.give_transitions(&kmer[1..])
            }
        } else {
            return self.give_last_transitions(kmer[0])
        }
        //we give 1,1,1,1 back if we can't find anything
    }

    //only in kmer_counting at the moment
    // merge actually absorbs
    fn merge(&mut self, other_trie: &KmerTrie) {
        if (self.order != other_trie.order) || (self.depth != other_trie.depth) { return; } //needs to be expanded for grammars
        for (key_mer, node_mer) in other_trie.children.iter() {
            match node_mer {
                KmerTrieNode::Trie(trie) => {
                    if let KmerTrieNode::Trie(my_trie) = self.children.entry(*key_mer).or_insert_with(|| new_leaf()) {
                        my_trie.merge(trie);
                    } //I
                    // self.children
                    // .entry(*key_mer)
                    // .and_modify(|node|node.merge(trie))
                    // .or_insert_with( || new_trie(self.order, self.depth + 4));
                },
                    //other trie loses values if i only have a leaf at this node
                KmerTrieNode::Leaf(leaf) => {
                    if let KmerTrieNode::Leaf(my_leaf) = self.children.entry(*key_mer).or_insert_with(|| new_leaf()) {
                        my_leaf.add_leaf(leaf.counts);
                    } //I have no answer for the situation where I have a branch at this point
                },
                KmerTrieNode::Empty => continue,
                //here we still only work when both tries are completly the same, 
                //with a grammar, or orther more comples structures, we need to look if we have a node here
            }
        }
    }

    fn collapse(&mut self) {
        //I collapse every node structure under me, so that my depth becomes full of leafs
        let extracted_tries: HashMap<u8, [usize; 4]> = self.children.extract_if(|_, n| n.is_trie())
            .map(|(k, node)| (k, node.give_first_base_transition()))
            .collect();
        for (key_mer, transition) in extracted_tries.iter() {
            if let KmerTrieNode::Leaf(leaf) = self.children.entry(*key_mer).or_insert_with(|| new_leaf()) {
                leaf.counts = *transition;
            }
        }
    }

    fn reduce_to_order(&mut self, k: u8) {
        self.order = self.depth + k;
        if k <= 5 {
            self.collapse();
            if k==5 { return; }
            self.depth = k-1; //last is always the leaf
            let mask: u8 = 0b1111_1100 << (8-k*2);
            let extracted_tries: HashMap<u8, KmerTrieNode> = self.children.extract_if(|k,_| (*k & mask) != *k).collect();
            for (key_mer, node_mer) in extracted_tries.iter() {
                if let KmerTrieNode::Leaf(leaf) = self.children.entry(*key_mer & mask).or_insert_with(|| new_leaf()) {
                    leaf.add_node(node_mer);
                }
            }
        } else {
            for node_mer in self.children.values.mut() {
                match node_mer {
                    KmerTrieNode::Trie(trie) => trie.reduce_to_order(k-4),
                    _ => continue,
                }
            }
        }
    }
}


//// TESTING FUNCTION FOR CORRECTNESS

// fn main() {
    
//     let mut letters = HashMap::new();
    
//     for ch in "a short treatise on fungi".chars() {
//         letters.entry(ch).and_modify(|counter| *counter += 1).or_insert(1);
//     }
//     let mut vec: Vec<u8> = vec![
//         b'A', b'A', b'A', b'A', b'A', b'A', b'A', 
//         b'3', b'B', b'c', b'A', b'O', b'A', b'i', 
//         b'G', b'C', b't', b'A', b'A', b'A', b'A', 
//         b'3', b'B', b'c', b'A', b'O', b'A', b'i', 
//         b'A', b'A', b'A', b't', b'A', b'c', b'C', 
//         b'3', b'B', b'c', b'T', b'O', b'A', b'i', 
//         b'A', b'A', b'A', b'A', b'A', b'g', b'A', 
//         b'3', b'B', b'c', b'A', b'O', b'g', b'i'
//     ];
//     let mut seed:u64 = 100;
    
//     let trie = rolling_kmer_end_hash(&mut vec[..], &mut seed, 5);

//     fn p1 (a: Vec<u8>) {
//         if let Ok(s) =  String::from_utf8(a.iter().map(|x| TO_CHAR[*x as usize]).collect()) {
//             print!("{}", s);
//         }
//     }
//     //println!("{}", trie.get_sum());
//     //println!("{}", trie.give_last_transitions(1u8)[0]);
//     print!("\n");
//     fn p(trie: &KmerTrie) {
//         for (key, node_mer) in trie.children.iter() {
//             println!("depth: {} kmer starts with",trie.depth);
//             p1(decompress_to_dna4(&[*key], 4));
//             match node_mer {
//                 KmerTrieNode::Trie(trie) => p(trie),
//                 KmerTrieNode::Leaf(leaf) => print!("leaf is here with a,c,g,t: {},{},{},{}", leaf.counts[0],leaf.counts[1],leaf.counts[2],leaf.counts[3]), //otherwise I over calculate every instance of leaf
//                 KmerTrieNode::Empty => continue,
//             }
//         }
//     }
// }