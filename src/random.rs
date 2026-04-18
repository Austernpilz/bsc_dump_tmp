use rand::distr::weighted::WeightedIndex;
use rand::rngs::StdRng;
use rand::prelude::*;
use rand::SeedableRng;
use std::fs::File;
use std::sync::mpsc; 
use std::thread;
use bio::io::fasta::{Writer, Record};
use std::io::BufWriter;
use std::io::Write;


use crate::io::{get_records, print_record, char_to_int, int_to_char, is_valid_path};
use crate::args::RandomArgs;



pub fn get_alphabet_and_distribution(
    input_path : &String, 
    alphabet : &Vec<char>, 
    weights : &Vec<usize>,
) -> (Vec<u8>, WeightedIndex<usize>) {

    let mut alphabet_u8 : Vec<u8> = Vec::new();
    let mut distribution : Vec<usize> = Vec::new();

    if !is_valid_path(input_path) {
        alphabet_u8 = alphabet.iter().map(char_to_int).collect();
        distribution = weights.clone();
    } else {
        let mut byte_count_array = [0usize; 256];

        for result in get_records(input_path.clone()) {
            let record = result.as_ref().expect("Error during fasta record parsing");
            for byte in record.seq() {
                let c = char_to_int(&mut int_to_char(byte));
                byte_count_array[c as usize] += 1;
            }
        }

        for (byte, &count) in byte_count_array.iter().enumerate() {
            let c = int_to_char(&(byte as u8));
            if count > 0 && !c.is_whitespace() {
                alphabet_u8.push(char_to_int(&c));
                distribution.push(count);
            }
        }
    }

    if alphabet_u8.len() != distribution.len() {
        eprintln!(
            "Error: Mismatch between alphabet and distribution lengths.\n\
            Alphabet has {} elements, but distribution has {} weights.",
            alphabet_u8.len(),
            distribution.len()
        );
        std::process::exit(1);
    }

    if alphabet_u8.is_empty() {
        eprintln!(
            "Error: Couldn't derive alphabet and distribution from path. \n\
            {}",
            input_path,
        );
        std::process::exit(1);
    }
    let dist = WeightedIndex::new(distribution).expect("Invalid weights");

    (alphabet_u8, dist)
} // end pub fn get_alphabet_and_distribution

pub fn generate_sequence(
    len: usize,
    alphabet: Vec<u8>,
    dist: WeightedIndex<usize>,
    seed: u64,
) -> Vec<u8> {

    let mut rng = StdRng::seed_from_u64(seed);
    let mut sequence = Vec::with_capacity(len);

    for _ in 0..len {
        let idx = dist.sample(&mut rng);
        sequence.push(alphabet[idx]);
    }

    sequence
} // end pub fn generate_sequence


pub fn create_random_fasta(args : &RandomArgs) {
    let file = File::create(&args.output).expect("Couldn't create file");
    let mut bufwriter = BufWriter::new(file);
    let mut fasta_writer = Writer::from_bufwriter(bufwriter);
    //fasta_writer.set_linewrap(Some(80));

    let (alphabet, dist) = get_alphabet_and_distribution(&args.input, &args.alphabet, &args.distribution);
    let seed = args.seed;

    for (nr, len) in args.lens.iter().enumerate() {

        let header = format!("random_generated_sequence_{}", nr);
        let desc = format!("length={}", len);
        let record = Record::with_attrs(&header, Some(&desc), &[]);
        fasta_writer.write_record(&record).expect("Header write failed");


        let (tx, rx) = mpsc::channel();
        let length = *len as u64;
        let chunk_size = length / args.threads;

        for i in 0..args.threads {
            let thread_tx = tx.clone();
            let alphabet_clone = alphabet.clone();
            let distribution_clone = dist.clone();
            let thread_seed = seed.wrapping_add(i as u64).wrapping_add(nr as u64);

            let start = i * chunk_size;
            let end = std::cmp::min(start + chunk_size, length);
            let chunk_length = end - start;

            thread::spawn(move || {
                let sequence = generate_sequence(chunk_length as usize, alphabet_clone, distribution_clone, thread_seed);
                thread_tx.send((i, sequence)).expect("Failed to send chunk");
            });
        }

        drop(tx);

        let mut received = std::collections::BTreeMap::new();
//        let inner_writer = writer.get_ref_mut();

        let mut next_to_write = 0;
       // let mut column_count = 0;

        while next_to_write < args.threads {

            if let Ok((i, sequence)) = rx.recv() {
                received.insert(i, sequence);
            }
            // write them in order, so the kmer stays the same
            while let Some(sequence) = received.remove(&next_to_write) {
                let _ = fasta_writer.write(&next_to_write.to_string(), None, &sequence);
                // for byte in sequence {
                //     writer.write_record(&[byte]);
                //     column_count += 1;
                //     if column_count == 80 {
                //         writer.write_record(b"\n");
                //         column_count = 0;
                //     }
                // }
                next_to_write += 1;
            }
        }
        //bufwriter.flush();
        let _ = fasta_writer.flush();
    }
}// end pub fn create_random_fasta
