use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct SimulateArgs {
    #[arg(short, long)]
    pub input: String,

    #[arg(short, long, default_value_t = String::from("seq.fasta"))]
    pub output: String,

    #[arg(short, long)]
    pub lens: Vec<usize>,

    #[arg(long, default_value_t = 3)]
    pub order: usize,

    #[arg(long, default_value_t = 42)]
    pub seed: u64,

    #[arg(short, long, default_value_t = false)]
    pub verbose: bool,
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct MutateArgs {
    #[arg(short, long)]
    pub input: String,

    #[arg(short, long, default_value_t = String::from("seq.fasta"))]
    pub output: String,

    #[arg(long, default_value_t = 42)]
    pub seed: u64,

    #[arg(short, long, default_value_t = 0.1)]
    pub error: f64,

    #[arg(short, long, default_value_t = false)]
    pub verbose: bool,

    #[arg(short, long, default_value_t = false)]
    pub debug: bool,
}


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct RandomArgs {
    #[arg(short, long, default_value_t = String::from(""))]
    pub input: String,

    #[arg(short, long, default_value_t = String::from("rand_seq.fasta"))]
    pub output: String,

    #[arg(short, long, num_args = 1.., required = true)]
    pub lens: Vec<usize>,

    #[arg(long, default_value_t = 42)]
    pub seed: u64,

    #[arg(short, long, default_value_t = false)]
    pub verbose: bool,

    // The alphabet to use, if --input is used, input will be filtered
    #[arg(short, long, value_delimiter = ',', num_args = 1.., default_values_t = vec!['A', 'C', 'G', 'T'])]
    pub alphabet: Vec<char>,

    // Weights for the alphabet, must match alphabet length.
    #[arg(short, long, value_delimiter = ',', num_args = 1.., default_values_t = vec![25, 25, 25, 25])]
    pub distribution: Vec<usize>,

    #[arg(long, default_value_t = 1)]
    pub threads: u64,
}

#[derive(Debug, Parser)] // requires `derive` feature
#[command(name = "markov_genome")]
#[command(about = "Markov chain sequence simulation", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    #[command(arg_required_else_help = true)]
    Simulate (SimulateArgs),
    Mutate (MutateArgs),
    Random (RandomArgs)
}