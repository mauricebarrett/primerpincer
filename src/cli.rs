use clap::Parser;

/// 🦀 PrimerPincer — A command-line tool for the rapid trimming of primers from long read amplicons
#[derive(Parser, Debug)]
#[command(
    author = "Maurice Barrett",
    version = "alpha",
    about = "A command-line tool for the rapid trimming of primers from long read amplicons",
    long_about = None
)]
pub struct Cli {
    /// Input FASTQ file
    #[arg(short = 'i', long = "input", value_name = "FILE", required = true)]
    pub input: String,

    /// Output FASTQ file
    #[arg(short = 'o', long = "output", value_name = "FILE", required = true)]
    pub output: String,

    /// Forward primer sequence (5' to 3' orientation)
    #[arg(
        short = 'f',
        long = "forward",
        value_name = "SEQUENCE",
        required = true
    )]
    pub forward_primer: String,

    /// Reverse primer sequence (5' to 3' orientation)
    #[arg(
        short = 'r',
        long = "reverse",
        value_name = "SEQUENCE",
        required = true
    )]
    pub reverse_primer: String,

    /// Maximum allowed mismatches in primer matching
    #[arg(
        short = 'm',
        long = "mismatches",
        value_name = "INT",
        default_value_t = 2
    )]
    pub max_mismatches: usize,

    /// Number of threads to use
    #[arg(short = 't', long = "threads", value_name = "INT", default_value_t = 4)]
    pub threads: usize,
}
