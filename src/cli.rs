use clap::{Parser, ValueEnum};

/// Algorithm selection for primer matching
#[derive(ValueEnum, Clone, Debug, Copy)]
pub enum Algorithm {
    /// Use standard Levenshtein edit distance calculation
    Levenshtein,
    /// Use Myers algorithm with IUPAC support
    Myers,
    /// Use local pairwise alignment with IUPAC-aware scoring
    Local,
    /// Use Sassy SIMD-accelerated search (fastest, requires AVX2/NEON)
    Sassy,
}

impl Default for Algorithm {
    fn default() -> Self {
        Algorithm::Sassy
    }
}

/// 🦀 PrimerPincer — A command-line tool for the rapid trimming of primers from long read amplicons
#[derive(Parser, Debug)]
#[command(
    author = "Maurice Barrett",
    version = "alpha",
    about = "PrimerPincer - a CLI primer trimming tool for long-read sequencing data",
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

    /// Algorithm selection for primer matching
    #[arg(
        short = 'a',
        long = "algorithm",
        value_enum,
        help = "Algorithm to use for primer matching: levenshtein, myers, local, or sassy",
        default_value_t = Algorithm::Sassy
    )]
    pub algorithm: Algorithm,

    /// Maximum allowed edit distance in primer matching
    #[arg(
        short = 'e',
        long = "edit-distance",
        value_name = "INT",
        help = "Maximum allowed edit distance in primer matching",
        default_value_t = 3
    )]
    pub edit_distance: usize,

    /// Maximum number of mismatches allowed in local alignment
    #[arg(
        short = 'm',
        long = "max-mismatch",
        value_name = "INT",
        help = "Maximum mismatches allowed when using the local alignment algorithm",
        default_value_t = 3
    )]
    pub max_mismatch: usize,

    #[arg(
        short = 'l',
        long = "search-length",
        value_name = "INT",
        help = "Length to search for primer at start and end of sequence",
        default_value_t = 100
    )]
    pub search_length: usize,

    /// Minimum overlap length (cutadapt -O)
    #[arg(
        short = 'O',
        long = "overlap",
        value_name = "MINLENGTH",
        help = "Minimum overlap length. Require MINLENGTH bases of the primer to match (like cutadapt -O)",
        default_value_t = 6
    )]
    pub overlap: usize,

    /// Number of threads to use
    #[arg(short = 't', long = "threads", value_name = "INT", default_value_t = 4)]
    pub threads: usize,
}
