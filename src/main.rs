mod amplicon_processor;
mod cli;
mod compression;
mod io;
mod preparing_input;
mod primer_search;
mod search_algos;
mod sinks;

use clap::Parser;
use cli::Cli;
use io::process_fastq;

fn main() -> anyhow::Result<()> {
    // Parse command-line arguments using Clap
    // Note: Clap automatically handles --help and -h flags
    // To programmatically print help, use: Cli::command().print_help()?;
    let args = Cli::parse();

    eprintln!("🦀 PrimerPincer — starting primer trimming");
    eprintln!("Input FASTQ: {}", args.input);
    eprintln!("Output FASTQ: {}", args.output);
    eprintln!("Forward primer: {}", args.forward_primer);
    eprintln!("Reverse primer: {}", args.reverse_primer);
    eprintln!("Algorithm: {:?}", args.algorithm);
    eprintln!("Maximum error rate: {:.1}%", args.error_rate * 100.0);
    eprintln!("Window size: {}", args.window_size);
    eprintln!("Minimum overlap: {}", args.overlap);
    eprintln!("Threads: {}", args.threads);
    eprintln!("Output compression: {:?}", args.compression);
    eprintln!(
        "Min length: {}",
        args.min_length
            .map(|v| v.to_string())
            .unwrap_or_else(|| "not set".to_string())
    );
    eprintln!(
        "Max length: {}",
        args.max_length
            .map(|v| v.to_string())
            .unwrap_or_else(|| "not set".to_string())
    );
    eprintln!(
        "Min average quality: {}",
        args.min_average_quality
            .map(|v| v.to_string())
            .unwrap_or_else(|| "not set".to_string())
    );
    eprintln!("Version: {}", env!("CARGO_PKG_VERSION"));

    // Process FASTQ file with primer trimming
    process_fastq(
        &args.input,
        &args.output,
        &args.forward_primer,
        &args.reverse_primer,
        args.window_size,
        args.algorithm,
        args.error_rate,
        args.overlap,
        args.compression,
        args.threads,
        args.min_length,
        args.max_length,
        args.min_average_quality,
    )?;

    eprintln!("✅ Primer trimming complete!");

    Ok(())
}
