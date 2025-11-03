use clap::Parser;
mod cli;
mod io;
mod primer_search;
mod primer_trim;

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
    eprintln!("Allowed edit distance: {}", args.edit_distance);
    eprintln!("Search length: {}", args.search_length);
    eprintln!("Minimum overlap: {}", args.overlap);
    eprintln!("Threads: {}", args.threads);

    // Process FASTQ file with primer trimming
    process_fastq(
        &args.input,
        &args.output,
        &args.forward_primer,
        &args.reverse_primer,
        args.search_length,
        args.algorithm,
        args.edit_distance,
        args.overlap,
        args.threads,
    )?;

    eprintln!("✅ Primer trimming complete!");

    Ok(())
}
