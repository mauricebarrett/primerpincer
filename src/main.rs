use clap::Parser;
mod cli;
mod io;
mod primer_search;
mod primer_trim;

use cli::Cli;
use io::process_fastq;
use primer_search::SearchMethod;

fn main() -> anyhow::Result<()> {
    // Parse command-line arguments using Clap
    let args = Cli::parse();

    // Parse algorithm string
    let algorithm = match args.algorithm.to_lowercase().as_str() {
        "alignment" | "align" | "cutadapt" | "cut" => SearchMethod::Alignment,
        "hamming" | "ham" | "minibar" | "min" => SearchMethod::Hamming,
        _ => {
            eprintln!("Error: Unknown algorithm '{}'. Use 'alignment' or 'hamming'", args.algorithm);
            std::process::exit(1);
        }
    };

    eprintln!("🦀 PrimerPincer — starting primer trimming");
    eprintln!("Input FASTQ: {}", args.input);
    eprintln!("Output FASTQ: {}", args.output);
    eprintln!("Forward primer: {}", args.forward_primer);
    eprintln!("Reverse primer: {}", args.reverse_primer);
    eprintln!("Allowed mismatches: {}", args.max_mismatches);
    eprintln!("Search length: {}", args.search_length);
    eprintln!("Algorithm: {:?}", algorithm);
    eprintln!("Threads: {}", args.threads);

    // Process FASTQ file with primer trimming
    process_fastq(
        &args.input,
        &args.output,
        &args.forward_primer,
        &args.reverse_primer,
        args.search_length,
        args.max_mismatches,
        args.threads,
        algorithm,
    )?;

    eprintln!("✅ Primer trimming complete!");

    Ok(())
}
