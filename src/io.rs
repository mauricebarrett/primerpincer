use crate::cli::Algorithm;
use crate::preparing_input::PrimerSet;
use crate::primer_trim::PrimerTrimmer;
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use paraseq::fastx;
use paraseq::prelude::*;
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::Path;

/// Process FASTQ file with parallel primer trimming
/// Handles both compressed (.gz) and uncompressed FASTQ files
pub fn process_fastq(
    input_path: &str,
    output_path: &str,
    forward_primer: &str,
    reverse_primer: &str,
    search_length: usize,
    algorithm: Algorithm,
    error_rate: f64,
    min_overlap: usize,
    threads: usize,
) -> anyhow::Result<()> {
    // Create output directory if it doesn't exist
    if let Some(parent) = Path::new(output_path).parent() {
        std::fs::create_dir_all(parent)?;
    }

    // Open output FASTQ file and wrap with gzip encoder if needed
    let output_file = File::create(output_path)?;
    let output: Box<dyn Write + Send> = if output_path.ends_with(".gz") {
        Box::new(GzEncoder::new(output_file, Compression::fast()))
    } else {
        Box::new(output_file)
    };

    let primers = PrimerSet::new(forward_primer.to_string(), reverse_primer.to_string());

    eprintln!("🎬 Starting FASTQ processing with {} threads", threads);

    // Create processor
    let mut processor = PrimerTrimmer::new(
        output,
        primers,
        search_length,
        algorithm,
        error_rate,
        min_overlap,
    )?;

    // Open input FASTQ file and handle decompression
    let input_file = File::open(input_path)?;

    // Process based on file format
    if input_path.ends_with(".gz") {
        let decoder = GzDecoder::new(input_file);
        // Use 256KB buffer for better decompression efficiency and fewer syscalls
        let buffered = BufReader::with_capacity(256 * 1024, decoder);
        let reader = fastx::Reader::new(buffered)?;
        reader.process_parallel(&mut processor, threads)?;
    } else {
        // Use 256KB buffer for better I/O throughput and fewer syscalls
        let buffered = BufReader::with_capacity(256 * 1024, input_file);
        let reader = fastx::Reader::new(buffered)?;
        reader.process_parallel(&mut processor, threads)?;
    }

    Ok(())
}
