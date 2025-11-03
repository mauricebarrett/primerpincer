use paraseq::prelude::*;
use paraseq::fastx;
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::Path;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use crate::primer_trim::PrimerTrimmer;

/// Process FASTQ file with parallel primer trimming
/// Handles both compressed (.gz) and uncompressed FASTQ files
pub fn process_fastq(
    input_path: &str,
    output_path: &str,
    forward_primer: &str,
    reverse_primer: &str,
    search_length: usize,
    max_mismatches: usize,
    min_overlap: usize,
    threads: usize,
) -> anyhow::Result<()> {
    // Create output directory if it doesn't exist
    if let Some(parent) = Path::new(output_path).parent() {
        std::fs::create_dir_all(parent)?;
    }
    
    // Open input FASTQ file and handle compression based on extension
    let input_file = File::open(input_path)?;
    
    // Open output FASTQ file - compress if extension is .gz
    let output_file = File::create(output_path)?;
    let output: Box<dyn Write + Send> = if output_path.ends_with(".gz") {
        // Gzip compress the output
        Box::new(GzEncoder::new(output_file, Compression::default()))
    } else {
        // Plain text output
        Box::new(output_file)
    };
    
    // Create processor
    let mut processor = PrimerTrimmer::new(
        output,
        forward_primer.to_string(),
        reverse_primer.to_string(),
        search_length,
        max_mismatches,
        min_overlap,
    )?;
    
    // Process in parallel - handle compressed and uncompressed files
    if input_path.ends_with(".gz") {
        // Gzip compressed file
        let decoder = GzDecoder::new(input_file);
        let reader = fastx::Reader::new(BufReader::new(decoder))?;
        reader.process_parallel(&mut processor, threads)?;
    } else {
        // Uncompressed file
        let reader = fastx::Reader::new(BufReader::new(input_file))?;
        reader.process_parallel(&mut processor, threads)?;
    }
    
    Ok(())
}
