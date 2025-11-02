use paraseq::prelude::*;
use paraseq::fastx;
use std::fs::File;
use crate::primer_trim::PrimerTrimmer;

/// Process FASTQ file with parallel primer trimming
pub fn process_fastq(
    input_path: &str,
    output_path: &str,
    forward_primer: &str,
    reverse_primer: &str,
    search_length: usize,
    max_mismatches: usize,
    threads: usize,
) -> anyhow::Result<()> {
    // Open input FASTQ file
    let input_file = File::open(input_path)?;
    let reader = fastx::Reader::new(input_file)?;
    
    // Open output FASTQ file
    let output_file = File::create(output_path)?;
    
    // Create processor
    let mut processor = PrimerTrimmer::new(
        Box::new(output_file),
        forward_primer.to_string(),
        reverse_primer.to_string(),
        search_length,
        max_mismatches,
    )?;
    
    // Process in parallel
    reader.process_parallel(&mut processor, threads)?;
    
    Ok(())
}
