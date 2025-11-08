use crate::cli::Algorithm;
use crate::preparing_input::PrimerSet;
use crate::primer_trim::PrimerTrimmer;
use crate::search_algos::MyersPatternCache;
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use paraseq::fastx;
use paraseq::prelude::*;
use std::fs::File;
use std::io::{BufReader, Write, Read};
use std::path::Path;
use std::time::Instant;

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

    let setup_start = Instant::now();
    let primers = PrimerSet::new(forward_primer, reverse_primer);

    // Build Myers matchers once if using Myers algorithm
    let myers_cache = if matches!(algorithm, Algorithm::Myers) {
        eprintln!("🔧 Pre-building Myers pattern matchers for 4 primer variants...");
        Some(MyersPatternCache::new(&primers))
    } else {
        None
    };

    eprintln!("🎬 Starting FASTQ processing with {} threads", threads);

    // Create processor with pre-built Myers cache
    let mut processor = PrimerTrimmer::new(
        output,
        primers,
        search_length,
        algorithm,
        error_rate,
        min_overlap,
        myers_cache,
        None,
    )?;

    // Peek at the file to detect compression using magic bytes
    let mut peek_buf = [0u8; 2];
    let mut temp_file = File::open(input_path)?;
    let bytes_read = temp_file.read(&mut peek_buf)?;
    
    // Detect compression format based on magic bytes
    // Gzip files start with 0x1f 0x8b
    let is_gzip = bytes_read >= 2 && peek_buf[0] == 0x1f && peek_buf[1] == 0x8b;
    
    let format_name = if is_gzip { "gzip" } else { "uncompressed" };
    eprintln!("📂 Detected format: {}", format_name);
    
    let setup_time = setup_start.elapsed();
    eprintln!("⏱️  Setup time: {:.2}s", setup_time.as_secs_f64());
    
    // Measure reading time separately (this includes decompression if needed)
    let read_start = Instant::now();
    let input_file = File::open(input_path)?;
    
    let processing_start = Instant::now();
    if is_gzip {
        let decoder = GzDecoder::new(input_file);
        let buffered = BufReader::with_capacity(256 * 1024, decoder);
        let reader = fastx::Reader::new(buffered)?;
        reader.process_parallel(&mut processor, threads)?;
    } else {
        // For uncompressed files, read directly
        let buffered = BufReader::with_capacity(256 * 1024, input_file);
        let reader = fastx::Reader::new(buffered)?;
        reader.process_parallel(&mut processor, threads)?;
    }
    
    let processing_time = processing_start.elapsed();
    let read_time = read_start.elapsed();
    
    eprintln!("⏱️  Total processing time: {:.2}s", processing_time.as_secs_f64());
    eprintln!("⏱️  Input reading + decompression time: {:.2}s", read_time.as_secs_f64());
    
    let read_time_s = read_time.as_secs_f64();
    
    // Print detailed timing statistics from the processor
    processor.print_timing_stats(read_time_s);

    Ok(())
}
