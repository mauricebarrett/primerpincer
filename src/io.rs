use crate::amplicon_processor::AmpliconRecordProcessor;
use crate::compression::OutputCompression;
use crate::preparing_input::{ExpandedPrimerSet, MyersPatternSet, PrimerSet};
use crate::search_algos::Algorithm;
use crate::sinks::{SizeFilterSink, WriterSink};
use niffler::send;
use paraseq::fastx;
use paraseq::prelude::*;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

// Buffer size matching deacon's optimized configuration
// 8MB buffers for both input and output for maximum throughput
const OUTPUT_BUFFER_SIZE: usize = 8 * 1024 * 1024; // 8MB buffer

/// Process FASTQ file with parallel primer trimming
/// Automatically handles compressed (.gz, .zst, .xz, .bz2) and uncompressed FASTQ files via niffler
#[allow(clippy::too_many_arguments)]
pub fn process_fastq(
    input_path: &str,
    output_path: &str,
    forward_primer: &str,
    reverse_primer: &str,
    search_length: usize,
    algorithm: Algorithm,
    error_rate: f64,
    min_overlap: usize,
    compression: OutputCompression,
    threads: usize,
    min_length: Option<usize>,
    max_length: Option<usize>,
) -> anyhow::Result<()> {
    // Create output directory if it doesn't exist
    if let Some(parent) = Path::new(output_path).parent() {
        std::fs::create_dir_all(parent)?;
    }

    // Open output FASTQ file with 8MB buffer and wrap with the requested compression codec.
    let output_file = File::create(output_path)?;
    let buffered = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, output_file);
    let output: Box<dyn Write + Send> = send::get_writer(
        Box::new(buffered),
        compression.to_format(),
        compression.default_level(),
    )?;

    let primers = PrimerSet::new(forward_primer, reverse_primer);

    // Build Myers matchers once if using Myers algorithm
    let myers_patterns = if matches!(algorithm, Algorithm::Myers) {
        eprintln!("🔧 Pre-building Myers pattern matchers for 4 primer variants...");
        Some(MyersPatternSet::new(&primers))
    } else {
        None
    };

    // Build expanded primer set once if using BNDM algorithm
    let expanded_primers = if matches!(algorithm, Algorithm::Bndm | Algorithm::Hamming) {
        eprintln!("🔧 Pre-expanding degenerate primers for BNDM...");
        Some(ExpandedPrimerSet::new(&primers))
    } else {
        None
    };

    // print the expanded primers if they are not none
    if let Some(ref expanded_primers) = expanded_primers {
        eprintln!("Expanded primers: {:?}", expanded_primers);
        eprintln!("Expanded primers forward: {:?}", expanded_primers.forward);
        eprintln!("Expanded primers reverse: {:?}", expanded_primers.reverse);
        eprintln!(
            "Expanded primers forward_rc: {:?}",
            expanded_primers.forward_rc
        );
        eprintln!(
            "Expanded primers reverse_rc: {:?}",
            expanded_primers.reverse_rc
        );
    }

    // Build sink chain: SizeFilter -> Writer, with shared QC counters
    let input_count = Arc::new(AtomicUsize::new(0));
    let trimmed_count = Arc::new(AtomicUsize::new(0));
    let written_count = Arc::new(AtomicUsize::new(0));
    let writer_sink = WriterSink::new(output, written_count.clone());
    let sink_chain = SizeFilterSink::new(writer_sink, min_length, max_length);

    // Create processor with pre-built caches and sink
    let mut processor = AmpliconRecordProcessor::new(
        sink_chain,
        primers,
        search_length,
        algorithm,
        error_rate,
        min_overlap,
        myers_patterns,
        expanded_primers,
        input_count.clone(),
        trimmed_count.clone(),
    )?;

    // Use niffler for automatic compression detection (gzip, zstd, xz, bzip2)
    let input_file = File::open(input_path)?;
    // Cast to Box<dyn Read + Send> for parallel processing compatibility
    let input_boxed: Box<dyn std::io::Read + Send> = Box::new(input_file);
    let (decompressed_reader, _format) = niffler::send::get_reader(input_boxed)?;
    let buffered = BufReader::with_capacity(OUTPUT_BUFFER_SIZE, decompressed_reader);
    let reader = fastx::Reader::new(buffered)?;
    reader.process_parallel(&mut processor, threads)?;

    // QC summary
    eprintln!(
        "QC — reads: input={}, post-trim={}, post-size={}",
        input_count.load(Ordering::Relaxed),
        trimmed_count.load(Ordering::Relaxed),
        written_count.load(Ordering::Relaxed)
    );

    Ok(())
}
