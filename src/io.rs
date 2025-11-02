use paraseq::prelude::*;
use paraseq::parallel::ParallelProcessor;
use paraseq::fastx;
use std::fs::File;
use std::io::Write;
use std::sync::{Arc, Mutex};
use crate::primer_search::SearchMethod;
use crate::primer_trim::PrimerMatcher;

/// Processor for trimming primers from FASTQ records
pub struct PrimerTrimmer {
    matcher: PrimerMatcher,
    local_output: String,
    global_output: Arc<Mutex<Box<dyn Write + Send>>>,
}

impl PrimerTrimmer {
    pub fn new(
        output: Box<dyn Write + Send>,
        forward_primer: String,
        reverse_primer: String,
        search_length: usize,
        max_mismatches: usize,
        algorithm: SearchMethod,
    ) -> Self {
        let matcher = PrimerMatcher::new(
            forward_primer,
            reverse_primer,
            search_length,
            max_mismatches,
            algorithm,
        );
        
        Self {
            matcher,
            local_output: String::new(),
            global_output: Arc::new(Mutex::new(output)),
        }
    }
}

impl<R: fastx::Record> ParallelProcessor<R> for PrimerTrimmer {
    fn process_record(&mut self, record: R) -> paraseq::parallel::Result<()> {
        // Get sequence and quality from record
        let seq = record.seq();
        let qual = record.qual();
        
        // Find primer match in first search_length bases
        if let Some((trim_pos, _which_primer, _orientation)) = self.matcher.find_primer_match(seq) {
            // Trim sequence and quality
            if let Some((trimmed_seq, trimmed_qual)) = self.matcher.trim_sequence(seq, qual, trim_pos) {
                // Write trimmed record manually
                use std::fmt::Write;
                writeln!(self.local_output, "@{}", record.id_str())?;
                writeln!(self.local_output, "{}", trimmed_seq)?;
                writeln!(self.local_output, "+")?;
                writeln!(self.local_output, "{}", trimmed_qual)?;
            }
        }
        // If no primer found, discard the read (don't write anything)
        
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::parallel::Result<()> {
        let mut global_out = self.global_output.lock().unwrap();
        global_out.write_all(self.local_output.as_bytes())?;
        global_out.flush()?;
        self.local_output.clear();
        Ok(())
    }
}

/// Process FASTQ file with parallel primer trimming
pub fn process_fastq(
    input_path: &str,
    output_path: &str,
    forward_primer: &str,
    reverse_primer: &str,
    search_length: usize,
    max_mismatches: usize,
    threads: usize,
    algorithm: SearchMethod,
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
        algorithm,
    );
    
    // Process in parallel
    reader.process_parallel(&mut processor, threads)?;
    
    Ok(())
}
