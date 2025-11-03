use crate::primer_search::{PrimerMatcher, PairedPrimerSearchResult, reverse_complement};
use paraseq::parallel::{ParallelProcessor, IntoProcessError};
use paraseq::Record;
use std::io::Write;
use std::sync::{Arc, Mutex};
use anyhow;

/// Trim sequence and quality based on search result
/// Handles reverse complementing if needed (after trimming)
/// Returns (trimmed_sequence, trimmed_quality) if trimming is successful, None otherwise
fn trim_sequence(
    seq: &str,
    qual: &str,
    result: &PairedPrimerSearchResult,
) -> Option<(String, String)> {
    if !result.found {
        return None;
    }
    
    // Trim from both ends first (using original read coordinates)
    // trim_start: position to start keeping (trim everything before this)
    // trim_end: position to stop keeping (trim everything from this position to end)
    if result.trim_end <= result.trim_start {
        // Invalid trimming coordinates
        return None;
    }
    
    // Ensure we have valid trim coordinates
    let trim_end = result.trim_end.min(seq.len());
    if trim_end <= result.trim_start {
        return None; // Invalid or empty result
    }
    
    // Trim the sequence and quality
    let trimmed_seq = &seq[result.trim_start..trim_end];
    let trimmed_qual = qual.get(result.trim_start..trim_end)?; // Returns None if out of bounds
    
    // Reverse complement the trimmed amplicon if needed
    if result.needs_reverse_complement {
        let seq_rc = reverse_complement(trimmed_seq);
        let qual_rc = trimmed_qual.chars().rev().collect::<String>();
        Some((seq_rc, qual_rc))
    } else {
        Some((trimmed_seq.to_string(), trimmed_qual.to_string()))
    }
}

/// Processor for trimming primers from FASTQ records
/// Implements ParallelProcessor for use with paraseq parallel processing
#[derive(Clone)]
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
        min_overlap: usize,
    ) -> anyhow::Result<Self> {
        let matcher = PrimerMatcher::new(
            forward_primer,
            reverse_primer,
            search_length,
            max_mismatches,
            min_overlap,
        )?;
        
        Ok(Self {
            matcher,
            local_output: String::new(),
            global_output: Arc::new(Mutex::new(output)),
        })
    }
}

impl<R: Record> ParallelProcessor<R> for PrimerTrimmer {
    fn process_record(&mut self, record: R) -> paraseq::parallel::Result<()> {
        // Convert sequence and quality to &str (keep bindings for lifetime)
        let seq_bytes = record.seq();
        let seq = std::str::from_utf8(&seq_bytes).map_err(|e| e.into_process_error())?;
        let qual_bytes = record.qual().ok_or_else(|| anyhow::anyhow!("Missing quality scores"))?;
        let qual = std::str::from_utf8(qual_bytes).map_err(|e| e.into_process_error())?;
        
        // Search for paired primers and trim if found
        let search_result = self.matcher.search_primers(seq);
        
        if let Some((trimmed_seq, trimmed_qual)) = trim_sequence(seq, qual, &search_result) {
            // Write trimmed record to local buffer
            use std::fmt::Write;
            writeln!(self.local_output, "@{}", record.id_str())
                .map_err(|e| e.into_process_error())?;
            writeln!(self.local_output, "{}", trimmed_seq)
                .map_err(|e| e.into_process_error())?;
            writeln!(self.local_output, "+")
                .map_err(|e| e.into_process_error())?;
            writeln!(self.local_output, "{}", trimmed_qual)
                .map_err(|e| e.into_process_error())?;
        }
        
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
