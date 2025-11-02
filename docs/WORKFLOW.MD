# Workflow of data through the code

## Entry point - main.rs

1. **CLI Argument Parsing** (`main()`)
   - Parses command-line arguments using Clap
   - Passes algorithm as string to `process_fastq()` (parsing deferred)

## I/O Layer - io.rs

2. **File Setup** (`process_fastq()`)
   - Opens input FASTQ file for reading using `paraseq::fastx::Reader`
   - Creates output FASTQ file for writing
   - Creates `PrimerTrimmer` processor (from `primer_trim.rs`) with configuration
   - Algorithm string passed through (not parsed here)
   - Initiates parallel processing via `reader.process_parallel()`

## Per-Record Processing - primer_trim.rs

3. **Parallel Processing Setup**
   - `PrimerTrimmer` implements `ParallelProcessor` trait from `paraseq`
   - Spawns worker threads (configurable via `--threads`)
   - Each thread processes records independently

4. **Record Processing** (`PrimerTrimmer::process_record()`)
   - Each FASTQ record is processed in parallel across threads
   - Extracts sequence and quality strings from record
   - Calls `PrimerMatcher::search_primers()` to find paired primers

## Primer Search - primer_search.rs

5. **Paired Primer Search** (`search_paired_primers()`)
   - Searches for primers using a two-scenario strategy:
   
   **Scenario 1** (normal orientation):
   - Forward primer searched at **start** (first `search_length` bases)
   - Reverse complement of reverse primer searched at **end** (last `search_length` bases)
   - If both found: returns `PairedPrimerSearchResult` with trim coordinates, `needs_reverse_complement = false`
   
   **Scenario 2** (reverse orientation - tried if Scenario 1 fails):
   - Reverse primer searched at **start** (first `search_length` bases)
   - Reverse complement of forward primer searched at **end** (last `search_length` bases)
   - If both found: returns `PairedPrimerSearchResult` with trim coordinates, `needs_reverse_complement = true`

6. **End Search Helper** (`find_primer_at_end()`)
   - Searches for primers in the last `search_length` bases
   - Extracts end region and searches using the same `PrimerSearcher` trait
   - Converts coordinates from end region to original read coordinates

7. **Primer Search Algorithms** (`primer_search.rs`)
   - **Alignment**: Semi-global alignment with edit distance (default 10% error rate)
   - **Hamming**: Hamming distance sliding window (exact mismatch count)
   - Implement `PrimerSearcher` trait with `find_primer()` method
   - Algorithm parsing happens in `PrimerMatcher::new()` (deferred from CLI)

## Trimming - primer_trim.rs

8. **Sequence Trimming** (`PrimerMatcher::trim_sequence()`)
   - If primers found (`result.found == true`):
     1. **Trim first**: Extract amplicon from original read using `trim_start` and `trim_end` coordinates
     2. **Reverse complement if needed**: If `needs_reverse_complement == true`, reverse complement the trimmed amplicon (not the whole read)
     3. Trim quality scores to match sequence trimming
     4. Validates trimmed sequence is not empty
   - Returns `Some((trimmed_seq, trimmed_qual))` if successful, `None` if invalid or empty

9. **Output Writing** (`PrimerTrimmer::process_record()`)
   - If primers found and sequence trimmed successfully:
     - Formats trimmed record as FASTQ (ID, sequence, '+', quality)
     - Writes to thread-local output buffer (`local_output`)
   - If no primers found:
     - Read is discarded (not written to output)

10. **Batch Completion** (`PrimerTrimmer::on_batch_complete()`)
   - Called periodically by `paraseq` after processing batches
   - Flushes thread-local output buffer to global output file
   - Thread-safe write using mutex (`global_output`)
   - Clears local buffer for next batch

## Data Flow Summary

```
Input FASTQ File
    ↓
[io.rs: process_fastq()]
    ├─ Open input file (fastx::Reader)
    └─ Create output file
        ↓
[primer_trim.rs: PrimerTrimmer::new()]
    └─ Creates PrimerMatcher (parses algorithm string here)
        ↓
[paraseq: reader.process_parallel()]
    ↓ (parallel processing across threads)
[primer_trim.rs: PrimerTrimmer::process_record()]
    ├─ Extract seq & qual from record
    └─ Call PrimerMatcher::search_primers()
        ↓
[primer_search.rs: search_paired_primers()]
    ├─ Scenario 1: Forward primer at start + Reverse_primer_rc at end
    │   └─ If found: return coordinates, needs_rc = false
    └─ Scenario 2: Reverse primer at start + Forward_primer_rc at end
        └─ If found: return coordinates, needs_rc = true
            ↓
[primer_search.rs: AlignmentSearcher OR HammingSearcher]
    └─ find_primer() searches in search_length window
        ↓
[primer_trim.rs: PrimerMatcher::trim_sequence()]
    ├─ Trim amplicon from original read (trim_start..trim_end)
    ├─ If needs_reverse_complement: Reverse complement trimmed amplicon
    └─ Trim quality scores
        ↓
Primers found & valid? → Yes → Format FASTQ
    ↓                    ↓
   No               Write to local buffer
    ↓                    ↓
Discard             [on_batch_complete()]
                          ↓
                   Flush to output file
                          ↓
                   Output FASTQ File
```

## Search Strategy Details

### Scenario 1: Normal Orientation
- Forward primer at 5' end (positions 0..search_length)
- Reverse complement of reverse primer at 3' end (positions read_len-search_length..read_len)
- Output: Amplicon between primers in original orientation

### Scenario 2: Reverse Orientation
- Reverse primer at 5' end (positions 0..search_length)
- Reverse complement of forward primer at 3' end (positions read_len-search_length..read_len)
- Output: Amplicon between primers, then reverse complemented

## Key Components

- **main.rs**: Entry point, CLI parsing, orchestrates execution
- **io.rs**: File I/O orchestration - opens files, sets up parallel processing, passes algorithm as string
- **primer_trim.rs**: Primer trimming module:
  - `PrimerMatcher`: Orchestrates search and trimming
    - Parses algorithm string (deferred from CLI)
    - Calls search functions from `primer_search.rs`
    - Handles trimming and reverse complementing
  - `PrimerTrimmer`: Parallel processor that implements `ParallelProcessor` trait
- **primer_search.rs**: Core search algorithms and paired search logic:
  - `reverse_complement()`: Utility function for DNA sequences
  - `search_paired_primers()`: Two-scenario paired primer search
  - `find_primer_at_end()`: Helper to search at 3' end
  - `AlignmentSearcher`: Semi-global alignment algorithm
  - `HammingSearcher`: Hamming distance algorithm
  - Both implement `PrimerSearcher` trait

## Parallelization

- Processing happens in parallel using `paraseq`
- Each thread processes records independently
- Local buffers accumulate results
- Batch completion synchronizes writes to output file
- Thread count configurable via `--threads` flag
