# 🦀 PrimerPincer

**Blazingly fast and accurate tool for trimming primers from FASTQ files derived from long-read sequencing (ONT & PacBio).**

PrimerPincer is a Rust-based command-line tool designed to detect and trim primers efficiently from single-end FASTQ reads, focusing on long-read sequencing data such as Oxford Nanopore and PacBio.

## Usage

Reads on stdin and writes to stdout.

```text
Usage: chopper [OPTIONS]

Options:
  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version

  -i, --input
          Path to input fastx file
          
  -o, --minlength <MINLENGTH>
          Sets a minimum read length
    
```


## Examples

```bash
pixi run cargo run -- \
 -i ./example_data/raw/ATCC-MSA1003-toy-example.fastq.gz \
 -o ./example_data/primerpincer/ATCC-MSA1003-toy-example.fastq.gz \
 -f "AGRGTTYGATYMTGGCTCAG" \
 -r "RGYTACCTTGTTACGACTT"  \
 -t 12

```