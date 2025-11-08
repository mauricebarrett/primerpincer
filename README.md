# 🦀 PrimerPincer 🦀

**Blazingly fast and accurate tool for trimming primers from FASTQ files derived from long-read sequencing (ONT & PacBio).**

[![Pixi Badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh)

PrimerPincer is a Rust-based command-line tool designed to efficiently detect and remove primers from single-end amplicon reads in FASTQ format, with a focus on long-read sequencing data such as PacBio or Oxford Nanopore platforms.

## Usage

```text
PrimerPincer - a CLI primer trimming tool for long-read sequencing data

Usage: primerpincer [OPTIONS] --input <FILE> --output <FILE> --forward <SEQUENCE> --reverse <SEQUENCE>

Options:
  -i, --input <FILE>
          Input FASTQ file

  -o, --output <FILE>
          Output FASTQ file

  -f, --forward <SEQUENCE>
          Forward primer sequence (5' to 3' orientation)

  -r, --reverse <SEQUENCE>
          Reverse primer sequence (5' to 3' orientation)

  -a, --algorithm <ALGORITHM>
          Algorithm to use for primer matching

          Possible values:
          - myers: Use Myers bit-parallel algorithm for approximate matching
          - sassy: Use Sassy SIMD-accelerated search (fastest, requires AVX2/NEON)
          - bndm:  Use BNDM for exact matching (fastest for short exact matches)
          
          [default: sassy]

  -e, --error-rate <FLOAT>
          Maximum error rate in primer matching (e.g., 0.15 for 15% errors)
          
          [default: 0.15]

  -l, --search-length <INT>
          Length to search for primer at start and end of sequence
          
          [default: 100]

  -O, --overlap <MINLENGTH>
          Minimum overlap length. Require MINLENGTH bases of the primer to match (like cutadapt -O)
          
          [default: 6]

  -t, --threads <INT>
          Number of threads to use
          
          [default: 4]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

## Examples

```bash
pixi run cargo run -- \
 -i ./example_data/raw/ATCC-MSA1003-toy-example.fastq.gz \
 -o ./example_data/primerpincer_proccesed/ATCC-MSA1003-toy-example.fastq.gz \
 -f "AGRGTTYGATYMTGGCTCAG" \
 -r "RGYTACCTTGTTACGACTT"  \
 -t 12 \
 -a sassy \
 -O 6 \
 -l 500

fqkit size -v ./example_data/primerpincer/ATCC-MSA1003-toy-example.fastq.gz 

```
