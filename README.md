# Bioinformatics Tools

A small Python project for simple DNA/RNA analysis tasks.

## Features

- GC content calculation
- Reverse complement generation
- DNA translation into protein sequence
- FASTA parsing and summary

## Installation

```bash
python -m pip install -e .
```

## Usage

```bash
bioinformatics-tools gc-content ATGC
bioinformatics-tools reverse-complement ATGC
bioinformatics-tools translate ATGGCCCTT
bioinformatics-tools fasta-summary examples/example.fasta
```

You can also run the module directly:

```bash
python -m bioinformatics_tools.cli gc-content ATGC
```

## Example FASTA

See the sample file in `examples/example.fasta`.
