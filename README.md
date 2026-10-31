# Genome Bioinformatics Analysis Tools

A lightweight Python toolkit for common DNA/RNA analysis tasks. This project is designed for quick sequence inspection, educational use, and simple bioinformatics workflows.

## Features

- GC content calculation
- Reverse complement generation
- DNA translation into protein sequence
- Sequence length and base counting
- Motif search across a sequence
- FASTA file parsing and summary

## Installation

```bash
python -m pip install -e .
```

## Quick Start

```bash
bioinformatics-tools gc-content ATGC
bioinformatics-tools reverse-complement ATGC
bioinformatics-tools translate ATGGCCCTT
bioinformatics-tools length ATGGCCCTT
bioinformatics-tools base-counts ATGGCCCTT
bioinformatics-tools motif-search ATGATGATG ATG
bioinformatics-tools fasta-summary examples/example.fasta
```

You can also invoke the module directly:

```bash
python -m bioinformatics_tools.cli gc-content ATGC
```

## Example Data

The repository includes a sample FASTA file in [examples/example.fasta](examples/example.fasta).

## Example Output

```bash
$ bioinformatics-tools gc-content ATGC
50.00%

$ bioinformatics-tools reverse-complement AUGC
GCAU

$ bioinformatics-tools translate ATGGCCCTT
MAL
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
