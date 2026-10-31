# Genome Bioinformatics Analysis Tools

A lightweight Python toolkit for common DNA/RNA analysis tasks, optimized for gene-length and moderate-size FASTA inputs such as sequences around 3000 bp. This project is designed for quick sequence inspection, educational use, and simple bioinformatics workflows.

## Project Scope

This tool is optimized for gene-length and moderate-size FASTA inputs, including common sequences around 3000 bp. For inputs in this range, analysis is fast and responsive in both CLI and browser-based workflows. Very large genome-scale files may still work, but performance depends on sequence size, the number of records, and the browser environment.

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

### CLI usage

```bash
bioinformatics-tools gc-content ATGC
bioinformatics-tools reverse-complement ATGC
bioinformatics-tools translate ATGGCCCTT
bioinformatics-tools length ATGGCCCTT
bioinformatics-tools base-counts ATGGCCCTT
bioinformatics-tools motif-search ATGATGATG ATG
bioinformatics-tools genome-visualize ATGCGATCGATCGATCG --output genome_plot.png --window-size 4
bioinformatics-tools chromosome-map ATGCGATCGATCGATCG --window-size 4 --width 30
bioinformatics-tools fasta-summary examples/example.fasta
```

### Web app

```bash
streamlit run app.py
```

Then open the local URL shown by Streamlit to analyze sequences using a browser interface with plots, motif search, and VCF variant inspection.

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
