"""Command-line interface for simple bioinformatics utilities."""

from __future__ import annotations

import argparse

from .core import (
    base_counts,
    find_motifs,
    gc_content,
    nucleotide_frequencies,
    read_fasta,
    reverse_complement,
    sequence_length,
    translate_dna,
    visualize_genome,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Simple bioinformatics tools for DNA/RNA analysis.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    gc_parser = subparsers.add_parser("gc-content", help="Calculate GC percentage for a DNA sequence.")
    gc_parser.add_argument("sequence", help="DNA or RNA sequence")

    rc_parser = subparsers.add_parser("reverse-complement", help="Compute the reverse complement of a DNA sequence.")
    rc_parser.add_argument("sequence", help="DNA or RNA sequence")
    rc_parser.add_argument("--rna", action="store_true", help="Treat the sequence as RNA")

    translate_parser = subparsers.add_parser("translate", help="Translate a DNA sequence into amino acids.")
    translate_parser.add_argument("sequence", help="DNA sequence")

    length_parser = subparsers.add_parser("length", help="Return the length of a sequence.")
    length_parser.add_argument("sequence", help="DNA or RNA sequence")

    counts_parser = subparsers.add_parser("base-counts", help="Count nucleotides in a sequence.")
    counts_parser.add_argument("sequence", help="DNA or RNA sequence")

    motif_parser = subparsers.add_parser("motif-search", help="Find all motif positions inside a sequence.")
    motif_parser.add_argument("sequence", help="DNA or RNA sequence")
    motif_parser.add_argument("motif", help="Motif to search for")

    stats_parser = subparsers.add_parser("stats", help="Display GC and nucleotide frequency summary.")
    stats_parser.add_argument("sequence", help="DNA or RNA sequence")

    fasta_parser = subparsers.add_parser("fasta-summary", help="Summarize the sequences in a FASTA file.")
    fasta_parser.add_argument("file_path", help="Path to a FASTA file")

    viz_parser = subparsers.add_parser("genome-visualize", help="Generate a GC-content visualization for a sequence.")
    viz_parser.add_argument("sequence", help="DNA or RNA sequence")
    viz_parser.add_argument("--output", default="genome_visualization.png", help="Path to save the visualization image")
    viz_parser.add_argument("--window-size", type=int, default=10, help="Window size for GC content calculation")

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "gc-content":
        print(f"{gc_content(args.sequence):.2f}%")
    elif args.command == "reverse-complement":
        print(reverse_complement(args.sequence, rna=args.rna))
    elif args.command == "translate":
        print(translate_dna(args.sequence))
    elif args.command == "length":
        print(sequence_length(args.sequence))
    elif args.command == "base-counts":
        print(base_counts(args.sequence))
    elif args.command == "motif-search":
        print(find_motifs(args.sequence, args.motif))
    elif args.command == "stats":
        seq = args.sequence
        print({
            "length": sequence_length(seq),
            "gc_content": round(gc_content(seq), 2),
            "base_counts": base_counts(seq),
            "nucleotide_frequencies": nucleotide_frequencies(seq),
        })
    elif args.command == "fasta-summary":
        records = read_fasta(args.file_path)
        for name, sequence in records.items():
            print(f">{name} len={len(sequence)} gc={gc_content(sequence):.2f}%")
    elif args.command == "genome-visualize":
        output_path = visualize_genome(args.sequence, output_path=args.output, window_size=args.window_size)
        print(output_path)
    else:
        parser.error(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()
