"""Command-line interface for simple bioinformatics utilities."""

from __future__ import annotations

import argparse

from .core import gc_content, read_fasta, reverse_complement, translate_dna


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Simple bioinformatics tools.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    gc_parser = subparsers.add_parser("gc-content", help="Calculate GC percentage for a DNA sequence.")
    gc_parser.add_argument("sequence", help="DNA or RNA sequence")

    rc_parser = subparsers.add_parser("reverse-complement", help="Compute the reverse complement of a DNA sequence.")
    rc_parser.add_argument("sequence", help="DNA or RNA sequence")
    rc_parser.add_argument("--rna", action="store_true", help="Treat the sequence as RNA")

    translate_parser = subparsers.add_parser("translate", help="Translate a DNA sequence into amino acids.")
    translate_parser.add_argument("sequence", help="DNA sequence")

    fasta_parser = subparsers.add_parser("fasta-summary", help="Summarize the sequences in a FASTA file.")
    fasta_parser.add_argument("file_path", help="Path to a FASTA file")

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
    elif args.command == "fasta-summary":
        records = read_fasta(args.file_path)
        for name, sequence in records.items():
            print(f">{name} len={len(sequence)} gc={gc_content(sequence):.2f}%")
    else:
        parser.error(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()
