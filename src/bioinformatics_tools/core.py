"""Core helpers for simple DNA and RNA analysis."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _normalize_sequence(sequence: str, *, allow_rna: bool = False) -> str:
    if sequence is None:
        raise ValueError("Sequence cannot be None.")

    normalized = "".join(ch for ch in str(sequence).strip().upper() if ch.isalpha())
    if not normalized:
        raise ValueError("Sequence must contain at least one nucleotide.")

    valid_bases = {"A", "C", "G", "T", "U"}
    invalid = sorted({base for base in normalized if base not in valid_bases})
    if invalid:
        raise ValueError(f"Sequence contains invalid nucleotides: {', '.join(invalid)}")

    if not allow_rna and "U" in normalized:
        raise ValueError("RNA sequences are not allowed in this context. Use rna=True or provide DNA data.")
    return normalized


def gc_content(sequence: str) -> float:
    seq = _normalize_sequence(sequence, allow_rna=True)
    gc_count = sum(1 for base in seq if base in {"G", "C"})
    return float((gc_count / len(seq)) * 100)


def base_counts(sequence: str) -> dict[str, int]:
    seq = _normalize_sequence(sequence, allow_rna=True)
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "U": 0}
    for base in seq:
        if base in counts:
            counts[base] += 1
    return counts


def sequence_length(sequence: str) -> int:
    return len(_normalize_sequence(sequence, allow_rna=True))


def nucleotide_frequencies(sequence: str) -> dict[str, float]:
    seq = _normalize_sequence(sequence, allow_rna=True)
    counts = base_counts(seq)
    total = sum(counts.values())
    if total == 0:
        return {base: 0.0 for base in counts}
    return {base: round((count / total) * 100, 2) for base, count in counts.items()}


def reverse_complement(sequence: str, rna: bool = False) -> str:
    seq = _normalize_sequence(sequence, allow_rna=True)
    use_rna = rna or "U" in seq
    complement_map = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "U": "A",
    }
    if use_rna:
        complement_map = {"A": "U", "U": "A", "C": "G", "G": "C"}

    return "".join(complement_map.get(base, base) for base in reversed(seq))


def translate_dna(sequence: str) -> str:
    dna = _normalize_sequence(sequence, allow_rna=True).replace("U", "T")
    protein: list[str] = []

    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        amino_acid = CODON_TABLE.get(codon, "X")
        if amino_acid == "*":
            break
        protein.append(amino_acid)

    return "".join(protein)


def find_motifs(sequence: str, motif: str) -> list[int]:
    seq = _normalize_sequence(sequence, allow_rna=True)
    motif_seq = _normalize_sequence(motif, allow_rna=True)
    if not motif_seq:
        return []

    positions: list[int] = []
    start = 0
    while True:
        index = seq.find(motif_seq, start)
        if index == -1:
            break
        positions.append(index)
        start = index + 1
    return positions


def visualize_genome(sequence: str, output_path: str = "genome_visualization.png", window_size: int = 10) -> str:
    seq = _normalize_sequence(sequence, allow_rna=True)
    if window_size <= 0:
        raise ValueError("window_size must be greater than zero.")

    scores: list[float] = []
    labels: list[str] = []
    for i in range(0, len(seq), window_size):
        segment = seq[i:i + window_size]
        if not segment:
            continue
        gc_percent = (sum(1 for base in segment if base in {"G", "C"}) / len(segment)) * 100
        scores.append(gc_percent)
        labels.append(f"{i + 1}-{i + len(segment)}")

    if not scores:
        raise ValueError("Sequence has no data to visualize.")

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12, 4))
    colors = ["#4CAF50" if value >= 50 else "#2196F3" for value in scores]
    ax.bar(labels, scores, color=colors, alpha=0.8)
    ax.axhline(50, color="black", linestyle="--", linewidth=1)
    ax.set_ylim(0, 100)
    ax.set_ylabel("GC %")
    ax.set_xlabel("Genome window")
    ax.set_title("Genome GC content visualization")
    ax.grid(axis="y", linestyle="--", alpha=0.3)
    plt.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)
    return str(output)


def chromosome_map(sequence: str, window_size: int = 10, width: int = 60) -> str:
    seq = _normalize_sequence(sequence, allow_rna=True)
    if window_size <= 0:
        raise ValueError("window_size must be greater than zero.")
    if width <= 0:
        raise ValueError("width must be greater than zero.")

    blocks: list[str] = []
    for i in range(0, len(seq), window_size):
        segment = seq[i:i + window_size]
        if not segment:
            continue
        gc_percent = (sum(1 for base in segment if base in {"G", "C"}) / len(segment)) * 100
        fullness = max(0, min(1, gc_percent / 100))
        filled = int(fullness * width)
        blocks.append("#" * filled + "." * (width - filled))

    if not blocks:
        raise ValueError("Sequence has no data to map.")

    header = "Chromosome map: " + "=" * 10 + "\n"
    return header + "\n".join(blocks)


def parse_vcf(vcf_text: str) -> list[dict[str, object]]:
    if vcf_text is None:
        return []

    records: list[dict[str, object]] = []
    for raw_line in str(vcf_text).splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]
        records.append({
            "CHROM": chrom,
            "POS": int(pos),
            "ID": vid,
            "REF": ref,
            "ALT": alt,
            "QUAL": qual,
            "FILTER": filt,
            "INFO": info,
        })
    return records


def read_fasta(path: str) -> dict[str, str]:
    records: dict[str, str] = {}
    current_name: str | None = None
    current_parts: list[str] = []

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records[current_name] = "".join(current_parts)
                current_name = line[1:].strip()
                current_parts = []
            else:
                current_parts.append(line)

    if current_name is not None:
        records[current_name] = "".join(current_parts)

    return records
