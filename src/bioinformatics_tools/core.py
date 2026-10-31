"""Core helpers for simple DNA and RNA analysis."""

from __future__ import annotations

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
    if not allow_rna and set(normalized) - {"A", "C", "G", "T", "U"}:
        invalid = sorted({base for base in normalized if base not in valid_bases})
        raise ValueError(f"Sequence contains invalid nucleotides: {', '.join(invalid)}")

    invalid = sorted({base for base in normalized if base not in valid_bases})
    if invalid:
        raise ValueError(f"Sequence contains invalid nucleotides: {', '.join(invalid)}")
    return normalized


def gc_content(sequence: str) -> float:
    seq = _normalize_sequence(sequence)
    gc_count = sum(1 for base in seq if base in {"G", "C"})
    return float((gc_count / len(seq)) * 100)


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
    dna = _normalize_sequence(sequence).replace("U", "T")
    protein: list[str] = []

    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        amino_acid = CODON_TABLE.get(codon, "X")
        if amino_acid == "*":
            break
        protein.append(amino_acid)

    return "".join(protein)


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
