from pathlib import Path

from bioinformatics_tools.core import (
    base_counts,
    chromosome_map,
    find_motifs,
    gc_content,
    parse_vcf,
    reverse_complement,
    translate_dna,
    visualize_genome,
)


def test_gc_content():
    assert gc_content("ATGC") == 50.0
    assert gc_content("AATT") == 0.0


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AUGC") == "GCAU"


def test_translate_dna():
    assert translate_dna("ATGGCCCTT") == "MAL"
    assert translate_dna("ATGAAATGA") == "MK"


def test_base_counts():
    counts = base_counts("AATGCC")
    assert counts == {"A": 2, "C": 2, "G": 1, "T": 1, "U": 0}


def test_find_motifs():
    assert find_motifs("ATGATGATG", "ATG") == [0, 3, 6]
    assert find_motifs("AAAA", "TT") == []


def test_visualize_genome(tmp_path):
    output_path = tmp_path / "genome_visualization.png"
    result = visualize_genome("ATGCGATCGATCGATCG", output_path=str(output_path), window_size=4)
    assert result == str(output_path)
    assert Path(result).exists()


def test_chromosome_map():
    map_text = chromosome_map("ATGCGATCGATCGATCG", window_size=4)
    assert isinstance(map_text, str)
    assert len(map_text) > 10
    assert "#" in map_text or "█" in map_text or "=" in map_text


def test_parse_vcf():
    vcf_text = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10\t.\tA\tT\t.\tPASS\tDP=10
chr1\t20\t.\tG\tC\t.\tPASS\tDP=15
"""
    variants = parse_vcf(vcf_text)
    assert len(variants) == 2
    assert variants[0]["CHROM"] == "chr1"
    assert variants[0]["POS"] == 10
    assert variants[1]["ALT"] == "C"
