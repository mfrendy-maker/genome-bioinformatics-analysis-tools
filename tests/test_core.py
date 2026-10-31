from bioinformatics_tools.core import (
    base_counts,
    find_motifs,
    gc_content,
    reverse_complement,
    translate_dna,
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
