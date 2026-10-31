from bioinformatics_tools.core import gc_content, reverse_complement, translate_dna


def test_gc_content():
    assert gc_content("ATGC") == 50.0
    assert gc_content("AATT") == 0.0


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AUGC") == "GCAU"


def test_translate_dna():
    assert translate_dna("ATGGCCCTT") == "MAL"
    assert translate_dna("ATGAAATGA") == "MK"
