from __future__ import annotations

from pathlib import Path

import streamlit as st

from bioinformatics_tools.core import (
    base_counts,
    chromosome_map,
    find_motifs,
    gc_content,
    nucleotide_frequencies,
    read_fasta,
    reverse_complement,
    sequence_length,
    translate_dna,
    visualize_genome,
)


st.set_page_config(page_title="Genome Bioinformatics Analysis", layout="wide")

st.title("Genome Bioinformatics Analysis Tools")
st.caption("Simple DNA/RNA analysis, genome visualization, and motif search")

uploaded = st.sidebar.file_uploader("Upload a FASTA file", type=["fa", "fasta", "txt"])
uploaded_sequence = ""
if uploaded is not None:
    content = uploaded.read().decode("utf-8", errors="ignore")
    records = read_fasta_from_text(content)
    if records:
        names = list(records.keys())
        selected_name = st.sidebar.selectbox("Select sequence", names)
        uploaded_sequence = records[selected_name]
        st.sidebar.success(f"Loaded {selected_name}")

sequence = st.text_area(
    "DNA / RNA sequence",
    value=uploaded_sequence or "ATGCGATCGATCGATCG",
    height=180,
)

if not sequence:
    st.warning("Please provide a nucleotide sequence.")
    st.stop()

window_size = st.slider("Visualization window size", min_value=2, max_value=50, value=4)

col1, col2, col3, col4 = st.columns(4)
with col1:
    st.metric("Length", sequence_length(sequence))
with col2:
    st.metric("GC %", f"{gc_content(sequence):.2f}%")
with col3:
    st.metric("A/T bases", base_counts(sequence)["A"] + base_counts(sequence)["T"])
with col4:
    st.metric("G/C bases", base_counts(sequence)["G"] + base_counts(sequence)["C"])

motif = st.text_input("Motif search", value="ATG")
positions = find_motifs(sequence, motif)
st.write(f"Motif positions for {motif!r}: {positions}")

st.subheader("Sequence summaries")
summary_cols = st.columns(2)
with summary_cols[0]:
    st.json({"base_counts": base_counts(sequence), "nucleotide_frequencies": nucleotide_frequencies(sequence)})
with summary_cols[1]:
    st.code(reverse_complement(sequence))

st.subheader("Protein translation")
prot = translate_dna(sequence)
st.code(prot if prot else "No valid protein sequence found")

st.subheader("Genome visualization")
plot_path = Path("output") / "streamlit_genome_plot.png"
plot_path.parent.mkdir(exist_ok=True)
visualize_genome(sequence, output_path=str(plot_path), window_size=window_size)
st.image(str(plot_path), caption="GC content across genome windows")

st.subheader("Chromosome map")
map_text = chromosome_map(sequence, window_size=window_size, width=50)
st.code(map_text)


def read_fasta_from_text(raw_text: str) -> dict[str, str]:
    records: dict[str, str] = {}
    current_name: str | None = None
    current_parts: list[str] = []

    for line in raw_text.splitlines():
        clean = line.strip()
        if not clean:
            continue
        if clean.startswith(">"):
            if current_name is not None:
                records[current_name] = "".join(current_parts)
            current_name = clean[1:].strip()
            current_parts = []
        else:
            current_parts.append(clean)

    if current_name is not None:
        records[current_name] = "".join(current_parts)

    return records
