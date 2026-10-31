from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from bioinformatics_tools.core import (
    base_counts,
    chromosome_map,
    find_motifs,
    gc_content,
    nucleotide_frequencies,
    parse_vcf,
    reverse_complement,
    sequence_length,
    translate_dna,
    visualize_genome,
)


st.set_page_config(page_title="Genome Browser", layout="wide")

st.title("Genome Browser & Analysis Dashboard")
st.caption("Upload a FASTA or VCF to inspect GC content, motifs, and variant calls")

uploaded = st.sidebar.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"])
uploaded_sequence = ""
if uploaded is not None:
    content = uploaded.read().decode("utf-8", errors="ignore")
    records = read_fasta_from_text(content)
    if records:
        names = list(records.keys())
        selected_name = st.sidebar.selectbox("Select sequence", names)
        uploaded_sequence = records[selected_name]
        st.sidebar.success(f"Loaded {selected_name}")

vcf_upload = st.sidebar.file_uploader("Upload VCF file", type=["vcf", "txt"])
vcf_variants = []
if vcf_upload is not None:
    vcf_text = vcf_upload.read().decode("utf-8", errors="ignore")
    vcf_variants = parse_vcf(vcf_text)
    st.sidebar.success(f"Loaded {len(vcf_variants)} variant(s)")

sequence = st.text_area(
    "DNA / RNA sequence",
    value=uploaded_sequence or "ATGCGATCGATCGATCG",
    height=180,
)

if not sequence:
    st.warning("Please provide a nucleotide sequence.")
    st.stop()

window_size = st.slider("Genome window size", min_value=2, max_value=50, value=4)

col1, col2, col3, col4 = st.columns(4)
counts = base_counts(sequence)
with col1:
    st.metric("Length", sequence_length(sequence))
with col2:
    st.metric("GC %", f"{gc_content(sequence):.2f}%")
with col3:
    st.metric("AT bases", counts["A"] + counts["T"])
with col4:
    st.metric("GC bases", counts["G"] + counts["C"])

motif = st.text_input("Motif search", value="ATG")
positions = find_motifs(sequence, motif)

overview, track, map_view, motifs, variants = st.tabs(["Overview", "GC Track", "Chromosome Map", "Motif Hits", "Variants"])

with overview:
    st.subheader("Sequence summary")
    st.json({
        "length": sequence_length(sequence),
        "gc_content": round(gc_content(sequence), 2),
        "base_counts": counts,
        "nucleotide_frequencies": nucleotide_frequencies(sequence),
        "reverse_complement": reverse_complement(sequence),
    })

    st.subheader("Protein translation")
    protein = translate_dna(sequence)
    st.code(protein if protein else "No valid protein translation found")

with track:
    st.subheader("GC content track")
    gc_values = []
    labels = []
    for i in range(0, len(sequence), window_size):
        segment = sequence[i:i + window_size]
        if not segment:
            continue
        gc_values.append((sum(1 for b in segment if b in {"G", "C"}) / len(segment)) * 100)
        labels.append(f"{i + 1}-{i + len(segment)}")

    if gc_values:
        chart_df = pd.DataFrame({"window": labels, "gc_percent": gc_values})
        st.line_chart(chart_df.set_index("window"))
    else:
        st.info("No GC track data available.")

    plot_path = Path("output") / "streamlit_genome_plot.png"
    plot_path.parent.mkdir(exist_ok=True)
    visualize_genome(sequence, output_path=str(plot_path), window_size=window_size)
    st.image(str(plot_path), caption="Genome GC content visualization")

with map_view:
    st.subheader("Chromosome-style map")
    st.code(chromosome_map(sequence, window_size=window_size, width=50))

with motifs:
    st.subheader("Motif positions")
    st.write(f"Motif {motif!r} found at positions: {positions}")
    if positions:
        highlight_text = []
        last = 0
        for pos in positions:
            highlight_text.append(sequence[last:pos])
            highlight_text.append(f"[{sequence[pos:pos + len(motif)]}]")
            last = pos + len(motif)
        highlight_text.append(sequence[last:])
        st.code("".join(highlight_text))

    motif_rows = [{"position": pos, "motif": sequence[pos:pos + len(motif)]} for pos in positions]
    st.dataframe(pd.DataFrame(motif_rows), use_container_width=True)

with variants:
    st.subheader("VCF variant overview")
    if not vcf_variants:
        st.info("Upload a VCF file to see variant details here.")
    else:
        variant_df = pd.DataFrame(vcf_variants)
        st.dataframe(variant_df, use_container_width=True)
        st.json({
            "total_variants": len(variant_df),
            "chromosomes": variant_df["CHROM"].value_counts().to_dict() if "CHROM" in variant_df.columns else {},
            "variant_types": {
                "SNVs": int((variant_df["REF"].str.len() == 1) & (variant_df["ALT"].str.len() == 1)).sum(),
                "indels": int((variant_df["REF"].str.len() != variant_df["ALT"].str.len())).sum(),
            },
        })


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
