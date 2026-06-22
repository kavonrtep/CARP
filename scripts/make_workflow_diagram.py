#!/usr/bin/env python3
"""Generate a high-level, publication-quality schematic of the pipeline.

The *structure* is derived live from the Snakefile (via ``snakemake
--rulegraph``) so the diagram never drifts from the actual workflow: rules
are grouped into a curated set of ~9 analysis stages, and the stage-to-stage
edges are collapsed from the real rule dependency graph.

Each stage node carries an SVG **tooltip** listing its member rules together
with the first line of each rule's Snakefile docstring — so the explanatory
comments are single-sourced in the Snakefile and stay current too.

A completeness check fails loudly if a rule exists in the Snakefile but is not
assigned to a stage here (so adding a rule forces a one-line update, instead of
silently dropping it from the figure).

Outputs ``figs/workflow_overview.svg`` and ``figs/workflow_overview.png``.

Usage::

    scripts/make_workflow_diagram.py
    SNAKEMAKE_BIN=/path/to/snakemake scripts/make_workflow_diagram.py \
        --snakefile Snakefile --outdir figs

Requires ``snakemake`` and Graphviz ``dot`` on PATH (or via SNAKEMAKE_BIN /
--dot).
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

# ── Curated rule → stage mapping ──────────────────────────────────────────
# Ordered top→bottom. Every Snakefile rule except `all` must appear here.
STAGES = [
    ("Input preparation", [
        "clean_genome_fasta", "index_genome_fasta", "calculate_seqlengths"]),
    ("DANTE protein-domain annotation", [
        "dante", "filter_dante"]),
    ("LTR retrotransposons (DANTE_LTR)", [
        "dante_ltr", "make_library_of_ltrs"]),
    ("TIR DNA transposons (DANTE_TIR)", [
        "dante_tir", "dante_tir_fallback", "merge_dante_tir_with_fallback",
        "make_tir_combined_library", "build_fallback_tir_library"]),
    ("LINEs (DANTE_LINE)", [
        "dante_line"]),
    ("Tandem repeats (TideCluster)", [
        "tidecluster_long", "tidecluster_short", "tidecluster_reannotate",
        "merge_tidecluster_default_and_short"]),
    ("Repeat library construction", [
        "make_subclass_2_library", "filter_ltr_rt_library",
        "concatenate_libraries", "reduce_library"]),
    ("RepeatMasker similarity annotation", [
        "repeatmasker", "subtract_satellites_from_rm", "merge_rm_and_dante"]),
    ("Unified annotation", [
        "make_unified_annotation", "validate_classifications"]),
    ("Density tracks, outputs & reports", [
        "make_summary_statistics_and_split_by_class", "make_bigwig_density",
        "calculate_bigwig_density", "make_unified_tandem_per_family_bigwig",
        "make_tidecluster_tandem_per_family_bigwig", "make_track_for_masking",
        "make_track_for_Ns", "add_top_level_outputs", "add_html_outputs",
        "make_summary_plots", "make_benchmark_report", "make_repeat_report"]),
]

# Stage-edges that must always be drawn (never removed by transitive
# reduction), so a primary dependency isn't hidden behind a secondary path.
# The DANTE protein-domain stage fans out to all three element callers; that
# fan-out is the primary flow, but TIR also depends on the LTR *library*
# (build_fallback_tir_library), which would otherwise let reduction drop the
# direct DANTE->TIR edge. Only edges that are genuine rule dependencies are
# pinned (others are silently ignored).
PINNED_EDGES = {
    ("DANTE protein-domain annotation", "LTR retrotransposons (DANTE_LTR)"),
    ("DANTE protein-domain annotation", "TIR DNA transposons (DANTE_TIR)"),
    ("DANTE protein-domain annotation", "LINEs (DANTE_LINE)"),
}

# Muted, print-friendly fill colours (one per stage, in order).
PALETTE = [
    "#E8EEF7", "#DCE9D5", "#F6E7CF", "#F3D9D5", "#E7Dced",
    "#D6ECEC", "#EDE3D0", "#F1D9E6", "#DfE7EE", "#E9E4DE",
]


def rule_to_stage():
    m = {}
    for stage, rules in STAGES:
        for r in rules:
            m[r] = stage
    return m


_SECTION_KW = re.compile(
    r"(?m)^\s+(input|output|params|log|benchmark|conda|threads|shell|run|"
    r"resources|priority|group|wildcard_constraints|message|shadow|retries):")


def parse_rule_docs(snakefile):
    """rule_name -> first docstring line (or '').

    Only a triple-quoted block that appears *before* the first rule-section
    keyword (input/output/shell/...) counts as the docstring — otherwise the
    rule's ``shell: \"\"\"...\"\"\"`` block would be mistaken for a docstring.
    """
    text = open(snakefile, encoding="utf-8").read()
    # Split into (name, body) blocks at each top-level `rule NAME:`.
    parts = re.split(r"(?m)^rule\s+(\w+):", text)
    docs = {}
    # parts = [pre, name1, body1, name2, body2, ...]
    for i in range(1, len(parts), 2):
        name, body = parts[i], parts[i + 1]
        kw = _SECTION_KW.search(body)
        head = body[:kw.start()] if kw else body
        m = re.search(r'"""(.*?)"""', head, re.S)
        desc = ""
        if m:
            for line in m.group(1).splitlines():
                if line.strip():
                    desc = line.strip()
                    break
        docs[name] = desc
    return docs


def get_rule_edges(snakefile, snakemake_bin):
    """Return a set of (producer_rule, consumer_rule) from `--rulegraph`.

    Uses a throwaway dummy genome + config so it works without real data.
    """
    with tempfile.TemporaryDirectory() as td:
        fasta = os.path.join(td, "dummy.fasta")
        with open(fasta, "w") as fh:
            fh.write(">chr1\n" + "ACGT" * 25 + "\n")
        cfg = os.path.join(td, "config.yaml")
        with open(cfg, "w") as fh:
            fh.write(f"genome_fasta: {fasta}\noutput_dir: {td}/out\n")
        out = subprocess.run(
            [snakemake_bin, "-s", snakefile, "--configfile", cfg,
             "--rulegraph", "--cores", "1"],
            capture_output=True, text=True)
        if out.returncode != 0:
            sys.exit("snakemake --rulegraph failed:\n" + out.stderr)
        dot = out.stdout

    id2name = dict(re.findall(r'(\d+)\[label = "([^"]+)"', dot))
    raw = [(a, b) for a, b in re.findall(r"(\d+) -> (\d+)", dot)]
    edges = {(id2name[a], id2name[b]) for a, b in raw}

    # Orient producer→consumer. snakemake points dependency→dependent, so the
    # final target `all` should be a pure sink. If instead `all` is a source,
    # flip every edge.
    if any(a == "all" for a, _ in edges) and not any(b == "all" for _, b in edges):
        edges = {(b, a) for a, b in edges}
    return edges


def transitive_reduction(edges, protected=frozenset()):
    """Drop edges implied by a longer path (DAG transitive reduction), so the
    overview shows immediate stage-to-stage flow instead of every skip edge.
    Edges in ``protected`` are never removed."""
    succ = {}
    for a, b in edges:
        succ.setdefault(a, set()).add(b)
    reduced = set(edges)
    for a, b in list(edges):
        if (a, b) in protected:
            continue
        stack = [s for s in succ.get(a, ()) if s != b]
        seen = set()
        while stack:
            x = stack.pop()
            if x == b:
                reduced.discard((a, b))
                break
            if x in seen:
                continue
            seen.add(x)
            stack.extend(succ.get(x, ()))
    return reduced


def dot_str(s):
    """Escape text for a DOT double-quoted string. Graphviz handles SVG/XML
    escaping of <,>,& itself, so we only escape backslash and double-quote.
    (Our literal ``\\n`` separators are added after this, as real newlines.)"""
    return s.replace("\\", "\\\\").replace('"', '\\"')


def build_dot(snakefile, snakemake_bin):
    docs = parse_rule_docs(snakefile)
    edges = get_rule_edges(snakefile, snakemake_bin)
    r2s = rule_to_stage()

    # Completeness check.
    all_rules = set(docs) - {"all"}
    unassigned = sorted(all_rules - set(r2s))
    if unassigned:
        sys.exit("Rules not assigned to a stage in STAGES (update "
                 "make_workflow_diagram.py): " + ", ".join(unassigned))

    stage_names = [s for s, _ in STAGES]
    stage_idx = {s: i for i, s in enumerate(stage_names)}

    # Collapse rule edges to stage edges (skip intra-stage and back edges).
    stage_edges = set()
    for a, b in edges:
        sa, sb = r2s.get(a), r2s.get(b)
        if sa and sb and sa != sb and stage_idx[sa] < stage_idx[sb]:
            stage_edges.add((sa, sb))
    protected = PINNED_EDGES & stage_edges
    stage_edges = transitive_reduction(stage_edges, protected)

    lines = [
        "digraph pipeline {",
        '  graph [rankdir=TB, fontname="Helvetica", labelloc="t", '
        'fontsize=20, label="Assembly Repeat Annotation Pipeline\\n'
        '(stage overview — hover any box for its rules)", nodesep=0.4, ranksep=0.55];',
        '  node [shape=box, style="rounded,filled", fontname="Helvetica", '
        'fontsize=13, margin="0.18,0.10", penwidth=1.2, color="#555555"];',
        '  edge [color="#777777", penwidth=1.3, arrowsize=0.8];',
    ]
    for i, (stage, rules) in enumerate(STAGES):
        tip_lines = [dot_str(f"{r} — {docs.get(r, '') or '(no description)'}") for r in rules]
        tooltip = "\\n".join(tip_lines)
        plural = "s" if len(rules) != 1 else ""
        label = f"{i + 1}. {dot_str(stage)}\\n({len(rules)} rule{plural})"
        lines.append(
            f'  "{stage}" [label="{label}", fillcolor="{PALETTE[i % len(PALETTE)]}", '
            f'tooltip="{tooltip}"];')
    for a, b in sorted(stage_edges, key=lambda e: (stage_idx[e[0]], stage_idx[e[1]])):
        lines.append(f'  "{a}" -> "{b}";')
    lines.append("}")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snakefile", default="Snakefile")
    ap.add_argument("--outdir", default="figs")
    ap.add_argument("--stem", default="workflow_overview")
    ap.add_argument("--snakemake", default=os.environ.get("SNAKEMAKE_BIN")
                    or shutil.which("snakemake") or "snakemake")
    ap.add_argument("--dot", default=shutil.which("dot") or "dot")
    args = ap.parse_args()

    dot_src = build_dot(args.snakefile, args.snakemake)
    os.makedirs(args.outdir, exist_ok=True)
    dot_path = os.path.join(args.outdir, args.stem + ".dot")
    with open(dot_path, "w") as fh:
        fh.write(dot_src)
    print("wrote", dot_path)

    for fmt in ("svg", "png"):
        outp = os.path.join(args.outdir, f"{args.stem}.{fmt}")
        r = subprocess.run([args.dot, f"-T{fmt}", dot_path, "-o", outp])
        if r.returncode != 0:
            sys.exit(f"dot failed for {fmt}")
        print("wrote", outp)


if __name__ == "__main__":
    main()
