#!/usr/bin/env python3
"""Generate the unified_multibatch determinism fixture.

Design (see README): 3 sequences declared 2 Mb each in the .fai so make_batches
splits them into 3 batches at --batch_size 1000000 (threads>1), while threads=1
takes the single-batch path. seq_te carries a TE-derived-TRC array (te_sat
trigger); seq_ovl carries overlapping tier-1 elements + overlapping DANTE domains
FAR from any te_sat (the bug's victims: pre-fix they are pre-disjoined only when
their batch shares the te_sat, i.e. only in the single batch). seq_pad pads to a
third batch.
"""
import os, sys

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)

CRM = "Class_I|LTR|Ty3/gypsy|chromovirus|CRM"

def w(name, lines):
    with open(os.path.join(OUT, name), "w") as f:
        f.write("##gff-version 3\n")
        for l in lines:
            f.write(l + "\n")

def gff(seq, src, typ, start, end, strand, attrs):
    return f"{seq}\t{src}\t{typ}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"

# ── DANTE_LTR.gff3 ────────────────────────────────────────────────────────────
ltr = []
# seq_te: 4 complete CRM elements tiling the TRC_1 array 100001-140000
for i in range(4):
    s = 100001 + i * 10000
    e = s + 9999
    ltr.append(gff("seq_te", "DANTE_LTR", "transposable_element", s, e, "+",
                   f"ID=TE_te_{i};Final_Classification={CRM};Name={CRM};Rank=DL"))
# seq_ovl: two OVERLAPPING complete CRM elements (the CFRUME023 shape),
# far from any te_sat. Overlap 105557-107030.
ltr.append(gff("seq_ovl", "DANTE_LTR", "transposable_element", 100000, 107030, "+",
               f"ID=TE_ovl_A;Final_Classification={CRM};Name={CRM};Rank=DL"))
ltr.append(gff("seq_ovl", "DANTE_LTR", "transposable_element", 105557, 115187, "+",
               f"ID=TE_ovl_B;Final_Classification={CRM};Name={CRM};Rank=DL"))
w("DANTE_LTR.gff3", ltr)

# ── DANTE_LTR_tandems.gff3 (no LTR_RT_TR containers) ──────────────────────────
w("DANTE_LTR_tandems.gff3", [])

# ── DANTE_TIR / DANTE_LINE (empty) ────────────────────────────────────────────
w("DANTE_TIR_combined.gff3", [])
w("DANTE_LINE.gff3", [])

# ── DANTE_filtered.gff3 (tier-2 domains) ──────────────────────────────────────
dante = []
# seq_te: one RT domain per 10 kb window so te_domain_rhythm occupancy = 1.0
for i in range(4):
    c = 103000 + i * 10000
    dante.append(gff("seq_te", "DANTE", "protein_domain", c, c + 400, "+",
                     f"Final_Classification={CRM};Name=RT"))
# seq_ovl: two DANTE domains that OVERLAP EACH OTHER but no higher tier, far from
# the tier-1 elements. This is the tier-2 twin of the tier-1 bug: pre-fix their
# strand flips between the single batch (which contains seq_te's te_sat, so
# trim_to_nonoverlap disjoins ALL of t2 -> these keep their source strand) and
# seq_ovl's own batch (no te_sat overlap -> trim short-circuited, leaving the
# overlap for resolve_within_tier -> strand '*'). The trim_to_nonoverlap fix
# (always decompose lower's internal overlaps) makes both batchings identical.
dante.append(gff("seq_ovl", "DANTE", "protein_domain", 120000, 121000, "+",
                 f"Final_Classification={CRM};Name=RT"))
dante.append(gff("seq_ovl", "DANTE", "protein_domain", 120500, 121500, "-",
                 f"Final_Classification={CRM};Name=INT"))
w("DANTE_filtered.gff3", dante)

# ── TideCluster clustering (default): TRC_1 array on seq_te ────────────────────
w("TideCluster_clustering.gff3", [
    gff("seq_te", "TideCluster", "repeat_region", 100001, 140000, "+", "Name=TRC_1"),
])
w("TideCluster_clustering_short.gff3", [])

# ── RM on TideCluster (empty) ─────────────────────────────────────────────────
w("RM_on_TideCluster_Library.gff3", [])

# ── RepeatMasker merged (tier 5): a hit on seq_pad so it isn't empty ──────────
w("RM_on_combined_library_plus_DANTE.gff3", [
    gff("seq_pad", "RepeatMasker", "transposable_element", 500000, 501000, "+",
        f"Name={CRM};classification={CRM}"),
])

# ── TideHunter residuals + raw (empty) ────────────────────────────────────────
for n in ("TideHunter_short_default.gff3", "TideHunter_short_short.gff3",
          "TideHunter_raw_default.gff3", "TideHunter_raw_short.gff3"):
    w(n, [])

# ── trc_table.tsv (default): TRC_1 monomer period 10000 ───────────────────────
with open(os.path.join(OUT, "trc_table_default.tsv"), "w") as f:
    f.write("TRC_ID\tmonomer_tarean\tmonomer_kite\n")
    f.write("TRC_1\t10000\t10000\n")
with open(os.path.join(OUT, "trc_table_short.tsv"), "w") as f:
    f.write("TRC_ID\tmonomer_tarean\tmonomer_kite\n")

# ── rDNA TSVs (header only) ───────────────────────────────────────────────────
for n in ("rdna_default.tsv", "rdna_short.tsv"):
    with open(os.path.join(OUT, n), "w") as f:
        f.write("TRC\trDNA_type\tcoverage\n")

# ── .fai: 3 sequences declared 2 Mb each (forces 3 batches) ───────────────────
with open(os.path.join(OUT, "genome.fai"), "w") as f:
    for i, s in enumerate(("seq_te", "seq_ovl", "seq_pad")):
        f.write(f"{s}\t2000000\t{i*2000001}\t80\t81\n")

print("wrote fixture to", OUT)
