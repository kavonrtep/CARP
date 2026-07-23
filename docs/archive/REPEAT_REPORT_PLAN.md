# Repeat Annotation Report — Implementation Plan

## 1. Overview

A single self-contained HTML file generated at the end of the pipeline.
All data embedded as JSON. Plotly.js embedded inline. No server, no CDN required at view time.

**Generator script:** `scripts/make_repeat_report.R`
**Output file:** `{output_dir}/repeat_annotation_report.html`
**Snakemake rule:** `make_repeat_report` (depends on all BigWig `.done` checkpoints + key GFF3 outputs)

---

## 2. Files to create / modify

| Action | File |
|--------|------|
| Create | `scripts/make_repeat_report.R` |
| Create | `data/README_plotly.md` — note about Plotly.js requirement |
| Modify | `Snakefile` — add `make_repeat_report` rule |
| Modify | `Snakefile` — add `repeat_annotation_report.html` to `rule all` |
| Modify | `Singularity` — add `wget` of `plotly.min.js` to `data/` at build time |

---

## 3. Configurable script parameters

All have defaults; passed as CLI arguments.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output_dir` | required | Pipeline output directory |
| `--bin_width` | `100000` | Density track bin size in bp. Must be ≥ native 10 kb resolution. |
| `--min_len_chart` | `500000` | Min sequence length to appear as own bar in per-sequence composition chart. Shorter sequences aggregated into "Other sequences (N)". |
| `--min_len_tracks` | `1000000` | Min sequence length to receive its own density track panel. |
| `--max_tracks` | `50` | Hard cap on number of density track sequences regardless of threshold. |
| `--top_sat_clusters` | `10` | Number of TideCluster satellite clusters shown in section 5. |

---

## 4. Data sources

| Source file | Section used | Notes |
|-------------|-------------|-------|
| `genome_cleaned.fasta.fai` | 1, 4a, 4b | Sequence names & lengths |
| `summary_statistics.csv` | 1, 2a, 2b | Similarity-based composition |
| `Repeat_Annotation_NoSat_split_by_class_bigwig/10k/*.bw` | 4a, 4b | 42 class BigWig files |
| `TideCluster/default/TideCluster_clustering_10k.bw` | 4b, 5 | Aggregated satellite density |
| `TideCluster/default/TideCluster_clustering_split_files_bigwig/10k/*.bw` | 5 | Per-cluster density |
| `DANTE_LTR/DANTE_LTR.gff3` | 3a | Count + bp per LTR lineage |
| `DANTE_TIR/TIR_classification_summary.txt` | 3b | Count per TIR superfamily |
| `DANTE_LINE/DANTE_LINE.gff3` | 3c | Count + bp (protein domain features) |

---

## 5. BigWig handling

### Native resolution
All similarity-based BigWig files are pre-computed at 10 kb windows by the pipeline.
TideCluster BigWigs are also at 10 kb.

### Re-binning strategy
The R script re-bins to `bin_width` using `IRanges::Views` + `viewMeans` on the Rle coverage:

```
1. import(bw_file, as="RleList")  → one Rle per sequence
2. For each sequence above threshold:
   a. Extract Rle for that sequence
   b. Create Views at bin_width intervals
   c. viewMeans(views) → density array (fraction 0-1 per bin)
3. Store as list: density[[class_label]][[seqname]] = numeric_vector
```

If `bin_width == 100000` and a 100 kb file exists (pattern `*_100k.bw`), use it directly.

### Track mapping (fixed set for section 4b)
Eight tracks, each mapped to one or more BigWig files.
When a track requires multiple files (e.g., TIR = 5 files), arrays are **summed element-wise**,
then clipped to [0, 1].

| Track label | BigWig file(s) |
|-------------|---------------|
| Satellites | `TideCluster_clustering_10k.bw` |
| Ty1/copia | `All_Ty1_Copia_10k.bw` |
| Ty3/gypsy | `All_Ty3_Gypsy_10k.bw` |
| LINE | `Class_I.LINE_10k.bw` |
| TIR (Class II) | `Class_II.Subclass_1.TIR.*.bw` + `Class_II.Subclass_2.Helitron_10k.bw` |
| Simple repeats | `Simple_repeat_10k.bw` or `Simple_repeats_10k.bw` (whichever exists) |
| Low complexity | `Low_complexity_10k.bw` |
| rDNA | `rDNA_10k.bw` |

Missing files are silently skipped (track omitted).

### Per-sequence coverage for section 4a
For each class BigWig and each sequence: `sum(score × width) / seq_length`.
Scores are fractions; multiplying by window width gives bp covered; dividing by seq_length gives fraction.

---

## 6. Classification hierarchy — sunburst

Each row in `summary_statistics.csv` represents **only elements unclassified further at that level**.
Parent = all path components except the last. Plotly sunburst with `branchvalues = "remainder"`:
parent node's displayed arc = its own CSV value; children displayed separately.

Build procedure:
```
1. For each row, split path by "/" to get id, label, parent
2. Top-level rows (no "/") have parent = "" (root)
3. Add synthetic root node id="" label="All repeats" value=0
4. Pass to Plotly sunburst as ids/labels/parents/values
```

Colour scheme: Class_I = blue family, Class_II = red/orange family,
Satellites/tandem = green, Low_complexity/Simple = grey tones, rDNA = purple, Unknown = light grey.

---

## 7. Report sections

### Section 1 — Overview (no Plotly)
Four summary cards (HTML only):
- Total genome size (Mb)
- Total repeat content % (sum of all similarity-based categories)
- Number of sequences (total) / sequences shown in charts (above `min_len_chart`)
- Total structure-annotated elements (DANTE_LTR + DANTE_TIR + DANTE_LINE domain clusters)

One-sentence annotation strategy note.

### Section 2 — Repeat composition (similarity-based)
**2a.** Plotly sunburst — interactive hierarchy drill-down.
  - `branchvalues = "remainder"`. Click zooms into subtree.
  - Hover shows: classification, bp, % genome.

**2b.** Sortable HTML table — full `summary_statistics.csv` content.
  - Columns: Classification path, bp, % of genome.
  - Indentation by depth (CSS padding-left based on "/" count).
  - Sorted by % descending by default.
  - Lightweight JS (no library) for click-to-sort.

### Section 3 — Structure-based annotation (DANTE)
Three sub-panels in a flex row (or stacked on small screens).

**3a. DANTE_LTR** — parse GFF3 for `transposable_element` features, extract `Final_Classification`.
  - Horizontal bar chart (count per lineage, sorted descending).
  - Small table: lineage, count, total Mb.
  - Link to `DANTE_LTR_report.html`.

**3b. DANTE_TIR** — parse `TIR_classification_summary.txt`.
  - Horizontal bar chart (count per superfamily).
  - Link to DANTE_TIR `report.html`.

**3c. DANTE_LINE** — parse GFF3, count unique `protein_domain` loci (by coordinate).
  - Single stats card: region count, total Mb.
  - Link to log file (no dedicated HTML report).

### Section 4 — Genomic distribution (Plotly)
**4a. Per-sequence repeat content** — horizontal stacked bar.
  - Sequences ≥ `min_len_chart` each get one bar; shorter → "Other sequences (N)".
  - Categories: Satellites, Ty1/copia, Ty3/gypsy, LINE, TIR, Simple, Low_complexity, rDNA, Unknown, Other.
  - Bars sorted by total repeat content (descending).
  - X axis: fraction of sequence (0–1). Hover: bp per class, sequence length.

**4b. Density profile tracks** — Plotly subplot figure.
  - One row per track class (up to 8), shared X axis.
  - Plotly `updatemenus` dropdown to switch displayed sequence.
    - All traces pre-loaded; dropdown sets `visible` array via `restyle`.
  - X axis: position in Mb. Y axis: density (fraction, 0–1).
  - Line + filled area (`fill = "tozeroy"`).
  - Caption shows sequence name, length, bin width used.
  - Sequences with length < `min_len_tracks` excluded; max `max_tracks` sequences in dropdown.

### Section 5 — Tandem repeats / Satellites (Plotly)
  - Table: top `top_sat_clusters` TideCluster clusters by total genome coverage.
    Columns: cluster ID, total Mb, % of genome.
  - Horizontal bar chart: cluster coverage per sequence (stacked), same threshold as 4a.
  - Link to `TideCluster_report.html`.

### Section 6 — Links to sub-reports
Cards linking to: `DANTE_LTR_report.html`, `TideCluster_report.html`, `DANTE_TIR/report.html`.
Relative paths (work when HTML files stay in the same output directory).

---

## 8. JavaScript / HTML architecture

- **Plotly.js**: loaded from `data/plotly.min.js` (embedded inline at build time).
  Fallback: `download.file()` at report generation; if that fails, CDN `<script src>` tag with warning.
- **Sort table**: ~30 lines of vanilla JS embedded in the HTML, no library.
- **CSS**: inline `<style>` block. Responsive flex grid. No external stylesheet.
- **Plotly charts**: each rendered with `Plotly.newPlot(divId, data, layout, {responsive: true})`.
  All trace data and layout objects emitted as inline `<script>` blocks by the R script using `jsonlite::toJSON()`.

### Plotly.js procurement (Singularity)
Add to `Singularity` `%post` section:
```
wget -q -O /opt/pipeline/data/plotly.min.js \
  https://cdn.plot.ly/plotly-2.35.2.min.js
```

---

## 9. Snakemake integration

### New rule `make_repeat_report`
```
input:
  fai    = genome_cleaned.fasta (seqkit creates .fai alongside)
  stats  = {output_dir}/summary_statistics.csv
  bw_rm  = {output_dir}/Repeat_Annotation_NoSat_split_by_class_bigwig/.done
  bw_tc  = {output_dir}/TideCluster/default/.bigwig_done
  ltr    = {output_dir}/DANTE_LTR/DANTE_LTR.gff3
  tir    = {output_dir}/DANTE_TIR/TIR_classification_summary.txt
  line   = {output_dir}/DANTE_LINE/DANTE_LINE.gff3
output:
  {output_dir}/repeat_annotation_report.html
log:   {output_dir}/logs/make_repeat_report.log / .err
benchmark: {output_dir}/benchmarks/make_repeat_report.tsv
conda: envs/tidecluster.yaml
shell: make_repeat_report.R --output_dir ... --bin_width ... etc.
```

Parameters passed from `config` where configurable; defaults used otherwise.
New output added to `rule all`.
Added to `BENCHMARKED_RULES` list.

---

## 10. R script internal structure

```
make_repeat_report.R
├── A. parse_args()                      optparse
├── B. load_data()
│   ├── load_genome_info()               read FAI → data.frame(seqname, length)
│   ├── load_composition()               read CSV
│   ├── load_dante_ltr_stats()           rtracklayer::import GFF3, filter TE features
│   ├── load_dante_tir_stats()           read TXT
│   ├── load_dante_line_stats()          rtracklayer::import GFF3, count domains
│   └── discover_bw_files()             list.files() → named path list
├── C. process_data()
│   ├── apply_seq_thresholds()          split seqs into chart_seqs / track_seqs / other
│   ├── build_sunburst_data()           parse "/" hierarchy → ids/labels/parents/values
│   ├── compute_per_seq_coverage()      import BW → sum(score×width)/seqlen per class
│   ├── build_density_arrays()          import BW → Views/viewMeans → nested list
│   └── build_satellite_coverage()      per-cluster BW → per-seq totals + density arrays
├── D. generate_json()
│   ├── json_sunburst()                 Plotly sunburst trace
│   ├── json_dante_bars()               3 horizontal bar subplots
│   ├── json_composition_bar()          stacked horizontal bar
│   ├── json_density_subplot()          subplots + updatemenus dropdown
│   └── json_satellite_bar()           stacked horizontal bar
├── E. render_html()
│   ├── load_plotly_js()               read file or download → inline string
│   ├── html_cards()                   overview section HTML
│   ├── html_table()                   sortable composition table HTML
│   └── assemble_and_write()           sprintf/paste HTML skeleton → writeLines
└── F. main()                          orchestrates A–E
```

Dependencies (all in `tidecluster` conda env):
`optparse`, `rtracklayer`, `IRanges`, `GenomicRanges`, `jsonlite`

`jsonlite` may need to be added to `envs/tidecluster.yaml` if not already present.

---

## 11. Implementation steps (ordered)

1. Check `jsonlite` availability; add to `envs/tidecluster.yaml` if missing.
2. Add `wget plotly.min.js` to `Singularity` `%post` section.
3. Write `scripts/make_repeat_report.R` in order: A → B → C → D → E → F.
4. Add `make_repeat_report` rule to `Snakefile`.
5. Add `repeat_annotation_report.html` to `rule all` inputs.
6. Add `"make_repeat_report"` to `BENCHMARKED_RULES`.
7. Test with `--dry-run`, then with the `tmp/output` test data.

---

## 12. Data volume estimate

| Component | Approx. JSON size |
|-----------|------------------|
| Plotly.js (embedded) | 3.5 MB |
| Density tracks: 30 seqs × 8 classes × ~200 bins | ~240 KB |
| Satellite tracks: 10 clusters × 30 seqs × ~200 bins | ~300 KB |
| Per-seq coverage matrix: 50 × 10 | < 5 KB |
| Sunburst + bar chart data | < 50 KB |
| HTML + CSS + JS | ~20 KB |
| **Total estimate** | **~4.2 MB** |

Well within the 20 MB target. Headroom available for larger or more complex genomes.
