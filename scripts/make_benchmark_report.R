#!/usr/bin/env Rscript
# make_benchmark_report.R
# Reads Snakemake benchmark TSV files and generates a standalone HTML report.
# Usage: make_benchmark_report.R <benchmark_dir> <output.html>

args <- commandArgs(trailingOnly = TRUE)
benchmark_dir <- args[1]
output_html   <- args[2]

# ── Read all TSV files ──────────────────────────────────────────────────────
tsv_files <- list.files(benchmark_dir, pattern = "\\.tsv$", full.names = TRUE)
if (length(tsv_files) == 0) stop("No benchmark TSV files found in: ", benchmark_dir)

read_benchmark <- function(f) {
  rule_name <- sub("\\.tsv$", "", basename(f))
  d <- tryCatch(
    read.table(f, header = TRUE, sep = "\t", check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d$rule <- rule_name
  d
}

df <- do.call(rbind, Filter(Negate(is.null), lapply(tsv_files, read_benchmark)))
colnames(df)[colnames(df) == "h:m:s"] <- "hms"

# Coerce numeric columns (Snakemake writes "-" when measurement unavailable)
for (col in c("s", "max_rss", "max_vms", "max_uss", "max_pss",
              "io_in", "io_out", "mean_load", "cpu_time")) {
  if (col %in% colnames(df)) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
}

# Sort by wall time descending for all plots
df <- df[order(-df$s, na.last = TRUE), ]

# ── SVG helper ─────────────────────────────────────────────────────────────
# Returns inline SVG string from a base-R plot function
make_svg <- function(plot_fn, width = 10, height = 5.5) {
  tmp <- tempfile(fileext = ".svg")
  on.exit(unlink(tmp))
  svg(tmp, width = width, height = height)
  tryCatch(plot_fn(), error = function(e) {
    plot.new(); text(0.5, 0.5, paste("Plot error:", e$message), cex = 0.8)
  })
  dev.off()
  paste(readLines(tmp, warn = FALSE), collapse = "\n")
}

bar_colors <- function(vals, hi = "#d73027", lo = "#4575b4") {
  ifelse(!is.na(vals) & vals == max(vals, na.rm = TRUE), hi, lo)
}

# ── Plot 1: wall-clock time ─────────────────────────────────────────────────
svg_time <- make_svg(function() {
  vals <- df$s / 60
  par(mar = c(11, 5, 3, 1))
  barplot(vals, names.arg = df$rule, las = 2, cex.names = 0.7,
          col = bar_colors(vals), border = NA,
          ylab = "Wall time (minutes)", main = "Wall clock time per rule",
          col.main = "#2c3e50")
  grid(nx = NA, ny = NULL, lty = "dotted", col = "#cccccc")
})

# ── Plot 2: peak RSS memory ─────────────────────────────────────────────────
svg_mem <- make_svg(function() {
  vals <- df$max_rss / 1024   # MB → GB
  par(mar = c(11, 5, 3, 1))
  barplot(vals, names.arg = df$rule, las = 2, cex.names = 0.7,
          col = bar_colors(vals), border = NA,
          ylab = "Peak RSS memory (GB)", main = "Peak memory usage per rule",
          col.main = "#2c3e50")
  grid(nx = NA, ny = NULL, lty = "dotted", col = "#cccccc")
})

# ── Plot 3: CPU time ────────────────────────────────────────────────────────
svg_cpu <- make_svg(function() {
  vals <- df$cpu_time / 60
  par(mar = c(11, 5, 3, 1))
  barplot(vals, names.arg = df$rule, las = 2, cex.names = 0.7,
          col = bar_colors(vals, hi = "#d73027", lo = "#91bfdb"), border = NA,
          ylab = "CPU time (minutes)", main = "CPU time per rule",
          col.main = "#2c3e50")
  grid(nx = NA, ny = NULL, lty = "dotted", col = "#cccccc")
})

# ── Plot 4: parallelisation efficiency (cpu_time / (s * mean_load * ncpus))
#    Simple proxy: cpu_time / s  (> 1 means parallel, ~1 serial)
svg_par <- make_svg(function() {
  vals <- ifelse(df$s > 0, df$cpu_time / df$s, NA)
  par(mar = c(11, 5, 3, 1))
  cols <- ifelse(!is.na(vals) & vals > 1, "#1a9850", "#4575b4")
  barplot(vals, names.arg = df$rule, las = 2, cex.names = 0.7,
          col = cols, border = NA,
          ylab = "CPU time / wall time", main = "Parallelisation factor (cpu / wall)",
          col.main = "#2c3e50")
  abline(h = 1, lty = 2, col = "#e74c3c")
  grid(nx = NA, ny = NULL, lty = "dotted", col = "#cccccc")
})

# ── Summary table ───────────────────────────────────────────────────────────
fmt <- function(x, digits = 2) ifelse(is.na(x), "-", formatC(x, digits = digits, format = "f"))

table_rows <- paste(apply(df, 1, function(r) {
  sprintf(
    '<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>',
    r["rule"],
    r["hms"],
    fmt(as.numeric(r["s"])),
    fmt(as.numeric(r["max_rss"]) / 1024),
    fmt(as.numeric(r["cpu_time"]) / 60),
    fmt(as.numeric(r["mean_load"])),
    fmt(as.numeric(r["io_in"])),
    fmt(as.numeric(r["io_out"]))
  )
}), collapse = "\n")

total_wall <- sum(df$s, na.rm = TRUE)
total_cpu  <- sum(df$cpu_time, na.rm = TRUE)
peak_mem   <- max(df$max_rss, na.rm = TRUE)

# ── Assemble HTML ───────────────────────────────────────────────────────────
html <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pipeline Benchmark Report</title>
<style>
  * { box-sizing: border-box; }
  body { font-family: "Segoe UI", Arial, sans-serif; max-width: 1200px;
         margin: 0 auto; padding: 24px; color: #333; background: #fafafa; }
  h1 { color: #2c3e50; margin-bottom: 4px; }
  .meta { color: #777; font-size: 0.9em; margin-bottom: 24px; }
  h2 { color: #2c3e50; border-bottom: 2px solid #e0e0e0;
       padding-bottom: 6px; margin-top: 40px; }
  .cards { display: flex; gap: 16px; flex-wrap: wrap; margin-bottom: 32px; }
  .card { background: white; border-radius: 8px; padding: 16px 24px;
          box-shadow: 0 1px 4px rgba(0,0,0,.1); min-width: 180px; }
  .card .label { font-size: 0.8em; color: #777; text-transform: uppercase;
                 letter-spacing: .05em; }
  .card .value { font-size: 1.6em; font-weight: 600; color: #2c3e50; }
  .plot-wrap { background: white; border-radius: 8px; padding: 12px;
               box-shadow: 0 1px 4px rgba(0,0,0,.1); margin-bottom: 24px; }
  .plot-wrap svg { width: 100%%; height: auto; display: block; }
  table { border-collapse: collapse; width: 100%%; background: white;
          border-radius: 8px; overflow: hidden;
          box-shadow: 0 1px 4px rgba(0,0,0,.1); font-size: 0.88em; }
  thead th { background: #2c3e50; color: white; padding: 10px 14px;
             text-align: left; white-space: nowrap; }
  tbody td { padding: 7px 14px; border-bottom: 1px solid #f0f0f0; }
  tbody tr:last-child td { border-bottom: none; }
  tbody tr:hover { background: #f0f4f8; }
  td:not(:first-child) { font-variant-numeric: tabular-nums; }
</style>
</head>
<body>
<h1>Pipeline Benchmark Report</h1>
<p class="meta">Generated: %s &nbsp;|&nbsp; Rules: %d</p>

<h2>Summary</h2>
<div class="cards">
  <div class="card"><div class="label">Total wall time</div>
    <div class="value">%s min</div></div>
  <div class="card"><div class="label">Total CPU time</div>
    <div class="value">%s min</div></div>
  <div class="card"><div class="label">Peak memory (RSS)</div>
    <div class="value">%s GB</div></div>
  <div class="card"><div class="label">Mean parallel factor</div>
    <div class="value">%.1f ×</div></div>
</div>

<h2>Wall Clock Time</h2>
<div class="plot-wrap">%s</div>

<h2>Peak Memory Usage</h2>
<div class="plot-wrap">%s</div>

<h2>CPU Time</h2>
<div class="plot-wrap">%s</div>

<h2>Parallelisation Factor (CPU time / wall time)</h2>
<div class="plot-wrap">%s</div>

<h2>Detailed Table</h2>
<table>
<thead><tr>
  <th>Rule</th><th>Duration (h:m:s)</th><th>Wall time (s)</th>
  <th>Max RSS (GB)</th><th>CPU time (min)</th><th>Mean CPU load</th>
  <th>I/O in (MB)</th><th>I/O out (MB)</th>
</tr></thead>
<tbody>%s</tbody>
</table>
</body>
</html>',
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  nrow(df),
  fmt(total_wall / 60),
  fmt(total_cpu  / 60),
  fmt(peak_mem   / 1024),
  mean(df$cpu_time / df$s, na.rm = TRUE),
  svg_time, svg_mem, svg_cpu, svg_par,
  table_rows
)

writeLines(html, output_html)
message("Benchmark report written to: ", output_html)
