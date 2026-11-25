#!/usr/bin/env Rscript
# HTML note about the DANTE extra protein database (make_repeat_report.R).
# The functions are taken from the script itself by parsing it (its main is
# unguarded, so it cannot be sourced), so there is no copy to keep in sync.
stopifnot(requireNamespace("jsonlite", quietly = TRUE))
REPO <- normalizePath(file.path(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])), ".."))
wanted <- c("load_provenance", "extra_db_summary_path", "load_extra_db_summary", "html_extra_db")
for (e in parse(file.path(REPO, "scripts", "make_repeat_report.R")))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% wanted) eval(e, globalenv())
stopifnot(all(vapply(wanted, exists, logical(1))))

d <- tempfile("extra_db_"); dir.create(file.path(d, "DANTE"), recursive = TRUE)
f <- file.path(d, "DANTE", "extra_database_summary.tsv")
write_summary <- function(rows) writeLines(c(
  "#source\tmy <proteins> & co.fa", "#sha256\tabc123", "#rexdb\tViridiplantae_v4.0",
  "#proteins\t3", "classification\tdomain\tproteins\tlevel", rows), f)
write_summary(c("Class_I|LTR|Ty1/copia\tRT\t1200\tinternal",
                "Class_I|LTR|Ty1/copia|Ale\tGAG\t1\tterminal",
                "Class_I|LTR|Ty1/copia|Ale\tRT\t2\tterminal"))

# content
h <- html_extra_db(load_extra_db_summary(f))
has <- function(x, s) grepl(s, x, fixed = TRUE)
stopifnot(length(h) == 1L, has(h, "Extra protein database used"),
          has(h, "REXdb Viridiplantae_v4.0 extended with 1,203 user proteins"),
          has(h, "my &lt;proteins&gt; &amp; co.fa"), !has(h, "<proteins>"),
          has(h, "1,200 of the 1,203 proteins have an internal class"),
          has(h, "GAG (1), RT (2)"), has(h, "abc123"))
# no internal classes -> no internal note; a single row still works
write_summary("Class_I|LINE\tRT\t1\tterminal")
h <- html_extra_db(load_extra_db_summary(f))
stopifnot(!has(h, "internal"), has(h, "1 user protein from"))

# when it is shown
prov <- function(json) writeLines(json, file.path(d, "run_provenance.json"))
stopifnot(identical(extra_db_summary_path(d, "x.tsv"), "x.tsv"))       # passed by the rule
stopifnot(is.null(extra_db_summary_path(d, "")))                       # no provenance
prov('{"config": {"genome_fasta": "g.fa"}}')
stopifnot(is.null(extra_db_summary_path(d, "")))                       # stale file, option unset
prov('{"config": {"dante_extra_database": ""}}')
stopifnot(is.null(extra_db_summary_path(d, "")))
prov('{"config": {"dante_extra_database": null}}')
stopifnot(is.null(extra_db_summary_path(d, "")))
prov('{"config": {"dante_extra_database": "p.fa"}}')
stopifnot(identical(extra_db_summary_path(d, ""), f))                  # re-rendered by hand
unlink(f)
stopifnot(is.null(extra_db_summary_path(d, "")))
cat("OK: extra protein database note\n")
