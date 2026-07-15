#!/usr/bin/env Rscript
# Regression test for the provenance-header prepend in finalise_output()
# (make_unified_annotation.R).
#
# The streaming prepend rewrites the finished unified GFF3 through a sibling temp
# file to insert `##pipeline-version / ##git-sha / ##run-started` after the
# `##gff-version` line, without pulling a multi-GB file into memory. Its on.exit
# cleanup MUST NOT call isOpen() on an already-closed connection: in R, isOpen()
# *raises* "invalid connection" on a closed/destroyed connection (it does NOT
# return FALSE). After the explicit close() on the success path, the on.exit
# isOpen() therefore crashed R at function exit — AFTER a clean write — so
# make_unified exited non-zero and snakemake failed the rule.
#
# This branch only runs when <output_dir>/run_provenance.json exists, i.e. every
# run_pipeline.py / container run (both CI fixtures) but NOT a bare `snakemake`
# invocation — which is exactly why it escaped local testing and turned both CI
# fixtures red on the 1.1.0 tag. The fix uses try(close(...), silent = TRUE) in
# on.exit instead of an isOpen() guard.
#
# prepend_provenance() is copied VERBATIM from finalise_output(); keep in sync.
# Calling it exercises the on.exit at function return, which is where the bug
# fired.

stopifnot(requireNamespace("jsonlite", quietly = TRUE))

# ---- copied verbatim from make_unified_annotation.R finalise_output() --------
prepend_provenance <- function(output_path) {
  prov_path <- file.path(dirname(output_path), "run_provenance.json")
  if (file.exists(prov_path) &&
      requireNamespace("jsonlite", quietly = TRUE)) {
    prov <- tryCatch(jsonlite::fromJSON(prov_path, simplifyVector = TRUE),
                      error = function(e) NULL)
    if (!is.null(prov)) {
      .or_unknown <- function(x) if (is.null(x) || length(x) == 0) "unknown" else as.character(x)
      header <- c(
        sprintf("##pipeline-version %s", .or_unknown(prov$pipeline_version)),
        sprintf("##git-sha %s", .or_unknown(prov$git_sha)),
        sprintf("##run-started %s", .or_unknown(prov$run_started))
      )
      tmp_out <- tempfile(tmpdir = dirname(output_path))
      con_in  <- file(output_path, "r")
      con_out <- file(tmp_out, "w")
      on.exit({ try(close(con_in),  silent = TRUE)
                try(close(con_out), silent = TRUE) }, add = TRUE)
      first <- readLines(con_in, n = 1L)
      writeLines(c(first, header), con_out)
      repeat {
        chunk <- readLines(con_in, n = 100000L)
        if (length(chunk) == 0L) break
        writeLines(chunk, con_out)
      }
      close(con_in); close(con_out)
      file.rename(tmp_out, output_path)
    }
  }
}
# -----------------------------------------------------------------------------

fresh_dir <- function(tag) {
  d <- tempfile(paste0("prov_", tag, "_")); dir.create(d); d
}

BODY <- c("chr1\tCARP\trepeat_region\t1\t100\t.\t+\t.\tID=x;Name=Foo",
          "chr1\tCARP\trepeat_region\t200\t300\t.\t-\t.\tID=y;Name=Bar")

test_success_path <- function() {
  d   <- fresh_dir("ok")
  gff <- file.path(d, "Repeat_Annotation_Unified.gff3")
  writeLines(c("##gff-version 3", BODY), gff)
  writeLines(paste0('{"schema_version":1,"pipeline_version":"9.9.9",',
                    '"git_sha":"deadbeef","run_started":"2026-07-15T00:00:00Z"}'),
             file.path(d, "run_provenance.json"))

  # This call is the regression trigger: the on.exit fires when the function
  # returns. The old isOpen() guard crashed HERE, after a fully-written file.
  err <- tryCatch({ prepend_provenance(gff); NULL },
                  error = function(e) conditionMessage(e))
  if (!is.null(err))
    stop("prepend_provenance errored on the success path (the 1.1.0 bug): ", err)

  out <- readLines(gff)
  stopifnot(identical(out[1], "##gff-version 3"))            # pragma stays first
  stopifnot(identical(out[2], "##pipeline-version 9.9.9"))
  stopifnot(identical(out[3], "##git-sha deadbeef"))
  stopifnot(identical(out[4], "##run-started 2026-07-15T00:00:00Z"))
  stopifnot(identical(out[5:6], BODY))                       # body intact, in order
  stopifnot(length(out) == 6L)                               # 1 pragma + 3 header + 2 body
  # temp file was renamed away, not left behind
  stopifnot(length(list.files(d, pattern = "^file")) == 0L)
  cat("  success path: header inserted after ##gff-version, body preserved, no on.exit crash\n")
}

test_missing_provenance_is_noop <- function() {
  d   <- fresh_dir("noprov")             # no run_provenance.json in this dir
  gff <- file.path(d, "Repeat_Annotation_Unified.gff3")
  orig <- c("##gff-version 3", BODY)
  writeLines(orig, gff)
  err <- tryCatch({ prepend_provenance(gff); NULL },
                  error = function(e) conditionMessage(e))
  if (!is.null(err)) stop("prepend_provenance errored when provenance absent: ", err)
  stopifnot(identical(readLines(gff), orig))                 # untouched
  cat("  missing provenance: file left unchanged (branch skipped cleanly)\n")
}

test_malformed_provenance_is_noop <- function() {
  d   <- fresh_dir("bad")
  gff <- file.path(d, "Repeat_Annotation_Unified.gff3")
  orig <- c("##gff-version 3", BODY)
  writeLines(orig, gff)
  writeLines("{ this is not valid json", file.path(d, "run_provenance.json"))
  err <- tryCatch({ prepend_provenance(gff); NULL },
                  error = function(e) conditionMessage(e))
  if (!is.null(err)) stop("prepend_provenance errored on malformed provenance: ", err)
  stopifnot(identical(readLines(gff), orig))                 # fromJSON->NULL, no prepend
  cat("  malformed provenance: parse failure tolerated, file left unchanged\n")
}

test_success_path()
test_missing_provenance_is_noop()
test_malformed_provenance_is_noop()
cat("test_provenance_prepend: PASSED\n")
