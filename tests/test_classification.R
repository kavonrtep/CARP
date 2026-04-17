#!/usr/bin/env Rscript
# Mirror of tests/test_classification.py.
# Runs scripts/classification.R against tests/classification_cases.tsv.
# Exit code 0 on full pass, 1 on any failure.

`%||%` <- function(a, b) if (is.null(a)) b else a

args_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(args_file) || !nzchar(args_file)) args_file <- "tests/test_classification.R"
REPO <- normalizePath(file.path(dirname(args_file), ".."), mustWork = FALSE)
source(file.path(REPO, "scripts", "classification.R"))

cases_path <- file.path(REPO, "tests", "classification_cases.tsv")
lines <- readLines(cases_path, warn = FALSE)

vocab <- load_vocabulary()

passed <- 0L
failures <- character(0)

for (i in seq_along(lines)) {
  line <- lines[i]
  if (!nzchar(line) || grepl("^\\s*#", line)) next
  if (startsWith(line, "source\t")) next
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(parts) < 3) next
  src <- if (nzchar(parts[1])) parts[1] else NULL
  raw <- parts[2]
  expected <- parts[3]
  notes <- if (length(parts) >= 4) parts[4] else ""

  if (identical(expected, "!!ERROR")) {
    result <- tryCatch(
      .canonicalise_one(raw, source = src, vocab = vocab, validate = TRUE),
      error = function(e) structure(conditionMessage(e), class = "expected_err")
    )
    if (inherits(result, "expected_err")) {
      passed <- passed + 1L
    } else {
      failures <- c(failures, sprintf(
        "L%d: expected error for source=%s raw='%s' but got '%s'",
        i, src %||% "NULL", raw, result))
    }
    next
  }

  result <- tryCatch(
    .canonicalise_one(raw, source = src, vocab = vocab, validate = TRUE),
    error = function(e) structure(conditionMessage(e), class = "unexpected_err")
  )
  if (inherits(result, "unexpected_err")) {
    failures <- c(failures, sprintf(
      "L%d: unexpected error for source=%s raw='%s': %s",
      i, src %||% "NULL", raw, as.character(result)))
    next
  }

  if (!identical(result, expected)) {
    failures <- c(failures, sprintf(
      "L%d: source=%s raw='%s'\n  expected: '%s'\n  got:      '%s'  (%s)",
      i, src %||% "NULL", raw, expected, result, notes))
  } else {
    passed <- passed + 1L
  }
}

total <- passed + length(failures)
cat(sprintf("%d/%d cases passed\n", passed, total), file = stderr())
if (length(failures) > 0) {
  for (f in failures) cat("FAIL: ", f, "\n", sep = "", file = stderr())
  quit(status = 1)
}
quit(status = 0)
