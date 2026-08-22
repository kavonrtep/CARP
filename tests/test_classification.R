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

# ── vectorised canonicalise(): memoisation must be transparent ───────────────
# canonicalise() resolves unique() values and maps back, because callers pass one
# value per annotated feature (millions) drawn from a few dozen distinct strings.
# These cases pin the properties that make that substitution safe: same answers in
# the same order however the duplicates fall, unnamed output, empty in / empty
# out, factors accepted, and NA still an error rather than a silently mapped value.
vec_fail <- function(msg) failures <<- c(failures, paste0("vectorised: ", msg))

raws <- c("Class_I/LTR/Ty1_copia/Ale", "Class_II/Subclass_1/TIR/hAT",
          "Class_I/LINE", "Class_I/LTR/Ty3_gypsy/chromovirus/Tekay")
x <- rep(raws, length.out = 500)[c(seq(2, 500, 2), seq(1, 499, 2))]   # shuffled dups
ref <- vapply(x, .canonicalise_one, character(1),
              source = NULL, vocab = vocab, validate = TRUE, USE.NAMES = FALSE)
got <- canonicalise(x, vocab = vocab)
if (!identical(got, ref)) vec_fail("result differs from the per-element reference")
if (!is.null(names(got)))  vec_fail("result must be unnamed")
if (!identical(canonicalise(character(0), vocab = vocab), character(0)))
  vec_fail("zero-length input must give character(0)")
if (!identical(canonicalise(factor(raws), vocab = vocab),
               canonicalise(raws, vocab = vocab)))
  vec_fail("factor input must match character input")
na_err <- tryCatch({ canonicalise(c(raws[1], NA_character_), vocab = vocab); FALSE },
                   error = function(e) TRUE)
if (!na_err) vec_fail("NA must still raise an error")
if (length(failures) == 0L) passed <- passed + 5L

total <- passed + length(failures)
cat(sprintf("%d/%d cases passed\n", passed, total), file = stderr())
if (length(failures) > 0) {
  for (f in failures) cat("FAIL: ", f, "\n", sep = "", file = stderr())
  quit(status = 1)
}
quit(status = 0)
