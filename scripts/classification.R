# ============================================================================
# CARP repeat-classification normalisation — R mirror of scripts/classification.py
# ============================================================================
# Reads classification_vocabulary.yaml (single source of truth) and provides:
#
#   load_vocabulary(path = NULL)          -> named list (cached)
#   canonicalise(x, source = NULL, ...)   -> character (vectorised)
#   is_canonical(x, ...)                  -> logical (vectorised)
#   strip_similarity_suffix(x)            -> character (vectorised)
#
# Behaviour must match classification.py exactly; tests/test_classification.py
# and tests/test_classification.R share the same tests/classification_cases.tsv.
# ============================================================================

suppressPackageStartupMessages({
  library(yaml)
})

.classification_env <- new.env(parent = emptyenv())
.classification_env$cache <- list()

.find_default_vocabulary <- function() {
  override <- Sys.getenv("CARP_VOCABULARY", unset = NA)
  if (!is.na(override) && nzchar(override)) return(override)
  this_file <- tryCatch({
    # when sourced: use the sys.call frame; when Rscript: use commandArgs
    cmd <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cmd, value = TRUE)
    if (length(file_arg)) sub("^--file=", "", file_arg[1])
    else normalizePath(sys.frames()[[1]]$ofile %||% "scripts/classification.R",
                       mustWork = FALSE)
  }, error = function(e) "scripts/classification.R")
  here <- normalizePath(dirname(this_file), mustWork = FALSE)
  candidates <- c(
    file.path(dirname(here), "classification_vocabulary.yaml"),
    file.path(here, "classification_vocabulary.yaml")
  )
  for (c in candidates) {
    if (file.exists(c)) return(c)
  }
  stop("Could not locate classification_vocabulary.yaml; set CARP_VOCABULARY.")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

load_vocabulary <- function(path = NULL) {
  p <- if (is.null(path)) .find_default_vocabulary() else path
  p <- normalizePath(p, mustWork = TRUE)
  cached <- .classification_env$cache[[p]]
  if (!is.null(cached)) return(cached)

  raw <- yaml::read_yaml(p)

  canonical <- unlist(raw$classifications %||% character(0), use.names = FALSE)
  leaf_aliases <- raw$leaf_aliases %||% list()
  special_classes <- raw$special_classes %||% list()
  aggregation_buckets <- names(raw$aggregation_buckets %||% list())
  tool_dialects <- raw$tool_dialects %||% list()

  # Build ordered TIR prefix pairs (longest first)
  tir_entries <- tool_dialects$DANTE_TIR$hierarchy_prefixes %||% list()
  tir_pairs <- lapply(tir_entries, function(e) {
    if (!is.list(e) || is.null(e$underscore) || is.null(e$canonical)) {
      stop("DANTE_TIR hierarchy_prefixes entries must be mappings with ",
           "'underscore' and 'canonical' keys")
    }
    c(underscore = e$underscore, canonical = e$canonical)
  })
  if (length(tir_pairs) > 0) {
    ord <- order(-nchar(vapply(tir_pairs, `[[`, character(1), "underscore")))
    tir_pairs <- tir_pairs[ord]
  }

  vocab <- list(
    canonical = canonical,
    canonical_set = setNames(rep(TRUE, length(canonical)), canonical),
    leaf_aliases = leaf_aliases,
    special_classes = special_classes,
    aggregation_buckets = aggregation_buckets,
    tool_dialects = tool_dialects,
    tir_prefixes = tir_pairs,
    sources = names(tool_dialects)
  )
  .classification_env$cache[[p]] <- vocab
  vocab
}

.TIDECLUSTER_SUFFIX_RE <- "\\s*\\([0-9.]+%25\\)\\s*$"

strip_similarity_suffix <- function(x) {
  trimws(sub(.TIDECLUSTER_SUFFIX_RE, "", as.character(x)), which = "right")
}

.is_canonical_one <- function(s, vocab) {
  hit <- vocab$canonical_set[s]
  if (!is.na(hit) && isTRUE(unname(hit))) return(TRUE)
  if (s %in% vocab$aggregation_buckets) return(TRUE)
  for (cls in names(vocab$special_classes)) {
    meta <- vocab$special_classes[[cls]]
    if (identical(s, cls)) return(TRUE)
    if (isTRUE(meta$accepts_subpath) &&
        startsWith(s, paste0(cls, "/")) &&
        nchar(s) > nchar(cls) + 1) return(TRUE)
  }
  FALSE
}

is_canonical <- function(x, vocab = NULL) {
  if (is.null(vocab)) vocab <- load_vocabulary()
  vapply(as.character(x), .is_canonical_one, logical(1), vocab = vocab, USE.NAMES = FALSE)
}

.canon_dante_one <- function(s, vocab) {
  parts <- strsplit(s, "|", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  mapped <- vapply(parts, function(p) {
    hit <- vocab$leaf_aliases[[p]]
    if (is.null(hit)) p else hit
  }, character(1), USE.NAMES = FALSE)
  paste(mapped, collapse = "/")
}

.canon_dante_tir_one <- function(s, vocab) {
  s <- trimws(s)
  # Longest prefix wins — vocab$tir_prefixes is sorted descending. Fully
  # canonical input matches no prefix and is returned unchanged.
  for (pair in vocab$tir_prefixes) {
    u <- pair[["underscore"]]
    if (startsWith(s, u)) {
      leaf <- substr(s, nchar(u) + 1L, nchar(s))
      return(paste0(pair[["canonical"]], leaf))
    }
  }
  s
}

.canonicalise_one <- function(s, source, vocab, validate) {
  if (is.na(s)) stop("classification is NA")
  s <- trimws(as.character(s))
  if (!nzchar(s)) stop("classification is empty")

  if (identical(source, "TideCluster") || grepl(.TIDECLUSTER_SUFFIX_RE, s))
    s <- strip_similarity_suffix(s)

  result <- if (identical(source, "DANTE") ||
                identical(source, "DANTE_LTR") ||
                identical(source, "DANTE_LINE")) {
    .canon_dante_one(s, vocab)
  } else if (identical(source, "DANTE_TIR")) {
    .canon_dante_tir_one(s, vocab)
  } else if (identical(source, "DANTE_TIR_FALLBACK") ||
             identical(source, "RepeatMasker") ||
             identical(source, "TideCluster") ||
             identical(source, "custom_library")) {
    s
  } else if (is.null(source) || identical(source, "")) {
    if (grepl("|", s, fixed = TRUE)) {
      .canon_dante_one(s, vocab)
    } else if (any(vapply(vocab$tir_prefixes,
                          function(pair) startsWith(s, pair[["underscore"]]),
                          logical(1)))) {
      .canon_dante_tir_one(s, vocab)
    } else {
      s
    }
  } else {
    stop(sprintf("Unknown classification source: %s", source))
  }

  if (isTRUE(validate) && !.is_canonical_one(result, vocab)) {
    stop(sprintf("'%s' (from raw '%s', source=%s) is not a canonical classification",
                 result, s, source %||% "NULL"))
  }
  result
}

# Vectorised — accepts character vector, returns character vector.
# source is scalar (one tool at a time).
canonicalise <- function(x, source = NULL, vocab = NULL, validate = TRUE) {
  if (is.null(vocab)) vocab <- load_vocabulary()
  vapply(as.character(x),
         .canonicalise_one,
         character(1),
         source = source, vocab = vocab, validate = validate,
         USE.NAMES = FALSE)
}

iter_canonical <- function(vocab = NULL) {
  if (is.null(vocab)) vocab <- load_vocabulary()
  sort(vocab$canonical)
}

# Validate a character vector; returns a data.frame of failures (raw, reason).
validate_values <- function(x, source = NULL, vocab = NULL) {
  if (is.null(vocab)) vocab <- load_vocabulary()
  uniq <- unique(as.character(x))
  errors <- list()
  for (v in uniq) {
    result <- tryCatch(
      .canonicalise_one(v, source = source, vocab = vocab, validate = TRUE),
      error = function(e) structure(conditionMessage(e), class = "error")
    )
    if (inherits(result, "error"))
      errors[[length(errors) + 1L]] <- data.frame(raw = v, reason = as.character(result))
  }
  if (length(errors) == 0) data.frame(raw = character(), reason = character())
  else do.call(rbind, errors)
}

# Exported getters — scripts should not inspect vocab directly.
get_canonical <- function(vocab = NULL) {
  if (is.null(vocab)) vocab <- load_vocabulary()
  vocab$canonical
}
