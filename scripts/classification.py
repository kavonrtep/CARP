#!/usr/bin/env python3
"""
CARP repeat-classification normalisation module.

Reads classification_vocabulary.yaml (the single source of truth) and
provides a programmatic API plus a CLI for:

  * canonicalising tool-native classification strings into slash form
  * validating classifications against the authoritative vocabulary
  * rewriting FASTA headers and GFF3 classification attributes in place

See CLASSIFICATION_REFACTOR_PLAN.md for the design rationale.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

import yaml

__all__ = [
    "UnknownClassification",
    "Vocabulary",
    "load_vocabulary",
    "canonicalise",
    "is_canonical",
    "strip_similarity_suffix",
    "iter_canonical",
]


class UnknownClassification(ValueError):
    """Raised when a classification string cannot be resolved against the vocabulary."""


_TIDECLUSTER_SUFFIX_RE = re.compile(r"\s*\([0-9.]+%25\)\s*$")


@dataclass
class Vocabulary:
    canonical: frozenset[str]
    leaf_aliases: dict[str, str]
    special_classes: dict[str, dict]
    aggregation_buckets: frozenset[str]
    tool_dialects: dict[str, dict]
    # Ordered longest-first so longest match wins.
    tir_prefixes: tuple[tuple[str, str], ...]
    sources: frozenset[str] = field(default_factory=frozenset)

    def is_canonical(self, s: str) -> bool:
        if s in self.canonical:
            return True
        if s in self.aggregation_buckets:
            return True
        for cls, meta in self.special_classes.items():
            if s == cls:
                return True
            if meta.get("accepts_subpath") and s.startswith(cls + "/") and len(s) > len(cls) + 1:
                return True
        return False


_VOCAB_CACHE: dict[Path, Vocabulary] = {}


def _find_default_vocabulary() -> Path:
    env_override = os.environ.get("CARP_VOCABULARY")
    if env_override:
        return Path(env_override)
    here = Path(__file__).resolve().parent
    for candidate in (here.parent / "classification_vocabulary.yaml", here / "classification_vocabulary.yaml"):
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "Could not locate classification_vocabulary.yaml; set CARP_VOCABULARY or pass --vocabulary."
    )


def load_vocabulary(path: str | Path | None = None) -> Vocabulary:
    p = Path(path) if path else _find_default_vocabulary()
    p = p.resolve()
    cached = _VOCAB_CACHE.get(p)
    if cached is not None:
        return cached

    with p.open() as fh:
        raw = yaml.safe_load(fh)

    tool_dialects = raw.get("tool_dialects", {}) or {}
    dante_tir = tool_dialects.get("DANTE_TIR", {}) or {}
    raw_prefixes = dante_tir.get("hierarchy_prefixes", []) or []
    prefix_pairs: list[tuple[str, str]] = []
    for entry in raw_prefixes:
        if not isinstance(entry, dict):
            raise ValueError(
                f"DANTE_TIR hierarchy_prefixes must be mappings with "
                f"'underscore' and 'canonical' keys, got: {entry!r}"
            )
        prefix_pairs.append((entry["underscore"], entry["canonical"]))
    prefix_pairs.sort(key=lambda p: len(p[0]), reverse=True)
    tir_prefixes = tuple(prefix_pairs)

    vocab = Vocabulary(
        canonical=frozenset(raw.get("classifications", []) or []),
        leaf_aliases=dict(raw.get("leaf_aliases", {}) or {}),
        special_classes=dict(raw.get("special_classes", {}) or {}),
        aggregation_buckets=frozenset((raw.get("aggregation_buckets", {}) or {}).keys()),
        tool_dialects=tool_dialects,
        tir_prefixes=tir_prefixes,
        sources=frozenset(tool_dialects.keys()),
    )
    _VOCAB_CACHE[p] = vocab
    return vocab


def strip_similarity_suffix(s: str) -> str:
    return _TIDECLUSTER_SUFFIX_RE.sub("", s).rstrip()


def _canon_dante(s: str, vocab: Vocabulary) -> str:
    parts = [p.strip() for p in s.split("|") if p.strip()]
    parts = [vocab.leaf_aliases.get(p, p) for p in parts]
    result = "/".join(parts)
    return result


def _canon_dante_tir(s: str, vocab: Vocabulary) -> str:
    s = s.strip()
    # Longest prefix wins — vocab.tir_prefixes is sorted descending. Matches
    # both the fully underscore-encoded GFF3 form (Class_II_Subclass_1_TIR_X)
    # and the mixed FASTA-header form (Class_II/Subclass_1/TIR_X). Fully
    # canonical input (Class_II/Subclass_1/TIR/X) matches no prefix → returned
    # unchanged.
    for underscore_prefix, canonical_prefix in vocab.tir_prefixes:
        if s.startswith(underscore_prefix):
            leaf = s[len(underscore_prefix):]
            return canonical_prefix + leaf
    return s


def canonicalise(s: str, source: str | None = None, vocab: Vocabulary | None = None, validate: bool = True) -> str:
    if vocab is None:
        vocab = load_vocabulary()
    if s is None:
        raise UnknownClassification("classification is None")

    s = s.strip()
    if not s:
        raise UnknownClassification("classification is empty")

    # TideCluster appends a similarity percentage — strip before parsing.
    if source == "TideCluster" or _TIDECLUSTER_SUFFIX_RE.search(s):
        s = strip_similarity_suffix(s)

    if source in ("DANTE", "DANTE_LTR", "DANTE_LINE"):
        result = _canon_dante(s, vocab)
    elif source == "DANTE_TIR":
        result = _canon_dante_tir(s, vocab)
    elif source in ("DANTE_TIR_FALLBACK", "RepeatMasker", "TideCluster", "custom_library"):
        result = s
    elif source is None:
        # Auto-detect by separator / known DANTE_TIR prefix.
        if "|" in s:
            result = _canon_dante(s, vocab)
        elif any(s.startswith(u) for u, _ in vocab.tir_prefixes):
            result = _canon_dante_tir(s, vocab)
        else:
            result = s
    else:
        raise ValueError(f"Unknown classification source: {source!r}")

    if validate and not vocab.is_canonical(result):
        raise UnknownClassification(
            f"{result!r} (from raw {s!r}, source={source}) is not a canonical classification"
        )
    return result


def is_canonical(s: str, vocab: Vocabulary | None = None) -> bool:
    if vocab is None:
        vocab = load_vocabulary()
    return vocab.is_canonical(s)


def iter_canonical(vocab: Vocabulary | None = None) -> Iterator[str]:
    if vocab is None:
        vocab = load_vocabulary()
    yield from sorted(vocab.canonical)


# ---------------------------------------------------------------------------
# FASTA / GFF3 helpers
# ---------------------------------------------------------------------------

_FASTA_HEADER_RE = re.compile(r"^>([^#\s]+)#(\S+)(.*)$")


def canonicalise_fasta_headers(
    infile: str | Path,
    outfile: str | Path,
    source: str | None = None,
    vocab: Vocabulary | None = None,
    validate: bool = True,
) -> int:
    """Rewrite #class suffix on every FASTA header. Returns number of headers rewritten."""
    if vocab is None:
        vocab = load_vocabulary()
    rewritten = 0
    with open(infile) as fin, open(outfile, "w") as fout:
        for line in fin:
            if not line.startswith(">"):
                fout.write(line)
                continue
            m = _FASTA_HEADER_RE.match(line.rstrip("\n"))
            if not m:
                fout.write(line)
                continue
            name, cls, tail = m.group(1), m.group(2), m.group(3)
            new_cls = canonicalise(cls, source=source, vocab=vocab, validate=validate)
            fout.write(f">{name}#{new_cls}{tail}\n")
            if new_cls != cls:
                rewritten += 1
    return rewritten


_GFF_ATTR_RE = re.compile(r"(?P<k>[^=;]+)=(?P<v>[^;]*)")


def _iter_gff_lines(path: str | Path) -> Iterator[str]:
    with open(path) as fh:
        for line in fh:
            yield line


def canonicalise_gff3_attribute(
    infile: str | Path,
    outfile: str | Path,
    attribute: str,
    source: str | None = None,
    vocab: Vocabulary | None = None,
    validate: bool = True,
) -> int:
    """Rewrite a single GFF3 attribute's value on every feature line."""
    if vocab is None:
        vocab = load_vocabulary()
    rewritten = 0
    needle = attribute + "="
    with open(outfile, "w") as fout:
        for line in _iter_gff_lines(infile):
            if line.startswith("#") or "\t" not in line:
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                fout.write(line)
                continue
            attrs = fields[8]
            if needle not in attrs:
                fout.write(line)
                continue
            new_parts = []
            changed = False
            for part in attrs.split(";"):
                if part.startswith(needle):
                    old = part[len(needle):]
                    new = canonicalise(old, source=source, vocab=vocab, validate=validate)
                    if new != old:
                        changed = True
                    new_parts.append(f"{attribute}={new}")
                else:
                    new_parts.append(part)
            fields[8] = ";".join(new_parts)
            fout.write("\t".join(fields) + "\n")
            if changed:
                rewritten += 1
    return rewritten


def collect_gff3_attribute(infile: str | Path, attribute: str) -> list[str]:
    values: list[str] = []
    needle = attribute + "="
    for line in _iter_gff_lines(infile):
        if line.startswith("#") or "\t" not in line:
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        for part in fields[8].split(";"):
            if part.startswith(needle):
                values.append(part[len(needle):])
                break
    return values


def collect_fasta_classes(infile: str | Path) -> list[str]:
    values: list[str] = []
    with open(infile) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            m = _FASTA_HEADER_RE.match(line.rstrip("\n"))
            if m:
                values.append(m.group(2))
    return values


def validate_values(
    values: Iterable[str],
    source: str | None,
    vocab: Vocabulary | None = None,
) -> list[tuple[str, str]]:
    """Return list of (raw, reason) for every value that fails canonicalisation."""
    if vocab is None:
        vocab = load_vocabulary()
    errors: list[tuple[str, str]] = []
    seen: set[str] = set()
    for v in values:
        if v in seen:
            continue
        seen.add(v)
        try:
            canonicalise(v, source=source, vocab=vocab, validate=True)
        except UnknownClassification as e:
            errors.append((v, str(e)))
        except ValueError as e:
            errors.append((v, str(e)))
    return errors


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--vocabulary", help="Override path to classification_vocabulary.yaml")
    p.add_argument("--source", help="Upstream tool dialect (DANTE / DANTE_TIR / RepeatMasker / ...)")


def _cmd_canonicalise(args: argparse.Namespace) -> int:
    vocab = load_vocabulary(args.vocabulary)
    stream = sys.stdin if args.input == "-" else open(args.input)
    try:
        for line in stream:
            s = line.rstrip("\n")
            if not s:
                print()
                continue
            try:
                out = canonicalise(s, source=args.source, vocab=vocab, validate=not args.no_validate)
            except (UnknownClassification, ValueError) as e:
                print(f"ERROR: {e}", file=sys.stderr)
                return 2
            print(out)
    finally:
        if stream is not sys.stdin:
            stream.close()
    return 0


def _cmd_canonicalise_fasta(args: argparse.Namespace) -> int:
    vocab = load_vocabulary(args.vocabulary)
    try:
        n = canonicalise_fasta_headers(args.input, args.output, source=args.source, vocab=vocab, validate=not args.no_validate)
    except (UnknownClassification, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    print(f"Rewrote {n} FASTA header(s)", file=sys.stderr)
    return 0


def _cmd_canonicalise_gff(args: argparse.Namespace) -> int:
    vocab = load_vocabulary(args.vocabulary)
    try:
        n = canonicalise_gff3_attribute(
            args.input, args.output, attribute=args.attribute,
            source=args.source, vocab=vocab, validate=not args.no_validate,
        )
    except (UnknownClassification, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    print(f"Rewrote {n} GFF3 feature line(s)", file=sys.stderr)
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    vocab = load_vocabulary(args.vocabulary)
    if args.mode == "fasta":
        values = collect_fasta_classes(args.input)
    elif args.mode == "gff3":
        if not args.attribute:
            print("--attribute is required for gff3 mode", file=sys.stderr)
            return 2
        values = collect_gff3_attribute(args.input, args.attribute)
    else:
        values = [line.rstrip("\n") for line in open(args.input) if line.strip()]
    errors = validate_values(values, source=args.source, vocab=vocab)
    if errors:
        for raw, reason in errors:
            print(f"INVALID\t{raw}\t{reason}", file=sys.stderr)
        print(f"FAILED: {len(errors)} invalid classification(s) in {args.input}", file=sys.stderr)
        return 1
    print(f"OK: all classifications in {args.input} are canonical", file=sys.stderr)
    return 0


def _cmd_list_canonical(args: argparse.Namespace) -> int:
    vocab = load_vocabulary(args.vocabulary)
    for name in iter_canonical(vocab):
        print(name)
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="CARP repeat-classification normaliser")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("canonicalise", help="Canonicalise one classification per line")
    _add_common(p)
    p.add_argument("--input", default="-", help="Input file (one classification per line), default stdin")
    p.add_argument("--no-validate", action="store_true")
    p.set_defaults(func=_cmd_canonicalise)

    p = sub.add_parser("canonicalise-fasta-headers", help="Rewrite #class suffix on FASTA headers")
    _add_common(p)
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--no-validate", action="store_true")
    p.set_defaults(func=_cmd_canonicalise_fasta)

    p = sub.add_parser("canonicalise-gff3-attribute", help="Rewrite one attribute on GFF3 lines")
    _add_common(p)
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--attribute", required=True, help="GFF3 attribute to rewrite (e.g. Classification, Name)")
    p.add_argument("--no-validate", action="store_true")
    p.set_defaults(func=_cmd_canonicalise_gff)

    p = sub.add_parser("validate", help="Fail if any classification in file is not canonical")
    _add_common(p)
    p.add_argument("input")
    p.add_argument("--mode", choices=("fasta", "gff3", "lines"), default="lines")
    p.add_argument("--attribute", help="GFF3 attribute to scan (required for mode=gff3)")
    p.set_defaults(func=_cmd_validate)

    p = sub.add_parser("list-canonical", help="Print every canonical classification")
    _add_common(p)
    p.set_defaults(func=_cmd_list_canonical)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
