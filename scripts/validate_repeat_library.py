#!/usr/bin/env python3
"""Validate a user repeat library for CARP and convert it to canonical form.

CARP takes two optional user libraries, both FASTA with RepeatMasker-style
headers ``>id#classification``:

* ``tandem_repeat_library`` (``--type tandem``) -- every class must be
  ``Satellite[/<name>]`` or an ``rDNA/...`` class. The user's word that the
  sequences are tandem repeats is trusted (a library entry is often a single
  monomer, so tandem structure is not re-checked); only the header is fixed:
  ``#CEN180`` -> ``#Satellite/CEN180``, ``#45S_rDNA`` -> ``#rDNA/45S_rDNA``.
* ``custom_library`` (``--type custom``) -- every class must be a canonical
  CARP classification (``classification_vocabulary.yaml``). Tool-native and
  RepeatMasker/Dfam names are translated where the meaning is unambiguous
  (``LTR/Gypsy`` -> ``Class_I/LTR/Ty3_gypsy``, ``Athila`` ->
  ``Class_I/LTR/Ty3_gypsy/non-chromovirus/OTA/Athila``); anything else is
  reported for a manual fix, never guessed.

Outputs (``-o PREFIX``, or each path set explicitly with ``--out-fasta`` /
``--out-tsv`` / ``--out-html``, as a Galaxy wrapper does): the converted
library FASTA, a TSV with one row per sequence, and a self-contained HTML
report with the formatting guideline.

The library FASTA is written only when every sequence is usable. With
``--drop-invalid`` the unusable sequences are left out instead (an opt-in:
dropping sequences silently changes what the user annotates with). Otherwise
no library is produced: with ``-o`` the file is absent; with an explicit
``--out-fasta`` it is created EMPTY, because a Galaxy tool must create every
declared output -- the exit status then says it is not a library.

Exit status: 0 = a library was written (valid as is, converted, or filtered;
the report says which), 2 = no library -- some sequences need a manual fix,
3 = unreadable input. Non-zero means "do not use the FASTA", which is also how
Galaxy reads it, so a successful conversion never shows as a failed job.
"""
from __future__ import annotations

import argparse
import html
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import yaml

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
import classification  # noqa: E402

sys.path.insert(0, str(_HERE.parent))
try:
    from version import __version__ as CARP_VERSION
except ImportError:  # container: version.py sits in /opt/pipeline
    sys.path.insert(0, "/opt/pipeline")
    try:
        from version import __version__ as CARP_VERSION  # type: ignore
    except ImportError:
        CARP_VERSION = "unknown"

# IUPAC nucleotide codes. Anything else (E, F, L, P, Q, *, digits, ...) means
# the sequence is not DNA -- most often a protein FASTA given by mistake.
_DNA_OK = set("ACGTRYSWKMBDHVN")
_GAP = set("-.")
_ID_BAD = re.compile(r"[^A-Za-z0-9._:\-+|()\[\]]")
_NAME_BAD = re.compile(r"[^A-Za-z0-9._\-+]")

STATUS_UNCHANGED = "unchanged"
STATUS_CONVERTED = "converted"
STATUS_ERROR = "error"
STATUS_DROPPED = "dropped"


@dataclass
class Record:
    index: int
    raw_header: str
    seq: str
    in_id: str = ""
    in_class: str | None = None
    out_id: str = ""
    out_class: str | None = None
    status: str = STATUS_UNCHANGED
    rule: str = ""
    changes: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    @property
    def out_header(self) -> str:
        return f"{self.out_id}#{self.out_class}" if self.out_class else self.out_id


# ---------------------------------------------------------------------------
# FASTA reading
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> tuple[list[Record], list[str]]:
    """Parse FASTA leniently; return records and file-level notes."""
    notes: list[str] = []
    data = path.read_bytes()
    if b"\r" in data:
        notes.append("Windows line endings (CRLF) were converted to Unix (LF).")
    text = data.decode("utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")
    records: list[Record] = []
    header: str | None = None
    chunks: list[str] = []
    preamble = False
    for line in text.split("\n"):
        if line.startswith(">"):
            if header is not None:
                records.append(Record(len(records) + 1, header, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        elif header is None:
            if line.strip():
                preamble = True
        else:
            chunks.append("".join(line.split()))
    if header is not None:
        records.append(Record(len(records) + 1, header, "".join(chunks)))
    if preamble:
        notes.append("Text before the first '>' header was ignored.")
    return records, notes


# ---------------------------------------------------------------------------
# Classification rules
# ---------------------------------------------------------------------------

def _tok(token: str) -> str:
    # A token bounded by anything that is not a letter, digit or '.'
    # ('_', '-', '/' and string ends count as boundaries).
    return r"(?<![A-Za-z0-9.])" + token + r"(?![A-Za-z0-9])"


# Order matters: subunits before units, and 45S before 5S ("45S" contains
# "5S"; the boundary rule also stops "45S" matching as "5S").
_RDNA_RULES: list[tuple[re.Pattern, str]] = [
    (re.compile(_tok(r"(18S|SSU)"), re.I), "rDNA/45S_rDNA/18S"),
    (re.compile(_tok(r"(25S|26S|28S|LSU)"), re.I), "rDNA/45S_rDNA/25S"),
    (re.compile(_tok(r"5\.8S"), re.I), "rDNA/45S_rDNA/5.8S"),
    (re.compile(_tok(r"ITS-?1"), re.I), "rDNA/45S_rDNA/ITS1"),
    (re.compile(_tok(r"ITS-?2"), re.I), "rDNA/45S_rDNA/ITS2"),
    (re.compile(_tok(r"IGS"), re.I), "rDNA/45S_rDNA/IGS"),
    (re.compile(_tok(r"(45S|35S|47S|48S)"), re.I), "rDNA/45S_rDNA"),
    (re.compile(_tok(r"5S"), re.I), "rDNA/5S_rDNA"),
    (re.compile(r"rDNA|rRNA|ribosomal", re.I), "rDNA"),
]


def rdna_class(s: str) -> str | None:
    for pat, target in _RDNA_RULES:
        if pat.search(s):
            return target
    return None


def sanitise_name(s: str) -> str:
    return _NAME_BAD.sub("_", s).strip("_") or "unnamed"


class Resolver:
    def __init__(self, vocab_path: Path | None):
        self.vocab_path = vocab_path or classification._find_default_vocabulary()
        self.vocab = classification.load_vocabulary(self.vocab_path)
        raw = yaml.safe_load(Path(self.vocab_path).read_text())
        self.aliases = [(re.compile(a["pattern"], re.I), a["pattern"], a["canonical"])
                        for a in raw.get("library_class_aliases") or []]
        for _, pat, target in self.aliases:
            if not self.vocab.is_canonical(target):
                raise ValueError(f"library_class_aliases target {target!r} "
                                 f"(pattern {pat!r}) is not canonical")
        self.canonical_sorted = sorted(self.vocab.canonical)
        self._lower = {c.lower(): c for c in self.canonical_sorted}

    def canonical(self, s: str) -> str | None:
        """Exact canonical form via the shared normaliser, else None."""
        try:
            return classification.canonicalise(s, source=None, vocab=self.vocab)
        except (classification.UnknownClassification, ValueError):
            return None

    def satellite_prefixed(self, s: str) -> str | None:
        m = re.fullmatch(r"satellites?(?:[/_\-](.+))?", s, re.I)
        if not m:
            return None
        return "Satellite" + (f"/{sanitise_name(m.group(1))}" if m.group(1) else "")

    def by_lineage(self, s: str) -> str | None:
        """Unique canonical class whose path ENDS with s (case-insensitive)."""
        t = s
        for alias, canon in self.vocab.leaf_aliases.items():
            t = t.replace(alias, canon)
        t = t.strip("/").lower()
        if t in self._lower:
            return self._lower[t]
        hits = [c for c in self.canonical_sorted if c.lower().endswith("/" + t)]
        return hits[0] if len(hits) == 1 else None

    def by_dictionary(self, s: str) -> tuple[str, str] | None:
        for rx, pat, target in self.aliases:
            if rx.fullmatch(s):
                return target, pat
        return None


def resolve_tandem(r: Record, R: Resolver) -> None:
    raw = r.in_class
    if raw is None:
        rd = rdna_class(r.in_id)
        if rd:
            r.out_class, r.rule = rd, "no class; rDNA name in the ID"
        else:
            r.out_class, r.rule = f"Satellite/{sanitise_name(r.in_id)}", "no class; ID used as satellite name"
        r.changes.append(f"class added: {r.out_class}")
        return
    s = raw
    if s.endswith("?"):
        s = s.rstrip("?")
        r.changes.append("RepeatMasker uncertainty mark '?' removed")
    c = R.canonical(s)
    if c is not None:
        if c.startswith(("Satellite", "rDNA")):
            r.out_class, r.rule = c, "already canonical"
            if c != raw:
                r.changes.append(f"class {raw} -> {c}")
            return
        r.errors.append(
            f"'{c}' is a valid CARP class but not a tandem-repeat class. A tandem "
            f"library may only hold Satellite or rDNA sequences; move this "
            f"sequence to custom_library.")
        return
    sat = R.satellite_prefixed(s)
    if sat:
        r.out_class, r.rule = sat, "Satellite spelling normalised"
    elif (rd := rdna_class(s)):
        r.out_class, r.rule = rd, "rDNA name recognised"
    elif (hit := R.by_dictionary(s)) and not hit[0].startswith(("Satellite", "rDNA")):
        r.errors.append(
            f"'{raw}' is a transposable-element class ({hit[0]}), not a tandem "
            f"repeat; move this sequence to custom_library.")
        return
    else:
        r.out_class, r.rule = f"Satellite/{sanitise_name(s)}", "Satellite/ prefix added"
    r.changes.append(f"class {raw} -> {r.out_class}")


def resolve_custom(r: Record, R: Resolver) -> None:
    raw = r.in_class
    if raw is None:
        r.errors.append("No classification. Add one after '#', e.g. "
                        f">{r.in_id}#Class_I/LTR/Ty1_copia/Ale")
        return
    s = raw
    if s.endswith("?"):
        s = s.rstrip("?")
        r.changes.append("RepeatMasker uncertainty mark '?' removed")
    c = R.canonical(s)
    rule = "already canonical"
    if c is None:
        c = R.satellite_prefixed(s)
        rule = "Satellite spelling normalised"
    if c is None:
        c = R.by_lineage(s)
        rule = "CARP class name matched"
    if c is None and (hit := R.by_dictionary(s)):
        c, pat = hit
        rule = f"RepeatMasker/Dfam name (rule '{pat}')"
    if c is None and (rd := rdna_class(s)):
        c, rule = rd, "rDNA name recognised"
    if c is None:
        r.errors.append(
            f"'{raw}' is not a CARP classification and has no unambiguous "
            f"translation. Replace it with a class from the guideline, e.g. "
            f"Class_I/LTR/Ty3_gypsy or Class_II/Subclass_1/TIR/hAT.")
        return
    r.out_class, r.rule = c, rule
    if c != raw:
        r.changes.append(f"class {raw} -> {c}")


# ---------------------------------------------------------------------------
# Per-record checks
# ---------------------------------------------------------------------------

def check_record(r: Record, kind: str, R: Resolver) -> None:
    head = r.raw_header
    parts = head.split(None, 1)
    token = parts[0] if parts else ""
    if len(parts) > 1:
        r.changes.append("description after the first space dropped "
                         "(RepeatMasker reads only the first word)")
    if "#" in token:
        r.in_id, cls = token.split("#", 1)
        r.in_class = cls or None
        if r.in_class and "#" in r.in_class:
            r.errors.append("More than one '#' in the header; use exactly one: >id#class")
    else:
        r.in_id = token
    if not r.in_id:
        r.errors.append("Empty sequence ID before '#'.")
    r.out_id = _ID_BAD.sub("_", r.in_id)
    if r.out_id != r.in_id:
        r.changes.append(f"ID characters replaced: {r.in_id} -> {r.out_id}")

    if not r.seq:
        r.errors.append("Empty sequence.")
    else:
        # Keep the user's letter case; only these edits change the sequence.
        if set(r.seq) & _GAP:
            r.seq = "".join(ch for ch in r.seq if ch not in _GAP)
            r.changes.append("alignment gap characters (- .) removed")
        if "U" in r.seq.upper():
            r.seq = r.seq.replace("U", "T").replace("u", "t")
            r.changes.append("RNA base U converted to T")
        bad = sorted(set(r.seq.upper()) - _DNA_OK)
        if bad:
            r.errors.append(f"Not a DNA sequence (contains {''.join(bad)}); "
                            f"a protein FASTA cannot be used as a repeat library.")

    (resolve_tandem if kind == "tandem" else resolve_custom)(r, R)

    if r.errors:
        r.status = STATUS_ERROR
    elif r.changes:
        r.status = STATUS_CONVERTED


def dedupe_ids(records: list[Record]) -> None:
    seen: dict[str, int] = {}
    for r in records:
        if r.status == STATUS_ERROR:
            continue
        n = seen.get(r.out_id, 0) + 1
        seen[r.out_id] = n
        if n > 1:
            new = f"{r.out_id}_dup{n}"
            while new in seen:
                n += 1
                new = f"{r.out_id}_dup{n}"
            seen[new] = 1
            r.changes.append(f"duplicate ID renamed: {r.out_id} -> {new}")
            r.out_id = new
            r.status = STATUS_CONVERTED


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

def _atomic_write(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def write_fasta(path: Path, records: list[Record]) -> None:
    out = []
    for r in records:
        if r.status in (STATUS_ERROR, STATUS_DROPPED):
            continue
        out.append(f">{r.out_header}\n")
        out.extend(r.seq[i:i + 60] + "\n" for i in range(0, len(r.seq), 60))
    _atomic_write(path, "".join(out))


def write_tsv(path: Path, records: list[Record]) -> None:
    cols = ["index", "input_id", "input_class", "output_id", "output_class",
            "status", "length_bp", "rule", "changes", "errors"]
    clean = lambda x: (x or "").replace("\t", " ").replace("\n", " ")
    rows = ["\t".join(cols)]
    for r in records:
        rows.append("\t".join(clean(str(v)) for v in (
            r.index, r.in_id, r.in_class or "", r.out_id if r.status != STATUS_ERROR else "",
            r.out_class if r.status != STATUS_ERROR else "", r.status, len(r.seq),
            r.rule, "; ".join(r.changes), "; ".join(r.errors))))
    _atomic_write(path, "\n".join(rows) + "\n")


_CSS = """
:root{--bg:#fff;--fg:#1d1d1f;--muted:#5f6368;--line:#e3e3e8;--card:#f7f7f9;
--ok:#1e7f43;--okbg:#e6f4ea;--conv:#1a5fb4;--convbg:#e7effa;--err:#b3261e;--errbg:#fbe9e7;
--drop:#8a5a00;--dropbg:#fff3d6;--code:#f1f1f4}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#17181a;--fg:#e8e8ea;
--muted:#a0a3a8;--line:#2e3034;--card:#1f2124;--ok:#6ccf8e;--okbg:#17301f;--conv:#8ab4f8;
--convbg:#172437;--err:#f28b82;--errbg:#3a1a18;--drop:#f2c56b;--dropbg:#352a12;--code:#26282c}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);
font:15px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:1100px;margin:0 auto;padding:24px 16px 48px}
h1{font-size:1.5em;margin:0 0 4px}h2{font-size:1.15em;margin:32px 0 8px}
.sub{color:var(--muted);margin:0 0 20px}
.verdict{border-radius:10px;padding:14px 18px;margin:0 0 20px;font-weight:600}
.v-ok{background:var(--okbg);color:var(--ok)}.v-conv{background:var(--convbg);color:var(--conv)}
.v-err{background:var(--errbg);color:var(--err)}
.verdict p{font-weight:400;color:var(--fg);margin:6px 0 0}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px}
.tile{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 14px}
.tile b{display:block;font-size:1.6em}.tile span{color:var(--muted);font-size:.9em}
.filters{margin:8px 0}.filters button{font:inherit;background:var(--card);color:var(--fg);
border:1px solid var(--line);border-radius:999px;padding:3px 12px;margin:0 6px 6px 0;cursor:pointer}
.filters button.on{border-color:var(--fg)}
.tw{overflow-x:auto;border:1px solid var(--line);border-radius:10px}
table{border-collapse:collapse;width:100%;font-size:.92em}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--line);vertical-align:top}
th{background:var(--card);position:sticky;top:0}tr:last-child td{border-bottom:0}
code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.92em;
background:var(--code);padding:1px 5px;border-radius:4px;word-break:break-all}
.pill{display:inline-block;border-radius:999px;padding:1px 9px;font-size:.85em;font-weight:600}
.unchanged{background:var(--okbg);color:var(--ok)}.converted{background:var(--convbg);color:var(--conv)}
.error{background:var(--errbg);color:var(--err)}.dropped{background:var(--dropbg);color:var(--drop)}
td ul{margin:0;padding-left:18px}.errtxt{color:var(--err)}.muted{color:var(--muted)}
details{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 14px;margin:8px 0}
summary{cursor:pointer;font-weight:600}.cols{columns:2 320px;font-size:.9em}
.cols code{display:inline-block;margin:1px 0}
"""

_JS = """
// Buttons are created here, not in the markup: without scripts the page still
// shows the full table, with no dead controls.
const fd=document.querySelector('.filters');
[['all','All'],['changed','Changed'],['error','Errors']].forEach(([f,l],i)=>{
 const b=document.createElement('button');b.dataset.f=f;b.textContent=l+' ('+fd.dataset[f]+')';
 if(!i)b.classList.add('on');fd.appendChild(b);});
document.querySelectorAll('.filters button').forEach(b=>b.addEventListener('click',()=>{
 document.querySelectorAll('.filters button').forEach(x=>x.classList.remove('on'));b.classList.add('on');
 const f=b.dataset.f;document.querySelectorAll('#recs tbody tr').forEach(tr=>{
 tr.style.display=(f==='all'||(f==='error'&&(tr.dataset.s==='error'||tr.dataset.s==='dropped'))||(f==='changed'&&tr.dataset.s!=='unchanged'))?'':'none';});}));
"""


def _guideline_html(kind: str, R: Resolver) -> str:
    e = html.escape
    if kind == "tandem":
        body = """
<p>Each header is <code>&gt;ID#CLASS</code>: a unique ID without spaces, one <code>#</code>, and a class.
Only two kinds of class are allowed in a tandem repeat library:</p>
<ul>
<li><code>Satellite</code> or <code>Satellite/&lt;name&gt;</code> for any satellite or other tandem repeat.
The name is free text (letters, digits, <code>. _ - +</code>), e.g. <code>&gt;X04322.1#Satellite/CEN180</code>.</li>
<li>An rDNA class: <code>rDNA/45S_rDNA</code>, <code>rDNA/5S_rDNA</code>, or a 45S subunit such as
<code>rDNA/45S_rDNA/18S</code>, <code>rDNA/45S_rDNA/ITS1</code>, <code>rDNA/45S_rDNA/IGS</code>.</li>
</ul>
<p>Transposable elements (LTR retrotransposons, DNA transposons, LINEs, ...) do not belong here,
even when they occur in tandem. Put them in <code>custom_library</code> with their TE class.</p>"""
    else:
        groups: dict[str, list[str]] = {}
        for c in R.canonical_sorted:
            groups.setdefault(c.split("/")[0], []).append(c)
        lists = "".join(
            f"<details><summary>{e(g)} ({len(v)})</summary><div class='cols'>"
            + "<br>".join(f"<code>{e(c)}</code>" for c in v) + "</div></details>"
            for g, v in sorted(groups.items()))
        body = f"""
<p>Each header is <code>&gt;ID#CLASS</code>: a unique ID without spaces, one <code>#</code>, and a
CARP classification written with <code>/</code> between levels, from the most general to the most
specific, e.g. <code>&gt;myCopia_1#Class_I/LTR/Ty1_copia/Ale</code>. Use the most specific level you
are sure of; <code>Class_I/LTR/Ty3_gypsy</code> is better than a wrong lineage.</p>
<p>Also accepted: <code>Satellite</code> or <code>Satellite/&lt;name&gt;</code>, <code>Simple_repeat</code>,
<code>Low_complexity</code>, <code>Unknown</code>. Common RepeatMasker/Dfam names
(<code>LTR/Gypsy</code>, <code>DNA/hAT-Ac</code>, <code>LINE/L1</code>, ...) and the three-letter codes of
Wicker et al. (2007) (<code>RLC</code>, <code>RLG</code>, <code>DTC</code>, <code>DHH</code>, <code>RXX</code>,
<code>DXX</code>, ...; a family name after the code, e.g. <code>RLC_Angela</code>, is allowed but only the
superfamily is used) are converted automatically;
classes that CARP does not model (tRNA, snRNA, ...) cannot be used.</p>
<p>All valid classes:</p>{lists}"""
    return body


def write_html(path: Path, records: list[Record], kind: str, src: Path,
               notes: list[str], R: Resolver, fasta_written: bool, out_fasta: Path) -> None:
    e = html.escape
    n = len(records)
    cnt = {s: sum(r.status == s for r in records)
           for s in (STATUS_UNCHANGED, STATUS_CONVERTED, STATUS_ERROR, STATUS_DROPPED)}
    lib = "tandem repeat library" if kind == "tandem" else "custom repeat library"
    cfg = "tandem_repeat_library" if kind == "tandem" else "custom_library"
    if cnt[STATUS_ERROR]:
        cls, title = "v-err", "Cannot be converted automatically"
        msg = (f"{cnt[STATUS_ERROR]} of {n} sequences need a manual fix (listed below as "
               f"<b>error</b>). No library file was written. Fix them in the input and run the "
               f"check again, or use <code>--drop-invalid</code> to write a library without them.")
    elif cnt[STATUS_CONVERTED] or cnt[STATUS_DROPPED] or notes:
        cls, title = "v-conv", "Converted &mdash; ready to use"
        msg = (f"{cnt[STATUS_CONVERTED]} of {n} sequences were changed"
               + (f" and {cnt[STATUS_DROPPED]} unusable ones left out" if cnt[STATUS_DROPPED] else "")
               + f". Use <code>{e(out_fasta.name)}</code> as <code>{cfg}</code> in your config.")
    else:
        cls, title = "v-ok", "Valid &mdash; no changes needed"
        msg = f"All {n} sequences already follow the CARP format. You can use the input file as is."
    note_html = ("<ul>" + "".join(f"<li>{e(x)}</li>" for x in notes) + "</ul>") if notes else ""

    rows = []
    for r in records:
        what = [f"<li>{e(c)}</li>" for c in r.changes] + \
               [f"<li class='errtxt'>{e(x)}</li>" for x in r.errors]
        out = (f"<code>&gt;{e(r.out_header)}</code>" if r.status not in (STATUS_ERROR, STATUS_DROPPED)
               else "<span class='muted'>&mdash;</span>")
        rows.append(
            f"<tr data-s='{r.status}'><td>{r.index}</td>"
            f"<td><code>&gt;{e(r.raw_header)}</code></td><td>{out}</td>"
            f"<td><span class='pill {r.status}'>{r.status}</span></td>"
            f"<td>{len(r.seq):,}</td>"
            f"<td>{('<ul>' + ''.join(what) + '</ul>') if what else '<span class=muted>none</span>'}"
            f"{('<div class=muted>' + e(r.rule) + '</div>') if r.rule and r.status == STATUS_CONVERTED else ''}</td></tr>")

    tiles = "".join(f"<div class='tile'><b>{v}</b><span>{k}</span></div>" for k, v in (
        ("sequences", n), ("unchanged", cnt[STATUS_UNCHANGED]), ("converted", cnt[STATUS_CONVERTED]),
        ("need a manual fix", cnt[STATUS_ERROR] + cnt[STATUS_DROPPED])))
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Library Check</title><style>{_CSS}</style></head><body><main>
<h1>CARP library check</h1>
<p class="sub">{e(src.name)} &middot; checked as a <b>{lib}</b> &middot; CARP {e(CARP_VERSION)}</p>
<div class="verdict {cls}">{title}<p>{msg}</p>{note_html}</div>
<div class="tiles">{tiles}</div>
<h2>Sequences</h2>
<div class="filters" data-all="{n}" data-changed="{cnt[STATUS_CONVERTED] + cnt[STATUS_ERROR] + cnt[STATUS_DROPPED]}"
 data-error="{cnt[STATUS_ERROR] + cnt[STATUS_DROPPED]}"></div>
<div class="tw"><table id="recs"><thead><tr><th>#</th><th>Input header</th><th>Output header</th>
<th>Status</th><th>Length (bp)</th><th>What was done</th></tr></thead>
<tbody>{''.join(rows)}</tbody></table></div>
<h2>Format guideline for a {lib}</h2>{_guideline_html(kind, R)}
</main><script>{_JS}</script></body></html>
"""
    _atomic_write(path, doc)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def run(input_path: Path, kind: str, out_fasta: Path, out_tsv: Path, out_html: Path,
        drop_invalid: bool = False, vocab_path: Path | None = None,
        empty_fasta_on_error: bool = False) -> int:
    R = Resolver(vocab_path)
    try:
        records, notes = read_fasta(input_path)
    except OSError as exc:
        print(f"ERROR: cannot read {input_path}: {exc}", file=sys.stderr)
        return 3
    if not records:
        print(f"ERROR: {input_path} contains no FASTA records.", file=sys.stderr)
        return 3
    for r in records:
        check_record(r, kind, R)
    dedupe_ids(records)

    n_err = sum(r.status == STATUS_ERROR for r in records)
    if n_err and drop_invalid:
        for r in records:
            if r.status == STATUS_ERROR:
                r.status = STATUS_DROPPED
    write_fasta_ok = n_err == 0 or drop_invalid

    # Self-check: everything written must pass the pipeline's own gate.
    for r in records:
        if r.status in (STATUS_UNCHANGED, STATUS_CONVERTED):
            assert r.out_class and R.vocab.is_canonical(r.out_class), r.out_header
            if kind == "tandem":
                assert r.out_class.startswith(("Satellite", "rDNA")), r.out_header

    for p in (out_fasta, out_tsv, out_html):
        p.parent.mkdir(parents=True, exist_ok=True)
    if write_fasta_ok:
        write_fasta(out_fasta, records)
    elif empty_fasta_on_error:
        _atomic_write(out_fasta, "")   # declared Galaxy output must exist
    elif out_fasta.exists():
        out_fasta.unlink()   # never leave a stale library from an earlier run
    write_tsv(out_tsv, records)
    write_html(out_html, records, kind, input_path, notes, R, write_fasta_ok, out_fasta)

    n = len(records)
    cnt = {s: sum(r.status == s for r in records)
           for s in (STATUS_UNCHANGED, STATUS_CONVERTED, STATUS_ERROR, STATUS_DROPPED)}
    print(f"{n} sequences: {cnt[STATUS_UNCHANGED]} unchanged, {cnt[STATUS_CONVERTED]} converted, "
          f"{cnt[STATUS_ERROR]} need a manual fix, {cnt[STATUS_DROPPED]} dropped.")
    for r in records:
        for x in r.errors:
            print(f"  ERROR  #{r.index} >{r.raw_header}: {x}")
    print(f"Report: {out_html} (and {out_tsv})")
    if n_err and not drop_invalid:
        print("No library written: fix the errors above, or rerun with --drop-invalid.")
        return 2
    print(f"Library: {out_fasta}")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Check a repeat library for CARP and convert it to the CARP format.",
        epilog="Exit status: 0 library written (valid, converted or filtered), "
               "2 no library - sequences need a manual fix, 3 unreadable input.")
    ap.add_argument("input", type=Path, help="library FASTA")
    ap.add_argument("-t", "--type", required=True, choices=["tandem", "custom"],
                    help="tandem = tandem_repeat_library, custom = custom_library")
    ap.add_argument("-o", "--output-prefix", type=Path,
                    help="writes PREFIX.fasta, PREFIX.tsv, PREFIX.html "
                         "[default: <input stem>.carp_<type>]")
    ap.add_argument("--out-fasta", type=Path, help="library FASTA path (overrides -o); "
                    "created empty when no library can be written")
    ap.add_argument("--out-tsv", type=Path, help="per-sequence TSV path (overrides -o)")
    ap.add_argument("--out-html", type=Path, help="HTML report path (overrides -o)")
    ap.add_argument("--drop-invalid", action="store_true",
                    help="write the library without the sequences that cannot be converted "
                         "(off by default: those sequences are then simply missing)")
    ap.add_argument("--vocabulary", type=Path, help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    prefix = a.output_prefix or a.input.with_name(f"{a.input.stem}.carp_{a.type}")
    side = lambda ext: prefix.with_name(prefix.name + ext)
    return run(a.input, a.type,
               a.out_fasta or side(".fasta"), a.out_tsv or side(".tsv"),
               a.out_html or side(".html"), a.drop_invalid, a.vocabulary,
               empty_fasta_on_error=a.out_fasta is not None)


if __name__ == "__main__":
    sys.exit(main())
