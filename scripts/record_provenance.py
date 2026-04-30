#!/usr/bin/env python3
"""Record per-run provenance for a CARP pipeline run.

A pipeline run produces ``<output_dir>/run_provenance.json``: a single
file describing exactly what was run, by whom, where, with which tool
versions. Three phases write to it:

  1. ``init`` — called by ``run_pipeline.py`` *before* invoking
     snakemake. Records pipeline version (from ``version.py``), git
     SHA + dirty flag, host, user, container detection, the user's
     config, and ``run_started`` timestamp. ``exit_status`` is set to
     ``"running"``. This phase is robust to snakemake startup
     failures: even if conda env solve dies, the file exists.
  2. ``extend-with-envs`` — called by the ``record_tool_versions``
     Snakemake rule (commit 4 of the versioning rollout) once envs
     have been created. Fills in the ``envs`` section with
     per-environment primary-package versions (hybrid: requested from
     envs/*.yaml and resolved from the actual conda env contents).
  3. ``finalise`` — called by ``run_pipeline.py`` after snakemake
     exits, with the snakemake exit code. Sets ``run_finished`` and
     ``exit_status`` (``"completed"`` or ``"failed"``).

The file is the canonical "what produced this output" record:
referenced from the HTML report footer, the Repeat_Annotation_Unified
GFF header, and 0INFO.txt on Metacentrum runs (commit 5).

Schema is at ``schema_version: 1``. Bump on incompatible changes.

Importable as a library (``init_provenance``, ``finalise_provenance``,
``query_env_versions``) from ``run_pipeline.py``; also a CLI for use
from Snakemake rules (``init`` / ``extend-with-envs`` / ``finalise``).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import getpass
import json
import os
import socket
import subprocess
import sys
from pathlib import Path

# Single source of truth for the pipeline version. The script lives in
# scripts/, version.py at the repo root, so step up one level.
_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent))
try:
    from version import __version__ as PIPELINE_VERSION
except ImportError:
    # Container-context fallback: scripts and version.py both copied to
    # /opt/pipeline/, but the relative parent walk above doesn't apply.
    sys.path.insert(0, "/opt/pipeline")
    from version import __version__ as PIPELINE_VERSION  # type: ignore


SCHEMA_VERSION = 1


def _utcnow_iso() -> str:
    """ISO 8601 UTC timestamp, second precision, ``Z`` suffix."""
    return (_dt.datetime.utcnow()
            .replace(microsecond=0).isoformat() + "Z")


def _git_sha() -> tuple[str, bool]:
    """Return ``(short_sha, dirty)``. Falls back to baked-in container
    artefact at ``/etc/carp_git_sha`` if .git is unavailable."""
    try:
        sha = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=_HERE.parent, stderr=subprocess.DEVNULL,
        ).decode().strip()
        # `git diff --quiet` exit 1 == dirty
        rc = subprocess.call(
            ["git", "diff", "--quiet"], cwd=_HERE.parent,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        return sha, (rc != 0)
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Inside the container .git is absent. A future Singularity
        # build can drop the SHA at /etc/carp_git_sha; until then,
        # report "container".
        try:
            return Path("/etc/carp_git_sha").read_text().strip(), False
        except FileNotFoundError:
            return "unknown", False


def _detect_container() -> str | None:
    """Return the container .sif path if we're running inside one, else
    ``None``. Singularity creates ``/.singularity.d/`` in any image."""
    if Path("/.singularity.d").exists():
        # The .sif path itself isn't readable from inside; record what
        # we know. Caller can pass an explicit override via
        # SINGULARITY_CONTAINER env var (Singularity sets it).
        return os.environ.get("SINGULARITY_CONTAINER",
                              "(running inside Singularity)")
    return None


def _filter_config(config: dict) -> dict:
    """Project the user's config down to keys this project recognises.
    Anything else is dropped — provenance is a record, not a dump."""
    keys = (
        "genome_fasta", "output_dir", "custom_library",
        "tandem_repeat_library", "repeatmasker_sensitivity",
        "reduce_library",
        "dante_tir_min_multiplicity",
        "dante_tir_fallback_min_alignments",
        "dante_tir_fallback_min_cluster_size",
        "include_dante_tir_fallback_in_library",
        "dante_tir_fallback_library_min_multiplicity",
    )
    return {k: config[k] for k in keys if k in config}


# ──────────────────────────────────────────────────────────────────────────
# Library API
# ──────────────────────────────────────────────────────────────────────────
def init_provenance(config_path: str | Path,
                    output_dir: str | Path) -> Path:
    """Phase-1 write. Returns the provenance JSON path."""
    import yaml
    config = yaml.safe_load(Path(config_path).read_text()) or {}
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sha, dirty = _git_sha()
    record = {
        "schema_version": SCHEMA_VERSION,
        "pipeline_version": PIPELINE_VERSION,
        "git_sha": sha,
        "git_dirty": dirty,
        "run_started": _utcnow_iso(),
        "run_finished": None,
        "exit_status": "running",
        "host": socket.gethostname(),
        "user": getpass.getuser(),
        "container_sif": _detect_container(),
        "config": _filter_config(config),
        "envs": {},
    }
    out = output_dir / "run_provenance.json"
    _atomic_write_json(out, record)
    return out


def finalise_provenance(output_dir: str | Path,
                        exit_status: str) -> Path:
    """Phase-3 write. ``exit_status`` ∈ {"completed", "failed"}."""
    out = Path(output_dir) / "run_provenance.json"
    if not out.exists():
        # init never ran or output_dir got wiped; emit a stub so a
        # consumer downstream still has *something*.
        record = {
            "schema_version": SCHEMA_VERSION,
            "pipeline_version": PIPELINE_VERSION,
            "exit_status": exit_status,
            "run_finished": _utcnow_iso(),
            "note": "init phase missing; this provenance record is partial.",
        }
    else:
        record = json.loads(out.read_text())
        record["run_finished"] = _utcnow_iso()
        record["exit_status"] = exit_status
    _atomic_write_json(out, record)
    return out


def _atomic_write_json(path: Path, obj: dict) -> None:
    """Write JSON via temp-file + rename so a crash mid-write doesn't
    leave a half-written record on disk."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=False))
    tmp.replace(path)


# ──────────────────────────────────────────────────────────────────────────
# Conda env query (used by commit 4 of the rollout — defined now so the
# library is complete)
# ──────────────────────────────────────────────────────────────────────────
def _parse_yaml_primaries(yaml_path: Path) -> dict[str, str | None]:
    """envs/<name>.yaml → {package_name: pinned_version_or_None}."""
    import yaml
    spec = yaml.safe_load(yaml_path.read_text()) or {}
    deps = spec.get("dependencies", []) or []
    out: dict[str, str | None] = {}
    for dep in deps:
        if not isinstance(dep, str):
            continue  # skip pip:- nested specs
        # "name=version", "name=version=build", or just "name"
        head = dep.split(" ", 1)[0]
        parts = head.split("=", 2)
        name = parts[0].strip()
        version = parts[1].strip() if len(parts) > 1 else None
        out[name] = version
    return out


def _resolve_env_dir(yaml_path: Path,
                     conda_prefix: Path) -> Path | None:
    """Find the ``.snakemake/conda/<hash>/`` (or container baked
    ``/opt/conda/envs/<name>/``) directory for a given envs/*.yaml."""
    # Snakemake-managed: the YAML is copied next to the env dir as
    # ``<hash>.yaml``. Match by content.
    target = yaml_path.read_bytes()
    if conda_prefix.exists():
        for candidate in conda_prefix.iterdir():
            if candidate.suffix == ".yaml" and candidate.read_bytes() == target:
                env_dir = candidate.with_suffix("")
                if env_dir.is_dir():
                    return env_dir
    # Container baked: env name == YAML basename without extension.
    container_dir = Path("/opt/conda/envs") / yaml_path.stem
    if container_dir.is_dir():
        return container_dir
    return None


def _conda_list(env_dir: Path) -> dict[str, str]:
    """Map ``package_name → installed_version`` for every package in
    the env. Reads ``conda-meta/*.json`` directly to avoid spawning
    conda (faster, no PATH dance)."""
    out: dict[str, str] = {}
    meta = env_dir / "conda-meta"
    if not meta.is_dir():
        return out
    for f in meta.glob("*.json"):
        try:
            doc = json.loads(f.read_text())
            name, ver = doc.get("name"), doc.get("version")
            if name and ver:
                out[name] = ver
        except (json.JSONDecodeError, OSError):
            continue
    return out


def query_env_versions(envs_dir: Path,
                       conda_prefix: Path) -> dict[str, dict]:
    """For every ``envs/*.yaml``, return:

    {env_name: {
        "yaml_path": "envs/foo.yaml",
        "conda_prefix": ".snakemake/conda/<hash>/" or "<query failed>",
        "primary_packages": {
            pkg: {"requested": "1.2.3" or None, "resolved": "1.2.3"} ...
        },
    }}
    """
    out: dict[str, dict] = {}
    for y in sorted(envs_dir.glob("*.yaml")):
        primaries = _parse_yaml_primaries(y)
        env_dir = _resolve_env_dir(y, conda_prefix)
        installed = _conda_list(env_dir) if env_dir else {}
        out[y.stem] = {
            "yaml_path": str(y),
            "conda_prefix": str(env_dir) if env_dir else "<not found>",
            "primary_packages": {
                pkg: {
                    "requested": req,
                    "resolved": installed.get(pkg, "<query failed>"),
                }
                for pkg, req in primaries.items()
            },
        }
    return out


def extend_with_envs(provenance_path: Path,
                     envs_dir: Path,
                     conda_prefix: Path) -> Path:
    """Phase-2 write. Replaces the ``envs:`` section."""
    record = json.loads(provenance_path.read_text())
    record["envs"] = query_env_versions(envs_dir, conda_prefix)
    _atomic_write_json(provenance_path, record)
    return provenance_path


# ──────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_init = sub.add_parser("init", help="Phase 1: write initial provenance.")
    p_init.add_argument("--config", required=True, type=Path)
    p_init.add_argument("--output-dir", required=True, type=Path)

    p_extend = sub.add_parser(
        "extend-with-envs",
        help="Phase 2: query conda envs and amend provenance JSON.")
    p_extend.add_argument("--provenance-json", required=True, type=Path)
    p_extend.add_argument("--envs-dir", required=True, type=Path)
    p_extend.add_argument("--conda-prefix", required=True, type=Path)

    p_final = sub.add_parser(
        "finalise", help="Phase 3: stamp run_finished + exit_status.")
    p_final.add_argument("--output-dir", required=True, type=Path)
    p_final.add_argument("--exit-status", required=True,
                         choices=("completed", "failed"))

    args = ap.parse_args()
    if args.cmd == "init":
        out = init_provenance(args.config, args.output_dir)
        print(f"phase 1: wrote {out}")
    elif args.cmd == "extend-with-envs":
        out = extend_with_envs(args.provenance_json,
                               args.envs_dir, args.conda_prefix)
        print(f"phase 2: amended envs in {out}")
    elif args.cmd == "finalise":
        out = finalise_provenance(args.output_dir, args.exit_status)
        print(f"phase 3: finalised {out} ({args.exit_status})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
