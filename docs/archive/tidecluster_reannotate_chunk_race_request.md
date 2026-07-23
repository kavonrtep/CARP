# TideCluster `tc_reannotate`: silent truncation on a RepeatMasker library-build race

## Summary

On the **first** `tc_reannotate` run in a freshly-built conda environment, the
chunked/pooled RepeatMasker runner can **silently drop a large fraction of the
genome**. Its parallel pool launches several RepeatMasker instances at once;
each one, on first use, tries to build RepeatMasker's shared
`Libraries/general/*.lib` blast databases via `makeblastdb`. Concurrent
instances **race on that build**, some error out and mask **nothing**, and
`tc_reannotate` **merges only the surviving chunks without failing** — emitting a
truncated GFF3 with a zero exit status.

## Evidence (TideCluster 1.16.2, A. thaliana T2T, 144 Mb, fresh env)

Same reduced dimer library + genome in every case (md5-identical inputs),
`--sensitivity default --chunk_size 50000000`:

| run | context | RM-on-TideCluster coverage |
|-----|---------|----------------------------|
| pipeline reannotate | fresh env, pool of 14, **first** RM use | **9.9 Mb** (≈42%) |
| standalone rerun | same env, **after** `is.lib` already built | **23.5 Mb** |

The truncated run's log contains, three times:

```
Building general libraries in: .../share/RepeatMasker/Libraries//general
RepeatMasker::createLib(): Error invoking .../bin//makeblastdb on file
    .../share/RepeatMasker/Libraries//general/is.lib.
```

Only ~2 of 5 genome chunks reached the masking stage (`identifying matches ...
batch N of M`); the other ~3 errored during library setup and contributed no
hits. `tc_reannotate` still exited 0 and wrote the ~9.9 Mb (2-chunk) result.
A rerun once `is.lib` exists processes all 5 chunks → 23.5 Mb. So ~58% of the
genome's tandem reannotation was silently lost on the first run.

Note: the command uses `-no_is`, but RepeatMasker still *builds* the general
libraries during setup, so `-no_is` does not avoid the race.

## Root cause (two independent defects, both in `tc_reannotate`)

1. **The race.** The chunked/pooled runner fans out concurrent RepeatMasker
   instances that each attempt to build the shared, per-installation
   `Libraries/general` blast DBs on first use. Nothing serialises that
   one-time setup.
2. **The silent truncation (the more dangerous one).** Per-chunk RepeatMasker
   failures are not detected: the merge step consumes whatever chunk outputs
   exist and the tool returns success, so downstream consumers receive a
   truncated annotation indistinguishable from a complete one.

## Requested fixes

1. **Build the RepeatMasker libraries once, before the pool.** Run a single
   serial RepeatMasker/`makeblastdb` warm-up (or an explicit library build) so
   the parallel chunks all reuse the prepared `Libraries/general`, eliminating
   the race.
2. **Fail on any chunk's RepeatMasker error.** Check each chunk's RepeatMasker
   exit status (and/or the presence of its expected `.out`), and make
   `tc_reannotate` exit **non-zero** if any chunk failed — never merge a partial
   set silently. Optionally retry a failed chunk once before failing.
3. **(Optional) Completeness assertion.** Before merging, assert that every
   input chunk produced an output; abort with a clear message otherwise.

## Impact / who hits it

Any first `tc_reannotate` run in a fresh environment with a chunked genome and
`--cpu > 1`. Deployments that pre-build the RepeatMasker libraries at
image-build time (e.g. a prepared container) avoid the *race*, but the
*silent-truncation-on-chunk-error* is a latent data-integrity risk regardless of
how the libraries got built.

CARP pins the TideCluster version and consumes `RM_on_TideCluster_Library.gff3`
directly (for satellite annotation and the superfamily-aware array recovery in
`scripts/tc_reannotate_sf_filter.py`), so a silent truncation here propagates
straight into the unified repeat annotation.
