# End-to-end check of the new LINE bounds: DO NOT SHIP AS-IS

`dante_line` run twice per genome on the exact inputs of a 1.7.0rc1 baseline —
same code, same masks, only `classification_vocabulary.yaml` swapped between the
rc1 bounds (5' 2000 / 3' 800 flat) and main's (5' 3400 / span 4500).

## The control behaves perfectly

**GCA_937616625.2** (baseline run-000230, 8 of 633 elements capped under rc1):

| | rc1 | new |
|---|---:|---:|
| total extension bp | 94,236 | 97,626 (+3.6%) |
| element length median | 2,163 | 2,163 |
| element length max | 5,020 | 5,020 |
| elements > 7000 bp | 0 | 0 |
| fraction of LINE bp inferred | 0.06 | 0.06 |

Exactly as designed: a small correction where the flat 3' cap was clipping,
nothing else moves.

## The sensitive genome does not

**GCA_973357735.1** (baseline run-000232, 494 of 550 capped under rc1):

| | rc1 | new |
|---|---:|---:|
| total extension bp | 1,362,872 | 2,381,085 (**+74.7%**) |
| element length median | 4,960 | **7,429** |
| element length p99 | 5,886 | 7,499 |
| elements > 7000 bp | 0 | **322 of 550** |
| fraction of LINE bp inferred | 0.53 | **0.66** |
| **median 5' extension** | **2,000** | **3,400** |

**The median 5' extension equals the cap exactly, under both settings.** For most
elements on this genome the flank inference wants more than it is allowed; the
cap is not a safety net, it is the only thing setting the element's length.
Raising it simply makes the elements longer. Two-thirds of LINE base pairs are
now inferred flank rather than observed domain, and 322 of 550 elements exceed
7 kb -- the top of the documented plant LINE range -- where rc1 produced none.

## Why the calibration missed this

The bounds came from 218 loci where BOTH ends were directly observed (poly-A
tail + intact TSD). Those are young elements, and on them the inference
converges: it picks a length below the cap. The genomes and elements where it
does NOT converge are exactly the ones excluded from that set -- I had already
found that on GCA_986270105 the layer extended 19.0% of all elements but 0 of 27
confirmed ones, and flagged the population as possibly unrepresentative. It is.

So the measurement was not wrong; it answered "how long is a LINE whose ends we
can see", which is not the same question as "how far should we extend when we
cannot see them".

## What this means

- **The span rule itself is sound** and the control shows it corrects a real
  defect harmlessly.
- **The bound values cannot be set from converged elements alone.** A cap that
  binds on the median element is doing load-bearing work, and its value is then
  a policy choice about un-converged cases, not a measurement.
- The natural fix is to make the two cases distinguishable rather than to pick a
  different constant: when the applied extension equals the cap, the boundary was
  never determined -- that is a different confidence class from an element whose
  inference stopped on its own. `Extension_capped` already records it per element;
  nothing downstream uses it.

## Recommendation

Do not ship the new bounds in 1.7.0. Either keep the rc1 values, or gate the
looser bound on the inference having converged. Ship the rest of main (the two
new no-op flags, the size guard, the docs) which is unaffected.
