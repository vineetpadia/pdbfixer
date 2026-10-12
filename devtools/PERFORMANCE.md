# Performance investigation

The changes remove avoidable Python work in sequence matching, hydrogen variant
preparation, and topology reconstruction. They preserve the repair force field,
minimization settings, and public API.

## Findings and changes

- `findMissingResidues()` incremented every remaining residue ID for each
  insertion code. A running offset now computes the same IDs in one pass.
  Alignment attempts also copied the entire remaining sequence; direct indexing
  of observed residues removes those allocations and skips wildcard gaps.
- `addMissingHydrogens()` called `residue.bonds()` for every applicable custom
  residue. OpenMM implements that iterator by scanning all topology bonds. One
  external-bond index now serves the whole call and is rebuilt after topology
  changes. Standard-only variant preparation does not need the index.
- `_addAtomsToTopology()` calculated every residue's center on each reconstruction
  even when no residues were missing. Those calculations now only run for missing
  residue placement, and existing centers are reused for loop endpoints.

## Measurements

Baseline: commit `5e658a4`. Environment: AMD Ryzen 7 3700X, Linux x86-64,
Python 3.13.11, OpenMM 8.4, NumPy 2.1.3, `OPENMM_CPU_THREADS=1`.
The table reports medians of seven samples after one warmup. Baseline and changed
versions ran sequentially. Cyclic garbage collection runs outside each timed
sample, preventing cleanup of a previous topology from distorting the next sample.

| Workload | Baseline | Changed | Ratio |
| --- | ---: | ---: | ---: |
| Unmatched sequence, 32,000 entries (synthetic) | 486.15 ms | 7.78 ms | 62.5x |
| Insertion-code chain, 8,000 residues (synthetic) | 827.01 ms | 3.07 ms | 269.3x |
| Hydrogen variant preparation, 2,000 custom residues (synthetic) | 248.10 ms | 3.80 ms | 65.3x |
| 4JSV topology reconstruction, missing-residue insertion disabled | 238.18 ms | 164.93 ms | 1.44x |
| 4JSV topology reconstruction, 250 missing residues inserted | 422.54 ms | 427.49 ms | 0.99x |
| `addMissingAtoms()` on pre-parsed, complete `test.pdb` | 3.41 ms | 1.53 ms | 2.22x |

The large synthetic gains are isolated stage measurements, not whole-repair
speedups. The 4JSV fixture contains 22,194 input atoms. Reconstruction timings
exclude file parsing and energy minimization. The complete-structure case measures
warm calls on a pre-parsed fixer. The loop-placement case is essentially unchanged.

An initial CPU profile of `addMissingAtoms(seed=1)` on 1BHL attributed 10.10 of
10.52 seconds (96%) to OpenMM dynamics and energy minimization, versus 0.14 seconds
to the two topology reconstructions. The changes therefore do not promise a large
speedup for difficult loop rebuilding. They improve Python preparation and its
scaling while retaining the numerical repair procedure.

Three fresh full-pipeline CPU samples (including parsing and detection) had medians
of 9.88 seconds before and 7.86 seconds after. The changed samples ranged from 6.80
to 10.65 seconds, so these variable simulation timings do not establish a
consistent whole-repair speedup. The stage measurements above are the evidence for
the performance improvements.

## Reproduce

Install the package and its dependencies in an environment, then run from the
changed repository root:

```bash
git worktree add --detach /tmp/pdbfixer-baseline 5e658a4
OPENMM_CPU_THREADS=1 PYTHONPATH=/tmp/pdbfixer-baseline python devtools/benchmark_performance.py --repeat 7 > /tmp/pdbfixer-before.json
OPENMM_CPU_THREADS=1 PYTHONPATH=. python devtools/benchmark_performance.py --repeat 7 > /tmp/pdbfixer-after.json
```

The script uses only bundled fixtures and synthetic inputs. Its JSON records the
resolved imported source path, dependency versions, and every timing sample.
The hydrogen benchmark supplies definitions locally and replaces OpenMM's
hydrogen-addition call with an output check to isolate variant preparation.
Add `--include-repair` to measure a fresh load, missing-residue/atom detection, and
seeded CPU repair of 1BHL, including minimization. This takes longer.

## Validation

The full suite passes: **54 tests**, including 30 new regression cases for sequence
matching and hydrogen variants. The suite includes remote downloads, nonstandard
residues, solvation, CLI calls, and reconstruction geometry.

Additional baseline comparisons verified:

- Identical missing-residue results for 2,000 randomized cases and all five bundled
  PDB fixtures.
- Identical atom/bond topology and coordinates for 12 reconstruction combinations:
  `test.pdb`, 1BHL, and 4JSV; with/without missing residues; both reconstruction modes.
- Internal and terminal gaps with Vec3/NumPy positions in nanometers/angstroms:
  nanometer results are identical; the largest angstrom conversion difference is
  `4.8e-14 nm`.

```bash
OPENMM_CPU_THREADS=1 python -m pytest -q pdbfixer --durations=10
```
