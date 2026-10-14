# Indexed-neighbor development benchmarks

This directory contains benchmark harness source and a compact record of the
13 September 2026 development measurements. It is not a confirmatory accuracy
study. See [BENCHMARK_REPORT.md](BENCHMARK_REPORT.md) and
[SCALABILITY_100K_1M.md](SCALABILITY_100K_1M.md) for the measured scope.

[BENCHMARK_ROWS.csv](BENCHMARK_ROWS.csv) preserves selected fields from all 98
attempted-run records, including failures and timeouts. [SUMMARY_RECEIPT.json](SUMMARY_RECEIPT.json)
records the archived source hash, cohort counts, and original manifest/raw-file
hashes. Full raw outputs, installed libraries, source archives, workstation
activation receipts, and old harness snapshots are not bundled.

## Reproducing a comparison

The frozen baseline is the blocked exact implementation, version 0.5.3; it
already used 256-query blocks. Supply that library and its build directory
explicitly. The build directory must contain `prime_0.5.3.tar.gz` and the
recorded `build_manifest.json`. The runner checks their fingerprints before
planning or executing. It never modifies the baseline library.

```sh
python3 benchmarks/scaling_20260913/run.py \
  --baseline-lib /path/to/frozen/baseline/library \
  --baseline-build /path/to/frozen/baseline/build \
  --new-lib /path/to/new/isolated/library \
  --output /path/to/new/plan-directory
```

Without `--execute`, this writes a plan only. A new output directory is always
required. Add `--n 100 --p 5 --modes screen,pme,neighbors --timeout 30
--memory-gib 1 --execute` for a bounded plumbing run. Dependencies must be
available to the selected R executable; use `--rscript` and `R_LIBS_USER` when
needed. `--new-build-manifest` can pin the new build's archive and installed files.

```sh
python3 benchmarks/scaling_20260913/validate.py /path/to/output
```

The harness supports macOS and Linux process-resource measurement. Package CI
separately checks Linux, macOS, and Windows; the benchmark harness does not claim
Windows process-monitor portability. `--help` lists all budgets and controls.

## Protocol and interpretation

- Independent uniform predictors; response `2*sin(pi*X1)+1.5*X2+X3^2` plus
  Gaussian noise with standard deviation 0.2. Coordinate-specific seeds keep
  active inputs fixed when nuisance dimensions are added.
- `screen` times `prime_screen`; `pme` fixes the active set to 1:3; `pipeline`
  decomposes screening and allocation with a harness-only selected-dimension
  budget (default 8). Its historical empty-screen sentinel is unavailable as a
  normalized allocation, rather than a valid zero allocation.
- Default package queries and neighbor search are all-query exact. Positive
  `--eps` changes neighbor-distance approximation; `--twin-r` separately uses
  representative queries while retaining all observations as neighbor candidates.
- `--threads` controls candidate workers; BLAS/OMP/MKL/Accelerate are set to one.
  Stage times, process wall time, and whole-process peak RSS are distinct fields.
- One repetition per archived condition is a development measurement, not a
  general runtime guarantee. Timeouts remain censored; undispatched work is
  counted in the cohort receipt and has no invented run row.
- ANN distance guarantees do not bound selection or allocation error. The
  one-million-row result evaluates a five-dimensional index at 16 queries,
  not full screening or PME. Use `neighbor_queries` for that executed scope.
- Exact PME remains exponential in retained dimension. Neighbor streaming does
  not bound the full projection's memory/runtime requirements.

Optional edge checks accept an explicit installed library:

```sh
Rscript --vanilla benchmarks/scaling_20260913/edge_validation/collinear_lambda_zero.R /path/to/library
python3 benchmarks/scaling_20260913/edge_validation/check_interrupt.py --library /path/to/library
```

These create their own bounded diagnostics. Generated outputs are ignored by Git;
they are not part of the archived measurements summarized here.
