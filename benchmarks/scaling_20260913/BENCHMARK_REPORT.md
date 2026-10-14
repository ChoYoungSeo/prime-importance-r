# PRIME indexed-neighbor development benchmark, 2026-09-13

The exact indexed implementation improves screening time and memory on this
paired development fixture. At n=50,000, p=25, all-query screening completes in
46.979 seconds with eight candidate threads; ANN eps=0.5 reduces this to 36.280
seconds. With the existing opt-in `twin_r=20`, the single-thread full pipeline
uses 2,500 query rows and finishes in 14.690 seconds with exact neighbors.
The default remains all-query exact. Earlier single-thread all-query p=25
screens exceeded their 120-second limit and remain preserved.

## Provenance and scope

- Frozen baseline: installed `prime 0.5.3` from Hilbert v11, already using
  256-query blocks. The comparison is blocked exact versus indexed exact/ANN.
- Measured development build: `prime 0.5.4.9000`; its source archive and cohort
  hashes are retained in [SUMMARY_RECEIPT.json](SUMMARY_RECEIPT.json).
- R 4.6.1, a fresh R process per operation, sequential runs. The primary paired
  baseline comparison uses one thread. A separately labeled production-scale
  screen uses eight PRIME candidate threads with BLAS/OMP/MKL/VECLIB at one.
- The primary comparison uses every observation as a query. A separately labeled
  `twin_r=20` pilot exercises the existing opt-in query approximation.
- The fixture has three active inputs with independent Uniform(-1,1) predictors,
  response 2 sin(pi X1) + 1.5 X2 + X3^2, and Gaussian noise sd=0.2.
- One deterministic paired replication per cell. Timings are development
  measurements, not replicated runtime estimates or a confirmatory accuracy study.
- The statistical buffer and other method settings are unchanged across these
  exact/ANN comparisons. This does not decide the separate buffer ablation.

## Screening

Times below use `screen_sec`, which measures the public `prime_screen` call.
RSS is the operating-system maximum for the entire fresh process, including R,
package loading, and data preparation. Approximate variants use Euclidean eps.
This first table uses one candidate thread throughout.

| n | p | Blocked baseline, s | Indexed exact, s | ANN 0.5, s | Exact speedup | Baseline/exact RSS, MiB |
|---:|---:|---:|---:|---:|---:|---:|
| 2,000 | 5 | 0.334 | 0.216 | 0.188 | 1.55x | 130.3 / 93.7 |
| 2,000 | 25 | 2.058 | 1.059 | 0.919 | 1.94x | 120.6 / 102.0 |
| 10,000 | 5 | 5.735 | 2.323 | 1.956 | 2.47x | 184.1 / 109.8 |
| 10,000 | 25 | 39.387 | 12.063 | 10.464 | 3.26x | 206.6 / 142.8 |
| 50,000 | 5 | Not attempted | 27.791 | 22.856 | Unavailable | Unavailable / 141.3 |
| 50,000 | 25 | Not attempted | Timed out at 120 s | Not attempted | Unavailable | Unavailable |

At n=50,000, p=25, ANN eps=0.1 also timed out at 120 seconds. The remaining
ANN eps=0.5 screen was not dispatched after the 300-second dispatch budget.
The sampled peak RSS on the two timed-out screens was 201.2 and 202.6 MiB,
respectively, far below the 4-GiB memory guard. These are time-censored outcomes;
they do not yield completed algorithm times or speed ratios. No n=50,000 old
baseline was run.

The separately authorized eight-candidate-thread all-query runs both complete:

| n | p | Query rows | Candidate threads | Neighbor rule | Screen, s | Whole process, s | OS peak RSS, MiB |
|---:|---:|---:|---:|---|---:|---:|---:|
| 50,000 | 25 | 50,000 | 8 | Indexed exact | 46.979 | 47.275 | 227.1 |
| 50,000 | 25 | 50,000 | 8 | ANN eps=0.5 | 36.280 | 36.565 | 232.8 |

Both candidate and final selected sets equal {1,2,3}. ANN screening is 1.295x
faster than indexed exact under the same eight-thread setting. No eight-thread
old-baseline run was made, so no old-baseline eight-thread speedup is claimed.

## Fixed-active-set PME

These runs hold the retained coordinates at 1:3, so all variants estimate the
same eight-coalition game. They isolate downstream engine cost from selection.
The stage timer is `pme_sec`.

| n | p | Baseline, s | Indexed exact, s | ANN 0.5, s | Baseline/exact RSS, MiB |
|---:|---:|---:|---:|---:|---:|
| 10,000 | 5 | 0.980 | 0.978 | 0.884 | 151.7 / 95.2 |
| 10,000 | 25 | 0.971 | 0.975 | 0.860 | 173.7 / 112.0 |
| 50,000 | 5 | Not attempted | 12.111 | 10.162 | Unavailable / 124.2 |
| 50,000 | 25 | Not attempted | 12.155 | 10.231 | Unavailable / 162.0 |

The main measured speed benefit is in screening. At n=10,000, fixed-s PME exact
time is similar to the baseline, with a substantially smaller memory footprint.

## Separate opt-in representative-query pipeline

The existing `twin_r=20` option at n=50,000, p=25 selects exactly 2,500 query
locations from standardized (X,y); all 50,000 observations remain eligible
neighbors. These runs use one candidate thread and the same seed/input arrays
as the saved all-query results. Default query selection is unchanged.

| Neighbor rule | Actual query count | Screen including query preparation, s | Subsequent PME, s | Whole pipeline process, s | OS peak RSS, MiB |
|---|---:|---:|---:|---:|---:|
| Indexed exact | 2,500 | 14.215 | 0.193 | 14.690 | 246.9 |
| ANN eps=0.5 | 2,500 | 13.023 | 0.182 | 13.546 | 247.1 |

The pipeline reuses the package's existing risk cache after screening. Its
subsequent PME stage should not be directly compared to standalone fixed-s PME
timing, which starts without that screen cache.

Both runs retain exactly {1,2,3} and share the same query indices. Their half-L1
allocation differences from the saved all-query indexed-exact fixed-{1,2,3} PME
reference are 0.000091365 for query approximation alone and 0.000053864 for query
approximation plus ANN. The latter smaller discrepancy may reflect cancellation
between two approximation errors; it does not establish that ANN is more accurate.
The ANN allocation difference from the same-query exact pipeline is 0.000135529.

## Fidelity and completeness

- Corrected n=100 smoke: 11/11 operations completed and validated.
- n=2,000: 30/30 operations completed and validated, including a decomposed full
  pipeline for p=5 and p=25.
- n=10,000: 22/22 operations completed and validated.
- n=50,000: 18 operations planned, 17 attempted, 15 completed; two screen
  timeouts and one screen not attempted remain explicitly recorded.
- Separate n=50,000, p=25 eight-candidate-thread all-query screen: 2/2 completed
  and validated, with no changes to candidate or selected sets.
- Separate n=50,000, p=25 `twin_r=20` pipeline: 2/2 completed and validated,
  including actual query counts/indices and all-query reference comparison.
- Across completed exact-baseline comparisons, candidate and selected sets are
  identical; the observed PME allocation half-L1 difference is exactly zero.
- Across completed ANN/exact screens at n>=2,000, candidate and selected sets
  are identical on this fixture. All retain inputs 1,2,3.
- Maximum ANN allocation half-L1 difference: 0.00144177 at n=2,000,
  0.00040659 at n=10,000, and 0.00010614 at n=50,000.
- Every completed native neighbor fixture passed the distance-bound and returned
  squared-distance checks; exact neighbor recall was one. ANN recall/error is
  retained per query in `neighbor_recall.csv`.
- The initial smoke v1 failed during package loading because `Rcpp` was absent
  from the default library path. Those 11 failed attempts are preserved; v2 uses
  the build's existing dependency library and succeeds.

## Archived evidence

The full output directories are retained in the original development archive and
are not bundled in this repository. [BENCHMARK_ROWS.csv](BENCHMARK_ROWS.csv) contains
an unchanged projection of their recorded metrics; [SUMMARY_RECEIPT.json](SUMMARY_RECEIPT.json)
records all ten cohorts' planned/attempted/completed counts and original file hashes.
The compact record includes dependency-load failures, timeouts, and the later
100k/million-row scope described in [SCALABILITY_100K_1M.md](SCALABILITY_100K_1M.md).
No current runtime claim is inferred solely from a shared development version label.

`operation_sec` additionally includes diagnostic serialization/bookkeeping;
`wall_sec` includes process startup and input preparation. Neither replaces the
stage-specific timers above. ANN distance tolerances do not by themselves bound
allocation error; the reported allocation changes are measured on these inputs.
