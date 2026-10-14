# Scope of the 100,000- and 1,000,000-row measurements

The archived 13 September 2026 development fixture has 25 predictors and three
active inputs. Each condition was measured once. The measured build is identified
by its archive hash in [SUMMARY_RECEIPT.json](SUMMARY_RECEIPT.json), not solely by
the reused development version label 0.5.4.9000. These are historical measurements,
not a new benchmark of the subsequent boundary and packaging repairs.

| Recorded operation | Queries / candidate workers | Recorded time | Whole-process peak RSS |
|---|---|---:|---:|
| 100k rows, ANN eps=0.5 screening | all / 8 | 103.709 s screen | 340.1 MiB |
| 100k rows, ANN eps=0.5 complete pipeline, twin_r=40 | 2,500 / 8 | 30.214 s process wall | 466.0 MiB |
| 1m rows, five-dimensional ANN index only, eps=0.5 | 16 / 1 | 0.171 s native index/query | 713.3 MiB |

The representative-query pipeline's screening and allocation stages took 29.571
and 0.265 seconds. It retained all 100,000 observations as possible neighbors;
no all-query PME reference was computed for that 100k condition.

The million-row check is not a full screening or PME run. Its 2.193-second
process wall time and 713.3-MiB RSS include generation of the 25-column input
fixture and direct-distance validation. They do not measure the index alone.
The executed query count is `neighbor_queries=16`; the generic `query_count`
field was initialized to n. The maximum kth-distance ratio, 1.0303, is below
the requested bound 1.5; mean neighbor recall was 0.9294.

At small retained dimension, K is 475 for 100k rows and 1,500 for one million
rows. All-query local fitting therefore visits about 47.5 million and 1.5 billion
neighbor rows per evaluated subset. Preprocessing copies, global ridge/SVD,
per-worker projected matrices/trees, and exact PME's exponential retained-size
cost remain. A fast complete million-row pipeline was not measured, and no
300,000-row/2,500-query receipt was found in the audited records.

The compact [BENCHMARK_ROWS.csv](BENCHMARK_ROWS.csv) and cohort hashes preserve
these scope distinctions. Full raw benchmark archives are not bundled.
