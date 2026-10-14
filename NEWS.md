# prime 0.5.5

- Update the reference manual and package citation for Cho and Joseph (2026).
- Clarify nearest-neighbor distance ties and plotting examples in the manual.

# prime 0.5.4.9001

- Compute ridge ranking and BIC on a centered unit-variance response so that
  numerical guards do not change screening when response units change. Ridge
  path RSS diagnostics remain in the original response units; BIC diagnostics
  omit the common response-scale constant.
- Compute paired comparison scales in response-variance units to avoid squaring
  extremely large or small raw risk differences; diagnostics retain their units.
- Normalize precomputed games by their positive full-set value and use
  log-space ratio potentials, returning allocations in the supplied game units.
- Apply the relative null tolerance only to interior subset values, preserve
  both endpoints, and remove positive-value snapping from projection. Setting
  `zero_tolerance = 0` now retains every representable positive projected value.
- Validate the null tolerance and reject unrepresentable normalized game ratios
  clearly. Zero-full games project to zero and return zero allocations.

# prime 0.5.4.9000

- Reject exact projection edge counts beyond 32-bit index capacity before
  subset evaluation or projection allocation, without introducing a RAM budget.
- Support a single predictor with the existing MASS ridge normalization and
  preserve screening under response-colliding or non-syntactic predictor names.
- Correct singleton-cap and stopping descriptions, and document the existing
  empty-screen allocation sentinel.
- Replace production forward all-row distance sweeps with subset-specific
  indexed search, fused with querywise local fits. No query-by-n distance or
  query-by-K neighbor array is materialized in the production risk path.
- Add explicit `neighbor_eps` (zero exact; positive Euclidean distance
  approximation) and `query_block_size` controls. All queries and the entire
  neighbor pool remain in use by default, with the existing K rule.
- Separate scalar risk caching from a selection cache retaining at most two
  query-contribution vectors; allocation no longer retains a vector for every subset.
- Native exact equal-distance ties use original row indices. This can differ from
  historical FNN tie ordering on duplicate or quantized data. Continuous-data
  and tied-data accuracy are verified separately.
- Vendor the unmodified BSD-2-Clause nanoflann 1.9.0 header and license.

# prime 0.5.3

* Replaced the all-query quadratic forward-distance cache with shared blocks
  of at most 256 queries. All queries, candidate neighborhoods, and statistical
  calibration rules remain unchanged; the distance cache now uses at most
  256 * n doubles rather than query_count * n doubles.
* Preserved per-candidate querywise SSE addition order and boundary-tie flags
  across blocks, including multithreaded candidate evaluation.

# prime 0.5.2

* Replaced the fixed ridge threshold with BIC selection along the ridge-ranked
  prefix path.
* Added paired querywise calibration to forward gains and sequential backward
  deletion losses.
* Preserved raw risk ordering when all forward candidate values lie below the
  zero-signal boundary; calibration still determines whether the leader enters.
* Constructed the final proportional allocation from noise-adjusted subset
  risks before monotone projection.
* Removed the obsolete `alpha` tuning argument.

# prime 0.4.0

* Initial release of the PRIME screening and allocation workflow.
* Added local-linear subset-value estimation, inclusive forward screening,
  dimension-adjusted backward reduction, and Data Twinning support.
* Added structured print, summary, and plot methods for fitted and screening
  objects.
