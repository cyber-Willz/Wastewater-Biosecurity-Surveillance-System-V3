# Changelog

## v0.4 — corrections + SHPINN (this archive)

Derived from the formal analysis in `docs/ww_biosec_theory.md`; measured effect in its Appendix C.

### Detection (`ww_detection`)
* **Fix:** the first observation of a `(site, analyte)` key was counted twice (`n = 2` after one sample). Warm-up is now exactly 7 prior observations.
* **Fix:** outbreaks inflated the variance they were scored against (self-masking). New `EwmaBaseline::ingest_robust` (winsorised EWMA, outlier-excluded Welford variance); `AnomalyDetector::observe_analyte` uses it with `clip_z = z_threshold`.
* **Fix:** `SewageNetwork` Rayleigh quotient normalised by `λ_max` (was `/2`, which only reached `[0, ½]`).
* **New:** `spectral_score_residual` (baseline-invariant network score on the mixing-scaled standardised residual), `residual_field`, `null_threshold` (deterministic empirical-null quantile), `laplacian()`, `sqrt_degrees()`.
* **New:** `AnomalyDetector::peek_z`, `with_spread_threshold` / `set_spread_threshold`, `DetectionSeverity::from_scores_with`.
* `spectral_score` (raw log₁₀ field) is kept for compatibility but documented as not baseline-invariant.
* 12 unit tests (exact Belize spectrum `{0, 1/21, 9/14, 1×5}`, score decomposition, warm-up count, robustness, severity table).

### Simulation (`ww_runner`)
* Network score computed from peeked z-scores (no raw-level input); upgrade threshold = null q95 of the score.
* Prints an SHPINN transport check per analyte category.

### New crate `ww_shpinn`
* Linear-feature physics-informed solver for `∂ₜc + v∂ₓc − D∂ₓₓc + kc = f` (exact least-squares, deterministic), data assimilation of measurements, analytic Gaussian reference, Theorem 7.1 constants (`kappa`, `stability_bound`).
* Hypergraph-regularised inverse source estimation (`prior_matrix`, `estimate_source`; Theorem 7.2 closed form).
* 7 unit tests, `examples/verify_forward.rs`.

### Docs
* README, doc comments (`lib.rs`, `analyte.rs`, `review.rs`, `log.rs`) corrected; `docs/ww_biosec_theory.md` added.

### Not changed
* The hypergraph is still built from catchments only; directed `site_flows_to` links are not used by the spectral model.
* `flow_liters` / `catchment_pop` are stored but do not enter detection.
* The audit log is an in-memory append-only vector (no hash chain / signatures).

### Measured effect on the shipped Belize scenario
| | v0.3 | v0.4 |
|---|---|---|
| alerts on pulse pairs / no-pulse pairs | 25 / 26 | 43 / 19 |
| injected pulses alerted (of 15) | 13 | 14 |
| spectral tier upgrades | 8 | 0 (calibrated threshold 0.966 never exceeded) |
