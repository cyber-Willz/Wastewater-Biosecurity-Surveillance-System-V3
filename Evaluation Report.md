# ww_biosec — Evaluation Report

Two evaluations, both run against the **exact shipped detection engine**
(`ww_detection`'s `EwmaBaseline`/`AnomalyDetector` + `SewageNetwork` spectral
hypergraph scoring), not a re-implementation:

1. **Real public wastewater data** — Public Health Scotland's national
   SARS-CoV-2 wastewater surveillance programme (May 2020–Feb 2022, 122 sites,
   N1 gene RT-qPCR; published with a *Scientific Data* paper,
   `BioRDM/COVID-Wastewater-Scotland` on GitHub).
2. **Monte Carlo** — 50,000 randomized outbreak scenarios (+ a 20,000-run
   independent-seed replicate) on the 8-site Belize network topology from
   `ww_runner::scenario`.

Code: new `ww_eval` crate added to the workspace (`ww_eval/src/lib.rs`,
`src/bin/real_data.rs`, `src/bin/monte_carlo.rs`). It calls
`ww_detection::{AnomalyDetector, SewageNetwork}` directly — same
`observe_analyte`/`peek_z`/`spectral_score_residual` calls as
`ww_runner::scenario::run_belize`.

**A note on CDC NWSS:** I attempted this first. The full-history CSV
(`data.cdc.gov` resource `j9g8-acpt`) is >30 MB raw and CDC's weekly
Wastewater Viral Activity Level dataset (`atcp-73re`) serializes to >50 MB —
both exceeded what I could fetch in this environment (`data.cdc.gov` isn't on
the sandbox's network allowlist, and the whole-dataset downloads exceed the
fetch tool's size limits; a scoped Socrata query would need a URL I can't
construct here). Scotland's dataset is public, national-scale, real, and
similarly structured, so I used it as the primary real-data test instead.

---

## 1. Real-data evaluation (Scotland)

**Setup:** top 12 sites by sample count, weekly-aggregated (mean of samples
that week), `log10(max(gc/L, 10))` (10 gc/L floor as an LOD proxy). Sites are
grouped into catchments by NHS Health Board — an administrative proxy for
shared sewer topology, since the dataset doesn't publish real sewer
connectivity. Two real catchment pairs emerge (Fife: Dunfermline+Levenmouth;
Greater Glasgow & Clyde: Dalmuir+Shieldhall); everything else is
background-only for the spectral score. Fed through `AnomalyDetector` with
the shipped SARS-CoV-2 profile (α=0.39, z_threshold=2.5), 89 weekly rounds.

**Result: only 5 AMBER+ alerts across 89 weeks × 12 sites** (~1,050
site-weeks). 3 fall inside a hand-labeled, widely-reported Alpha/Delta/Omicron
wave window; 2 don't (plausibly local outbreaks the national wave label is
too coarse for — wastewater is inherently local). **The Omicron wave —
the largest surge in the whole period — produced zero alerts anywhere.**
I pulled the full weekly z-score trace (not just alerting weeks) to see why:
at every one of the 12 sites, the peak z-score during the Omicron window
(Dec 2021–Jan 2022) stayed under 1.0, nowhere near the 2.5 threshold.

The Shieldhall (Glasgow) trace above shows the pattern clearly: a sharp,
detected spike in Aug 2020 (z=2.9), then the EWMA baseline chases every
subsequent rise — including Omicron in Dec 2021 — closely enough that z never
crosses threshold.

**Root cause — first-pass hypothesis, tested and refined below:**
`decay_rate_k` and the derived α=`clamp(1−e^{−k}, 0.10, 0.40)` are calibrated
assuming **daily** sampling — the doc comment says "in one observation period
(1 day) a fraction `1−e^{−k}` decays." Real-world sampling here is
**weekly**, which I suspected was the dominant cause. Testing it directly
(next section) confirms it's a real, secondary effect — but not the
dominant one. The dominant one turns out to be unbounded variance
accumulation across repeated waves; see below.

Full per-alert and per-week z-score data: `real_data_alerts_scotland.csv`,
`real_data_full_trace_scotland.csv`.

### Follow-up: is it really α, or something else?

I tested the α-cadence hypothesis directly rather than leave it as a guess:
reran `real_data` with α swept from 0.05 to 0.39, everything else identical.

| α | total alerts (89 weeks × 12 sites) | Omicron wave detected? |
|---|---|---|
| 0.05 | 33 | **No** |
| 0.10 | 12 | No |
| 0.15 | 8 | No |
| 0.20 | 7 | No |
| 0.25 | 7 | No |
| 0.30 | 7 | No |
| 0.39 (shipped) | 5 | No |

Lower α does increase overall sensitivity monotonically (more alerts,
consistent with the cadence argument), and roughly doubles the peak z-score
reached during the Omicron window at every site (e.g. Levenmouth 0.62→1.80,
Seafield 0.77→1.61) — so the α-cadence effect is real and non-trivial. But
**even at α=0.05, no site ever reaches z=2.5 during Omicron.** α alone
doesn't explain it. I pulled the implied standard deviation over time
(`(obs − ewma) / z`, backed out from the trace) for Shieldhall to check for a
second effect:

| week | conc | ewma | z | implied σ |
|---|---|---|---|---|
| 2020-07-20 | 2.56 | 2.79 | -0.93 | **0.26** |
| 2020-09-07 | 4.07 | 3.24 | 2.17 | 0.39 |
| 2020-11-16 | 4.26 | 4.65 | -0.48 | 0.82 |
| 2021-06-28 | 5.14 | 4.72 | 0.60 | 0.71 |
| 2021-12-20 | 5.58 | 5.03 | 0.73 | **0.75** |

**The standard deviation roughly triples over the series' first four months
(0.26→0.82) and then plateaus around 0.7–0.8 for the rest of the 18-month
run.** That's the dominant effect, not α: `EwmaBaseline`'s Welford variance
accumulates over the *entire* history by design (only single-observation
outliers beyond `1.5·clip_z·σ` are excluded; a whole *wave's* worth of
moderately-elevated weeks, each individually under that exclusion bar, still
feeds the running variance). Once the Aug 2020 and Alpha (Jan 2021) waves are
baked into "normal" spread, the z-score denominator has tripled — so a wave
of *the same absolute size* as the Aug 2020 spike (which scored z=2.9) would
now only score roughly z≈1. Omicron's absolute rise is comparable to
earlier waves', but by Dec 2021 it's competing against a σ three times
larger than what caught the first spike. **This is the real root cause**,
and it's structural: any purely-cumulative variance estimator gets
progressively less sensitive to repeated events of the same relative size,
which is a poor fit for a system meant to catch every wave, not just the
first. α-cadence is a real, compounding secondary effect, not the main one —
I was too quick to name it as "the" cause in the first pass.

**Better fix:** bound the variance estimator's memory — an exponentially-
weighted (EWMA) variance instead of unbounded Welford, or a rolling window
(e.g. trailing 52 weeks), so a σ from 18 months ago stops permanently
depressing this month's z-scores. This is a more targeted and higher-value
fix than the α-cadence one I originally proposed; do both, but this one
matters more.

---

## 2. Monte Carlo (50,000 + 20,000 independent-seed runs)

**Setup:** the real 8-site Belize network + catchments from `scenario.rs`,
rebuilt via `SewageNetwork::build` (same Laplacian, `null_threshold` spread
threshold computed once and reused). Per run: one random outbreak (site,
peak day ∈ [9,16), amplitude ∈ [0.8, 3.5) log₁₀ units, width σ ∈ [1,3) days),
with 50% probability of a correlated "spreading" companion pulse (attenuated
0.4–0.7×, lagged 1–3 days) at another site in the same catchment. 22
simulated days (7-day warm-up) through the identical statistical + spectral
pipeline `run_belize` uses. Deterministic splitmix64 PRNG, independently
seeded per run — reproducible, and a second batch with a different base seed
reproduces the same statistics (sensitivity 0.9317 vs 0.9294; FA rate 0.0240
vs 0.0240), so results aren't an artifact of one seed.

| Metric | 50k-run result |
|---|---|
| Overall sensitivity (pulse detected AMBER+ while present) | **93.1%** |
| Median detection lead time vs. nominal peak day | **−2.4 days** (i.e. flagged ~2.4 days *before* the peak, on the rising edge) |
| Null-condition false-alarm rate (non-pulse site-days) | **2.40%** per site-day |
| Spectral tier-upgrade rate | **7.4%** of runs |

**Sensitivity scales with effect size as expected** (good — the detector
isn't just noise-triggered):

| Amplitude (log₁₀ units) | n | Detection rate | Median lead |
|---|---|---|---|
| 0.8–1.5 (weak) | 12,704 | 81.4% | −1.8 d |
| 1.5–2.2 (moderate) | 12,960 | 94.7% | −2.3 d |
| 2.2–2.9 (strong) | 13,100 | 98.0% | −2.6 d |
| 2.9–3.5 (severe) | 11,236 | 99.0% | −2.7 d |

**Wide pulses are detected less reliably than narrow ones** (87.6% at
σ=2.3–3.0d vs 97.1% at σ=1.0–1.7d for the same amplitude range) — a slow
ramp is harder for a fast-adapting EWMA to see as anomalous, echoing the
real-data finding above at a smaller scale.

**A genuine finding about the spectral "network-spread" upgrade:** its own
name suggests it should fire *more* on multi-site "spreading" scenarios. It
does the opposite — **9.4% upgrade rate on isolated single-site pulses vs.
1.9% on spreading ones.** This is actually consistent with the score's true
definition once you read the spectral.rs doc comment carefully: it's a
**high-frequency / roughness** statistic — high when concentrations diverge
sharply *between adjacent sites* (a still-localized anomaly), low when a
catchment moves together (which a correlated "spreading" pulse across
sites in the same catchment does, almost by construction). So the
score is working as mathematically defined, but "network-spread upgrade"
undersells what it detects — it's closer to a **"sharp, spatially isolated
anomaly" upgrade**, which is close to the opposite of what "network-spread"
suggests to a reader. Also worth noting: this evaluation finds it firing
7.4% of the time under broad randomized conditions, vs. 0% on the one fixed
scenario the CHANGELOG reports for v0.4 — the calibrated q95 threshold isn't
*never* reachable, just rarely, and the fixed scenario happened to land
below it every time.

Full per-run data (site, scenario type, amplitude, σ, detected, lead time,
severity, spectral-upgrade flag) for all 50,000 runs:
`monte_carlo_runs_50000.csv`.

---
