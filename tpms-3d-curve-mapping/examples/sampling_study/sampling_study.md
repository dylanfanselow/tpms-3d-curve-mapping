# Sampling sensitivity study (gyroid, t = 0.40, cell 5 mm)

Each table changes one setting and keeps the others at their defaults (nSeeds = 48, nGrid = 64, stepScale = 1). Cells show the deviation from the finest setting of that factor, with the 95% CI half-width in brackets for curve metrics; **bold** = fails the rule for that factor. Tolerances: 3% relative, 0.02 absolute for 0-1 / -1..+1 metrics, 10 points for % metrics.

- **nSeeds** must make every CI narrow enough (precision) and agree with the finest run.
- **nGrid** also reshuffles the seeds, so differences within sampling noise pass.
- **stepScale** keeps the same seeds, so any difference is a real discretisation error.

**Recommended (cheapest setting that passes):** nSeeds = 160 (precision target NOT reached at 160; about 438 seeds needed, driven by "Zero-bending lines: direction memory"; see the table below), nGrid = 32, stepScale = 2.

### Seeds needed per metric (95% CI half-width -> tolerance)

| Metric | CI half-width at 160 seeds | Tolerance | Seeds needed |
|---|---|---|---|
| Zero-bending lines: direction memory | 4.96% | 3% | 438 |
| Bending lines: twist pitch | 4.26% | 3% | 323 |
| Min-bending lines (k2): share forming closed loops | 7.81 pts | 10 pts | 98 |
| Bending lines: sideways turning radius | 1.78% | 3% | 57 |
| Handedness index | 0.00979 | 0.02 | 39 |
| Max-bending lines (k1): straightness over 1 cell | 0.00856 | 0.02 | 30 |
| Zero-bending lines: straightness over 1 cell | 0.00331 | 0.02 | 5 |
| Max-bending lines (k1): direction memory | 0.108% | 3% | 1 |
| Min-bending lines (k2): straightness over 1 cell | 0.00144 | 0.02 | 1 |
| Min-bending lines (k2): direction memory | 0.026% | 3% | 1 |

Run time grows roughly linearly with seeds (see the run-time row below). Metrics needing many seeds are medians over individual curves; their CIs in gyroid_summary.md show the precision actually achieved, so a wide interval there is the signal to add seeds.

## nSeeds

| Metric | Unit | Finest value | 12 | 24 | 48 | 96 | 160 |
|---|---|---|---|---|---|---|---|
| Typical surface bending radius | mm | 1.132 | +0.0% | +0.0% | +0.0% | +0.0% | +0.0% |
| Downward-facing overhang (> 45 deg) | % | 11.59 | +0 | +0 | +0 | +0 | +0 |
| Max-bending lines (k1): straightness over 1 cell | 0-1 | 0.2624 | **+0.0161 (+-0.024)** | +0.00891 (+-0.0094) | +0.00224 (+-0.011) | -0.00868 (+-0.0084) | +0 |
| Max-bending lines (k1): direction memory | mm | 1.126 | +0.4% (+-0.4%) | +0.3% (+-0.3%) | +0.0% (+-0.2%) | -0.1% (+-0.2%) | +0.0% |
| Max-bending lines (k1): share forming closed loops | % | 0 | +0 | +0 | +0 | +0 | +0 |
| Min-bending lines (k2): straightness over 1 cell | 0-1 | 0.5525 | -0.00314 (+-0.0021) | -0.0012 (+-0.0023) | +0.000427 (+-0.0013) | -0.000114 (+-0.0015) | +0 |
| Min-bending lines (k2): direction memory | mm | 1.582 | -0.0% (+-0.1%) | +0.0% (+-0.1%) | +0.0% (+-0.0%) | -0.0% (+-0.0%) | +0.0% |
| Min-bending lines (k2): share forming closed loops | % | 65.62 | **-40.6 (+-21)** | **-24 (+-19)** | **+1.04 (+-11)** | +2.08 (+-8.9) | +0 |
| Bending lines: sideways turning radius | mm | 5.423 | **+4.5% (+-4.1%)** | **+2.8% (+-4.0%)** | **-3.6% (+-3.4%)** | **-2.2% (+-3.1%)** | +0.0% |
| Bending lines: twist pitch | mm | 8.183 | **-2.6% (+-13.9%)** | **+6.4% (+-9.8%)** | **+2.1% (+-4.9%)** | **-2.5% (+-3.2%)** | **+0.0%** |
| Handedness index | -1..+1 | 0.1158 | **-0.0153 (+-0.036)** | **-0.0103 (+-0.026)** | -0.00101 (+-0.018) | +0.00962 (+-0.012) | +0 |
| Zero-bending lines: straightness over 1 cell | 0-1 | 0.933 | -0.00138 (+-0.0095) | -0.000867 (+-0.0067) | -0.000275 (+-0.0046) | -0.000275 (+-0.0028) | +0 |
| Zero-bending lines: direction memory | mm | 7.579 | **-6.1% (+-17.2%)** | **-1.3% (+-12.0%)** | **-0.4% (+-9.2%)** | **-0.4% (+-6.4%)** | **+0.0%** |
| Zero-bending families imbalance | -1..+1 | 0.1811 | +0 | +0 | +0 | +0 | +0 |
| Max-bending direction near vertical (< 30 deg) | % | 10.78 | +0 | +0 | +0 | +0 | +0 |
| *run time (s)* | | | 8 | 9 | 12 | 18 | 24 |

## nGrid

| Metric | Unit | Finest value | 32 | 48 | 64 | 96 |
|---|---|---|---|---|---|---|
| Typical surface bending radius | mm | 1.132 | -0.1% | -0.1% | -0.0% | +0.0% |
| Downward-facing overhang (> 45 deg) | % | 11.6 | +0.117 | +0.0302 | -0.0132 | +0 |
| Max-bending lines (k1): straightness over 1 cell | 0-1 | 0.2609 | +0.0105 (+-0.021) | +0.00632 (+-0.016) | +0.00373 (+-0.01) | +0 |
| Max-bending lines (k1): direction memory | mm | 1.126 | +0.1% (+-0.2%) | +0.1% (+-0.2%) | -0.0% (+-0.2%) | +0.0% |
| Max-bending lines (k1): share forming closed loops | % | 0 | +0 | +0 | +0 | +0 |
| Min-bending lines (k2): straightness over 1 cell | 0-1 | 0.5526 | +0.000936 (+-0.0019) | +0.000747 (+-0.0014) | +0.000256 (+-0.0011) | +0 |
| Min-bending lines (k2): direction memory | mm | 1.581 | +0.1% (+-0.0%) | +0.1% (+-0.0%) | +0.1% (+-0.0%) | +0.0% |
| Min-bending lines (k2): share forming closed loops | % | 64.58 | +2.08 (+-12) | +8.33 (+-11) | +2.08 (+-14) | +0 |
| Bending lines: sideways turning radius | mm | 5.415 | -0.4% (+-2.3%) | +0.0% (+-3.5%) | -3.5% (+-3.7%) | +0.0% |
| Bending lines: twist pitch | mm | 8.263 | +3.3% (+-7.3%) | +3.3% (+-5.9%) | +1.1% (+-4.3%) | +0.0% |
| Handedness index | -1..+1 | 0.1116 | -0.00974 (+-0.018) | -0.00883 (+-0.021) | +0.00309 (+-0.016) | +0 |
| Zero-bending lines: straightness over 1 cell | 0-1 | 0.9331 | +0.00148 (+-0.0044) | +0.00211 (+-0.0045) | -0.000389 (+-0.0045) | +0 |
| Zero-bending lines: direction memory | mm | 7.378 | +1.7% (+-9.9%) | +3.8% (+-9.6%) | +2.4% (+-8.7%) | +0.0% |
| Zero-bending families imbalance | -1..+1 | 0.1782 | -0.0072 | -0.00296 | +0.00293 | +0 |
| Max-bending direction near vertical (< 30 deg) | % | 10.83 | +0.21 | +0.0741 | -0.0515 | +0 |
| *run time (s)* | | | 12 | 12 | 12 | 13 |

## stepScale

| Metric | Unit | Finest value | 2 | 1 | 0.5 |
|---|---|---|---|---|---|
| Typical surface bending radius | mm | 1.132 | +0.0% | +0.0% | +0.0% |
| Downward-facing overhang (> 45 deg) | % | 11.59 | +0 | +0 | +0 |
| Max-bending lines (k1): straightness over 1 cell | 0-1 | 0.2647 | -6.86e-05 (+-0.01) | -5.56e-05 (+-0.01) | +0 |
| Max-bending lines (k1): direction memory | mm | 1.126 | +0.0% (+-0.2%) | +0.0% (+-0.2%) | +0.0% |
| Max-bending lines (k1): share forming closed loops | % | 0 | +0 | +0 | +0 |
| Min-bending lines (k2): straightness over 1 cell | 0-1 | 0.5529 | -3.19e-05 (+-0.0015) | -2.5e-05 (+-0.0012) | +0 |
| Min-bending lines (k2): direction memory | mm | 1.58 | +0.1% (+-0.1%) | +0.1% (+-0.0%) | +0.0% |
| Min-bending lines (k2): share forming closed loops | % | 66.67 | +0 (+-12) | +0 (+-12) | +0 |
| Bending lines: sideways turning radius | mm | 5.244 | -0.6% (+-4.9%) | -0.3% (+-3.2%) | +0.0% |
| Bending lines: twist pitch | mm | 8.367 | -0.5% (+-4.5%) | -0.2% (+-4.8%) | +0.0% |
| Handedness index | -1..+1 | 0.1146 | +0.00175 (+-0.019) | +9.79e-05 (+-0.018) | +0 |
| Zero-bending lines: straightness over 1 cell | 0-1 | 0.9327 | +0.00012 (+-0.0047) | +6.62e-05 (+-0.0046) | +0 |
| Zero-bending lines: direction memory | mm | 7.552 | +0.0% (+-8.6%) | -0.0% (+-9.1%) | +0.0% |
| Zero-bending families imbalance | -1..+1 | 0.1811 | +0 | +0 | +0 |
| Max-bending direction near vertical (< 30 deg) | % | 10.78 | +0 | +0 | +0 |
| *run time (s)* | | | 9 | 13 | 19 |

