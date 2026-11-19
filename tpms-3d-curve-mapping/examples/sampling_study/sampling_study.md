# Sampling sensitivity study (gyroid, t = 0.40, cell 5 mm)

Each table changes one setting and keeps the others at their defaults (nSeeds = 48, nGrid = 64, stepScale = 1). Cells show the deviation from the finest setting of that factor, with the 95% CI half-width in brackets for curve metrics; **bold** = fails the rule for that factor. Tolerances: 3% relative, 0.02 absolute for 0-1 / -1..+1 metrics, 10 points for % metrics.

- **nSeeds** must make every CI narrow enough (precision) and agree with the finest run.
- **nGrid** also reshuffles the seeds, so differences within sampling noise pass.
- **stepScale** keeps the same seeds, so any difference is a real discretisation error.

**Recommended (cheapest setting that passes):** nSeeds = 48, nGrid = 32, stepScale = 2.

### Seeds needed per metric (95% CI half-width -> tolerance)

| Metric | CI half-width at 160 seeds | Tolerance | Seeds needed |
|---|---|---|---|
| Handedness index | 0.00933 | 0.02 | 35 |
| Max-bending lines (k1): straightness over 1 cell | 0.00716 | 0.02 | 21 |
| Min-bending lines (k2): straightness over 1 cell | 0.00139 | 0.02 | 1 |

Run time grows roughly linearly with seeds (see the run-time row below). Metrics needing many seeds are medians over individual curves; their CIs in gyroid_summary.md show the precision actually achieved, so a wide interval there is the signal to add seeds.

## nSeeds

| Metric | Unit | Finest value | 12 | 24 | 48 | 96 | 160 |
|---|---|---|---|---|---|---|---|
| Typical surface bending radius | mm | 1.132 | +0.0% | +0.0% | +0.0% | +0.0% | +0.0% |
| Downward-facing overhang (> 45 deg) | % | 11.59 | +0 | +0 | +0 | +0 | +0 |
| Curvature-change length | mm | 8.61 | +0.0% | +0.0% | +0.0% | +0.0% | +0.0% |
| Bending lines: sideways turning radius | mm | 4.444 | +0.0% | +0.0% | +0.0% | +0.0% | +0.0% |
| Max-bending lines (k1): straightness over 1 cell | 0-1 | 0.2624 | **+0.0161 (+-0.024)** | +0.00891 (+-0.0084) | +0.00224 (+-0.01) | -0.00868 (+-0.0087) | +0 |
| Min-bending lines (k2): straightness over 1 cell | 0-1 | 0.5525 | -0.00314 (+-0.0021) | -0.0012 (+-0.0023) | +0.000427 (+-0.0013) | -0.000114 (+-0.0015) | +0 |
| Handedness index | -1..+1 | 0.1158 | **-0.0153 (+-0.034)** | **-0.0103 (+-0.026)** | -0.00101 (+-0.017) | +0.00962 (+-0.013) | +0 |
| *run time (s)* | | | 5 | 7 | 8 | 13 | 18 |

## nGrid

| Metric | Unit | Finest value | 32 | 48 | 64 | 96 |
|---|---|---|---|---|---|---|
| Typical surface bending radius | mm | 1.132 | -0.1% | -0.1% | -0.0% | +0.0% |
| Downward-facing overhang (> 45 deg) | % | 11.6 | +0.117 | +0.0302 | -0.0132 | +0 |
| Curvature-change length | mm | 8.638 | +2.4% | +0.4% | -0.3% | +0.0% |
| Bending lines: sideways turning radius | mm | 4.459 | +2.1% | +0.3% | -0.3% | +0.0% |
| Max-bending lines (k1): straightness over 1 cell | 0-1 | 0.2609 | +0.0105 (+-0.019) | +0.00632 (+-0.016) | +0.00373 (+-0.01) | +0 |
| Min-bending lines (k2): straightness over 1 cell | 0-1 | 0.5526 | +0.000936 (+-0.0021) | +0.000747 (+-0.0015) | +0.000256 (+-0.0013) | +0 |
| Handedness index | -1..+1 | 0.1116 | -0.00974 (+-0.02) | -0.00883 (+-0.022) | +0.00309 (+-0.017) | +0 |
| *run time (s)* | | | 9 | 9 | 9 | 9 |

## stepScale

| Metric | Unit | Finest value | 2 | 1 | 0.5 |
|---|---|---|---|---|---|
| Typical surface bending radius | mm | 1.132 | +0.0% | +0.0% | +0.0% |
| Downward-facing overhang (> 45 deg) | % | 11.59 | +0 | +0 | +0 |
| Curvature-change length | mm | 8.61 | +0.0% | +0.0% | +0.0% |
| Bending lines: sideways turning radius | mm | 4.444 | +0.0% | +0.0% | +0.0% |
| Max-bending lines (k1): straightness over 1 cell | 0-1 | 0.2647 | -6.86e-05 (+-0.0095) | -5.56e-05 (+-0.01) | +0 |
| Min-bending lines (k2): straightness over 1 cell | 0-1 | 0.5529 | -3.19e-05 (+-0.0015) | -2.5e-05 (+-0.0013) | +0 |
| Handedness index | -1..+1 | 0.1146 | +0.00175 (+-0.017) | +9.79e-05 (+-0.017) | +0 |
| *run time (s)* | | | 6 | 9 | 14 |

