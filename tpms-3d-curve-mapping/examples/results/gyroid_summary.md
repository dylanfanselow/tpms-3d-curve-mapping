# Gyroid curve-mapping summary

Lattice: **sheet** gyroid, cell size **5 mm**. Build direction [0  0  1]. Generated 2026-11-16 05:48.

**t** is the level-set value of the gyroid formula: raising it thickens the sheet walls (or, for a network lattice, adds material to one side). Each column below is one t value.

Numbers in brackets are 95% confidence intervals from resampling the traced curves (narrow = trustworthy), except for wall thickness and bending radius, where they give the 5-95% spread over the surface. "> x" means longer than the longest distance checked.

## What t does to the part

| Metric | Unit | t = 0.00 | t = 0.20 | t = 0.40 | t = 0.60 | t = 0.80 | Trend |
|---|---|---|---|---|---|---|---|
| Wall thickness (sheet) | mm | 0 | 0.209 [0.188, 0.224] | 0.422 [0.38, 0.452] | 0.645 [0.579, 0.691] | 0.886 [0.793, 0.948] | rises from 0 to 0.886 |
| Solid volume fraction | % | 0 | 12.9 | 25.8 | 39.0 | 52.2 | rises from 0 to 52.2 |
| Surface area per bounding volume | 1/mm | 1.24 | 1.23 | 1.2 | 1.15 | 1.07 | falls 13% |
| Surface-to-solid-volume ratio (S/V) | 1/mm | n/a | 9.54 | 4.66 | 2.95 | 2.06 | falls 78% |
| Typical surface bending radius | mm | 1.18 [0.949, 3.6] | 1.17 [0.939, 3.54] | 1.13 [0.908, 3.31] | 1.07 [0.855, 2.99] | 0.979 [0.772, 2.57] | falls 17% |
| Downward-facing overhang (> 45 deg) | % | 11.6 | 11.6 | 11.6 | 11.6 | 11.7 | steady |

**What the metrics mean**

- **Wall thickness (sheet)**: Mean wall thickness; the range covers 90% of the wall area.
- **Solid volume fraction**: Share of the bounding volume filled with material.
- **Surface area per bounding volume**: Surface area packed into each unit of part volume (higher = more surface to finish, coat or exchange through).
- **Surface-to-solid-volume ratio (S/V)**: Surface area per unit of solid material (same quantity as the STL tool's S/V).
- **Typical surface bending radius**: Radius of the sphere that bends like the surface does (1/curvedness). Point-based, same information as the STL tool.
- **Downward-facing overhang (> 45 deg)**: Share of surface facing downward more steeply than the overhang limit (typically the roughest in AM).

## How the surface curves behave

| Metric | Unit | t = 0.00 | t = 0.20 | t = 0.40 | t = 0.60 | t = 0.80 | Trend |
|---|---|---|---|---|---|---|---|
| Curvature-change length | mm | 11.4 | 10.5 | 8.61 | 6.82 | 5.32 | falls 53% |
| Bending lines: sideways turning radius | mm | 4.63 | 4.56 | 4.44 | 4.24 | 3.94 | falls 15% |
| Max-bending lines (k1): straightness over 1 cell | 0-1 | 0.441 [0.434, 0.449] | 0.358 [0.352, 0.373] | 0.265 [0.252, 0.272] | 0.141 [0.128, 0.163] | 0.0426 [0.0369, 0.0506] | falls 90% |
| Min-bending lines (k2): straightness over 1 cell | 0-1 | 0.444 [0.435, 0.449] | 0.505 [0.499, 0.507] | 0.553 [0.552, 0.554] | 0.6 [0.599, 0.601] | 0.641 [0.641, 0.642] | rises 45% |
| Handedness index | -1..+1 | 0.001 [-0.0207, 0.0249] | 0.067 [0.0487, 0.086] | 0.115 [0.0961, 0.13] | 0.147 [0.136, 0.159] | 0.17 [0.165, 0.174] | rises from 0.001 to 0.17 |

**What the metrics mean**

- **Curvature-change length**: How far you move along the surface before its bending changes by about its own typical amount (typical curvature / typical mean-curvature gradient). Short = rapidly varying curvature.
- **Bending lines: sideways turning radius**: Typical radius of the sideways turns the lines of curvature make within the surface (small = sharp steering). Evaluated over the whole surface, so no sampling noise.
- **Max-bending lines (k1): straightness over 1 cell**: 1 = runs straight, 0 = comes back to where it started. Net progress per distance travelled.
- **Min-bending lines (k2): straightness over 1 cell**: 1 = runs straight, 0 = comes back to where it started. Net progress per distance travelled.
- **Handedness index**: Share of line length twisting right minus share twisting left. 0 = balanced (mirror-symmetric); +-1 = all one way. Point-curvature statistics cannot see this.

## Electroplating (for ranking designs, not absolute prediction)

| Metric | Unit | t = 0.00 | t = 0.20 | t = 0.40 | t = 0.60 | t = 0.80 | Trend |
|---|---|---|---|---|---|---|---|
| Wagner number, geometry part | 1/m | 2.02 | 1.77 | 1.55 | 1.33 | 1.11 | falls 45% |
| Lattice Wagner number at 0.35 mA/cm2 | - | 0.196 | 0.172 | 0.15 | 0.128 | 0.108 | falls 45% |
| Max current density for Wa >= 1 | mA/cm2 | 0.0685 | 0.0601 | 0.0524 | 0.0449 | 0.0376 | falls 45% |
| Flow needed for Da^-1 = 440 | m/s | 0.0395 | 0.0392 | 0.0383 | 0.0367 | 0.0343 | falls 13% |
| Flow needed for Da^-1 = 2280 | m/s | 0.205 | 0.203 | 0.199 | 0.19 | 0.178 | falls 13% |

**What the metrics mean**

- **Wagner number, geometry part**: Geometry-only uniformity factor eps/(a_v (N_field L)^2): compare designs without bath data. Higher = more even coating.
- **Lattice Wagner number at 0.35 mA/cm2**: Coating uniformity. High = interior plates nearly as thick as the outside; low = current crowds onto the faces nearest the anode (paper: Octet 0.86 ~2x, IsoTruss 0.24 ~4.6x top-to-bottom).
- **Max current density for Wa >= 1**: Plating "speed limit": the highest average current that still meets the target uniformity.
- **Flow needed for Da^-1 = 440**: Mean electrolyte velocity through the lattice needed to reach Da^-1 = 440 (strut-lattice thresholds for >= 90% of bulk concentration range from 440, Octet, to 2280, simple cubic; the TPMS threshold is not yet calibrated).
- **Flow needed for Da^-1 = 2280**: Mean electrolyte velocity through the lattice needed to reach Da^-1 = 2280 (strut-lattice thresholds for >= 90% of bulk concentration range from 440, Octet, to 2280, simple cubic; the TPMS threshold is not yet calibrated).

**Source and caveats.** Lattice Wagner and inverse Damkohler numbers after Clemens et al., *ACS Appl. Eng. Mater.* 2025, 3(3), 600-612, doi:10.1021/acsaenm.4c00561. The equations are reconstructed from the paper's variable definitions: the scaling is right but absolute values may be off by a constant (about 1.5x). Use them to **rank designs**, not as absolute predictions. Wa predicts uniformity only when Da^-1 is high (enough flow); otherwise depletion dominates. The paper's Da^-1 thresholds are for strut lattices; TPMS thresholds still need calibration.

**Assumptions.** Tafel kinetics; uniform electrolyte conductivity; perfectly conductive substrate; convection-dominated ion supply. For a sheet lattice both faces are wetted.

**Process inputs used** (placeholders from the paper unless you changed `cfg.plating`): T = 298 K, alpha_c = 0.485, kappa = 6.4 S/m, i_avg = 0.35 mA/cm2, n = 2, C_o = 10 mol/m^3, v = not set (Da^-1 not computed), N_field = 4 cells (toward the anode), N_flow = 8 cells (along the flow).

## Choosing t for a target wall thickness (sheet)

| t | Wall thickness (mm) | 5-95% spread | Thickness / cell | Solid fraction |
|---|---|---|---|---|
| 0.1 | 0.104 | 0.0938-0.112 | 0.021 | 6% |
| 0.2 | 0.209 | 0.188-0.224 | 0.042 | 13% |
| 0.3 | 0.315 | 0.283-0.337 | 0.063 | 19% |
| 0.4 | 0.422 | 0.38-0.452 | 0.084 | 26% |
| 0.5 | 0.532 | 0.478-0.57 | 0.106 | 32% |
| 0.6 | 0.645 | 0.579-0.691 | 0.129 | 39% |
| 0.7 | 0.763 | 0.684-0.817 | 0.153 | 46% |
| 0.8 | 0.886 | 0.793-0.948 | 0.177 | 52% |
| 0.9 | 1.02 | 0.908-1.09 | 0.203 | 59% |
| 1.0 | 1.16 | 1.03-1.24 | 0.232 | 66% |
| 1.1 | 1.31 | 1.16-1.41 | 0.263 | 73% |
| 1.2 | 1.49 | 1.32-1.6 | 0.299 | 80% |
| 1.3 | 1.72 | 1.49-1.84 | 0.344 | 88% |
| 1.4 | 2.05 (unreliable: sheet near breakup) | 1.71-2.26 | 0.411 | 96% |

Rule of thumb: thickness = 0.21 x t x cell size for t up to about 0.5; above that the exact value grows faster. Use `tpms_t_for_thickness(target, cellSize)` for an exact inverse. The full table is in `gyroid_t_vs_wall_thickness.csv`.

## Electroplating numbers vs t (sheet lattice, ranking only)

| t | Porosity | a_v (1/m) | Wa_geom (1/m) | Wa | i_max for Wa >= 1 (mA/cm2) | Flow for Da^-1 = 440 (m/s) | Flow for Da^-1 = 2280 (m/s) |
|---|---|---|---|---|---|---|---|
| 0.0 | 100% | 1237 | 2.02 | 0.196 | 0.0685 | 0.0395 | 0.205 |
| 0.1 | 94% | 1235 | 1.89 | 0.183 | 0.0642 | 0.0394 | 0.204 |
| 0.2 | 87% | 1228 | 1.77 | 0.172 | 0.0601 | 0.0392 | 0.203 |
| 0.3 | 81% | 1217 | 1.66 | 0.161 | 0.0562 | 0.0388 | 0.201 |
| 0.4 | 74% | 1200 | 1.55 | 0.15 | 0.0524 | 0.0383 | 0.199 |
| 0.5 | 68% | 1179 | 1.44 | 0.139 | 0.0487 | 0.0376 | 0.195 |
| 0.6 | 61% | 1151 | 1.33 | 0.128 | 0.0449 | 0.0367 | 0.19 |
| 0.7 | 54% | 1117 | 1.22 | 0.118 | 0.0413 | 0.0356 | 0.185 |
| 0.8 | 48% | 1075 | 1.11 | 0.108 | 0.0376 | 0.0343 | 0.178 |
| 0.9 | 41% | 1023 | 1 | 0.0969 | 0.0339 | 0.0327 | 0.169 |
| 1.0 | 34% | 960 | 0.886 | 0.0858 | 0.03 | 0.0306 | 0.159 |
| 1.1 | 27% | 881 | 0.767 | 0.0742 | 0.026 | 0.0281 | 0.146 |
| 1.2 | 20% | 778 | 0.635 | 0.0615 | 0.0215 | 0.0248 | 0.129 |
| 1.3 | 12% | 636 | 0.484 | 0.0469 | 0.0164 | 0.0203 | 0.105 |
| 1.4 | 4% | 387 | 0.288 | 0.0279 | 0.010 | 0.0124 | 0.0641 |

The full table is in `gyroid_electroplating.csv`.
