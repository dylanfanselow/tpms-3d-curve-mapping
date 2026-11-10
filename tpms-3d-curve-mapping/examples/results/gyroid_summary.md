# Gyroid curve-mapping summary

Lattice: **sheet** gyroid, cell size **5 mm**. Build direction [0  0  1]. Generated 2026-10-20 00:45.

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
| Max-bending lines (k1): straightness over 1 cell | 0-1 | 0.44 [0.435, 0.446] | 0.353 [0.349, 0.364] | 0.254 [0.247, 0.265] | 0.143 [0.133, 0.157] | 0.047 [0.0377, 0.0503] | falls 89% |
| Max-bending lines (k1): direction memory | mm | 1.35 [1.35, 1.36] | 1.24 [1.24, 1.24] | 1.12 [1.12, 1.13] | 1 [1, 1.01] | 0.872 [0.87, 0.875] | falls 36% |
| Max-bending lines (k1): share forming closed loops | % | 0 | 0 | 0 | 0 | 0 | steady |
| Min-bending lines (k2): straightness over 1 cell | 0-1 | 0.441 [0.436, 0.444] | 0.505 [0.5, 0.506] | 0.552 [0.551, 0.554] | 0.6 [0.599, 0.6] | 0.641 [0.641, 0.642] | rises 45% |
| Min-bending lines (k2): direction memory | mm | 1.35 [1.35, 1.36] | 1.47 [1.47, 1.47] | 1.58 [1.58, 1.58] | 1.7 [1.7, 1.7] | 1.83 [1.82, 1.83] | rises 35% |
| Min-bending lines (k2): share forming closed loops | % | 0 | 0 | 67.7 [57.3, 76.0] | 97.9 [94.8, 100] | 100 | rises from 0 to 100 |
| Bending lines: sideways turning radius | mm | 5.45 [5.34, 5.59] | 5.42 [5.26, 5.53] | 5.3 [5.13, 5.48] | 5.09 [5.02, 5.19] | 5 [4.95, 5.05] | falls 8% |
| Bending lines: twist pitch | mm | 8.63 [8.27, 9.14] | 8.4 [7.95, 9.13] | 7.98 [7.71, 8.34] | 8.35 [8.19, 8.44] | 7.66 [7.49, 7.97] | falls 11% |
| Handedness index | -1..+1 | 0.002 [-0.0133, 0.0148] | 0.0725 [0.0595, 0.0853] | 0.125 [0.113, 0.138] | 0.147 [0.139, 0.153] | 0.169 [0.166, 0.172] | rises from 0.002 to 0.169 |
| Zero-bending lines: straightness over 1 cell | 0-1 | 0.932 [0.932, 0.933] | 0.933 [0.932, 0.934] | 0.933 [0.93, 0.935] | 0.93 [0.928, 0.933] | 0.94 [0.938, 0.941] | steady |
| Zero-bending lines: direction memory | mm | 6.41 [6.38, 6.44] | 6.75 [6.49, 7.03] | 7.55 [7.1, 8.03] | 9.35 [8.65, 9.89] | > 9.99 | rises 56% |
| Zero-bending families imbalance | -1..+1 | 0 | 0.0917 | 0.181 | 0.248 | 0.292 | rises from 0 to 0.292 |
| Max-bending direction near vertical (< 30 deg) | % | 9.14 | 10.1 | 10.8 | 11.4 | 11.5 | rises 26% |

**What the metrics mean**

- **Max-bending lines (k1): straightness over 1 cell**: 1 = runs straight, 0 = comes back to where it started. Net progress per distance travelled.
- **Max-bending lines (k1): direction memory**: How far a line runs before it has "forgotten" its starting direction.
- **Max-bending lines (k1): share forming closed loops**: Lines that close on themselves, e.g. rings around struts.
- **Min-bending lines (k2): straightness over 1 cell**: 1 = runs straight, 0 = comes back to where it started. Net progress per distance travelled.
- **Min-bending lines (k2): direction memory**: How far a line runs before it has "forgotten" its starting direction.
- **Min-bending lines (k2): share forming closed loops**: Lines that close on themselves, e.g. rings around struts.
- **Bending lines: sideways turning radius**: Typical radius of the sideways turns the lines make within the surface (small = sharp steering).
- **Bending lines: twist pitch**: Distance along a line for one full corkscrew twist (small = strongly twisted).
- **Handedness index**: Share of line length twisting right minus share twisting left. 0 = balanced (mirror-symmetric); +-1 = all one way. Point-curvature statistics cannot see this.
- **Zero-bending lines: straightness over 1 cell**: Asymptotic curves: paths along which the surface does not bend. Near 1 = almost straight.
- **Zero-bending lines: direction memory**: How far a zero-bending line runs before it has "forgotten" its starting direction ("> x" = longer than the longest lag checked).
- **Zero-bending families imbalance**: Difference between the two mirror-image families of zero-bending lines; 0 when the part is mirror-balanced.
- **Max-bending direction near vertical (< 30 deg)**: Share of surface whose sharpest bend runs close to the build direction (13.4% if directions were random).

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
