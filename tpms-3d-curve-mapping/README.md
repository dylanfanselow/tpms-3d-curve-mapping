# TPMS 3D Curve Mapping

MATLAB tools for **3D curve mapping on triply periodic minimal surface (TPMS) lattices**. The tools trace and measure curve networks across the surface, computed directly from the analytic level-set equation rather than from meshed STL files. The current release implements the **gyroid**; other TPMS types are planned.

![Lines of curvature and asymptotic curves on the gyroid at t = 0 and t = 0.8](docs/images/preview_curve_networks.png)

## What it does

Triply periodic minimal surfaces (TPMS) such as the gyroid are widely used in additively manufactured lattices. They are usually characterised point by point, using curvature values at each spot on a meshed surface. This project adds information that per-point curvature cannot provide: **how the surface curves are organised across the structure**.

It works from the surface's level-set equation. For the gyroid, the surface implemented so far, that is

```
sin x cos y + sin y cos z + sin z cos x = t
```

It traces two families of curves on the surface:

- **Lines of curvature** follow the directions of maximum and minimum bending.
- **Asymptotic curves** follow the directions along which the surface does not bend.

From these it reports plain-language metrics, for example:

| Metric | What it tells you |
|---|---|
| Straightness | How far a curve gets compared with how far it travels (0 to 1) |
| Direction memory | How far a curve runs before it "forgets" its starting direction |
| Closed-loop share | How many curves close on themselves, e.g. as rings around struts |
| Turning radius and twist pitch | How sharply the curves steer and corkscrew, in mm |
| Handedness index | Whether twisting is balanced or biased left/right; point curvature can't see this |
| Wall thickness, solid fraction, S/V | What the level-set value t physically does to the part |

Each curve-based number comes with a 95% confidence interval, so it is clear how trustworthy it is.

## Example results

For a sheet gyroid with a 5 mm cell, as the level-set value t rises from 0 to 0.8 (walls thicken from 0 to 0.89 mm):

- Max-bending lines lose almost all straightness (0.44 → 0.05) as they wrap around the thickening struts.
- Min-bending lines go from never closing to 100% closed loops.
- The handedness index rises from 0.00 (balanced) to 0.17, revealing the gyroid's chirality.

The full output from this run is in [`examples/results/gyroid_summary.md`](examples/results/gyroid_summary.md).

## Quick start

1. Open `run_gyroid_curve_mapping.m` in MATLAB.
2. In the configuration block at the top, set `cfg.cellSize` (your cell size, default 5 mm) and `cfg.latticeType` (`'sheet'` or `'network'`).
3. Run it. Results are written to `results/gyroid_curve_results/`. Start with `gyroid_summary.md`.

Other entry points:

- `run_sampling_study.m` shows how the results depend on sampling settings and recommends values.
- `tpms_t_for_thickness(0.4, 5)` returns the level-set value t that gives 0.4 mm walls in a 5 mm cell. Run `addpath src` first if you haven't run a script in this session.

### Requirements

- MATLAB R2016b or later.
- The Symbolic Math Toolbox is optional. Without it, set `cfg.derivMethod = 'manual'`, which uses hand-derived derivatives with identical results.
- The code is also tested in GNU Octave 8.4 with the `symbolic` package; figure generation there is untested.

## Outputs

| File | Contents |
|---|---|
| `gyroid_summary.md` / `.csv` | Plain-language metrics per t, with meanings and 95% confidence intervals |
| `gyroid_t_vs_wall_thickness.csv` | t → wall thickness, solid fraction, surface-to-volume ratio |
| `gyroid_curve_descriptors_PCA.csv` | Descriptors carrying genuinely new information, one row per t, ready for PCA |
| `gyroid_curve_descriptors_long.csv` | Every descriptor, tagged by role (new, redundant, reference, …) |
| `gyroid_curves_percurve.csv` | One row per traced curve |
| `*.png` | Curve networks, distributions, trends |

## Validation

The run starts with 26 self-checks. They include:

- symbolic vs hand-derived derivatives;
- closed-form vs eigen-solved curvature;
- analytic vs finite-difference geodesic curvature (second-order convergence confirmed);
- known differential-geometry identities;
- curves staying on the surface to about 1e-15;
- cubic and inversion symmetry of the gyroid.

## Repository layout

```
run_gyroid_curve_mapping.m   main analysis
run_sampling_study.m         sampling sensitivity study
src/                         all functions (tpms_*.m)
examples/                    outputs from a test run
docs/TECHNICAL.md            full method notes, conventions, limitations
```

## Documentation

[`docs/TECHNICAL.md`](docs/TECHNICAL.md) covers:

- the mathematics and sign conventions;
- which descriptors are genuinely new and which restate point curvature;
- the tortuosity normalisation;
- the wall-thickness method;
- the sampling-study results;
- known limitations.

## Status

This is a research tool. The 3D curve-mapping method is general to any TPMS defined by a level-set equation, but this release implements and validates it for the **gyroid only**. The code that derives the analytic fields accepts any level-set expression, so adding other TPMS types (Diamond, Schwarz P, Neovius, IWP) and blends is mainly a matter of adding their equations and re-running the validation. File names such as `run_gyroid_curve_mapping.m` and `gyroid_summary.md` refer to the gyroid run specifically.

## Citation

If you use this in academic work, please cite it using [`CITATION.cff`](CITATION.cff). On GitHub, the "Cite this repository" button generates the citation for you.

## License

Copyright (C) 2026 YOUR NAME

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [GNU General Public License](LICENSE) for more details.
