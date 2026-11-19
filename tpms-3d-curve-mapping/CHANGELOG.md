# Changelog

## 0.6.0

### Added
- **Non-cubic add-on.** Lattices whose repeating cell is not a cube, still written as trig level sets:
  - Lonsdaleite (hexagonal, orthohexagonal box);
  - Schwarz tP, tD and gyroid tG (tetragonal, c/a = 1.5).

  They come from `src/tpms_library_noncubic.m` and are merged by `src/tpms_library_all.m`. They are flagged in the new `cell` CSV column and starred (*) in the report and heatmap. `cfg.includeNonCubic` switches them off.
- `src/tpms_period3.m`: per-axis box period, used wherever a 2π cube was assumed.

### Changed
- Per-axis periods in sampling, seeding, closure detection, straightness, the genus translation search, wall-thickness and plating grids, and the duplicate check (same-box pairs only, box-preserving axis swaps). For cubic cells the arithmetic is unchanged.
  - The gyroid pipeline outputs are byte-identical to the previous version.
  - A cubic library run matches exactly on all metric columns.
- `examples/results/` refreshed. The previous files predated the v0.5 two-face wetted-area change.

### Fixed
- Corrected the explanation for Q*: its t = 0 level set is singular. It was not sampled in the wrong cell.

## 0.5.0

### Added
- `src/tpms_library.m`: 56 entries, 54 distinct TPMS and related level-set lattices.
  - Equations checked against the Fisher et al. 2023 catalog (*Data in Brief*; its MATLAB files), Bernhard 2026 (arXiv:2604.07160, Tables 6–7), nTop, Gandy et al. 2001, Wohlgemuth et al. 2001 and Vasile et al. 2024.
  - Conflicting published variants are kept as separate entries.
- `run_tpms_library.m`: maps the core metrics across the library at t = 0 and ranks the most extreme lattices.
  - Resumable, with one saved file per lattice.
  - Writes `tpms_library_report.md`, `tpms_library_metrics.csv` and a heatmap.
- `src/tpms_library_check.m`: per-surface checks.
  - Genus per primitive cell, from Gauss–Bonnet and translation symmetries (π/4 grid), with automatic grid refinement.
  - Side-swapping translations.
  - Minimality, smoothness, cubic symmetry and curvedness-change length.
- `src/tpms_library_duplicates.m`: proves which equations describe the same lattice, run automatically by `run_tpms_library`. Proven duplicates are marked `duplicateOf` and skipped.
- Sources added: Schwarz & Gompper (arXiv:cond-mat/9903406), MiniSurf (Hsieh & Valdevit 2020) and Lord & Mackay 2003 (genus values).
- `src/tpms_library_plating.m`: lattice Wagner and inverse Damköhler numbers for every library lattice at a matched solid fraction (default 16%). Shown in the report, CSV and heatmap.
- `examples/tpms_library/`: report and metrics from a full run.

### Changed
- `tpms_wall_thickness`: the sheet wetted area now adds the areas of F = +t and F = −t, which differ for unbalanced surfaces. Gyroid results are unchanged.
- `tpms_build_fields` accepts any library key, library entry or `struct('key', ..., 'expr', ...)`. The gyroid results are unchanged.

### Found
- The widely used nodal Lidinoid equation gives a genus-17 surface per primitive cell (true Lidinoid: genus 3).
- Axolotl's Double Primitive form and the Tubular Diamond as published break into closed pieces.
- D's alternative form and Double Diamond 2 duplicated earlier entries. Four duplicate candidates were rejected during the literature search.
- corax's C(Y) has genus 21, not C(Y)'s 13. The Fisher et al. form matches.

## 0.4.0

### Changed
- New lean core set as the default (`cfg.descriptorSet = 'core'`). Five curvature descriptors are the PCA inputs and the summary's curvature section:
  - mean-curvature gradient;
  - pooled geodesic curvature of the lines of curvature;
  - straightness of k1-lines and of k2-lines;
  - handedness.
- All other descriptors are still computed and saved with role `EXTENDED`. `'full'` restores every descriptor as a PCA input.
- The turning radius in the summary now comes from the geodesic-curvature surface field, so it has no sampling noise. The curve-sampled value remains in full mode.
- Default `nSeeds` lowered from 96 to 48. The sampling study shows at most 35 seeds are needed for the core set; tracing time halves.
- Example outputs (`examples/results/`, `examples/sampling_study/`) regenerated with v0.4.0.

### Added
- `field_kg_LoC`: pooled geodesic curvature of both line-of-curvature families.
- `handedness_LoC`: the handedness index as a PCA scalar.
- "Curvature-change length" in the summary: typical curvature ÷ typical mean-curvature gradient, in mm.
- `cfg.coreBuildAlignment` to add build-direction alignment to the core set.
- An early error for an invalid `cfg.descriptorSet`.

## 0.3.0

### Added
- Electroplating descriptors after Clemens et al., *ACS Appl. Eng. Mater.* 2025, 3(3), 600–612. These use reconstructed equations and are for ranking designs only.
  - Lattice Wagner number (coating uniformity) and inverse Damköhler number (depletion and dendrite risk).
  - Geometry-only parts (Wa_geom, Da_geom) for ranking without bath data.
  - Maximum current density for a target Wa, and the electrolyte flow needed for target Da⁻¹ values.
  - New `gyroid_electroplating.csv`, an electroplating section in `gyroid_summary.md`, and the `cfg.plating` settings with placeholder bath values from the paper.
  - `src/tpms_plating_numbers.m` takes porosity, wetted area and cell size, so it also works for strut and STL-based lattices.
- Validation check 10: self-consistency of the electroplating formulas.

### Changed
- Example outputs in `examples/results/` regenerated with v0.3.0, including `gyroid_electroplating.csv`.

## 0.2.0

### Added
- Plain-language summary (`gyroid_summary.md` / `.csv`): about 20 metrics per level set, each with a meaning, a trend and a 95% bootstrap confidence interval.
- Wall-thickness table (`gyroid_t_vs_wall_thickness.csv`):
  - exact sheet thickness along the mid-surface normal;
  - solid fraction (sheet and network), surface-area density and S/V;
  - `tpms_t_for_thickness` for the inverse (target thickness → t).
- Straightness (chord/window, bounded 0..1) for windows of ¼, ½, 1 and 2 cells.
- Direction memory (persistence length) per curve family.
- `run_sampling_study.m`:
  - varies the seed count, surface grid and step size one at a time;
  - recommends settings;
  - estimates the seeds each metric needs.
- Downward-facing overhang fraction and build-direction alignment of the principal directions.

### Changed
- Straightness replaces tortuosity as the PCA input. Tortuosity is unbounded when curves loop and stays in the long table.
- Default `nSeeds` raised from 48 to 96. Default `cellSize` is 5 (mm).
- Figures are forced to a light style with larger, black text, so they read correctly in MATLAB's dark theme.

### Fixed
- Closure detection now uses the seed-to-segment distance. Previously coarse steps missed closed loops.
- The handedness index uses the signed share of twisting length. The length-weighted mean torsion of a line of curvature telescopes to about 0, so it cannot measure handedness.
- MATLAB compatibility: invalid struct field name `asym-` in the validation.

## 0.1.0

### Added
- Analytic gyroid fields (Symbolic Toolbox via `matlabFunction`, plus a hand-derived cross-check).
- Lines of curvature and asymptotic curves:
  - traced by vectorised adaptive RK4 with Newton projection onto the surface;
  - with umbilic and parabolic-line handling.
- Darboux and Frenet descriptors. Geodesic curvature is computed analytically via Codazzi (lines of curvature) and Liouville (asymptotic curves).
- A role tag on every descriptor (PCA / CHECK / REDUNDANT / MIXED / REFERENCE / CONTEXT) separating new information from restated point curvature.
- 26-check validation suite.
