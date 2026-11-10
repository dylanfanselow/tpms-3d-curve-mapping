# Changelog

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
