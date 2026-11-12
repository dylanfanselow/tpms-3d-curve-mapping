# Technical documentation

Full method notes for TPMS 3D curve mapping (Method 1, surface curve networks), currently implemented for the gyroid. For a quick start see the [README](../README.md).

Lines of curvature (k1, k2 families) and asymptotic curves (asym+, asym−) traced directly on the analytic gyroid level set

  sin x cos y + sin y cos z + sin z cos x = t

Every quantity is evaluated from analytic derivatives of F generated with the Symbolic Math Toolbox (`matlabFunction`). No STL is used anywhere. A marching-cubes mesh supplies quadrature points and area weights only. Each vertex is Newton-projected onto the exact surface, and all geometry is evaluated analytically at the projected point.

## Version 0.2 additions

- **A plain-language summary** (`gyroid_summary.md`, start here). It has about 20 numbers per level set, each with a one-line explanation and a 95% confidence interval. See *Summary for people* below.
- **A table mapping t to wall thickness** (`gyroid_t_vs_wall_thickness.csv`). Sheet thickness is measured exactly along the mid-surface normal, alongside solid fraction (sheet and network), surface-area density and S/V. `tpms_t_for_thickness(target, cellSize)` gives the inverse.
- **Straightness replaces tortuosity as the PCA input.**
  - Straightness = chord/window, bounded between 0 and 1. Tortuosity (window/chord) blows up when curves loop: medians of 25 and 95th percentiles of 77 at t = 0.8.
  - Straightness is computed for windows of ¼, ½, 1 and 2 cells.
  - Tortuosity stays in the long table.
- **Direction memory (persistence length)** of each curve family: the arc length over which a curve "forgets" its direction.
- **`run_sampling_study.m`** varies the seed count, surface grid and step size one at a time, then recommends settings and estimates how many seeds each metric needs.
- **Fixes:**
  - Closure detection now uses the distance from the seed to each step segment. Previously coarse steps missed closures, cutting the k2-loop share by 23 points.
  - Handedness is now measured correctly (see below).
- **Defaults:** `cfg.cellSize = 5` (mm) and `cfg.nSeeds = 96`.

## Running

Open `run_gyroid_curve_mapping.m` in the repository root and run it. The script adds `src/` to the path. All settings are in the configuration block at the top. **Set `cfg.cellSize` to your cell size** (default 5 mm) and `cfg.latticeType` to `'sheet'` or `'network'`.

Run `run_sampling_study.m` to see how the numbers depend on sampling settings. It takes about 4 minutes and needs no Symbolic Toolbox.

Requirements are R2016b or later (implicit expansion). The Symbolic Math Toolbox is needed only for `cfg.derivMethod = 'symbolic'`. The `'manual'` option uses hand-derived gyroid derivatives, needs no toolbox, and gives identical results; validation check 1 compares the two. The default run (5 level sets, 96 seeds per family, 4 families) takes about 2.5 minutes in Octave, including validation and the thickness table. MATLAB should be faster.

Internally coordinates are in radians, so one cell is [0, 2π]³. With `cfg.cellSize` set, every output is converted to physical units: lengths ×λ, curvatures ÷λ, |∇H| ÷λ², and so on, where λ = cellSize/2π. With `cfg.cellSize = []`, the summary gives lengths as fractions of the cell.

## Files

The two `run_*.m` scripts are in the repository root; every `tpms_*.m` function is in `src/`.

| File | Purpose |
|---|---|
| `run_gyroid_curve_mapping.m` | Driver: config, validation, t-sweep, output |
| `tpms_build_fields.m` | Symbolic → numeric F, ∇F, Hess F, 3rd derivatives (plus the hand-derived cross-check) |
| `tpms_geometry.m` | k1, k2, H, K, principal and asymptotic directions, curvature derivatives, Codazzi/Liouville geodesic curvature |
| `tpms_project.m` | Newton projection onto F = t |
| `tpms_trace_curves.m` | Vectorised adaptive RK4 tracer for line fields, with stop rules |
| `tpms_curve_descriptors.m` | Darboux/Frenet quantities along curves, windowed tortuosity, analytic-vs-numeric checks |
| `tpms_surface_sample.m`, `tpms_seed_points.m` | Quadrature points and periodic farthest-point seeding |
| `tpms_analyze_level.m` | Builds the descriptor catalogue for one t, each descriptor tagged with a role |
| `tpms_validate.m` | 26 self-checks (see below) |
| `tpms_write_results.m`, `tpms_plot_results.m` | CSV/MAT output and figures |
| `tpms_wstats.m`, `tpms_default_trace_opts.m` | Weighted statistics (same set as the STL tool) and tracer options |
| `tpms_summary_metrics.m` | The plain-language metrics, with bootstrap CIs over curves |
| `tpms_write_summary.m`, `tpms_plot_summary.m` | Summary Markdown/CSV and trend plus wall-thickness figures |
| `tpms_wall_thickness.m`, `tpms_t_for_thickness.m` | t → wall thickness, solid fraction, S/V; and the inverse |
| `tpms_persistence.m` | Direction memory from the tangent autocorrelation |
| `run_sampling_study.m`, `tpms_metric_deviation.m` | Sampling sensitivity study |

Outputs go to `cfg.outDir`, which defaults to `results/gyroid_curve_results/` in the repository (git-ignored). Example outputs from a test run are in `examples/`. Start with `gyroid_summary.md`. The rest:

- `gyroid_summary.csv` holds the same table as the Markdown summary. `gyroid_summary_detail.csv` adds confidence intervals.
- `gyroid_t_vs_wall_thickness.csv` is the thickness lookup table.

- `…_PCA.csv` has one row per t and holds only descriptors with role PCA.
- `…_long.csv` holds every descriptor with its role and a note.
- `…_checks.csv` holds numerical checks per t and family.
- `…_percurve.csv` has one row per curve: length, stop reasons, closure, lattice translation, closure gap.
- The `.mat` file holds the curves and surface samples.

## Summary for people

Each metric in `gyroid_summary.md` has a plain name, a unit, one value per t, a trend ("rises 45%") and a one-line meaning. Lengths are in your cell-size unit. Brackets show 95% confidence intervals from resampling the traced curves (a narrow interval means the number is trustworthy). For wall thickness and bending radius, the brackets instead give the 5–95% spread over the surface.

| Metric | Plain meaning | Source |
|---|---|---|
| Wall thickness, solid fraction, S/V, surface per volume | What t physically does to the part | Exact geometry, no curves |
| Typical surface bending radius | Radius of the sphere that bends like the surface does | Point-based (same as the STL tool), for reference |
| Downward-facing overhang | Share of surface facing down more steeply than 45° | Normal field |
| Straightness over 1 cell (k1, k2, zero-bending) | Net progress per distance travelled, 0..1 | Traced curves |
| Direction memory | How far a line runs before forgetting its direction | Traced curves |
| Share forming closed loops | Lines that close on themselves (e.g. rings round struts) | Traced curves |
| Sideways turning radius | Radius of the turns the lines make within the surface | Codazzi geodesic curvature |
| Twist pitch | Distance for one full corkscrew twist | Frenet torsion |
| Handedness index | Twisting right minus twisting left, −1..+1 | Torsion sign |
| Zero-bending families imbalance | Mirror imbalance of the two asymptotic families | Liouville geodesic curvature |
| Max-bending direction near vertical | Share of surface whose sharpest bend runs within 30° of the build direction (13.4% if random) | Direction field |

**About the handedness index.** For a line of curvature, torsion is τ = −dφ/ds, so its length-weighted mean telescopes to (φ_start − φ_end)/length. That is about 0 for any long curve, whatever the geometry. The index therefore uses the share of length twisting each way instead. It is 0.002 at t = 0, with a confidence interval straddling zero, as the inversion symmetry requires. It then rises steadily to 0.17 at t = 0.8.

## Tortuosity: scaling and normalisation

- **Cell size.** Tortuosity and straightness are ratios of two lengths. With the window defined in cell lengths (`cfg.windowsCells`), they are identical for any cell size, so nothing needs rescaling. Only a window fixed in mm would change with cell size.
- **Window length matters most.** For k1-lines at t = 0, tortuosity goes 1.06 → 1.23 → 2.30 → 5.48 for windows of ¼, ½, 1 and 2 cells. Any single tortuosity number is tied to its window, so the window is always stated in the output.
- **Step size doesn't matter.** Halving the tracing step changes tortuosity by less than 0.2%.
- **Use straightness, not tortuosity.** Straightness = chord/window, which is 1/tortuosity, bounded 0..1 and well-behaved when curves loop.
- **Direction memory is a complementary measure.** It is a single length per family that needs no window choice.

## Wall thickness vs t (5 mm cell, sheet gyroid)

| t | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 | 1.2 |
|---|---|---|---|---|---|---|
| Mean wall (mm) | 0.209 | 0.422 | 0.645 | 0.886 | 1.16 | 1.49 |
| Exact ÷ first-order 2t/\|∇F\| | 1.00 | 1.01 | 1.03 | 1.06 | 1.11 | 1.20 |
| Solid fraction (sheet) | 13% | 26% | 39% | 52% | 66% | 80% |
| Solid fraction (network) | 56% | 63% | 70% | 76% | 83% | 90% |

Rule of thumb: wall ≈ 0.21 × t × cell size up to t ≈ 0.5, growing faster beyond that. Within the middle 90% of the wall area, thickness varies by about ±10% around the mean.

Above t ≈ 1.35 the normal rays start missing the opposite face because the sheet is close to breaking apart. Those rows are flagged `reliable = 0`.

The thickness is the distance between the two faces measured along the mid-surface normal. A maximal-inscribed-sphere thickness (as used by some CT software) would read slightly smaller where the wall is curved.

## Sampling study results (t = 0.4)

- **Grid and step size:** the coarsest tested settings (nGrid 32, step ×2) already pass. Results don't depend on them.
- **Number of seeds:** this matters only for metrics that are medians over individual curves. Seeds needed for a ±3% interval (±10 points for loop shares, ±0.02 for 0–1 metrics):

  | Metric | Seeds needed |
  |---|---|
  | Straightness, k1/k2 direction memory, handedness | ≤ 40 |
  | Turning radius, k2 loop share | about 60–100 |
  | Twist pitch | about 320 |
  | Zero-bending direction memory | about 440 |

  The default of 96 gives about ±4–6% on those last two, which their confidence intervals in the summary show.

## Electroplating descriptors (v0.3)

Lattice Wagner and inverse Damköhler numbers after Clemens et al., *ACS Appl. Eng. Mater.* 2025, 3(3), 600–612, doi:10.1021/acsaenm.4c00561. The equation images were unavailable, so the equations below are reconstructed from the paper's variable definitions. The scaling is right, but absolute values may be off by a constant (about 1.5× when checked). **Use them to rank designs, not as absolute predictions.**

**Lattice Wagner number (coating uniformity).** The ratio of surface-reaction resistance to electrolyte resistance:

  Wa = (R T / (α_c F)) · κ · ε / (i_avg · a_v · (N_field L)²)

- **High Wa:** current spreads evenly.
- **Low Wa:** current crowds onto the outer faces nearest the anode. In the paper, Octet had Wa = 0.86 at 0.35 mA/cm² with about 2× top-to-bottom thickness variation; IsoTruss had 0.24 at 1.4 mA/cm² with about 4.6×.

**Lattice inverse Damköhler number (depletion and dendrite risk).** The rate at which flow delivers metal ions versus the rate plating consumes them:

  Da⁻¹ ∝ n F v C_o / (i_avg · a_v · N_flow L)

The exact prefactor is uncertain and probably includes a factor like n/α_c. The paper's thresholds for keeping the surface concentration at ≥ 90% of bulk, for strut lattices with ε = 0.84, N_field = 4 and N_flow = 8, were:

| Lattice | Da⁻¹ threshold |
|---|---|
| Simple cubic | 2280 |
| FCC | 920 |
| BCC | 800 |
| IsoTruss | 600 |
| Octet | 440 |

TPMS thresholds are unknown and need calibration. Wa predicts uniformity only when Da⁻¹ is high; otherwise depletion dominates.

**What the code computes.** In `tpms_plating_numbers.m`, first the geometry-only parts, so designs can be ranked without bath data:

- Wa_geom = ε / (a_v (N_field L)²), in 1/m
- Da_geom = 1 / (a_v N_flow L), dimensionless

Then, given process inputs:

- Wa = Wa_geom · (R T/(α_c F)) · κ / i_avg
- Da⁻¹ = Da_geom · n F v C_o / i_avg
- **Speed limit:** i_max = (R T/(α_c F)) · κ · Wa_geom / Wa_target
- **Required flow:** v_min = Da⁻¹_target · i_avg / (Da_geom n F C_o), for each target in `cfg.plating.DaInvTarget` (default 440 and 2280, spanning the strut thresholds)

**Geometry inputs.**

- **Sheet lattice:** ε = 1 − solid fraction and a_v = 2A / cell volume, because both channel systems are open to the electrolyte.
- **Network lattice:** ε = 1 − solid fraction and a_v = A / cell volume.

A and the solid fractions come from `tpms_wall_thickness`. **N_field** is the number of cells toward the anode; if the part is plated from several sides, use the depth from the nearest exposed face. **N_flow** is the number of cells along the flow direction.

**Process inputs.** These are placeholders from the paper; replace them with your own bath values in `cfg.plating`:

| Input | Value |
|---|---|
| T | 298 K |
| α_c | 0.485 |
| κ | 6.4 S/m |
| n | 2 |
| C_o | 10 mol/m³ |
| i_avg | 0.35 mA/cm² |
| v | not set |

**Assumptions:** Tafel kinetics; uniform electrolyte conductivity, with no tortuosity correction beyond ε; a perfectly conductive substrate; convection-dominated ion supply.

**Gyroid results** (5 mm cell, N_field = 4, N_flow = 8, placeholder bath):

| t | Porosity (sheet) | a_v sheet (1/m) | Wa_geom (1/m) | Wa at 0.35 mA/cm² | i_max for Wa = 1 (mA/cm²) | Flow for Da⁻¹ = 440 / 2280, sheet (m/s) |
|---|---|---|---|---|---|---|
| 0.2 | 87% | 1228 | 1.77 | 0.172 | 0.060 | 0.039 / 0.20 |
| 0.4 | 74% | 1200 | 1.55 | 0.150 | 0.052 | 0.038 / 0.20 |
| 0.8 | 48% | 1075 | 1.11 | 0.108 | 0.038 | 0.034 / 0.18 |
| 1.2 | 20% | 778 | 0.64 | 0.062 | 0.022 | 0.025 / 0.13 |

What these show:

- **Thicker walls plate less evenly.** Wa falls with t, because porosity drops faster than surface area.
- **Sheet and network gyroids at the same t have identical Wa.** The network's void is exactly one of the sheet's two channels: half the porosity and half the wetted area, so ε/a_v is unchanged. The sheet needs twice the flow for the same Da⁻¹, because twice as much surface consumes ions.
- **Absolute comparisons with the paper's strut lattices are not possible here,** because their cell sizes and surface areas aren't in the reconstructed equations. Compare designs within the same model.

## Conventions

- **Normal.** n = ∇F/|∇F| (fixed global orientation).
- **Shape operator.** S = −P·HessF·P/|∇F|, and κ_n(u) = ⟨γ″, n⟩ = −uᵀ·HessF·u/|∇F|. Under this convention a sphere F = |x|² − R² has k1 = k2 = −1/R.
- **Principal frame.** k1 ≥ k2, and d2 = n × d1, so (d1, d2, n) is right-handed.
- **Asymptotic directions.** a± = cos θ_A·d1 ± sin θ_A·d2, with tan θ_A = √(k1/(−k2)). These exist only where K < 0.
- **Label invariance.** Because d2 = n × d1, the labels k1, k2, asym+ and asym− are geometric: they do not depend on the arbitrary sign of d1. The asym+ family has τ_g = −√(−K) and asym− has τ_g = +√(−K).
- **Signs.** Signs of k1, k2 and H flip if the normal orientation is flipped. K, |κ_g| and the Frenet torsion τ do not.

## What is new, and what is not

The handoff warned that curve statistics can quietly restate k1 and k2. Every descriptor therefore carries a role in the output.

| Role | Descriptors | Why |
|---|---|---|
| **PCA** (new) | `field_gradH`, `field_gradK` | 3rd-derivative information: how curvature varies across the surface |
| **PCA** (new) | `field_kg_k1lines`, `field_kg_k2lines` | Geodesic curvature of the lines of curvature, via Codazzi: κ_g = e₂(k1)/(k1−k2) and e₁(k2)/(k1−k2) |
| **PCA** (new) | `field_kg_asymP`, `field_kg_asymM` | Geodesic curvature of the asymptotic curves, via Liouville: κ_g = θ′ + cos θ·κ_g1 + sin θ·κ_g2 |
| **PCA** (new) | `field_d1_build`, `field_d2_build`, `A4_d1` | Direction fields relative to the build direction, and cubic anisotropy of d1. The STL tool discards directions. |
| **PCA** (new) | `curve_k1_tau`, `curve_k2_tau` | Frenet torsion of the lines of curvature (τ = −dφ/ds, φ = atan2(κ_g, κ_n)). Signed, so it is sensitive to chirality. |
| **PCA** (new, non-local) | `curve_*_straight`, `*_persistence`, `*_frac_singular_ends`, `*_mean_length_rel`, `curve_asym*_build` | Need the whole curve, not a point. `curve_*_straight` replaces tortuosity as the PCA input. |
| CONTEXT | `curve_*_tort`, `curve_*_straight_W*` | Tortuosity (unbounded) and straightness at the other windows |
| CHECK | `curve_*_kg` | The κ_g field re-sampled along curves. Length-weighted, so biased by curve density (medians about 10–30% below the area-weighted field). |
| MIXED | `curve_k*_kappa` | √(k_i² + κ_g²) mixes a redundant part with a new part |
| REDUNDANT | `curve_asym*_tau` | τ = ∓√(−K) (Beltrami–Enneper) |
| REFERENCE | `ref_k1`, `ref_k2`, `ref_H`, `ref_K` | Point-based, identical in meaning to the STL tool. Kept for the analytic-vs-STL comparison (build step 2). |

The redundant identities are checked numerically but not exported as descriptors: κ_n = k_i and τ_g = 0 on lines of curvature, κ_n = 0 on asymptotic curves.

Most of the new information turns out to be fields computable at every surface point, which matches the handoff's own conclusion in the roughness discussion. Tracing is still required for three things:

- the non-local descriptors (straightness, direction memory, closure, termination at singularities, curve length);
- Frenet torsion without 4th derivatives;
- independent verification of the analytic κ_g.

## Validation (all pass; see `tpms_validate.m`)

Check 10 (added in v0.3) verifies that the electroplating numbers match their formulas, that Wa equals the target at i_max and Da⁻¹ equals the target at v_min, and that scaling the cell size by 2 with a_v halved halves Wa_geom and leaves Da_geom unchanged.


1. **Derivatives.** Symbolic and hand-derived derivatives agree exactly.
2. **Curvatures.** Eigen-solve K and H agree with the closed-form level-set formulas to about 1e-15.
3. **Curvature derivatives.** The Hellmann–Feynman derivatives of k1 and k2 agree with finite differences on the surface to about 1e-9 relative.
4. **Geodesic curvature.** Analytic κ_g (Codazzi/Liouville) matches finite differences of the traced polylines, with observed convergence order 2.00, 2.00, 2.00 and 1.99 for the four families. Both `hmax` and `maxTurnDeg` must be refined, because on lines of curvature the turn limit sets the step.
5. **Redundancy identities.** These hold to about 1e-15.
6. **Surface drift.** Curves stay on the surface to about 1e-15 (up to 5e-14 at t = 0.8).
7. **Area.** Area per cell converges with the sampling grid: 3.04·a² at t = 0.3.
8. **Cubic symmetry.** The area-weighted orientation tensor of d1 equals I/3 to 5e-5.
9. **Inversion symmetry.** Under x → −x, F = t maps to F = −t with k1 ↔ −k2, k1-lines ↔ k2-lines and asym+ ↔ asym−. The statistics match to about 1e-3.
   - Practical consequence for **sheet gyroids**: the two faces F = +t and F = −t are congruent with the families swapped, so only one face needs to be analysed.

## First results (default run, t = 0 … 0.8)

- **Chirality shows up.** At t = 0 the median signed torsion is −0.096 for k1-lines and +0.101 for k2-lines, which is equal and opposite as inversion symmetry requires. As t grows, both shift positive (k1: +0.20, k2: +0.25 at t = 0.8), because one labyrinth, which is chiral, thickens. The handedness index summarises this: 0.002 → 0.17 across the sweep. The *mean* signed torsion is about 0 by construction (it telescopes), so the effect lives in the shape of the distribution, not its mean. Point-based k1/k2 statistics are mirror-invariant and cannot see this.
- **The asymptotic families split with t.** Median κ_g is identical at t = 0 (0.170). By t = 0.8 it is 0.140 for asym+ and 0.255 for asym−.
- **Topology of the curve network changes.** k2-lines go from never closing (t ≤ 0.2) to all closing (t = 0.8, loops of length about 13.6 rad). At t = 0.6, asym+ curves become lattice-periodic along ⟨011⟩ (length about 9.7 vs. 8.89 for the straight diagonal). At high t, k1-lines wrap around the thickening struts without quite closing: straightness over one cell falls from 0.44 to 0.05, and their direction memory falls from 1.35 to 0.87 mm in a 5 mm cell.
- **|∇_s H|** rises from 0.047 to 0.122 across the sweep: t = 0 is closest to minimal, where H ≈ 0 everywhere.

## Limitations and open items

- **Test environment.**
  - The code was tested end to end in GNU Octave 8.4 with the symbolic package 3.2.1 (SymPy backend), including the `matlabFunction` path. The syntax is kept MATLAB-compatible, but it has not been run in MATLAB itself.
  - `tpms_plot_results.m` could not be executed here because Octave's font renderer was broken in the test container. It is wrapped in try/catch so a plotting error cannot stop a run.
  - The preview image was rendered separately from the exported `.mat` file.
- **Curve descriptors depend on tracing settings.** Curve-level descriptors (torsion, straightness, termination and closure fractions) depend on `nSeeds`, `Lmax` and the windows, and they are length-weighted rather than area-weighted. Keep these settings fixed across any dataset used for PCA. `run_sampling_study.m` quantifies the seed dependence.
- **Possible next improvements** (not built yet):
  - *Evenly spaced curve seeding* (Jobard–Lefer style), so curves cover the surface uniformly. This would remove the curve-density bias in curve-sampled statistics.
  - *A field version of Frenet torsion* (4th derivatives, area-weighted). Twist pitch would then need no seeds at all.
  - *Surface maps of the summary metrics* (e.g. where on the cell the curves are straightest or most twisted).
  - *A node-vs-strut breakdown* of every metric.
  - *Inscribed-sphere thickness*, to match CT-software conventions.
- **Singular neighbourhoods.** Near umbilics (lines of curvature) and parabolic lines (asymptotic curves), κ_g diverges by construction. Nodes within 10% of the median k1−k2, or 5% of the median |K|, are excluded from statistics. The excluded fraction (1–5%) is written to `…_checks.csv`.
- **Closure detection is tolerance-based.** A curve counts as closed if it returns within 0.05 rad (modulo the lattice) with the same tangent. The closure gap is reported per curve.
- **Signed torsion statistics.** CV is not meaningful for signed torsion when its mean is near zero.
- **PCA.** A gyroid-only t-sweep gives few PCA rows. The intended use is to join `…_PCA.csv` to the STL tool's per-STL table on matching t.
- **STL comparison not yet done.** `lattice_curvature_heatmap_v10.m` was not available, so the analytic-vs-STL comparison has not been run. The `ref_*` columns use the same statistic set (mean, std, CV, p5, median, p95, RMS, area-weighted) so they can be compared directly once the STL output is available.
- **Scope.** The code only builds the gyroid. The symbolic path takes any expression, so adding other TPMS types is mostly adding a `case` in `tpms_build_fields.m` plus re-running validation. Blends would need `w` as an extra symbolic parameter.
