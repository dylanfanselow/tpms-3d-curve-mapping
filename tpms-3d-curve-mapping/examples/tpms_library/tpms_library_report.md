# TPMS library: core curvature metrics

58 lattices at t = 0 (literature mid-surface), unit cell 5 mm, sheet interpretation. Lengths in mm. Generated 2026-11-16 05:55.

## 1. Sanity checks

The genus of the surface in one box comes from the curvature (Gauss-Bonnet). Dividing out the translations that map the surface onto itself while keeping its sides gives the genus per primitive cell, the value quoted in the literature. A mismatch means the equation describes a different surface. *Minimality* is median |H| / median curvedness: near 0 for a good approximation of a minimal surface; entries labelled minimal but above 0.35 are flagged. *Sides congruent* means a translation swaps the two sides (labyrinths). Notes on individual entries are in `src/tpms_library.m`.

| Lattice | Kind | Genus / box | Translations | Genus / primitive cell | Literature | Minimality | Smoothness | Sides congruent | Status |
|---|---|---|---|---|---|---|---|---|---|
| Schwarz Primitive (P) | minimal | 3 | 1 | 3 (matches) | 3 | 0.106 | 0.763 | yes | ok |
| Schwarz Diamond (D) | minimal | 9 | 4 | 3 (matches) | 3 | 0.022 | 0.961 | yes | ok |
| Schoen Gyroid (G) | minimal | 5 | 2 | 3 (matches) | 3 | 0.016 | 0.931 | no | ok |
| Neovius (C(P)) | minimal | 9 | 1 | 9 (matches) | 9 | 0.121 | 0.356 | yes | ok |
| Schoen I-WP | minimal | 7 | 2 | 4 (matches) | 4 | 0.229 | 0.770 | no | ok |
| Schoen F-RD (4-term form) | minimal | 21 | 4 | 6 (matches) | 6 | 0.251 | 0.704 | no | ok |
| Schoen F-RD (8-term form) | minimal | 21 | 4 | 6 (matches) | 6 | 0.244 | 0.799 | no | ok |
| Fischer-Koch S | minimal | 21 | 2 | 11 (matches) | 11 | 0.070 | 0.692 | no | ok |
| Lidinoid (nodal approximation) | minimal | 33 | 2 | 17 (fine grid) **MISMATCH** | 3 | 0.599 **far from minimal** | 0.320 | no | ok |
| Split P | minimal | 17 | 2 | 9 | ? | 0.141 | 0.533 | no | ok |
| Split P (corax variant) | minimal | 10 | 1 | 10 | ? | 0.390 **far from minimal** | 0.228 | no | ok |
| Schoen O,C-TO (OCTO) | minimal | 10 | 1 | 10 | ? | 0.498 **far from minimal** | 0.285 | no | ok |
| OCTO (Axolotl form) | minimal | 10 | 1 | 10 | ? | 0.247 | 0.468 | no | ok |
| Complementary D (C(D)) | minimal | 121.47 | 4 | 31.12 (not integer: inconclusive) (fine grid) **MISMATCH** | 19 | 0.241 | 0.282 | yes | ok |
| Complementary S (C(S)) | minimal | 65 | 2 | 33 (fine grid) | ? | 0.174 | 0.289 | no | ok |
| Fischer-Koch Y | minimal | 13 | 1 | 13 | ? | 0.361 **far from minimal** | 0.210 | yes | ok |
| Fischer-Koch C(Y) (corax form) | minimal | 21 | 1 | 21 (fine grid) **MISMATCH** | 13 | 0.138 | 0.131 | yes | ok |
| PMY (+/-Y) | minimal | 21 | 1 | 21 | ? | 0.252 | 0.420 | yes | ok |
| Complementary PMY | minimal | 13 | 1 | 13 | ? | 0.074 | 0.868 | yes | ok |
| G Prime (corax form) | minimal | 25 | 2 | 13 (fine grid) | ? | 0.413 **far from minimal** | 0.371 | no | ok |
| G Prime 2 | minimal | 25 | 2 | 13 | ? | 0.075 | 0.848 | no | ok |
| D Prime (corax form) | minimal | 25.24 | 4 | 7 (not integer: inconclusive) (fine grid) | ? | 0.468 **far from minimal** (nearly singular: treat t = 0 metrics with caution) | 0.049 | no | ok |
| Karcher K | minimal | 12 | 1 | 12 | ? | 0.184 | 0.587 | no | ok |
| I2-Y** | minimal | 9 | 2 | 5 | ? | 0.359 **far from minimal** | 0.715 | no | ok |
| Complementary I2-Y** | minimal | 25 | 2 | 13 | ? | 0.256 | 0.522 | no | ok |
| Complementary G (C(G)) | minimal | 29 | 2 | 15 | ? | 0.316 | 0.744 | no | ok |
| Slotted P | minimal | 11 | 1 | 11 | ? | 0.178 | 0.612 | no | ok |
| Fischer-Koch (I-WP-like form) | minimal | 7 | 2 | 4 | ? | 0.218 | 0.594 | no | ok |
| W surface | minimal | n/a | 1 | n/a | ? | 0.668 | 0.000 | | singular point on surface; excluded from rankings |
| F surface | minimal | n/a | 4 | n/a | ? | 0.852 | 0.001 | | singular point on surface; excluded from rankings |
| P + C(P) (corax form) | non-minimal | 12 | 1 | 12 | ? | 0.288 | 0.625 | no | ok |
| Double Primitive | non-minimal | 5 | 2 | 3 | ? | 0.501 | 0.391 | no | ok |
| Double Primitive (Axolotl form) | non-minimal | -1 | 2 | 0 **splits into closed pieces; excluded from rankings** | ? | 0.695 | 0.893 | no | ok |
| Double Gyroid | non-minimal | 9 | 2 | 5 | ? | 0.381 | 0.834 | no | ok |
| Double Diamond 1 | non-minimal | 3 | 1 | 3 | ? | 0.373 | 0.828 | no | ok |
| Tubular Primitive | non-minimal | 3 | 1 | 3 | ? | 0.671 | 0.368 | no | ok |
| Tubular Gyroid | non-minimal | 5 | 2 | 3 | ? | 0.622 | 0.811 | no | ok |
| Tubular Diamond | non-minimal | -7 | 4 | -1 **splits into closed pieces; excluded from rankings** | ? | 0.952 | 0.811 | no | ok |
| Bionic Bone 1 | non-minimal | 5 | 2 | 3 | ? | 0.089 | 0.927 | no | ok |
| Fischer-Koch C(Y) | minimal | 13 | 1 | 13 (matches) | 13 | 0.137 | 0.587 | yes | ok |
| G Prime (G') | minimal | 25 | 2 | 13 | ? | 0.042 | 0.531 | no | ok |
| D Prime (D') | minimal | 33.10 | 4 | 9 (not integer: inconclusive) (fine grid) | ? | 0.333 | 0.443 | no | ok |
| P + C(P) | non-minimal | 12 | 1 | 12 | ? | 0.202 | 0.721 | no | ok |
| Q* surface | minimal | n/a | 1 | n/a | ? | 0.173 | 0.000 | | singular point on surface; excluded from rankings |
| S3 (Vasile et al.) | non-minimal | 5 | 2 | 3 | ? | 0.401 | 0.138 | no | ok |
| S4 (Vasile et al.) | non-minimal | 5 | 2 | 3 | ? | 0.400 | 0.141 | no | ok |
| S7 (Vasile et al.) | non-minimal | 2 | 1 | 2 | ? | 0.500 (nearly singular: treat t = 0 metrics with caution) | 0.017 | no | ok |
| S8 (Vasile et al.) | non-minimal | 10.47 | 2 | 5.74 (not integer: inconclusive) (fine grid) | ? | 0.457 (nearly singular: treat t = 0 metrics with caution) | 0.033 | yes | ok |
| S9 (Vasile et al.) | non-minimal | 5 | 2 | 3 | ? | 0.307 | 0.741 | yes | ok |
| C(P) (Schwarz & Gompper nodal form) | minimal | 9 | 1 | 9 (matches) | 9 | 0.262 | 0.669 | yes | ok |
| Double Gyroid (MiniSurf, offset -0.95) | non-minimal | 9 | 2 | 5 | ? | 0.457 | 0.856 | no | ok |
| FRP (Axolotl) | minimal | 17 | 4 | 5 | ? | 0.324 | 0.335 | no | ok |
| KP (Axolotl) | minimal | 6 | 1 | 6 | ? | 0.701 **far from minimal** | 0.272 | no | ok |
| Bionic Bone 2 | non-minimal | 5 | 2 | 3 | ? | 0.501 | 0.787 | no | ok |
| Lonsdaleite (hexagonal) * | minimal | 9 | 2 | 5 | ? | 0.161 | 0.647 | no | ok |
| Schwarz tP (c/a = 1.5) * | minimal | 3 | 1 | 3 (matches) | 3 | 0.226 | 0.548 | yes | ok |
| Schwarz tD (c/a = 1.5) * | minimal | 9 | 4 | 3 (matches) | 3 | 0.078 | 0.670 | yes | ok |
| Gyroid tG (c/a = 1.5) * | minimal | 5 | 2 | 3 (matches) | 3 | 0.185 | 0.660 | no | ok |

13 of 16 entries with a literature genus match it.

### Duplicate check

Every pair of equations is tested for describing the same lattice (same surface up to rotation, reflection, translation, scale or doubled frequency; `tpms_library_duplicates`).

- `schwarz_d` and `schwarz_d_alt` are the same lattice (frequency x1, axes [3 2 1], signs [1 1 1], shift [0.25 0.25 0.25]*pi, scale -1.41).
- `double_d_1` and `double_d_2` are the same lattice (frequency x1, axes [3 2 1], signs [1 1 1], shift [0.5 0.5 0.5]*pi, scale 1).

Not analysed because they duplicate an earlier entry: Schwarz Diamond (D), alternative form (= schwarz_d); Double Diamond 2 (= double_d_1).

## 2. Core metrics

Sorted from most to least extreme: mean |z| over curvedness-change length, turning radius and the two straightnesses. Brackets: 95% CI from resampling curves. The mean-curvature-based curvature-change length (the core metric) is shown in grey italics because at t = 0 it mostly reflects the nodal-approximation error of minimal surfaces.

| Lattice | Curvedness-change length (mm) | Turning radius (mm) | k1 straightness | k2 straightness | Handedness | *Curvature-change length (mm)* | Bending radius (mm) | Turning / bending | Extremeness |
|---|---|---|---|---|---|---|---|---|---|

\* = non-cubic cell (tetragonal or hexagonal; `cell` column in the CSV, add-on `tpms_library_noncubic`). Their straightness window is one length a along x, and the cubic-symmetry check does not apply.

| Tubular Primitive | 1.38 | 13.8 | 0.0648 [0.055, 0.0858] | 0.828 [0.827, 0.83] | -0.00679 [-0.0171, 0.00212] | *1.84* | 0.383 | 35.9 | 2.88 |
| Double Diamond 1 | 1.55 | 5.63 | 0.342 [0.253, 0.384] | 0.836 [0.834, 0.837] | 0.00351 [-0.00424, 0.0115] | *5.9* | 1.45 | 3.88 | 1.81 |
| Gyroid tG (c/a = 1.5) * | 1.11 | 4.95 | 0.563 [0.525, 0.569] | 0.557 [0.53, 0.563] | 0.00333 [-0.0134, 0.019] | *3.37* | 1.5 | 3.29 | 1.41 |
| Schwarz tP (c/a = 1.5) * | 0.779 | 4.96 | 0.652 [0.625, 0.695] | 0.653 [0.626, 0.693] | 0.000813 [-0.000717, 0.00233] | *1.49* | 2 | 2.48 | 1.37 |
| S7 (Vasile et al.) | 0.37 | 3.15 | 0.675 [0.645, 0.706] | 0.66 [0.642, 0.687] | -0.00134 [-0.0074, 0.00374] | *0.494* | 1.56 | 2.02 | 1.27 |
| Lonsdaleite (hexagonal) * | 1.09 | 4.1 | 0.503 [0.489, 0.519] | 0.533 [0.525, 0.544] | 0.00307 [-0.00496, 0.00956] | *1.59* | 1.46 | 2.82 | 1.15 |
| Schoen Gyroid (G) | 1.14 | 4.6 | 0.431 [0.426, 0.441] | 0.435 [0.428, 0.44] | 0.00429 [-0.0157, 0.023] | *11.3* | 1.18 | 3.9 | 1.05 |
| Schwarz Primitive (P) | 0.792 | 4.48 | 0.551 [0.501, 0.575] | 0.544 [0.512, 0.579] | -0.00314 [-0.00991, 0.000894] | *1.73* | 1.61 | 2.79 | 1.05 |
| Double Gyroid | 0.836 | 4.05 | 0.0636 [0.0617, 0.0651] | 0.656 [0.655, 0.656] | -0.00409 [-0.038, 0.0287] | *1.59* | 0.939 | 4.31 | 0.989 |
| I2-Y** | 0.597 | 3.13 | 0.652 [0.652, 0.653] | 0.0708 [0.0686, 0.0728] | 0.0017 [-0.0246, 0.0278] | *0.878* | 0.984 | 3.18 | 0.981 |
| Split P (corax variant) | 1.13 | 3.12 | 0.144 [0.11, 0.17] | 0.635 [0.626, 0.661] | 0.00892 [-0.0148, 0.0357] | *3.57* | 1.02 | 3.05 | 0.967 |
| Tubular Gyroid | 0.391 | 2.33 | 0.0895 [0.0869, 0.0918] | 0.758 [0.758, 0.759] | -0.162 [-0.169, -0.154] | *0.644* | 0.396 | 5.88 | 0.961 |
| Bionic Bone 1 | 1.07 | 4.57 | 0.358 [0.353, 0.367] | 0.506 [0.501, 0.511] | 0.0807 [0.0591, 0.0956] | *4.12* | 1.17 | 3.91 | 0.947 |
| Complementary D (C(D)) | 0.104 | 0.696 | 0.209 [0.2, 0.275] | 0.202 [0.193, 0.268] | -0.00197 [-0.00581, 0.00175] | *0.237* | 1.04 | 0.668 | 0.922 |
| KP (Axolotl) | 0.343 | 1.12 | 0.199 [0.108, 0.289] | 0.0223 [0.0141, 0.0803] | 0.000328 [-0.00706, 0.00962] | *0.384* | 0.621 | 1.8 | 0.875 |
| P + C(P) (corax form) | 0.364 | 1.24 | 0.0897 [0.0502, 0.12] | 0.201 [0.17, 0.362] | -0.00157 [-0.0102, 0.00453] | *0.416* | 0.741 | 1.68 | 0.822 |
| Double Gyroid (MiniSurf, offset -0.95) | 0.818 | 3.62 | 0.163 [0.159, 0.174] | 0.69 [0.69, 0.69] | 0.00713 [-0.0179, 0.0324] | *1.61* | 0.823 | 4.4 | 0.8 |
| C(P) (Schwarz & Gompper nodal form) | 0.399 | 2.25 | 0.104 [0.08, 0.133] | 0.0823 [0.0698, 0.126] | 0.00985 [-0.00116, 0.0251] | *0.722* | 1.06 | 2.13 | 0.767 |
| Double Primitive | 0.724 | 3.14 | 0.158 [0.135, 0.181] | 0.778 [0.774, 0.781] | -0.000965 [-0.00284, 0.00124] | *1.46* | 1.17 | 2.69 | 0.761 |
| Schoen F-RD (4-term form) | 0.355 | 1.54 | 0.15 [0.094, 0.191] | 0.159 [0.139, 0.19] | -0.00432 [-0.00854, -0.000563] | *0.773* | 0.812 | 1.9 | 0.741 |
| Lidinoid (nodal approximation) | 0.356 | 1.12 | 0.31 [0.254, 0.315] | 0.139 [0.134, 0.153] | 0.00903 [-0.0746, 0.0809] | *0.659* | 1.19 | 0.937 | 0.73 |
| OCTO (Axolotl form) | 0.259 | 1.52 | 0.197 [0.134, 0.243] | 0.19 [0.15, 0.214] | 0.00221 [-0.00164, 0.00726] | *0.561* | 1.55 | 0.976 | 0.718 |
| S4 (Vasile et al.) | 0.556 | 1.91 | 0.519 [0.488, 0.555] | 0.461 [0.437, 0.486] | 0.0207 [-0.0119, 0.0536] | *0.718* | 1.35 | 1.41 | 0.706 |
| Complementary S (C(S)) | 0.316 | 1.15 | 0.202 [0.181, 0.214] | 0.204 [0.196, 0.216] | 0.0015 [-0.00657, 0.00876] | *0.808* | 0.575 | 2 | 0.694 |
| Complementary PMY | 0.88 | 2.75 | 0.178 [0.042, 0.305] | 0.0371 [0.0331, 0.3] | 0.00245 [-0.00455, 0.00899] | *2.81* | 0.843 | 3.26 | 0.662 |
| D Prime (D') | 0.363 | 1.39 | 0.237 [0.224, 0.237] | 0.123 [0.0984, 0.153] | 0.00116 [-0.00295, 0.0049] | *0.831* | 0.804 | 1.72 | 0.657 |
| S9 (Vasile et al.) | 0.812 | 2.24 | 0.461 [0.436, 0.477] | 0.446 [0.432, 0.471] | -0.00141 [-0.00267, 2.5e-05] | *1.38* | 1.43 | 1.57 | 0.655 |
| Bionic Bone 2 | 0.478 | 3.01 | 0.198 [0.198, 0.199] | 0.707 [0.707, 0.708] | 0.0916 [0.0806, 0.104] | *0.717* | 0.774 | 3.89 | 0.645 |
| D Prime (corax form) | 0.361 | 1.36 | 0.203 [0.187, 0.217] | 0.191 [0.159, 0.217] | 0.0279 [-0.018, 0.0744] | *0.54* | 0.862 | 1.57 | 0.642 |
| Schwarz Diamond (D) | 0.874 | 3.76 | 0.17 [0.142, 0.21] | 0.187 [0.139, 0.207] | -0.0021 [-0.00422, -0.000144] | *4.64* | 0.905 | 4.16 | 0.64 |
| Fischer-Koch S | 0.518 | 1.72 | 0.132 [0.13, 0.133] | 0.132 [0.128, 0.133] | 0.00685 [-0.0251, 0.041] | *1.44* | 0.817 | 2.1 | 0.64 |
| S3 (Vasile et al.) | 0.753 | 2.7 | 0.401 [0.366, 0.475] | 0.602 [0.578, 0.61] | 0.1 [0.0691, 0.127] | *1.38* | 1.83 | 1.48 | 0.62 |
| Fischer-Koch C(Y) | 0.841 | 2.65 | 0.14 [0.0988, 0.177] | 0.134 [0.0922, 0.216] | 0.0911 [0.0751, 0.112] | *1.84* | 0.937 | 2.83 | 0.594 |
| Schoen F-RD (8-term form) | 0.432 | 1.69 | 0.187 [0.164, 0.193] | 0.175 [0.0859, 0.197] | 0.00105 [-0.00399, 0.0073] | *0.844* | 0.779 | 2.17 | 0.581 |
| PMY (+/-Y) | 0.263 | 1.26 | 0.282 [0.271, 0.33] | 0.294 [0.276, 0.341] | 0.0113 [-0.0051, 0.0254] | *0.62* | 1.23 | 1.02 | 0.58 |
| Complementary I2-Y** | 0.463 | 1.64 | 0.168 [0.162, 0.19] | 0.189 [0.188, 0.193] | 0.0101 [-0.000424, 0.0226] | *1.19* | 0.814 | 2.01 | 0.577 |
| Neovius (C(P)) | 0.51 | 2.15 | 0.145 [0.119, 0.178] | 0.137 [0.113, 0.176] | 0.00211 [-0.00212, 0.00923] | *1.26* | 1.33 | 1.61 | 0.565 |
| Fischer-Koch C(Y) (corax form) | 0.524 | 1.82 | 0.153 [0.141, 0.164] | 0.157 [0.134, 0.171] | 0.134 [0.124, 0.144] | *1.39* | 1.01 | 1.79 | 0.563 |
| Karcher K | 0.451 | 2.08 | 0.211 [0.148, 0.403] | 0.103 [0.0748, 0.128] | 0.00176 [-0.00179, 0.00711] | *0.927* | 0.996 | 2.08 | 0.556 |
| Schoen O,C-TO (OCTO) | 0.473 | 1.64 | 0.168 [0.143, 0.216] | 0.457 [0.0541, 0.484] | 7.69e-05 [-0.00348, 0.00303] | *1.12* | 1.83 | 0.897 | 0.536 |
| S8 (Vasile et al.) | 0.359 | 1.66 | 0.362 [0.331, 0.402] | 0.348 [0.293, 0.385] | -0.0062 [-0.0148, 0.00283] | *0.62* | 1.13 | 1.47 | 0.534 |
| G Prime (corax form) | 0.421 | 1.56 | 0.197 [0.0561, 0.279] | 0.235 [0.197, 0.239] | -0.00208 [-0.0363, 0.0292] | *0.74* | 0.995 | 1.56 | 0.529 |
| G Prime 2 | 0.507 | 2.23 | 0.161 [0.144, 0.193] | 0.166 [0.156, 0.175] | -0.00167 [-0.0137, 0.00996] | *1.89* | 0.688 | 3.23 | 0.501 |
| FRP (Axolotl) | 0.431 | 1.57 | 0.252 [0.216, 0.325] | 0.208 [0.194, 0.221] | 0.112 [0.0747, 0.143] | *0.78* | 1.02 | 1.54 | 0.461 |
| Slotted P | 0.61 | 1.9 | 0.227 [0.215, 0.238] | 0.0695 [0.0483, 0.133] | -0.00467 [-0.0161, 0.00936] | *1.39* | 1.06 | 1.8 | 0.456 |
| G Prime (G') | 0.587 | 2.34 | 0.0813 [0.024, 0.127] | 0.25 [0.248, 0.252] | 0.00567 [-0.0324, 0.0458] | *2.88* | 0.804 | 2.91 | 0.455 |
| Fischer-Koch (I-WP-like form) | 0.496 | 2.4 | 0.271 [0.267, 0.301] | 0.108 [0.0306, 0.17] | -0.000787 [-0.00315, 0.00127] | *0.88* | 1.01 | 2.37 | 0.42 |
| Fischer-Koch Y | 0.775 | 3.51 | 0.217 [0.192, 0.231] | 0.218 [0.195, 0.228] | 0.0589 [0.0396, 0.0782] | *0.99* | 1.04 | 3.37 | 0.419 |
| Schwarz tD (c/a = 1.5) * | 0.802 | 3.8 | 0.334 [0.319, 0.339] | 0.331 [0.32, 0.336] | -0.00202 [-0.00701, 0.000987] | *3.06* | 1.2 | 3.16 | 0.413 |
| Complementary G (C(G)) | 0.448 | 1.74 | 0.243 [0.145, 0.28] | 0.241 [0.163, 0.269] | 0.0032 [-0.00648, 0.0105] | *0.723* | 0.579 | 3 | 0.405 |
| Split P | 0.575 | 2.48 | 0.0844 [0.06, 0.112] | 0.323 [0.306, 0.331] | -0.00355 [-0.0495, 0.0351] | *1.72* | 0.936 | 2.65 | 0.364 |
| P + C(P) | 0.631 | 1.57 | 0.185 [0.135, 0.205] | 0.325 [0.218, 0.36] | 0.000841 [-0.00228, 0.00346] | *1.12* | 0.852 | 1.84 | 0.283 |
| Schoen I-WP | 0.659 | 2.32 | 0.211 [0.181, 0.251] | 0.252 [0.209, 0.305] | 0.00115 [-0.00107, 0.00354] | *1.26* | 1.15 | 2.02 | 0.248 |
| W surface | 0.0514 | n/a | n/a | n/a | n/a | *n/a* | n/a | n/a | n/a |
| F surface | 0.13 | n/a | n/a | n/a | n/a | *n/a* | n/a | n/a | n/a |
| Double Primitive (Axolotl form) | 0.541 | 2.21 | n/a | 0.332 [0.19, 0.377] | -0.00151 [-0.00644, 0.0044] | *0.509* | 1.31 | 1.69 | n/a |
| Tubular Diamond | 0.469 | 0.702 | 0.118 [0.108, 0.122] | 0.156 [0.105, 0.191] | 0.000305 [-0.00362, 0.00401] | *0.506* | 0.528 | 1.33 | n/a |
| Q* surface | 1.02 | n/a | n/a | n/a | n/a | *n/a* | n/a | n/a | n/a |

## 3. Who is most extreme?

- **Curvedness-change length**: highest Double Diamond 1 (1.55), Tubular Primitive (1.38), Schoen Gyroid (G) (1.14); lowest Complementary D (C(D)) (0.104), OCTO (Axolotl form) (0.259), PMY (+/-Y) (0.263).
- **Turning radius**: highest Tubular Primitive (13.8), Double Diamond 1 (5.63), Schwarz tP (c/a = 1.5) (4.96); lowest Complementary D (C(D)) (0.696), Lidinoid (nodal approximation) (1.12), KP (Axolotl) (1.12).
- **k1 straightness**: highest S7 (Vasile et al.) (0.675), I2-Y** (0.652), Schwarz tP (c/a = 1.5) (0.652); lowest Double Gyroid (0.0636), Tubular Primitive (0.0648), G Prime (G') (0.0813).
- **k2 straightness**: highest Double Diamond 1 (0.836), Tubular Primitive (0.828), Double Primitive (0.778); lowest KP (Axolotl) (0.0223), Complementary PMY (0.0371), Slotted P (0.0695).
- **Chirality** (|handedness| >= 0.02 and its CI excludes 0): Fischer-Koch Y (0.0589), Fischer-Koch C(Y) (corax form) (0.134), Tubular Gyroid (-0.162), Bionic Bone 1 (0.0807), Fischer-Koch C(Y) (0.0911), S3 (Vasile et al.) (0.1), FRP (Axolotl) (0.112), Bionic Bone 2 (0.0916). A non-zero value at t = 0 means the structure itself is single-handed.
- **Sharpest steering relative to its own bending** (turning / bending radius, low = sharp): Complementary D (C(D)) (0.668), Schoen O,C-TO (OCTO) (0.897), Lidinoid (nodal approximation) (0.937). Gentlest: Tubular Primitive (35.9), Tubular Gyroid (5.88), Double Gyroid (MiniSurf, offset -0.95) (4.4).
- **Most extreme overall**: Tubular Primitive (2.88), Double Diamond 1 (1.81), Gyroid tG (c/a = 1.5) (1.41).

## 4. Electroplating at matched solid fraction (ranking only)

At t = 0 a sheet has zero wall thickness, so the lattices are compared at the same solid fraction, **16%** (porosity 0.84, as the reference strut lattices of Clemens et al. 2025), sheet interpretation. Lattice Wagner number Wa at 0.35 mA/cm2 (higher = more even coating) and inverse Damkohler number Da^-1 at a reference flow of 0.05 m/s (higher = less ion depletion), N_field = 4, N_flow = 8, placeholder bath values from the paper. The equations are reconstructed from the paper (absolute values may be off by ~1.5x): use them to rank lattices.

**At equal porosity both numbers scale as 1/(wetted surface per volume)**, so they rank the lattices in the same order: lattices with less surface per volume plate more evenly *and* need less flow. Their absolute values still answer different questions (uniformity vs. depletion).

| Lattice | t | Wall thickness (mm) | Wetted area a_v (1/m) | Wa | Da^-1 at 0.05 m/s | Max current for Wa >= 1 (mA/cm2) | Flow for Da^-1 = 440 (m/s) |
|---|---|---|---|---|---|---|---|
| Tubular Primitive | 2.004 | n/a | 354.5 | 0.573 | 1943.9 | 0.201 | 0.0113 |
| Tubular Gyroid | 3.594 | n/a | 472.5 | 0.43 | 1458.5 | 0.151 | 0.0151 |
| Bionic Bone 2 | 2.584 | 0.383 | 828.3 | 0.245 | 832.0 | 0.0859 | 0.0264 |
| Schwarz tP (c/a = 1.5) | 0.280 | 0.38 | 836.3 | 0.243 | 824.1 | 0.0851 | 0.0267 |
| Schwarz Primitive (P) | 0.282 | 0.344 | 929.0 | 0.219 | 741.9 | 0.0766 | 0.0297 |
| Lonsdaleite (hexagonal) | 0.124 | 0.293 | 1086.4 | 0.187 | 634.4 | 0.0655 | 0.0347 |
| Gyroid tG (c/a = 1.5) | 0.248 | 0.289 | 1099.4 | 0.185 | 626.8 | 0.0647 | 0.0351 |
| Bionic Bone 1 | 4.955 | 0.26 | 1213.9 | 0.167 | 567.7 | 0.0586 | 0.0388 |
| Schoen Gyroid (G) | 0.248 | 0.26 | 1222.9 | 0.166 | 563.5 | 0.0582 | 0.039 |
| Double Diamond 1 | 0.266 | 0.255 | 1257.7 | 0.162 | 548.0 | 0.0566 | 0.0401 |
| S3 (Vasile et al.) | 0.365 | 0.251 | 1269.2 | 0.16 | 543.0 | 0.0561 | 0.0405 |
| S9 (Vasile et al.) | 0.609 | 0.241 | 1325.8 | 0.153 | 519.8 | 0.0537 | 0.0423 |
| Double Primitive | 0.099 | 0.24 | 1352.0 | 0.15 | 509.8 | 0.0526 | 0.0432 |
| Schwarz tD (c/a = 1.5) | 0.196 | 0.233 | 1364.1 | 0.149 | 505.2 | 0.0522 | 0.0435 |
| Neovius (C(P)) | 0.362 | 0.23 | 1367.0 | 0.149 | 504.2 | 0.0521 | 0.0436 |
| Schoen I-WP | 0.604 | 0.225 | 1405.1 | 0.145 | 490.5 | 0.0506 | 0.0449 |
| S4 (Vasile et al.) | 0.321 | 0.223 | 1408.0 | 0.144 | 489.5 | 0.0505 | 0.0449 |
| Schoen O,C-TO (OCTO) | 0.074 | 0.214 | 1426.6 | 0.143 | 483.1 | 0.0499 | 0.0455 |
| OCTO (Axolotl form) | 0.557 | 0.224 | 1436.1 | 0.142 | 479.9 | 0.0495 | 0.0458 |
| P + C(P) (corax form) | 0.047 | 0.218 | 1463.6 | 0.139 | 470.9 | 0.0486 | 0.0467 |
| P + C(P) | 0.067 | 0.214 | 1470.2 | 0.138 | 468.8 | 0.0484 | 0.0469 |
| C(P) (Schwarz & Gompper nodal form) | 0.233 | 0.214 | 1478.4 | 0.138 | 466.2 | 0.0481 | 0.0472 |
| Schwarz Diamond (D) | 0.195 | 0.209 | 1517.8 | 0.134 | 454.1 | 0.0469 | 0.0485 |
| Fischer-Koch (I-WP-like form) | 0.407 | 0.201 | 1539.2 | 0.132 | 447.7 | 0.0462 | 0.0491 |
| Karcher K | 0.153 | 0.205 | 1551.2 | 0.131 | 444.3 | 0.0459 | 0.0495 |
| S7 (Vasile et al.) | 0.270 | 0.193 | 1582.2 | 0.128 | 435.6 | 0.045 | 0.0505 |
| Slotted P | 0.799 | 0.182 | 1740.3 | 0.117 | 396.0 | 0.0409 | 0.0556 |
| Complementary PMY | 0.320 | 0.178 | 1777.8 | 0.114 | 387.7 | 0.04 | 0.0568 |
| Fischer-Koch C(Y) | 0.362 | 0.171 | 1854.2 | 0.11 | 371.7 | 0.0384 | 0.0592 |
| S8 (Vasile et al.) | 0.838 | 0.167 | 1856.5 | 0.11 | 371.2 | 0.0383 | 0.0593 |
| Double Gyroid (MiniSurf, offset -0.95) | 0.671 | 0.169 | 1872.5 | 0.109 | 368.0 | 0.038 | 0.0598 |
| KP (Axolotl) | 0.129 | 0.166 | 1888.7 | 0.108 | 364.9 | 0.0377 | 0.0603 |
| Schoen F-RD (8-term form) | 0.608 | 0.165 | 1903.4 | 0.107 | 362.1 | 0.0374 | 0.0608 |
| Schoen F-RD (4-term form) | 0.369 | 0.164 | 1912.2 | 0.106 | 360.4 | 0.0372 | 0.061 |
| FRP (Axolotl) | 0.531 | 0.165 | 1934.9 | 0.105 | 356.2 | 0.0368 | 0.0618 |
| Split P (corax variant) | 0.194 | 0.157 | 1950.1 | 0.104 | 353.4 | 0.0365 | 0.0623 |
| Fischer-Koch C(Y) (corax form) | 0.258 | 0.15 | 1969.9 | 0.103 | 349.9 | 0.0361 | 0.0629 |
| PMY (+/-Y) | 0.189 | 0.156 | 2015.3 | 0.101 | 342.0 | 0.0353 | 0.0643 |
| Split P | 0.221 | 0.157 | 2022.3 | 0.101 | 340.8 | 0.0352 | 0.0646 |
| Double Gyroid | 0.581 | 0.152 | 2100.4 | 0.0968 | 328.1 | 0.0339 | 0.067 |
| Fischer-Koch S | 0.162 | 0.147 | 2148.9 | 0.0946 | 320.7 | 0.0331 | 0.0686 |
| I2-Y** | 0.438 | 0.148 | 2151.6 | 0.0945 | 320.3 | 0.0331 | 0.0687 |
| G Prime (G') | 0.136 | 0.145 | 2174.1 | 0.0935 | 317.0 | 0.0327 | 0.0694 |
| G Prime 2 | 0.846 | 0.145 | 2197.3 | 0.0925 | 313.6 | 0.0324 | 0.0701 |
| Fischer-Koch Y | 0.155 | 0.141 | 2218.0 | 0.0917 | 310.7 | 0.0321 | 0.0708 |
| D Prime (corax form) | 0.102 | 0.144 | 2236.0 | 0.0909 | 308.2 | 0.0318 | 0.0714 |
| Complementary I2-Y** | 0.381 | 0.139 | 2278.5 | 0.0892 | 302.5 | 0.0312 | 0.0727 |
| Complementary G (C(G)) | 0.813 | 0.14 | 2286.2 | 0.0889 | 301.5 | 0.0311 | 0.073 |
| G Prime (corax form) | 0.222 | 0.132 | 2434.5 | 0.0835 | 283.1 | 0.0292 | 0.0777 |
| Lidinoid (nodal approximation) | 0.103 | 0.131 | 2443.8 | 0.0832 | 282.0 | 0.0291 | 0.078 |
| D Prime (D') | 0.093 | 0.129 | 2460.6 | 0.0826 | 280.1 | 0.0289 | 0.0785 |
| Complementary S (C(S)) | 0.414 | 0.102 | 3144.4 | 0.0647 | 219.2 | 0.0226 | 0.1 |
| Complementary D (C(D)) | 0.145 | 0.0916 | 3487.6 | 0.0583 | 197.6 | 0.0204 | 0.111 |

- **Most even coating / least flow needed**: Tubular Primitive (0.573), Tubular Gyroid (0.43), Bionic Bone 2 (0.245).
- **Least even coating / most flow needed**: Complementary D (C(D)) (0.0583), Complementary S (C(S)) (0.0647), D Prime (D') (0.0826).

**Reading the results.** All lattices share the same cell size, so surfaces with more surface per cell (higher-frequency equations such as F-RD or Split P) naturally have smaller radii and lower straightness. The turning / bending ratio removes that scale effect. Handedness is ~0 at t = 0 for centrosymmetric surfaces, so at t = 0 it flags structurally chiral lattices (such as the tubular gyroid, which thickens only one of the two gyroid channels); for the others it becomes informative at t ~= 0.

**Sources.** Equations were checked against published collections (Bernhard 2026, arXiv:2604.07160, Tables 6-7, including the "corax" collection; nTop documentation; Gandy et al. 2001; Wohlgemuth et al. 2001). Where sources disagree, each variant is a separate entry. Community-collected formulas are secondary sources: the checks in section 1 show which ones behave like the surfaces they are named after.
