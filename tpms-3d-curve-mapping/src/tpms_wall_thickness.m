function TT = tpms_wall_thickness(fld, tGrid, cfg)
%TPMS_WALL_THICKNESS  How the level-set value t maps to physical geometry.
%
%   TT = tpms_wall_thickness(fld, tGrid, cfg)    struct array, one entry per t
%
%   SHEET gyroid  (solid = -t < F < t):
%     wall thickness measured exactly along the mid-surface normal: from each
%     point x0 on F = 0, Newton-solve F(x0 + s*n) = +t and = -t; thickness =
%     s(+t) - s(-t). Area-weighted over the mid-surface. The first-order
%     estimate 2t/|gradF| is reported alongside for comparison.
%   NETWORK gyroid (solid = F < t): relative density only (no single wall).
%
%   Fields (lengths in mm when cfg.cellSize is set, else radians; *_cell are
%   fractions of the cell length and are always unitless):
%     t, thick_mean, thick_p5, thick_p95, thick_min, thick_max, thick_cv,
%     thick_firstorder_mean, thick_mean_cell, thick_fail_frac,
%     vf_sheet, vf_network       solid volume fractions (0..1)
%     area_per_cell              area of ONE surface F = t in one cell
%     area_density_sheet         2*area / cell volume     (1/length)
%     area_density_network       area / cell volume
%     SV_sheet, SV_network       surface / SOLID volume   (1/length)
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.

if isempty(cfg.cellSize), lam = 1; else, lam = cfg.cellSize/(2*pi); end
if isfield(cfg, 'vfGrid'), nv = cfg.vfGrid; else, nv = 128; end
L = 2*pi;

% mid-surface for thickness
S0 = tpms_surface_sample(fld, 0, cfg.nGrid);
g0 = fld.grad(S0.v); gn0 = sqrt(sum(g0.^2, 2)); n0 = g0 ./ gn0;
w0 = S0.w;

% voxel-centre grid for volume fractions (F evaluated once)
c  = ((1:nv) - 0.5) * L / nv;
[X, Y, Z] = ndgrid(c, c, c);
Fv = fld.F0([X(:) Y(:) Z(:)]);
Vc = L^3;

TT = struct([]);
for k = 1:numel(tGrid)
    t = tGrid(k);
    r.t = t;
    if t > 0
        [sp, okp] = rayHit(fld, S0.v, n0, gn0,  t);
        [sm, okm] = rayHit(fld, S0.v, n0, gn0, -t);
        ok = okp & okm;
        th = sp - sm;
        st = tpms_wstats(th(ok), w0(ok));
        fo = tpms_wstats(2*t ./ gn0, w0);
        r.thick_mean = st.mean*lam;  r.thick_p5 = st.p5*lam;  r.thick_p95 = st.p95*lam;
        r.thick_min = min(th(ok))*lam; r.thick_max = max(th(ok))*lam; r.thick_cv = st.CV;
        r.thick_firstorder_mean = fo.mean*lam;
        r.thick_mean_cell = st.mean/L;
        r.thick_fail_frac = sum(w0(~ok))/sum(w0);
    else
        r.thick_mean = 0; r.thick_p5 = 0; r.thick_p95 = 0; r.thick_min = 0; r.thick_max = 0;
        r.thick_cv = NaN; r.thick_firstorder_mean = 0; r.thick_mean_cell = 0; r.thick_fail_frac = 0;
    end
    r.vf_sheet   = mean(abs(Fv) < t);
    r.vf_network = mean(Fv < t);
    try
        A = tpms_surface_sample(fld, t, cfg.nGrid).A;
    catch
        A = 0;                                  % no surface at this t
    end
    r.area_per_cell        = A*lam^2;
    r.area_density_sheet   = 2*A/Vc/lam;
    r.area_density_network = A/Vc/lam;
    r.SV_sheet   = ifpos(r.vf_sheet,   2*A/(r.vf_sheet*Vc)/lam);
    r.SV_network = ifpos(r.vf_network, A/(r.vf_network*Vc)/lam);
    if isempty(TT), TT = r; else, TT(end+1) = r; end %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function [s, ok] = rayHit(fld, x0, n, gn, target)
% Newton along the fixed line x0 + s*n for F = target, from the first-order guess
s = target ./ gn;
for it = 1:40
    x  = x0 + s.*n;
    f  = fld.F0(x) - target;
    df = sum(fld.grad(x).*n, 2);
    step = f ./ df;
    step(~isfinite(step)) = 0;
    s = s - max(min(step, 0.5), -0.5);         % damped
    if max(abs(f)) < 1e-11, break; end
end
f  = fld.F0(x0 + s.*n) - target;
ok = abs(f) < 1e-9 & sign(s) == sign(target) & abs(s) < 3*abs(target)./gn + 0.5;
end

function v = ifpos(d, val)
if d > 0, v = val; else, v = NaN; end
end
