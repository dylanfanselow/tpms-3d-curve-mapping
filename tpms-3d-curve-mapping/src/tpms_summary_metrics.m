function M = tpms_summary_metrics(R, thick, cfg)
%TPMS_SUMMARY_METRICS  A short list of plain-language numbers for one level set.
%
%   M = tpms_summary_metrics(R, thick, cfg)
%     R      output of tpms_analyze_level
%     thick  matching row of tpms_wall_thickness (or [] to skip structure metrics)
%
%   M is a struct array with fields
%     key, label, unit, value, lo, hi, rangeType ('CI95' | 'P5-P95' | ''), meaning, group
%
%   Curve-based numbers carry a 95% bootstrap confidence interval obtained by
%   resampling whole curves (cfg.nBoot resamples). It shows how much the number
%   would move if a different set of curves had been traced: a wide interval
%   means more seeds are needed (see run_sampling_study.m).
%   Lengths are in the unit of cfg.cellSize (e.g. mm); with no cell size they
%   are given as fractions of the cell length.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.

if isempty(cfg.cellSize)
    lam = 1/(2*pi); uL = 'cell lengths';            % report lengths as fraction of cell
else
    lam = R.lambda; uL = cfg.lengthUnit;
end
if isfield(cfg, 'nBoot'), B = cfg.nBoot; else, B = 400; end
if isfield(cfg, 'latticeType'), lt = cfg.latticeType; else, lt = 'sheet'; end
M = struct('key',{},'label',{},'unit',{},'value',{},'lo',{},'hi',{},'rangeType',{},'meaning',{},'group',{});
    function add(key, label, unit, v, lo, hi, rt, meaning, group)
        M(end+1) = struct('key',key,'label',label,'unit',unit,'value',v,'lo',lo,'hi',hi, ...
                          'rangeType',rt,'meaning',meaning,'group',group);
    end
sc = @(name) R.scalars(strcmp({R.scalars.name}, name)).value;

% ---------------- structure (what t does to the part) ----------------------
if ~isempty(thick)
    if isempty(cfg.cellSize), f = 1/(2*pi); else, f = 1; end   % thick already in mm if cellSize set
    if strcmp(lt, 'sheet')
        add('wall_thickness', 'Wall thickness (sheet)', uL, thick.thick_mean*f, thick.thick_p5*f, ...
            thick.thick_p95*f, 'P5-P95', 'Mean wall thickness; the range covers 90% of the wall area.', 'structure');
        add('solid_fraction', 'Solid volume fraction', '%', 100*thick.vf_sheet, NaN, NaN, '', ...
            'Share of the bounding volume filled with material.', 'structure');
        ad = thick.area_density_sheet; sv = thick.SV_sheet;
    else
        add('solid_fraction', 'Solid volume fraction', '%', 100*thick.vf_network, NaN, NaN, '', ...
            'Share of the bounding volume filled with material.', 'structure');
        ad = thick.area_density_network; sv = thick.SV_network;
    end
    if isempty(cfg.cellSize), ad = ad*2*pi; sv = sv*2*pi; uA = '1/cell length';
    else, uA = ['1/' cfg.lengthUnit]; end
    add('surface_per_volume', 'Surface area per bounding volume', uA, ad, NaN, NaN, '', ...
        'Surface area packed into each unit of part volume (higher = more surface to finish, coat or exchange through).', 'structure');
    add('surface_per_solid', 'Surface-to-solid-volume ratio (S/V)', uA, sv, NaN, NaN, '', ...
        'Surface area per unit of solid material (same quantity as the STL tool''s S/V).', 'structure');
end
Cv = sqrt((R.G.k1.^2 + R.G.k2.^2)/2);
sC = tpms_wstats(Cv, R.S.w);
add('surface_bend_radius', 'Typical surface bending radius', uL, lam/sC.median, lam/sC.p95, lam/sC.p5, 'P5-P95', ...
    'Radius of the sphere that bends like the surface does (1/curvedness). Point-based, same information as the STL tool.', 'structure');
dn = sc('areafrac_normal_down'); up = sc('areafrac_normal_up');
if strcmp(lt, 'sheet'), ov = (dn + up)/2; else, ov = dn; end
add('downskin', sprintf('Downward-facing overhang (> %g deg)', cfgval(cfg,'overhangDeg',45)), '%', 100*ov, NaN, NaN, '', ...
    'Share of surface facing downward more steeply than the overhang limit (typically the roughest in AM).', 'structure');

% ---------------- lines of curvature ----------------------------------------
pk1 = getPC(R, 'k1'); pk2 = getPC(R, 'k2');
for fam = {'k1','k2'}
    key = fam{1}; pc = getPC(R, key);
    if isempty(pc), continue; end
    if strcmp(key,'k1'), nice = 'Max-bending lines (k1)'; else, nice = 'Min-bending lines (k2)'; end
    v = [pc.straight1];
    [m, lo, hi] = bootStat(@(i) nmed(v(i)), numel(v), B);
    add([key '_straightness'], [nice ': straightness over 1 cell'], '0-1', m, lo, hi, 'CI95', ...
        '1 = runs straight, 0 = comes back to where it started. Net progress per distance travelled.', 'curves');
    A = R.acc.(key);
    [m, lo, hi] = bootStat(@(i) tpms_persistence(A.lags, sum(A.sum(i,:),1), sum(A.cnt(i,:),1)), size(A.sum,1), B);
    [m, lo, hi, rt] = censor(m, lo, hi, A.lags(end));
    add([key '_memory'], [nice ': direction memory'], uL, m*lam, lo*lam, hi*lam, rt, ...
        'How far a line runs before it has "forgotten" its starting direction.', 'curves');
    cl = double([pc.closed]);
    [m, lo, hi] = bootStat(@(i) 100*mean(cl(i)), numel(cl), B);
    add([key '_loops'], [nice ': share forming closed loops'], '%', m, lo, hi, 'CI95', ...
        'Lines that close on themselves, e.g. rings around struts.', 'curves');
end
pLoC = [pk1, pk2];
if ~isempty(pLoC)
    kgv = [pLoC.median_abs_kg]; tav = [pLoC.median_abs_tau];
    [m, lo, hi] = bootStat(@(i) 1/nmed(kgv(i)), numel(kgv), B);
    add('turning_radius', 'Bending lines: sideways turning radius', uL, m*lam, lo*lam, hi*lam, 'CI95', ...
        'Typical radius of the sideways turns the lines make within the surface (small = sharp steering).', 'curves');
    [m, lo, hi] = bootStat(@(i) 2*pi/nmed(tav(i)), numel(tav), B);
    add('twist_pitch', 'Bending lines: twist pitch', uL, m*lam, lo*lam, hi*lam, 'CI95', ...
        'Distance along a line for one full corkscrew twist (small = strongly twisted).', 'curves');
    % NOTE: the length-weighted MEAN torsion of a line of curvature telescopes
    % (tau = -dphi/ds) and is ~0 for any long curve, so handedness is measured
    % as (length twisting +) - (length twisting -), normalised.
    ws = [pLoC.tau_wsign]; wt = [pLoC.wsum];
    [m, lo, hi] = bootStat(@(i) sum(ws(i))/sum(wt(i)), numel(ws), B);
    add('handedness', 'Handedness index', '-1..+1', m, lo, hi, 'CI95', ...
        ['Share of line length twisting right minus share twisting left. 0 = balanced (mirror-symmetric); ' ...
         '+-1 = all one way. Point-curvature statistics cannot see this.'], 'curves');
end

% ---------------- asymptotic curves ----------------------------------------
pa = [getPC(R,'asymP'), getPC(R,'asymM')];
if ~isempty(pa)
    v = [pa.straight1];
    [m, lo, hi] = bootStat(@(i) nmed(v(i)), numel(v), B);
    add('asym_straightness', 'Zero-bending lines: straightness over 1 cell', '0-1', m, lo, hi, 'CI95', ...
        'Asymptotic curves: paths along which the surface does not bend. Near 1 = almost straight.', 'curves');
    Ap = R.acc.asymP; Am = R.acc.asymM;
    As = [Ap.sum; Am.sum]; Ac = [Ap.cnt; Am.cnt];
    [m, lo, hi] = bootStat(@(i) tpms_persistence(Ap.lags, sum(As(i,:),1), sum(Ac(i,:),1)), size(As,1), B);
    [m, lo, hi, rt] = censor(m, lo, hi, Ap.lags(end));
    add('asym_memory', 'Zero-bending lines: direction memory', uL, m*lam, lo*lam, hi*lam, rt, ...
        'How far a zero-bending line runs before it has "forgotten" its starting direction ("> x" = longer than the longest lag checked).', 'curves');
    dP = R.desc(strcmp({R.desc.name}, 'field_kg_asymP')); dM = R.desc(strcmp({R.desc.name}, 'field_kg_asymM'));
    sP = tpms_wstats(dP.values, dP.weights); sM = tpms_wstats(dM.values, dM.weights);
    add('asym_imbalance', 'Zero-bending families imbalance', '-1..+1', ...
        (sM.median - sP.median)/(sM.median + sP.median), NaN, NaN, '', ...
        'Difference between the two mirror-image families of zero-bending lines; 0 when the part is mirror-balanced.', 'curves');
end

% ---------------- build direction --------------------------------------------
add('bend_dir_vertical', 'Max-bending direction near vertical (< 30 deg)', '%', ...
    100*sc('areafrac_d1_near_build'), NaN, NaN, '', ...
    'Share of surface whose sharpest bend runs close to the build direction (13.4% if directions were random).', 'curves');
end

% -------------------------------------------------------------------------
function pc = getPC(R, key)
if isfield(R.perCurve, key), pc = R.perCurve.(key); else, pc = []; end
if isempty(pc), pc = []; end
end

function [m, lo, hi] = bootStat(fun, n, B)
% statistic on all curves + percentile bootstrap CI over curves
if n == 0, m = NaN; lo = NaN; hi = NaN; return; end
m = fun(1:n);
if n < 3 || B < 1, lo = NaN; hi = NaN; return; end
bs = zeros(B,1);
for b = 1:B, bs(b) = fun(randi(n, n, 1)); end
bs = sort(bs(isfinite(bs)));
if numel(bs) < 0.8*B, lo = NaN; hi = NaN; return; end
lo = bs(max(1, round(0.025*numel(bs)))); hi = bs(max(1, round(0.975*numel(bs))));
end

function [m, lo, hi, rt] = censor(m, lo, hi, maxLag)
% persistence longer than the largest lag -> report as a lower bound
rt = 'CI95';
if ~isfinite(m), m = maxLag; lo = NaN; hi = NaN; rt = 'LOWER'; end
end

function v = cfgval(cfg, f, d)
if isfield(cfg, f), v = cfg.(f); else, v = d; end
end

function m = nmed(x)
x = x(isfinite(x));
if isempty(x), m = NaN; else, m = median(x); end
end
