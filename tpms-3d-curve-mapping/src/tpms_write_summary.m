function files = tpms_write_summary(Ms, tList, thickFine, cfg)
%TPMS_WRITE_SUMMARY  Human-readable summary of the curve mapping.
%
%   Ms        cell array, Ms{k} = tpms_summary_metrics(...) for tList(k)
%   thickFine tpms_wall_thickness on a fine t grid (lookup table)
%
%   Writes to cfg.outDir:
%     gyroid_summary.md               tables with plain-language meanings (open in any viewer)
%     gyroid_summary.csv              metric x t table (one number per cell, opens in Excel)
%     gyroid_summary_detail.csv       same with 95% CI / ranges, long format
%     gyroid_t_vs_wall_thickness.csv  t -> wall thickness, volume fraction, S/V lookup
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
keys = {Ms{1}.key};
nT = numel(tList);
if isempty(cfg.cellSize), uL = 'cell lengths'; else, uL = cfg.lengthUnit; end

% ---------------- wide CSV --------------------------------------------------
files.csv = fullfile(cfg.outDir, 'gyroid_summary.csv');
fid = fopen(files.csv, 'w');
fprintf(fid, 'group,metric,unit');
fprintf(fid, ',t=%.2f', tList); fprintf(fid, ',trend,meaning\n');
for k = 1:numel(keys)
    [m0, v] = rowVals(Ms, keys{k});
    fprintf(fid, '%s,"%s",%s', m0.group, m0.label, m0.unit);
    fprintf(fid, ',%.4g', v);
    fprintf(fid, ',"%s","%s"\n', trendText(v), m0.meaning);
end
fclose(fid);

% ---------------- detail CSV ------------------------------------------------
files.detail = fullfile(cfg.outDir, 'gyroid_summary_detail.csv');
fid = fopen(files.detail, 'w');
fprintf(fid, 't,group,key,metric,unit,value,low,high,range_type,meaning\n');
for j = 1:nT
    for m = Ms{j}
        fprintf(fid, '%.4g,%s,%s,"%s",%s,%.6g,%.6g,%.6g,%s,"%s"\n', tList(j), m.group, m.key, ...
            m.label, m.unit, m.value, m.lo, m.hi, m.rangeType, m.meaning);
    end
end
fclose(fid);

% ---------------- thickness lookup CSV --------------------------------------
files.thickness = fullfile(cfg.outDir, 'gyroid_t_vs_wall_thickness.csv');
fid = fopen(files.thickness, 'w');
fprintf(fid, ['t,wall_thickness_mean_%s,wall_thickness_p5_%s,wall_thickness_p95_%s,wall_thickness_min_%s,' ...
    'wall_thickness_max_%s,wall_thickness_firstorder_%s,wall_thickness_over_cell,thickness_cv,' ...
    'solid_fraction_sheet,solid_fraction_network,area_per_cell_%s2,SV_sheet_per_%s,SV_network_per_%s,reliable\n'], ...
    uL, uL, uL, uL, uL, uL, uL, uL, uL);
for r = thickFine
    rel = r.thick_fail_frac < 0.01;
    fprintf(fid, '%.4g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.4g,%.5g,%.5g,%.5g,%.5g,%.5g,%d\n', r.t, ...
        r.thick_mean, r.thick_p5, r.thick_p95, r.thick_min, r.thick_max, r.thick_firstorder_mean, ...
        r.thick_mean_cell, r.thick_cv, r.vf_sheet, r.vf_network, r.area_per_cell, r.SV_sheet, r.SV_network, rel);
end
fclose(fid);

% ---------------- Markdown report -------------------------------------------
files.md = fullfile(cfg.outDir, 'gyroid_summary.md');
fid = fopen(files.md, 'w');
fprintf(fid, '# Gyroid curve-mapping summary\n\n');
if isempty(cfg.cellSize)
    fprintf(fid, 'Lattice: **%s** gyroid. Lengths are fractions of the cell length (no cell size set). ', cfg.latticeType);
else
    fprintf(fid, 'Lattice: **%s** gyroid, cell size **%g %s**. ', cfg.latticeType, cfg.cellSize, cfg.lengthUnit);
end
fprintf(fid, 'Build direction [%s]. Generated %s.\n\n', num2str(cfg.buildDir), datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf(fid, ['**t** is the level-set value of the gyroid formula: raising it thickens the sheet walls ' ...
    '(or, for a network lattice, adds material to one side). Each column below is one t value.\n\n']);
fprintf(fid, ['Numbers in brackets are 95%% confidence intervals from resampling the traced curves ' ...
    '(narrow = trustworthy), except for wall thickness and bending radius, where they give the 5-95%% ' ...
    'spread over the surface. "> x" means longer than the longest distance checked.\n\n']);
groups = {'structure', 'What t does to the part'; 'curves', 'How the surface curves behave'};
for g = 1:size(groups,1)
    sel = find(strcmp({Ms{1}.group}, groups{g,1}));
    if isempty(sel), continue; end
    fprintf(fid, '## %s\n\n| Metric | Unit |', groups{g,2});
    fprintf(fid, ' t = %.2f |', tList); fprintf(fid, ' Trend |\n|---|---|');
    fprintf(fid, '%s---|\n', repmat('---|', 1, nT));
    for k = sel
        [m0, v, lo, hi] = rowVals(Ms, keys{k});
        fprintf(fid, '| %s | %s |', m0.label, m0.unit);
        for j = 1:nT, fprintf(fid, ' %s |', cellText(v(j), lo(j), hi(j), rtAt(Ms, keys{k}, j))); end
        fprintf(fid, ' %s |\n', trendText(v));
    end
    fprintf(fid, '\n**What the metrics mean**\n\n');
    for k = sel, fprintf(fid, '- **%s**: %s\n', Ms{1}(k).label, Ms{1}(k).meaning); end
    fprintf(fid, '\n');
end
fprintf(fid, '## Choosing t for a target wall thickness (sheet)\n\n');
fprintf(fid, '| t | Wall thickness (%s) | 5-95%% spread | Thickness / cell | Solid fraction |\n|---|---|---|---|---|\n', uL);
for r = thickFine
    if abs(r.t/0.1 - round(r.t/0.1)) > 1e-9 || r.t == 0, continue; end
    if r.thick_fail_frac < 0.01, flag = ''; else, flag = ' (unreliable: sheet near breakup)'; end
    fprintf(fid, '| %.1f | %.3g%s | %.3g-%.3g | %.3f | %.0f%% |\n', r.t, r.thick_mean, flag, ...
        r.thick_p5, r.thick_p95, r.thick_mean_cell, 100*r.vf_sheet);
end
fprintf(fid, ['\nRule of thumb: thickness = %.2f x t x cell size for t up to about 0.5; above that the exact ' ...
    'value grows faster. Use `tpms_t_for_thickness(target, cellSize)` for an exact inverse. The full table ' ...
    'is in `gyroid_t_vs_wall_thickness.csv`.\n'], ruleOfThumb(thickFine));
fclose(fid);
end

% -------------------------------------------------------------------------
function [m0, v, lo, hi] = rowVals(Ms, key)
nT = numel(Ms); v = nan(1,nT); lo = v; hi = v; m0 = [];
for j = 1:nT
    i = find(strcmp({Ms{j}.key}, key), 1);
    if isempty(i), continue; end
    m = Ms{j}(i); if isempty(m0), m0 = m; end
    v(j) = m.value; lo(j) = m.lo; hi(j) = m.hi;
end
end

function s = cellText(v, lo, hi, rt)
if ~isfinite(v), s = 'n/a'; return; end
if strcmp(rt, 'LOWER'), s = ['> ' fmt(v)]; return; end
s = fmt(v);
if isfinite(lo) && isfinite(hi) && hi > lo, s = sprintf('%s [%s, %s]', s, fmt(lo), fmt(hi)); end
end

function rt = rtAt(Ms, key, j)
i = find(strcmp({Ms{j}.key}, key), 1);
if isempty(i), rt = ''; else, rt = Ms{j}(i).rangeType; end
end

function s = fmt(v)
if abs(v) < 5e-4, s = '0';
elseif abs(v) >= 100, s = sprintf('%.0f', v);
elseif abs(v) >= 10, s = sprintf('%.1f', v);
elseif abs(v) < 0.01, s = sprintf('%.3f', v);
else, s = sprintf('%.3g', v); end
end

function s = trendText(v)
ok = find(isfinite(v));
if numel(ok) < 2, s = ''; return; end
a = v(ok(1)); b = v(ok(end));
d = b - a; scale = max(abs([a b]));
if scale == 0 || abs(d) < 0.05*scale, s = 'steady'; return; end
if sign(a) == sign(b) && abs(a) >= 0.1*abs(b)
    p = 100*d/abs(a);
    if d > 0, s = sprintf('rises %.0f%%', p); else, s = sprintf('falls %.0f%%', -p); end
else
    if d > 0, s = sprintf('rises from %s to %s', fmt(a), fmt(b)); else, s = sprintf('falls from %s to %s', fmt(a), fmt(b)); end
end
end

function c = ruleOfThumb(TT)
sel = [TT.t] > 0 & [TT.t] <= 0.5;
c = mean([TT(sel).thick_mean_cell] ./ [TT(sel).t]);
end
