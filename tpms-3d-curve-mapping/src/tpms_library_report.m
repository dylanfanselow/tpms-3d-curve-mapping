function files = tpms_library_report(cfg)
%TPMS_LIBRARY_REPORT  Combine the per-lattice results of run_tpms_library.
%
%   files = tpms_library_report(cfg)
%
%   Reads results/tpms_library/per_lattice/*.mat and writes
%     tpms_library_metrics.csv   one row per lattice
%     tpms_library_report.md     checks, metric table, rankings
%     tpms_library_heatmap.png   (if plotting works)
%
%   "Extremeness" of a lattice = mean |z-score| over four metrics:
%   curvedness-change length, turning radius and the two straightnesses.
%   The core curvature-change length (from the MEAN-curvature gradient) is
%   shown but not scored: for minimal surfaces H = 0, so at t = 0 it mostly
%   measures the error of the nodal approximation. Handedness is shown but not
%   scored: at t = 0 it is ~0 for every centrosymmetric surface, so it acts as
%   a chirality flag rather than a ranking metric.
%   Two scale-free ratios (turning radius and curvature-change length divided
%   by the typical bending radius) separate "small because the surface is
%   finer" from "steers sharply relative to its own bending".

perDir = fullfile(cfg.outDir, 'per_lattice');
d = dir(fullfile(perDir, '*.mat'));
rows = {};
for i = 1:numel(d)
    s = load(fullfile(perDir, d(i).name));
    rows{end+1} = s.row; %#ok<AGROW>
end
if isempty(rows), error('tpms_library_report:none', 'No results in %s.', perDir); end
lib = tpms_library_all(true);
[~, order] = ismember({lib.key}, cellfun(@(r) r.key, rows, 'UniformOutput', false));
rows = rows(order(order > 0));                      % library order
libK = lib(order > 0);
for i = 1:numel(rows)                               % names, kinds and literature values
    rows{i}.name = libK(i).name;                    % always from the current library
    rows{i}.kind = libK(i).kind;
    rows{i}.genusRef = libK(i).genusRef;
    rows{i}.cell = libK(i).cell;
    if strcmp(libK(i).cell, 'cubic'), rows{i}.star = ''; else, rows{i}.star = ' *'; end
end
anyNC = any(cellfun(@(r) ~isempty(r.star), rows));

mk  = {'curv_change_length', 'turning_radius', 'k1_straightness', 'k2_straightness', 'handedness'};
uL = cfg.lengthUnit;
n = numel(rows);
V = nan(n, numel(mk)); LO = V; HI = V; bend = nan(n,1);
for i = 1:n
    M = rows{i}.M;
    if isempty(M), continue; end
    for j = 1:numel(mk)
        q = M(strcmp({M.key}, mk{j}));
        if ~isempty(q), V(i,j) = q.value; LO(i,j) = q.lo; HI(i,j) = q.hi; end
    end
    q = M(strcmp({M.key}, 'surface_bend_radius')); if ~isempty(q), bend(i) = q.value; end
end
ccl = nan(n,1);                                        % curvedness-change length
for i = 1:n
    if ~isempty(rows{i}.diag) && isfield(rows{i}.diag, 'curvednessChangeLength')
        ccl(i) = rows{i}.diag.curvednessChangeLength * cfg.cellSize/(2*pi);
    end
end
pl = nan(n, 8);                                       % t, porosity, a_v, thickness, Wa, DaInv, i_max, v_min(1)
for i = 1:n
    if isfield(rows{i}, 'plating') && ~isempty(rows{i}.plating)
        P = rows{i}.plating;
        pl(i,:) = [P.t, P.porosity, P.a_v, P.thickness, P.Wa, P.DaInv, P.i_max/10, P.v_min(1)];
    end
end
relTurn = V(:,2) ./ bend;
relChange = V(:,1) ./ bend;

% z-scores and extremeness over analysed lattices
valid = false(n,1);                                     % analysed and not broken
for i = 1:n
    D = rows{i}.diag;
    valid(i) = strcmp(rows{i}.status, 'ok') && ~isempty(D) && D.genusRaw >= 0.5;
end
S4 = [ccl, V(:,2:4)];                                   % scored metrics
S4(~valid,:) = NaN;
s4N = {'Curvedness-change length', 'Turning radius', 'k1 straightness', 'k2 straightness'};
ok = all(isfinite(S4), 2);
Z = nan(size(S4));
for j = 1:size(S4,2)
    col = S4(ok,j); sd = std(col);
    if sd > 0, Z(ok,j) = (S4(ok,j) - mean(col)) / sd; end
end
extreme = nan(n,1);
extreme(ok) = mean(abs(Z(ok,:)), 2);

% ---------------- CSV ----------------------------------------------------------
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
files.csv = fullfile(cfg.outDir, 'tpms_library_metrics.csv');
fid = fopen(files.csv, 'w');
fprintf(fid, 'key,name,kind,symmetry,status,genus_ref_per_primitive,genus_box,genus_raw,translations,side_swapping_translations,genus_per_primitive,minimality,smoothness,cubic_dev,area_per_cell,half_period_axes');
for j = 1:numel(mk), fprintf(fid, ',%s,%s_ci_low,%s_ci_high', mk{j}, mk{j}, mk{j}); end
fprintf(fid, [',curvedness_change_length_%s,surface_bend_radius_%s,turning_radius_rel,curv_change_length_rel,extremeness,' ...
    'plating_t,plating_porosity,plating_a_v_per_m,plating_wall_thickness_%s,Wa,DaInv_at_ref_flow,i_max_mA_per_cm2,' ...
    'flow_for_DaInv_target1_m_per_s,kg_check_relerr,drift_max,seconds,cell\n'], uL, uL, uL);
for i = 1:n
    r = rows{i}; D = r.diag;
    if isempty(D), D = struct('genus',NaN,'genusRaw',NaN,'nTrans',NaN,'nSwap',NaN,'genusPrim',NaN,'minimality',NaN, ...
                              'smoothness',NaN,'cubicDev',NaN,'areaPerCell',NaN,'halfPeriodAxes',''); end
    fprintf(fid, '%s,"%s",%s,%s,"%s",%g,%g,%.4f,%g,%g,%.4g,%.4g,%.4g,%.3g,%.4g,%s', r.key, r.name, r.kind, r.symmetry, ...
        r.status, r.genusRef, D.genus, D.genusRaw, D.nTrans, D.nSwap, D.genusPrim, D.minimality, D.smoothness, ...
        D.cubicDev, D.areaPerCell, D.halfPeriodAxes);
    for j = 1:numel(mk), fprintf(fid, ',%.6g,%.6g,%.6g', V(i,j), LO(i,j), HI(i,j)); end
    [kgE, dr] = checkVals(r);
    fprintf(fid, ',%.6g,%.6g,%.6g,%.6g,%.4g', ccl(i), bend(i), relTurn(i), relChange(i), extreme(i));
    fprintf(fid, ',%.5g', pl(i,:));
    fprintf(fid, ',%.3g,%.3g,%.1f,"%s"\n', kgE, dr, r.seconds, r.cell);
end
fclose(fid);

% ---------------- Markdown -------------------------------------------------------
files.md = fullfile(cfg.outDir, 'tpms_library_report.md');
fid = fopen(files.md, 'w');
fprintf(fid, '# TPMS library: core curvature metrics\n\n');
fprintf(fid, ['%d lattices at t = %g (literature mid-surface), unit cell %g %s, sheet interpretation. ' ...
    'Lengths in %s. Generated %s.\n\n'], n, cfg.t, cfg.cellSize, uL, uL, datestr(now, 'yyyy-mm-dd HH:MM'));

fprintf(fid, '## 1. Sanity checks\n\n');
fprintf(fid, ['The genus of the surface in one box comes from the curvature (Gauss-Bonnet). Dividing out ' ...
    'the translations that map the surface onto itself while keeping its sides gives the genus per primitive ' ...
    'cell, the value quoted in the literature. A mismatch means the equation describes a different surface. ' ...
    '*Minimality* is median |H| / median curvedness: near 0 for a good approximation of a minimal surface; ' ...
    'entries labelled minimal but above 0.35 are flagged. *Sides congruent* means a translation swaps the two ' ...
    'sides (labyrinths). Notes on individual entries are in `src/tpms_library.m`.\n\n']);
fprintf(fid, '| Lattice | Kind | Genus / box | Translations | Genus / primitive cell | Literature | Minimality | Smoothness | Sides congruent | Status |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|\n');
nMatch = 0; nRef = 0;
for i = 1:n
    r = rows{i}; D = r.diag;
    if isempty(D)
        fprintf(fid, '| %s | %s | | | | %s | | | | %s |\n', r.name, r.kind, gtxt(r.genusRef), r.status); continue
    end
    gp = 1 + (D.genusRaw - 1)/D.nTrans; flag = ''; mflag = '';
    if ~strcmp(r.status, 'ok')
        fprintf(fid, '| %s | %s | n/a | %d | n/a | %s | %.3f | %.3f | | %s; excluded from rankings |\n', ...
            r.name, r.kind, D.nTrans, gtxt(r.genusRef), D.minimality, D.smoothness, r.status);
        continue
    end
    if D.genusRaw < 0.5, flag = ' **splits into closed pieces; excluded from rankings**'; end
    if abs(D.genusRaw - D.genus) > 0.1, flag = [flag ' (not integer: inconclusive)']; end
    if isfield(D, 'refined') && D.refined, flag = [flag ' (fine grid)']; end
    if isfinite(r.genusRef)
        nRef = nRef + 1;
        if abs(gp - r.genusRef) < 0.1, flag = [flag ' (matches)']; nMatch = nMatch + 1;
        else, flag = [flag ' **MISMATCH**']; end
    end
    if strcmp(r.kind, 'minimal') && D.minimality > 0.35, mflag = ' **far from minimal**'; end
    if D.smoothness < 0.05, mflag = [mflag ' (nearly singular: treat t = 0 metrics with caution)']; end
    if D.nSwap > 0, sc = 'yes'; else, sc = 'no'; end
    fprintf(fid, '| %s%s | %s | %s | %d | %s%s | %s | %.3f%s | %.3f | %s | %s |\n', r.name, r.star, r.kind, gnum(D.genusRaw), ...
        D.nTrans, gnum(gp), flag, gtxt(r.genusRef), D.minimality, mflag, D.smoothness, sc, r.status);
end
fprintf(fid, '\n%d of %d entries with a literature genus match it.\n', nMatch, nRef);
fprintf(fid, ['\n### Duplicate check\n\nEvery pair of equations is tested for describing the same lattice ' ...
    '(same surface up to rotation, reflection, translation, scale or doubled frequency; `tpms_library_duplicates`).\n\n']);
fd = fullfile(cfg.outDir, 'duplicates.mat');
if exist(fd, 'file')
    sd = load(fd);
    if isempty(sd.dups), fprintf(fid, 'No duplicates found.\n');
    else
        for q = 1:numel(sd.dups)
            fprintf(fid, '- `%s` and `%s` are the same lattice (%s).\n', sd.dups(q).a, sd.dups(q).b, sd.dups(q).how);
        end
    end
end
dl = lib(~cellfun(@isempty, {lib.duplicateOf}));
if ~isempty(dl)
    fprintf(fid, '\nNot analysed because they duplicate an earlier entry: %s.\n', ...
        strjoin(arrayfun(@(e) sprintf('%s (= %s)', e.name, e.duplicateOf), dl, 'UniformOutput', false), '; '));
end
fprintf(fid, '\n## 2. Core metrics\n\n');
fprintf(fid, ['Sorted from most to least extreme: mean |z| over curvedness-change length, turning radius and the ' ...
    'two straightnesses. Brackets: 95%% CI from resampling curves. The mean-curvature-based curvature-change ' ...
    'length (the core metric) is shown in grey italics because at t = 0 it mostly reflects the nodal-approximation ' ...
    'error of minimal surfaces.\n\n']);
fprintf(fid, '| Lattice | Curvedness-change length (%s) | Turning radius (%s) | k1 straightness | k2 straightness | Handedness | *Curvature-change length (%s)* | Bending radius (%s) | Turning / bending | Extremeness |\n', uL, uL, uL, uL);
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|\n');
if anyNC
    fprintf(fid, ['\n\\* = non-cubic cell (tetragonal or hexagonal; `cell` column in the CSV, add-on ' ...
        '`tpms_library_noncubic`). Their straightness window is one length a along x, and the cubic-symmetry ' ...
        'check does not apply.\n\n']);
end
[~, srt] = sort(extreme, 'descend');
srt = [srt(isfinite(extreme(srt))); srt(~isfinite(extreme(srt)))];
for i = srt.'
    fprintf(fid, '| %s%s | %s | %s | %s | %s | %s | *%s* | %s | %s | %s |\n', rows{i}.name, rows{i}.star, f3(ccl(i)), f3(V(i,2)), ...
        ci(V(i,3), LO(i,3), HI(i,3)), ci(V(i,4), LO(i,4), HI(i,4)), ci(V(i,5), LO(i,5), HI(i,5)), ...
        f3(V(i,1)), f3(bend(i)), f3(relTurn(i)), f3(extreme(i)));
end

fprintf(fid, '\n## 3. Who is most extreme?\n\n');
for j = 1:4
    col = S4(:,j); fin = find(isfinite(col));
    if numel(fin) < 2, continue; end
    [~, o] = sort(col(fin), 'descend'); hiI = fin(o(1:min(3,end))); loI = fin(o(end:-1:max(1,end-2)));
    fprintf(fid, '- **%s**: highest %s; lowest %s.\n', s4N{j}, lst(rows, hiI, col), lst(rows, loI, col));
end
hcol = V(:,5); chir = find(valid & isfinite(hcol) & abs(hcol) >= 0.02 & (LO(:,5) > 0 | HI(:,5) < 0));
if isempty(chir)
    fprintf(fid, '- **Chirality**: no lattice has a handedness index clearly different from 0.\n');
else
    fprintf(fid, ['- **Chirality** (|handedness| >= 0.02 and its CI excludes 0): %s. A non-zero value at t = 0 means the ' ...
        'structure itself is single-handed.\n'], lst(rows, chir, hcol));
end
fin = find(valid & isfinite(relTurn));
if numel(fin) >= 2
    [~, o] = sort(relTurn(fin)); 
    fprintf(fid, ['- **Sharpest steering relative to its own bending** (turning / bending radius, low = sharp): %s. ' ...
        'Gentlest: %s.\n'], lst(rows, fin(o(1:min(3,end))), relTurn), lst(rows, fin(o(end:-1:max(1,end-2))), relTurn));
end
if any(ok)
    fprintf(fid, '- **Most extreme overall**: %s.\n', lst(rows, srt(1:min(3, sum(ok))), extreme));
end
if any(isfinite(pl(:,5)) & valid)
    pc = cfg.plating;
    fprintf(fid, '\n## 4. Electroplating at matched solid fraction (ranking only)\n\n');
    fprintf(fid, ['At t = 0 a sheet has zero wall thickness, so the lattices are compared at the same solid fraction, ' ...
        '**%.0f%%** (porosity %.2f, as the reference strut lattices of Clemens et al. 2025), %s interpretation. ' ...
        'Lattice Wagner number Wa at %.3g mA/cm2 (higher = more even coating) and inverse Damkohler number ' ...
        'Da^-1 at a reference flow of %g m/s (higher = less ion depletion), N_field = %g, N_flow = %g, ' ...
        'placeholder bath values from the paper. The equations are reconstructed from the paper (absolute values ' ...
        'may be off by ~1.5x): use them to rank lattices.\n\n'], 100*cfg.platingSolidFraction, ...
        1 - cfg.platingSolidFraction, cfg.latticeType, pc.i_avg/10, pc.v, pc.N_field, pc.N_flow);
    fprintf(fid, ['**At equal porosity both numbers scale as 1/(wetted surface per volume)**, so they rank the ' ...
        'lattices in the same order: lattices with less surface per volume plate more evenly *and* need less ' ...
        'flow. Their absolute values still answer different questions (uniformity vs. depletion).\n\n']);
    fprintf(fid, '| Lattice | t | Wall thickness (%s) | Wetted area a_v (1/m) | Wa | Da^-1 at %g m/s | Max current for Wa >= %g (mA/cm2) | Flow for Da^-1 = %g (m/s) |\n', ...
        uL, pc.v, pc.WaTarget, pc.DaInvTarget(1));
    fprintf(fid, '|---|---|---|---|---|---|---|---|\n');
    pv = find(valid & isfinite(pl(:,5)));
    [~, o] = sort(pl(pv,5), 'descend');
    for i = pv(o).'
        fprintf(fid, '| %s | %.3f | %s | %s | %s | %s | %s | %s |\n', rows{i}.name, pl(i,1), f3(pl(i,4)), ...
            f3(pl(i,3)), f3(pl(i,5)), f3(pl(i,6)), f3(pl(i,7)), f3(pl(i,8)));
    end
    k = numel(pv);
    fprintf(fid, '\n- **Most even coating / least flow needed**: %s.\n', lst(rows, pv(o(1:min(3,k))), pl(:,5)));
    fprintf(fid, '- **Least even coating / most flow needed**: %s.\n', lst(rows, pv(o(end:-1:max(1,k-2))), pl(:,5)));
end
fprintf(fid, ['\n**Reading the results.** All lattices share the same cell size, so surfaces with more ' ...
    'surface per cell (higher-frequency equations such as F-RD or Split P) naturally have smaller radii and ' ...
    'lower straightness. The turning / bending ratio removes that scale effect. Handedness is ~0 at t = 0 ' ...
    'for centrosymmetric surfaces, so at t = 0 it flags structurally chiral lattices (such as the tubular gyroid, ' ...
    'which thickens only one of the two gyroid channels); for the others it becomes informative at t ~= 0.\n\n']);
fprintf(fid, ['**Sources.** Equations were checked against published collections (Bernhard 2026, arXiv:2604.07160, ' ...
    'Tables 6-7, including the "corax" collection; nTop documentation; Gandy et al. 2001; Wohlgemuth et al. 2001). ' ...
    'Where sources disagree, each variant is a separate entry. Community-collected formulas are secondary ' ...
    'sources: the checks in section 1 show which ones behave like the surfaces they are named after.\n']);
fclose(fid);

% ---------------- heatmap -----------------------------------------------------------
if isfield(cfg, 'makePlots') && cfg.makePlots
    try
        fg = figure('Visible', cfg.figVisible, 'Color', 'w', 'Position', [50 50 1300 max(650, 22*n)], 'InvertHardcopy', 'off');
        try, theme(fg, 'light'); catch, end
        ax = axes('Parent', fg);
        ZH = zcol(V(:,5), valid); ZW = zcol(pl(:,5), valid); ZD = zcol(pl(:,6), valid);
        Zmap = [Z, ZH, ZW, ZD];
        mapN = [s4N, {'Handedness', 'Wagner Wa', 'Damkohler Da^{-1}'}];
        imagesc(ax, Zmap(srt, :)); colormap(ax, bwr(256)); caxis(ax, [-2.5 2.5]);
        cb = colorbar(ax); set(cb, 'Color', 'k', 'FontSize', 12); ylabel(cb, 'z-score', 'Color', 'k', 'FontSize', 13);
        names = cellfun(@(r) [r.name r.star], rows, 'UniformOutput', false);
        set(ax, 'YTick', 1:n, 'YTickLabel', names(srt), 'XTick', 1:numel(mapN), 'XTickLabel', mapN, ...
            'XTickLabelRotation', 20, 'FontSize', 12, 'XColor', 'k', 'YColor', 'k', 'TickLabelInterpreter', 'none');
        title(ax, sprintf('Curvature metrics at t = %g (z-scores; rows: most to least extreme)', cfg.t), 'Color', 'k', 'FontSize', 15);
        print(fg, fullfile(cfg.outDir, 'tpms_library_heatmap.png'), '-dpng', '-r130');
        if strcmp(cfg.figVisible, 'off'), close(fg); end
    catch err
        warning('Heatmap failed: %s', err.message);
    end
end
end

% -------------------------------------------------------------------------
function [kgE, dr] = checkVals(r)
kgE = NaN; dr = NaN;
if isstruct(r.checks) && isfield(r.checks, 'k1')
    kgE = r.checks.k1.kg_relerr_median; dr = r.checks.k1.drift_max;
end
end
function s = gnum(g)
if abs(g - round(g)) < 0.1, s = sprintf('%d', round(g)); else, s = sprintf('%.2f', g); end
end
function s = gtxt(g)
if isfinite(g), s = sprintf('%d', g); else, s = '?'; end
end
function s = f3(v)
if ~isfinite(v), s = 'n/a'; elseif abs(v) >= 10, s = sprintf('%.1f', v); else, s = sprintf('%.3g', v); end
end
function s = ci(v, lo, hi)
s = f3(v);
if isfinite(lo) && isfinite(hi), s = sprintf('%s [%s, %s]', s, f3(lo), f3(hi)); end
end
function s = lst(rows, idx, col)
parts = cell(1, numel(idx));
for k = 1:numel(idx), parts{k} = sprintf('%s (%s)', rows{idx(k)}.name, f3(col(idx(k)))); end
s = strjoin(parts, ', ');
end
function z = zcol(v, valid)
z = nan(size(v)); ok = valid & isfinite(v);
if sum(ok) > 1 && std(v(ok)) > 0, z(ok) = (v(ok) - mean(v(ok))) / std(v(ok)); end
end
function c = bwr(m)
h = floor(m/2); u = linspace(0, 1, h)';
c = [[u u ones(h,1)]; [ones(m-h,1) flipud(linspace(0,1,m-h)') flipud(linspace(0,1,m-h)')]];
end
