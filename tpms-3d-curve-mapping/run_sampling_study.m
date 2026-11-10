%% RUN_SAMPLING_STUDY
% How much do the summary numbers depend on the sampling settings?
%
% Varies ONE setting at a time around the defaults, at a single level set:
%   nSeeds     number of traced curves per family     (statistical sampling)
%   nGrid      surface quadrature grid per cell edge  (field statistics, seed pool)
%   stepScale  tracing step size and turn limit       (curve discretisation)
% and compares every summary metric with the finest setting of that factor.
%
% Output (cfg.outDir = results/gyroid_sampling_study in the repository):
%   sampling_study.csv   long table: factor, setting, metric, value, CI, deviation
%   sampling_study.md    tables + recommended settings
%   sampling_study.png   deviation from the finest setting vs setting
%
% Deviation is reported in the metric's own units for bounded metrics
% (straightness 0-1, handedness -1..+1, percentages) and as a relative % for
% everything else. Tolerances: study.tolRel (%), study.tolAbs (0-1 metrics),
% study.tolPct (percentage points).
%
% Recommendation rules (cheapest setting that passes):
%   nSeeds     PRECISION: every 95% CI half-width <= tolerance, and the value
%              agrees with the finest run within tolerance + both CI half-widths
%   nGrid      changes the seed set too, so differences are allowed up to
%              tolerance + both CI half-widths (sampling noise is not bias)
%   stepScale  same seeds in every run, so any difference is discretisation
%              error and must be within the tolerance itself
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.

clear; clc;
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(fullfile(thisDir, 'src'));

%% ---------------------------- configuration ----------------------------
study.t          = 0.4;
study.nSeeds     = [12 24 48 96 160];
study.nGrid      = [32 48 64 96];
study.stepScale  = [2 1 0.5];          % multiplies hmax and maxTurnDeg
study.tolRel     = 3;                  % % for unbounded metrics (lengths, ratios)
study.tolAbs     = 0.02;               % for 0-1 and -1..+1 metrics
study.tolPct     = 10;                 % percentage points for % metrics (loop shares are binomial: noisy)

% base settings (same as run_gyroid_curve_mapping defaults)
cfg.cellSize = 5; cfg.lengthUnit = 'mm'; cfg.latticeType = 'sheet'; cfg.buildDir = [0 0 1];
cfg.nGrid = 64; cfg.nSeeds = 48; cfg.families = {'k1','k2','asym+','asym-'};
cfg.h0 = 0.05; cfg.hmin = 1e-5; cfg.hmax = 0.08; cfg.maxTurnDeg = 4; cfg.Lmax = 3*2*pi;
cfg.maxSteps = 12000; cfg.closeTol = 0.05; cfg.tortWindow = 2*pi; cfg.windowsCells = [0.25 0.5 1 2];
cfg.acStep = 0.05; cfg.acMaxLag = 4*pi;
cfg.umbStopRel = 0.02; cfg.umbExclRel = 0.10; cfg.parStopRel = 0.01; cfg.parExclRel = 0.05;
cfg.nBoot = 300; cfg.overhangDeg = 45;
cfg.outDir = fullfile(thisDir, 'results', 'gyroid_sampling_study'); cfg.figVisible = 'off'; cfg.makePlots = true;

rng(1);
fld = tpms_build_fields('gyroid', 'manual');   % identical numbers to 'symbolic', no toolbox needed
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end

%% ---------------------------- runs --------------------------------------
factors = {'nSeeds', study.nSeeds; 'nGrid', study.nGrid; 'stepScale', study.stepScale};
runs = struct('factor', {}, 'setting', {}, 'M', {}, 'seconds', {});
for f = 1:size(factors,1)
    for v = factors{f,2}
        c = cfg;
        switch factors{f,1}
            case 'nSeeds',    c.nSeeds = v;
            case 'nGrid',     c.nGrid = v;
            case 'stepScale', c.hmax = cfg.hmax*v; c.h0 = cfg.h0*v; c.maxTurnDeg = cfg.maxTurnDeg*v;
        end
        fprintf('%-9s = %-5g ... ', factors{f,1}, v); tic;
        R = tpms_analyze_level(fld, study.t, c);
        M = tpms_summary_metrics(R, [], c);
        el = toc; fprintf('%.1f s\n', el);
        runs(end+1) = struct('factor', factors{f,1}, 'setting', v, 'M', M, 'seconds', el); %#ok<SAGROW>
    end
end

%% ---------------------------- compare with finest ----------------------
finest = struct('nSeeds', max(study.nSeeds), 'nGrid', max(study.nGrid), 'stepScale', min(study.stepScale));
keys = {runs(1).M.key};
rows = struct('factor',{},'setting',{},'key',{},'label',{},'unit',{},'value',{},'lo',{},'hi',{}, ...
              'dev',{},'devUnit',{},'hw',{},'ok',{},'seconds',{});
for r = 1:numel(runs)
    ref = runs([runs.setting] == finest.(runs(r).factor) & strcmp({runs.factor}, runs(r).factor));
    for k = 1:numel(keys)
        m = runs(r).M(strcmp({runs(r).M.key}, keys{k}));
        m0 = ref.M(strcmp({ref.M.key}, keys{k}));
        [dev, du, okT, hw] = tpms_metric_deviation(m, m0, study);
        [~, ~, ~, hw0] = tpms_metric_deviation(m0, m0, study);
        switch du, case 'abs', tol = study.tolAbs; case 'points', tol = study.tolPct; otherwise, tol = study.tolRel; end
        agree = ~isfinite(dev) && okT || abs(dev) <= tol + hw + hw0;
        switch runs(r).factor
            case 'nSeeds',    ok = agree && hw <= tol;
            case 'nGrid',     ok = agree;
            case 'stepScale', ok = okT;
        end
        rows(end+1) = struct('factor', runs(r).factor, 'setting', runs(r).setting, 'key', keys{k}, ...
            'label', m.label, 'unit', m.unit, 'value', m.value, 'lo', m.lo, 'hi', m.hi, ...
            'dev', dev, 'devUnit', du, 'hw', hw, 'ok', ok, 'seconds', runs(r).seconds); %#ok<SAGROW>
    end
end

%% ---------------------------- recommendation ---------------------------
% seeds needed per metric for its CI half-width to reach the tolerance
% (half-width ~ 1/sqrt(nSeeds)), estimated from the finest nSeeds run
nMax = max(study.nSeeds);
need = struct('label', {}, 'hw', {}, 'tol', {}, 'unit', {}, 'n', {});
for k = 1:numel(keys)
    q = rows(strcmp({rows.factor}, 'nSeeds') & [rows.setting] == nMax & strcmp({rows.key}, keys{k}));
    if q.hw <= 0, continue; end
    switch q.devUnit, case 'abs', tol = study.tolAbs; case 'points', tol = study.tolPct; otherwise, tol = study.tolRel; end
    need(end+1) = struct('label', q.label, 'hw', q.hw, 'tol', tol, 'unit', q.devUnit, ...
                         'n', ceil(nMax*(q.hw/tol)^2)); %#ok<SAGROW>
end
rec = struct(); recNote = struct();
for f = 1:size(factors,1)
    fn = factors{f,1}; vals = factors{f,2};
    if strcmp(fn, 'stepScale'), order = sort(vals, 'descend'); else, order = sort(vals, 'ascend'); end
    rec.(fn) = NaN; recNote.(fn) = '';
    for v = order                               % cheapest first
        sel = strcmp({rows.factor}, fn) & [rows.setting] == v;
        if all([rows(sel).ok]), rec.(fn) = v; break; end
    end
    if ~isfinite(rec.(fn))
        rec.(fn) = order(end);
        if strcmp(fn, 'nSeeds') && ~isempty(need)
            [nn, j] = max([need.n]);
            recNote.(fn) = sprintf([' (precision target NOT reached at %d; about %d seeds needed, driven by ' ...
                '"%s"; see the table below)'], nMax, nn, need(j).label);
        else
            recNote.(fn) = ' (finest tested; tolerance not reached)';
        end
    end
end

%% ---------------------------- output -------------------------------------
fid = fopen(fullfile(cfg.outDir, 'sampling_study.csv'), 'w');
fprintf(fid, 't,factor,setting,key,metric,unit,value,ci_low,ci_high,deviation_from_finest,ci_halfwidth,deviation_unit,passes,run_seconds\n');
for r = rows
    fprintf(fid, '%.3g,%s,%g,%s,"%s",%s,%.6g,%.6g,%.6g,%.4g,%.4g,%s,%d,%.1f\n', study.t, r.factor, r.setting, ...
        r.key, r.label, r.unit, r.value, r.lo, r.hi, r.dev, r.hw, r.devUnit, r.ok, r.seconds);
end
fclose(fid);

fid = fopen(fullfile(cfg.outDir, 'sampling_study.md'), 'w');
fprintf(fid, '# Sampling sensitivity study (gyroid, t = %.2f, cell %g %s)\n\n', study.t, cfg.cellSize, cfg.lengthUnit);
fprintf(fid, ['Each table changes one setting and keeps the others at their defaults (nSeeds = %d, ' ...
    'nGrid = %d, stepScale = 1). Cells show the deviation from the finest setting of that factor, ' ...
    'with the 95%% CI half-width in brackets for curve metrics; **bold** = fails the rule for that factor. ' ...
    'Tolerances: %.0f%% relative, %.2f absolute for 0-1 / -1..+1 metrics, %.0f points for %% metrics.\n\n' ...
    '- **nSeeds** must make every CI narrow enough (precision) and agree with the finest run.\n' ...
    '- **nGrid** also reshuffles the seeds, so differences within sampling noise pass.\n' ...
    '- **stepScale** keeps the same seeds, so any difference is a real discretisation error.\n\n'], ...
    cfg.nSeeds, cfg.nGrid, study.tolRel, study.tolAbs, study.tolPct);
fprintf(fid, '**Recommended (cheapest setting that passes):** nSeeds = %g%s, nGrid = %g%s, stepScale = %g%s.\n\n', ...
    rec.nSeeds, recNote.nSeeds, rec.nGrid, recNote.nGrid, rec.stepScale, recNote.stepScale);
if ~isempty(need)
    fprintf(fid, '### Seeds needed per metric (95%% CI half-width -> tolerance)\n\n');
    fprintf(fid, '| Metric | CI half-width at %d seeds | Tolerance | Seeds needed |\n|---|---|---|---|\n', nMax);
    [~, o] = sort([need.n], 'descend');
    for j = o
        u = need(j).unit; if strcmp(u, 'abs'), u = ''; elseif strcmp(u, 'points'), u = ' pts'; end
        fprintf(fid, '| %s | %.3g%s | %.3g%s | %d |\n', need(j).label, need(j).hw, u, need(j).tol, u, need(j).n);
    end
    fprintf(fid, ['\nRun time grows roughly linearly with seeds (see the run-time row below). Metrics needing ' ...
        'many seeds are medians over individual curves; their CIs in gyroid_summary.md show the precision ' ...
        'actually achieved, so a wide interval there is the signal to add seeds.\n\n']);
end
for f = 1:size(factors,1)
    fn = factors{f,1}; vals = factors{f,2};
    fprintf(fid, '## %s\n\n| Metric | Unit | Finest value |', fn);
    for v = vals, fprintf(fid, ' %g |', v); end
    fprintf(fid, '\n|---|---|---|%s\n', repmat('---|', 1, numel(vals)));
    for k = 1:numel(keys)
        sel = find(strcmp({rows.factor}, fn) & strcmp({rows.key}, keys{k}));
        ref = rows(sel([rows(sel).setting] == finest.(fn)));
        fprintf(fid, '| %s | %s | %s |', ref.label, ref.unit, num2str(ref.value, 4));
        for v = vals
            q = rows(sel([rows(sel).setting] == v));
            if strcmp(q.devUnit, '%'), txt = sprintf('%+.1f%%', q.dev); else, txt = sprintf('%+.3g', q.dev); end
            if ~isfinite(q.dev), txt = 'n/a'; end
            if q.hw > 0 && q.setting ~= finest.(fn)
                if strcmp(q.devUnit, '%'), txt = sprintf('%s (+-%.1f%%)', txt, q.hw); else, txt = sprintf('%s (+-%.2g)', txt, q.hw); end
            end
            if ~q.ok, txt = ['**' txt '**']; end
            fprintf(fid, ' %s |', txt);
        end
        fprintf(fid, '\n');
    end
    secs = arrayfun(@(v) rows(find(strcmp({rows.factor}, fn) & [rows.setting] == v, 1)).seconds, vals);
    fprintf(fid, '| *run time (s)* | | |%s\n\n', sprintf(' %.0f |', secs));
end
fclose(fid);

% plot: |deviation| vs setting, one panel per factor
if cfg.makePlots
    try
        fg = figure('Visible', cfg.figVisible, 'Color', 'w', 'Position', [50 50 1600 500], 'InvertHardcopy', 'off');
        try, theme(fg, 'light'); catch, end
        for f = 1:size(factors,1)
            fn = factors{f,1}; vals = factors{f,2};
            ax = subplot(1, 3, f, 'Parent', fg); hold(ax, 'on');
            for k = 1:numel(keys)
                sel = find(strcmp({rows.factor}, fn) & strcmp({rows.key}, keys{k}));
                d = zeros(size(vals));
                for q = 1:numel(vals)
                    rr = rows(sel([rows(sel).setting] == vals(q)));
                    d(q) = tpms_metric_deviation(rr, [], study);
                end
                plot(ax, vals, d, 'o-', 'LineWidth', 1.3);
            end
            plot(ax, vals, ones(size(vals)), 'k--', 'LineWidth', 1.5);
            set(ax, 'XScale', 'log', 'YScale', 'log', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'FontSize', 13, 'Box', 'on');
            xlabel(ax, fn, 'Color', 'k', 'FontSize', 14);
            ylabel(ax, 'deviation / tolerance', 'Color', 'k', 'FontSize', 14);
            title(ax, sprintf('%s (dashed = tolerance)', fn), 'Color', 'k', 'FontSize', 15);
        end
        print(fg, fullfile(cfg.outDir, 'sampling_study.png'), '-dpng', '-r130');
        close(fg);
    catch err
        warning('Plotting failed: %s', err.message);
    end
end
fprintf('\nRecommended: nSeeds = %g%s\n             nGrid = %g%s, stepScale = %g%s\nSee %s\n', ...
    rec.nSeeds, recNote.nSeeds, rec.nGrid, recNote.nGrid, rec.stepScale, recNote.stepScale, ...
    fullfile(cfg.outDir, 'sampling_study.md'));
