%% RUN_TPMS_LIBRARY
% Map the core curvature metrics across every surface in tpms_library and
% rank which lattices have the most extreme values.
%
% Every lattice is analysed at t = 0 (its literature mid-surface) with the
% same unit-cell size, using the lean core metric set. Each surface is first
% sanity-checked: genus from Gauss-Bonnet (compared with the literature value
% where known), closeness to a minimal surface, smoothness, cubic symmetry.
%
% Resumable: each finished lattice is saved to results/tpms_library/per_lattice/,
% and a rerun skips those (cfg.resume). cfg.timeBudget stops starting new
% lattices after that many seconds, so long runs can be done in chunks.
%
% Electroplating (lattice Wagner and inverse Damkohler numbers) is compared at a
% matched solid fraction, cfg.platingSolidFraction, not at t = 0.
%
% Non-cubic lattices (tetragonal, hexagonal) come from the add-on
% tpms_library_noncubic; they are flagged in the 'cell' column and starred (*)
% in the report and heatmap. cfg.includeNonCubic = false gives the cubic-only run.
%
% Output (results/tpms_library):
%   tpms_library_report.md     START HERE: checks, metric table, rankings
%   tpms_library_metrics.csv   one row per lattice, every number
%   tpms_library_heatmap.png   z-scores of the core metrics (MATLAB figures)

clear; clc;
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(fullfile(thisDir, 'src'));

%% ---------------------------- configuration ----------------------------
cfg.keys         = {};             % {} = every entry of tpms_library; or e.g. {'gyroid','iwp'}
cfg.t            = 0;              % level-set value (0 = literature mid-surface)
cfg.cellSize     = 5;              % unit-cell size used for every lattice
cfg.lengthUnit   = 'mm';
cfg.latticeType  = 'sheet';
cfg.buildDir     = [0 0 1];

cfg.nGrid        = 96;             % finer than the single-lattice default: some surfaces
                                   % have frequency-2 terms (finer features)
cfg.genusRefineGrid = 160;         % if the genus is not an integer, recompute it on this grid
cfg.nSeeds       = 48;
cfg.families     = {'k1','k2'};    % core set needs only the lines of curvature
cfg.descriptorSet = 'core';
cfg.coreBuildAlignment = false;

cfg.h0 = 0.05; cfg.hmin = 1e-5; cfg.hmax = 0.08; cfg.maxTurnDeg = 4;
cfg.Lmax = 3*2*pi; cfg.maxSteps = 12000; cfg.closeTol = 0.05;
cfg.tortWindow = 2*pi; cfg.windowsCells = [0.25 0.5 1 2]; cfg.acStep = 0.05; cfg.acMaxLag = 4*pi;
cfg.umbStopRel = 0.02; cfg.umbExclRel = 0.10; cfg.parStopRel = 0.01; cfg.parExclRel = 0.05;
cfg.nBoot = 300; cfg.overhangDeg = 45;
% electroplating numbers, compared at a MATCHED SOLID FRACTION (at t = 0 a sheet
% has zero wall thickness). Placeholder bath values from Clemens et al. (2025):
% RANK designs only. Set cfg.plating = [] to skip.
cfg.platingSolidFraction = 0.16;   % 16% solid = porosity 0.84, as the paper's strut lattices
cfg.vfGrid = 96;
cfg.plating = struct('N_field', 4, 'N_flow', 8, 'T', 298, 'alpha_c', 0.485, 'kappa', 6.4, ...
    'i_avg', 3.5, 'n', 2, 'C_o', 10, 'v', 0.05, 'WaTarget', 1, 'DaInvTarget', [440 2280]);
                                   % v = 0.05 m/s is a REFERENCE flow for comparing Da^-1

cfg.resume       = true;           % skip lattices already saved
cfg.timeBudget   = Inf;            % seconds; stop starting new lattices after this
cfg.outDir       = fullfile(thisDir, 'results', 'tpms_library');
cfg.figVisible   = 'off';
cfg.makePlots    = true;

cfg.checkDuplicates = true;        % prove which equations describe the same lattice
cfg.includeNonCubic = true;        % add-on: tetragonal / hexagonal lattices from
                                   % tpms_library_noncubic (starred * in the report)

rng(1);
perDir = fullfile(cfg.outDir, 'per_lattice');
if ~exist(perDir, 'dir'), mkdir(perDir); end

%% ---------------------------- run ---------------------------------------
lib = tpms_library_all(cfg.includeNonCubic);
if ~isempty(cfg.keys), lib = lib(ismember({lib.key}, lower(cfg.keys))); end
if cfg.checkDuplicates
    fprintf('Checking for duplicate lattices ... '); td = tic;
    dups = tpms_library_duplicates(lib);
    save(fullfile(cfg.outDir, 'duplicates.mat'), 'dups', '-v7');
    fprintf('%d duplicate pair(s) (%.1f s)\n', numel(dups), toc(td));
    for q = 1:numel(dups), fprintf('  %s == %s (%s)\n', dups(q).a, dups(q).b, dups(q).how); end
end
t0 = tic;
for k = 1:numel(lib)
    E = lib(k);
    f = fullfile(perDir, [E.key '.mat']);
    if ~isempty(E.duplicateOf)
        fprintf('%-16s duplicate of %s: skipped\n', E.key, E.duplicateOf);
        if exist(f, 'file'), delete(f); end
        continue
    end
    if cfg.resume && exist(f, 'file')
        s = load(f); row = s.row;
        needPlating = ~isempty(cfg.plating) && strcmp(row.status, 'ok') && ...
                      (~isfield(row, 'plating') || isempty(row.plating));
        if ~needPlating, fprintf('%-16s already done\n', E.key); continue; end
        if toc(t0) > cfg.timeBudget, fprintf('Time budget reached; rerun to continue.\n'); break; end
        fprintf('%-16s adding plating ... ', E.key); tl = tic;
        try
            row.plating = tpms_library_plating(tpms_build_fields(E, 'symbolic'), cfg);
        catch err
            row.plating = []; fprintf('plating failed: %s ', err.message);
        end
        save(f, 'row', '-v7'); fprintf('(%.1f s)\n', toc(tl));
        continue
    end
    if toc(t0) > cfg.timeBudget, fprintf('Time budget reached; rerun to continue.\n'); break; end
    fprintf('%-16s ... ', E.key); tl = tic;
    row = struct('key', E.key, 'name', E.name, 'kind', E.kind, 'symmetry', E.symmetry, ...
                 'genusRef', E.genusRef, 'expr', E.expr, 'ref', E.ref, 'status', 'ok', ...
                 'diag', [], 'M', [], 'checks', [], 'plating', [], 'seconds', NaN);
    try
        fld = tpms_build_fields(E, 'symbolic');
        S = tpms_surface_sample(fld, cfg.t, cfg.nGrid);
        G = tpms_geometry(fld, S.v, true);
        row.diag = tpms_library_check(fld, S, G);
        row.diag.refined = false;
        if abs(row.diag.genusRaw - row.diag.genus) > 0.1 && row.diag.smoothness >= 1e-3
            Sf = tpms_surface_sample(fld, cfg.t, cfg.genusRefineGrid);
            Df = tpms_library_check(fld, Sf, tpms_geometry(fld, Sf.v, true));
            row.diag.genusRaw = Df.genusRaw; row.diag.genus = Df.genus;
            row.diag.genusPrim = Df.genusPrim; row.diag.refined = true;
        end
        if row.diag.smoothness < 1e-3
            row.status = 'singular point on surface';
        else
            R = tpms_analyze_level(fld, cfg.t, cfg);
            row.M = tpms_summary_metrics(R, [], cfg);
            row.checks = R.checks;
            if ~isempty(cfg.plating) && row.diag.genusRaw >= 0.5
                row.plating = tpms_library_plating(fld, cfg);
            end
        end
    catch err
        row.status = ['failed: ' err.message];
    end
    row.seconds = toc(tl);
    save(f, 'row', '-v7');
    fprintf('%s (%.1f s)\n', row.status, row.seconds);
end

%% ---------------------------- report ------------------------------------
files = tpms_library_report(cfg);
fprintf('\nSee %s\n', files.md);
