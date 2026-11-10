%% RUN_GYROID_CURVE_MAPPING
% Method 1 (surface curve network) on the analytic gyroid level set
%     sin x cos y + sin y cos z + sin z cos x = t
%
% Traces lines of curvature (k1, k2) and asymptotic curves (asym+, asym-)
% directly on the level set, with every geometric quantity evaluated from the
% analytic derivatives of F (Symbolic Toolbox -> matlabFunction). No STL.
%
% Output (cfg.outDir = results/gyroid_curve_results in the repository):
%   gyroid_summary.md                  START HERE: plain-language numbers per t, with meanings
%   gyroid_summary.csv / _detail.csv   the same as tables (detail has 95% CIs)
%   gyroid_t_vs_wall_thickness.csv     t -> wall thickness, solid fraction, S/V lookup
%   gyroid_curve_descriptors_PCA.csv   one row per t, PCA-role descriptors only
%   gyroid_curve_descriptors_long.csv  every descriptor, every role, with notes
%   gyroid_curve_checks.csv            numerical checks per t and family
%   gyroid_curves_percurve.csv         one row per traced curve
%   gyroid_curve_results.mat           curves, surface samples, metadata
%   *.png                              curve networks and distributions
%
% See README.md for conventions and for which descriptors are genuinely new.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.

clear; clc;

% Put the repository's src/ folder (all tpms_*.m functions) on the path, so the
% script works no matter what MATLAB's current folder is.
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir), thisDir = pwd; end
addpath(fullfile(thisDir, 'src'));
if exist('tpms_build_fields', 'file') ~= 2
    error(['tpms_build_fields.m not found. Keep this script in the repository root, ' ...
           'next to the src/ folder that holds the tpms_*.m functions.']);
end

%% ---------------------------- configuration ----------------------------
cfg.derivMethod  = 'symbolic';     % 'symbolic' (Symbolic Toolbox) | 'manual' (hand-derived)
cfg.tList        = 0:0.2:0.8;      % level-set values to analyse
cfg.cellSize     = 5;              % <-- SET THIS: physical cell length; [] -> radian units (cell = 2*pi)
cfg.lengthUnit   = 'mm';           % unit of cellSize, used in labels
cfg.latticeType  = 'sheet';        % 'sheet' (solid -t<F<t) or 'network' (solid F<t): summary only
cfg.buildDir     = [0 0 1];        % build / print direction for alignment descriptors

cfg.nGrid        = 64;             % marching-cubes grid for quadrature points (per cell edge)
cfg.nSeeds       = 96;             % seeds per curve family (see run_sampling_study: ~100 for most summary
                                   % numbers to +-3%; twist pitch / zero-bending memory need ~300-450)
cfg.families     = {'k1','k2','asym+','asym-'};

% tracer (all lengths in radians; one cell = 2*pi)
cfg.h0           = 0.05;
cfg.hmin         = 1e-5;
cfg.hmax         = 0.08;
cfg.maxTurnDeg   = 4;              % max tangent rotation per accepted step
cfg.Lmax         = 3*2*pi;         % max arc length in EACH direction from the seed
cfg.maxSteps     = 6000;
cfg.closeTol     = 0.05;           % closure distance (modulo lattice)
cfg.tortWindow   = 2*pi;           % main window for straightness/tortuosity = one cell length
cfg.windowsCells = [0.25 0.5 1 2]; % straightness profile windows, in cell lengths (cell-size independent)
cfg.acStep       = 0.05;           % arc-length grid for the direction-memory autocorrelation
cfg.acMaxLag     = 4*pi;           % longest lag (2 cells)

% singularity handling, relative to the median over the surface of
% (k1-k2) for umbilics and |K| for parabolic lines
cfg.umbStopRel   = 0.02;           % stop a line of curvature when k1-k2 < this * median
cfg.umbExclRel   = 0.10;           % exclude from statistics when k1-k2 < this * median
cfg.parStopRel   = 0.01;           % stop an asymptotic curve when K > -this * median|K|
cfg.parExclRel   = 0.05;

% summary for people
cfg.nBoot          = 400;          % bootstrap resamples for 95% confidence intervals
cfg.overhangDeg    = 45;           % downward-facing overhang threshold
cfg.thicknessTGrid = 0:0.05:1.4;   % fine t grid for the t -> wall thickness table
cfg.vfGrid         = 128;          % voxel grid for volume fractions

% validation
cfg.validate      = true;
cfg.validation_t  = 0.3;           % non-zero so the inversion-symmetry check is active
cfg.validationTol.stopOnFail = false;

% output
cfg.outDir       = fullfile(thisDir, 'results', 'gyroid_curve_results');
cfg.makePlots    = true;
cfg.figVisible   = 'off';

rng(1);   % seeding is deterministic anyway; this fixes the validation sample points

%% ---------------------------- fields ------------------------------------
fprintf('Building analytic fields (%s)...\n', cfg.derivMethod);
tic; fld = tpms_build_fields('gyroid', cfg.derivMethod); fprintf('  done in %.1f s\n', toc);

%% ---------------------------- validation --------------------------------
if cfg.validate
    V = tpms_validate(fld, cfg);
end

%% ---------------------------- analysis ----------------------------------
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
results = [];
for it = 1:numel(cfg.tList)
    t = cfg.tList(it);
    fprintf('Level set t = %.3f ... ', t); tic;
    R = tpms_analyze_level(fld, t, cfg);
    fprintf('%.1f s  (area/cell %.3f, K>0 area frac %.3f)\n', toc, ...
        R.scalars(strcmp({R.scalars.name}, 'area_per_cell')).value, ...
        R.scalars(strcmp({R.scalars.name}, 'areafrac_Kpos')).value);
    if cfg.makePlots
        try
            tpms_plot_results(R, cfg);
        catch err
            warning('Plotting failed for t = %.3f: %s', t, err.message);
        end
    end
    if isempty(results), results = R; else, results(end+1) = R; end %#ok<SAGROW>
end

%% ---------------------------- output ------------------------------------
files = tpms_write_results(results, cfg);

%% ---------------------------- summary for people -------------------------
fprintf('Wall thickness / volume fraction table ... '); tic;
tAll  = unique(round([cfg.thicknessTGrid(:); cfg.tList(:)]*1e9)/1e9).';
thick = tpms_wall_thickness(fld, tAll, cfg);
fprintf('%.1f s\n', toc);
Ms = cell(1, numel(results));
for k = 1:numel(results)
    ik = find(abs([thick.t] - results(k).t) < 1e-8, 1);
    Ms{k} = tpms_summary_metrics(results(k), thick(ik), cfg);
end
thickFine = thick(ismember(round([thick.t]*1e9), round(cfg.thicknessTGrid*1e9)));
sfiles = tpms_write_summary(Ms, [results.t], thickFine, cfg);
if cfg.makePlots
    try
        tpms_plot_summary(Ms, [results.t], thickFine, cfg);
    catch err
        warning('Summary plotting failed: %s', err.message);
    end
end
fprintf('\nWrote:\n'); disp(files); disp(sfiles);
fprintf('Start with %s\n', sfiles.md);

%% ---------------------------- summary -----------------------------------
fprintf('Median of PCA-role curve descriptors by level set:\n');
fprintf('%-26s', 'descriptor'); fprintf('%10.2f', cfg.tList); fprintf('\n');
names = {results(1).desc(strcmp({results(1).desc.role}, 'PCA')).name};
for k = 1:numel(names)
    fprintf('%-26s', names{k});
    for r = 1:numel(results)
        i = find(strcmp({results(r).desc.name}, names{k}), 1);
        if isempty(i), fprintf('%10s', '-'); continue; end
        st = tpms_wstats(results(r).desc(i).values, results(r).desc(i).weights);
        fprintf('%10.4f', st.median);
    end
    fprintf('\n');
end
