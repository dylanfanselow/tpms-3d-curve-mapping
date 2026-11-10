function files = tpms_write_results(results, cfg)
%TPMS_WRITE_RESULTS  Write CSV tables (long, PCA-wide, checks, per-curve) and a MAT file.
%   Plain fprintf CSV so it works without the table class; column layout
%   mirrors the STL tool's statistic set (mean std CV p5 median p95 RMS).
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
statN = {'mean','std','CV','p5','median','p95','RMS'};
files = struct();

% ---- long table ------------------------------------------------------------
files.long = fullfile(cfg.outDir, 'gyroid_curve_descriptors_long.csv');
fid = fopen(files.long, 'w');
fprintf(fid, 't,name,kind,source,family,role,inv_length_power,%s,n,note\n', strjoin(statN, ','));
for r = 1:numel(results)
    R = results(r);
    for d = 1:numel(R.desc)
        D = R.desc(d); st = tpms_wstats(D.values, D.weights);
        fprintf(fid, '%.6g,%s,dist,%s,%s,%s,%d', R.t, D.name, D.source, D.family, D.role, D.dim);
        for q = 1:7, fprintf(fid, ',%.10g', st.(statN{q})); end
        fprintf(fid, ',%d,"%s"\n', st.n, D.note);
    end
    for s = 1:numel(R.scalars)
        Sc = R.scalars(s);
        fprintf(fid, '%.6g,%s,scalar,,,%s,,%.10g,,,,,,,1,"%s"\n', R.t, Sc.name, Sc.role, Sc.value, Sc.note);
    end
end
fclose(fid);

% ---- wide PCA-ready table (role == 'PCA' only) --------------------------------
files.pca = fullfile(cfg.outDir, 'gyroid_curve_descriptors_PCA.csv');
hdr = {'t'}; rows = zeros(numel(results), 0);
for r = 1:numel(results)
    R = results(r); vals = R.t; h = {'t'};
    for d = 1:numel(R.desc)
        if ~strcmp(R.desc(d).role, 'PCA'), continue; end
        st = tpms_wstats(R.desc(d).values, R.desc(d).weights);
        for q = 1:7
            h{end+1} = [R.desc(d).name '_' statN{q}]; vals(end+1) = st.(statN{q}); %#ok<AGROW>
        end
    end
    for s = 1:numel(R.scalars)
        if ~strcmp(R.scalars(s).role, 'PCA'), continue; end
        h{end+1} = R.scalars(s).name; vals(end+1) = R.scalars(s).value; %#ok<AGROW>
    end
    if r == 1, hdr = h; rows = vals;
    else
        % align columns by name (a family can be missing for some t)
        row = nan(1, numel(hdr));
        for c = 1:numel(h)
            j = find(strcmp(hdr, h{c}), 1);
            if isempty(j), hdr{end+1} = h{c}; rows(:,end+1) = NaN; row(end+1) = vals(c); %#ok<AGROW>
            else, row(j) = vals(c); end
        end
        rows(end+1, 1:numel(row)) = row; %#ok<AGROW>
    end
end
fid = fopen(files.pca, 'w');
fprintf(fid, '%s\n', strjoin(hdr, ','));
for r = 1:size(rows,1)
    fprintf(fid, '%s\n', strjoin(arrayfun(@(v) sprintf('%.10g', v), rows(r,:), 'UniformOutput', false), ','));
end
fclose(fid);

% ---- numerical checks per level and family ------------------------------------
files.checks = fullfile(cfg.outDir, 'gyroid_curve_checks.csv');
fid = fopen(files.checks, 'w');
fprintf(fid, 't,family,kg_relerr_median,kn_resid_max,tg_resid_max,drift_max,frac_nodes_excluded\n');
for r = 1:numel(results)
    fn = fieldnames(results(r).checks);
    for k = 1:numel(fn)
        c = results(r).checks.(fn{k});
        fprintf(fid, '%.6g,%s,%.4g,%.4g,%.4g,%.4g,%.4g\n', results(r).t, fn{k}, c.kg_relerr_median, ...
            c.kn_resid_max, c.tg_resid_max, c.drift_max, c.frac_nodes_excluded);
    end
end
fclose(fid);

% ---- per-curve table --------------------------------------------------------------
files.perCurve = fullfile(cfg.outDir, 'gyroid_curves_percurve.csv');
fid = fopen(files.perCurve, 'w');
fprintf(fid, 't,family,curve,length,statusStart,statusEnd,closed,tx,ty,tz,closeGap,straightness_1cell,mean_abs_kg,mean_kappa,mean_tau,nNodes\n');
for r = 1:numel(results)
    lam = results(r).lambda;
    fn = fieldnames(results(r).perCurve);
    for k = 1:numel(fn)
        pc = results(r).perCurve.(fn{k});
        for m = 1:numel(pc)
            p = pc(m);
            fprintf(fid, '%.6g,%s,%d,%.6g,%s,%s,%d,%g,%g,%g,%.4g,%.4g,%.6g,%.6g,%.6g,%d\n', results(r).t, fn{k}, m, ...
                p.length*lam, p.statusStart, p.statusEnd, p.closed, p.translation + 0, p.closeGap*lam, p.straight1, ...
                p.mean_abs_kg/lam, p.mean_kappa/lam, p.mean_tau/lam, p.nNodes);
        end
    end
end
fclose(fid);

% ---- MAT file (curves in radian coordinates, cell = [0, 2*pi]^3) ----------------------
files.mat = fullfile(cfg.outDir, 'gyroid_curve_results.mat');
out = struct();
for r = 1:numel(results)
    out(r).t = results(r).t;
    out(r).curves = results(r).curves;
    out(r).checks = results(r).checks;
    out(r).desc = rmfield(results(r).desc, {'values','weights'});
    out(r).scalars = results(r).scalars;
    out(r).surface_v = results(r).S.v;
    out(r).surface_f = results(r).S.f;
end
cfgSaved = cfg; %#ok<NASGU>
save(files.mat, 'out', 'cfgSaved', '-v7');
end
