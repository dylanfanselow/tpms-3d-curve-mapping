function R = tpms_analyze_level(fld, t, cfg)
%TPMS_ANALYZE_LEVEL  Method-1 analysis of one gyroid level set F0 = t.
%
%   R = tpms_analyze_level(fld, t, cfg)
%
%   Every descriptor is tagged with a ROLE so the redundancy question from the
%   handoff is answered in the output itself:
%     'PCA'        new information (not a function of k1, k2 at a point);
%                  written to the wide PCA-ready table
%     'CHECK'      duplicate of another PCA descriptor sampled differently
%                  (e.g. kappa_g along traced curves vs the kappa_g field)
%     'REDUNDANT'  provably a restatement of k1, k2 / K (identities)
%     'MIXED'      combines new and redundant parts (Frenet kappa of lines of curvature)
%     'REFERENCE'  point-based k1, k2, H, K -- same as the STL tool, kept for
%                  the analytic-vs-STL comparison (build step 2)
%     'CONTEXT'    geometry context (areas, fractions)
%     'EXTENDED'   valid descriptor left out of the lean core set (cfg.descriptorSet =
%                  'core', the default); still computed and written to the long table
%
%   CORE SET (cfg.descriptorSet = 'core'): the PCA inputs are
%     field_gradH        mean-curvature gradient (surface field)
%     field_kg_LoC       geodesic curvature of both line-of-curvature families, pooled
%     curve_k1_straight  straightness of k1-lines over one cell
%     curve_k2_straight  straightness of k2-lines over one cell
%     handedness_LoC     signed share of twisting length (chirality)
%   plus field_d1_build if cfg.coreBuildAlignment is true.
%   cfg.descriptorSet = 'full' restores every PCA-role descriptor.
%
%   Units: coordinates are radians (cell = 2*pi) unless cfg.cellSize is given,
%   in which case every descriptor is converted to physical units of that length.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.

S   = tpms_surface_sample(fld, t, cfg.nGrid);
G   = tpms_geometry(fld, S.v, true);
opt = tpms_default_trace_opts(cfg, G);
b   = cfg.buildDir(:).' / norm(cfg.buildDir);
opt.buildDir = b;
if isempty(cfg.cellSize), lam = 1; else, lam = cfg.cellSize/(2*pi); end

desc = struct('name',{},'source',{},'family',{},'role',{},'dim',{},'values',{},'weights',{},'note',{});
scal = struct('name',{},'role',{},'value',{},'note',{});
    function addD(name, src, fam, role, dim, v, w, note)
        desc(end+1) = struct('name',name,'source',src,'family',fam,'role',role, ...
            'dim',dim,'values',v(:)*lam^(-dim),'weights',w(:),'note',note);
    end
    function addS(name, role, value, note)
        scal(end+1) = struct('name',name,'role',role,'value',value,'note',note);
    end

okU = G.gap >= opt.umbExcl;           % away from umbilics
okA = G.K   <= -opt.parExcl;          % asymptotic directions well defined
wU  = S.w .* okU;
wA  = S.w .* okA .* okU;

% ---- surface fields (area-weighted; no tracing needed) ----------------------
addD('field_gradH', 'field', '', 'PCA', 2, vnorm(G.gradS_H), S.w, '|surface gradient of H| (3rd derivatives of F)');
addD('field_gradK', 'field', '', 'PCA', 3, vnorm(G.gradS_K), S.w, '|surface gradient of K|');
addD('field_kg_k1lines', 'field', 'k1', 'PCA', 1, abs(G.kg1), wU, '|geodesic curvature| of k1-line field (Codazzi)');
addD('field_kg_k2lines', 'field', 'k2', 'PCA', 1, abs(G.kg2), wU, '|geodesic curvature| of k2-line field (Codazzi)');
addD('field_kg_asymP', 'field', 'asym+', 'PCA', 1, abs(G.kgAp), wA, '|geodesic curvature| of asym+ field (Liouville)');
addD('field_kg_asymM', 'field', 'asym-', 'PCA', 1, abs(G.kgAm), wA, '|geodesic curvature| of asym- field (Liouville)');
addD('field_d1_build', 'field', 'k1', 'PCA', 0, abs(G.d1*b.'), wU, '|cos| angle between d1 and build direction');
addD('field_d2_build', 'field', 'k2', 'PCA', 0, abs(G.d2*b.'), wU, '|cos| angle between d2 and build direction');
addD('field_kg_LoC', 'field', 'k1+k2', 'PCA', 1, [abs(G.kg1); abs(G.kg2)], [wU; wU], ...
     '|geodesic curvature| of both line-of-curvature families, pooled (Codazzi)');
addD('ref_k1', 'field', '', 'REFERENCE', 1, G.k1, S.w, 'point-based (same as STL tool)');
addD('ref_k2', 'field', '', 'REFERENCE', 1, G.k2, S.w, 'point-based (same as STL tool)');
addD('ref_H',  'field', '', 'REFERENCE', 1, G.H,  S.w, 'point-based (same as STL tool)');
addD('ref_K',  'field', '', 'REFERENCE', 2, G.K,  S.w, 'point-based (same as STL tool)');

wn = wU/sum(wU);
A4 = sum(wn .* sum(G.d1.^4, 2));
addS('A4_d1', 'PCA', A4, 'cubic anisotropy <d1x^4+d1y^4+d1z^4>; isotropic = 0.6');
addS('area_per_cell', 'CONTEXT', S.A*lam^2, 'surface area in one unit cell');
addS('areafrac_Kpos', 'CONTEXT', sum(S.w.*(G.K > 0))/S.A, 'no asymptotic curves here');
addS('areafrac_near_umbilic', 'CONTEXT', sum(S.w.*~okU)/S.A, 'excluded from kappa_g statistics');
if isfield(cfg, 'overhangDeg'), ovh = cfg.overhangDeg; else, ovh = 45; end
nb = G.n*b.';
addS('areafrac_normal_down', 'CONTEXT', sum(S.w.*(nb < -cosd(ovh)))/S.A, ...
     sprintf('n.b < -cos(%g deg), n = gradF/|gradF| (outward normal of the F<t solid)', ovh));
addS('areafrac_normal_up', 'CONTEXT', sum(S.w.*(nb > cosd(ovh)))/S.A, ...
     sprintf('n.b > cos(%g deg)', ovh));
addS('areafrac_d1_near_build', 'CONTEXT', sum(wU.*(abs(G.d1*b.') > cosd(30)))/sum(wU), ...
     'max-bending direction within 30 deg of build direction; isotropic = 0.134');

% ---- traced curves ---------------------------------------------------------
fams = cfg.families;
R.curves = struct(); R.checks = struct(); R.perCurve = struct(); R.acc = struct();
for f = 1:numel(fams)
    fam = fams{f}; key = famKey(fam);
    isLoC = any(strcmp(fam, {'k1','k2'}));
    if isLoC
        allowed = G.gap > 2*opt.umbExcl;
    else
        allowed = G.gap > 2*opt.umbExcl & G.K < -2*opt.parExcl;
    end
    si = tpms_seed_points(S.v, cfg.nSeeds, allowed, [], fld.period);
    if isempty(si)
        addS([key '_nCurves'], 'CONTEXT', 0, 'no valid seeds (field undefined on this surface)');
        continue
    end
    C = tpms_trace_curves(fld, t, S.v(si,:), fam, opt);
    [D, pc] = tpms_curve_descriptors(fld, t, C, opt);
    nd = D.nodes;

    iDef = find(abs(D.windowsUsed - opt.tortWindow) < 1e-9, 1);
    for q = 1:numel(D.windowsUsed)
        Wc = D.windowsUsed(q)/(2*pi);
        if q == iDef
            addD(['curve_' key '_straight'], 'curve', fam, 'PCA', 0, D.straight{q}, ones(size(D.straight{q})), ...
                 sprintf('straightness chord/W over W = %.3g cell (bounded 0..1)', Wc));
        else
            addD(sprintf('curve_%s_straight_W%s', key, strrep(sprintf('%g', Wc), '.', 'p')), 'curve', fam, ...
                 'CONTEXT', 0, D.straight{q}, ones(size(D.straight{q})), sprintf('straightness over W = %.3g cell', Wc));
        end
    end
    addD(['curve_' key '_tort'], 'curve', fam, 'CONTEXT', 0, D.windows, ones(size(D.windows)), ...
         sprintf('tortuosity W/chord, W = %.3g cell; unbounded when curves loop -> use _straight', opt.tortWindow/(2*pi)));
    addD(['curve_' key '_kg'], 'curve', fam, 'CHECK', 1, abs(nd.kg), nd.w, 'field_kg_* re-sampled along curves; length weighting is not area weighting (curve density varies), so medians sit ~10-30% lower; diagnostic only');
    if isLoC
        addD(['curve_' key '_tau'], 'curve', fam, 'PCA', 1, nd.tau, nd.w, 'Frenet torsion (signed; chirality-sensitive)');
        addD(['curve_' key '_kappa'], 'curve', fam, 'MIXED', 1, nd.kappa, nd.w, 'sqrt(k_i^2 + kappa_g^2)');
    else
        addD(['curve_' key '_tau'], 'curve', fam, 'REDUNDANT', 1, nd.tau, nd.w, 'equals -/+ sqrt(-K) (Beltrami-Enneper)');
        addD(['curve_' key '_build'], 'curve', fam, 'PCA', 0, nd.bz, nd.w, '|cos| angle between curve tangent and build direction');
    end
    ends = [{C.statusStart}, {C.statusEnd}];
    sing = ismember(ends, {'umbilic','parabolic','step_underflow'});
    addS([key '_frac_singular_ends'], 'PCA', mean(sing), 'fraction of curve ends at a field singularity');
    addS([key '_mean_length_rel'], 'PCA', mean([C.length])/(2*opt.Lmax), 'mean curve length / maximum possible');
    addS([key '_frac_closed'], 'CONTEXT', mean([C.closed]), 'closed within 2*Lmax');
    addS([key '_persistence'], 'PCA', D.persistence*lam, 'direction memory: arc length where <T(s).T(s+lag)> = 1/e');
    addS([key '_nCurves'], 'CONTEXT', numel(C), '');

    k = D.check;
    R.checks.(key) = struct('kg_relerr_median', median(abs(k.kg_ana-k.kg_num))/median(abs(k.kg_ana)), ...
        'kn_resid_max', max(k.kn_resid), 'tg_resid_max', max(k.tg_resid), ...
        'drift_max', max(k.drift), 'frac_nodes_excluded', k.nExcluded/max(k.nNodes,1));
    R.curves.(key)   = C;
    R.perCurve.(key) = pc;
    R.acc.(key)      = struct('lags', D.acLags, 'sum', D.acSum, 'cnt', D.acCnt);
end

% handedness of the lines of curvature, as a single PCA value
pcs = [];
for fk = {'k1', 'k2'}
    if isfield(R.perCurve, fk{1}) && ~isempty(R.perCurve.(fk{1})), pcs = [pcs, R.perCurve.(fk{1})]; end %#ok<AGROW>
end
if ~isempty(pcs)
    addS('handedness_LoC', 'PCA', sum([pcs.tau_wsign])/sum([pcs.wsum]), ...
         'share of line-of-curvature length twisting + minus share twisting - (-1..+1)');
end

% ---- core / full descriptor set ------------------------------------------------
if isfield(cfg, 'descriptorSet') && ~isempty(cfg.descriptorSet), ds = lower(cfg.descriptorSet); else, ds = 'core'; end
core = {'field_gradH', 'field_kg_LoC', 'curve_k1_straight', 'curve_k2_straight', 'handedness_LoC'};
if isfield(cfg, 'coreBuildAlignment') && cfg.coreBuildAlignment, core{end+1} = 'field_d1_build'; end
switch ds
    case 'core'
        for i = 1:numel(desc)
            if strcmp(desc(i).role, 'PCA') && ~any(strcmp(desc(i).name, core)), desc(i).role = 'EXTENDED'; end
        end
        for i = 1:numel(scal)
            if strcmp(scal(i).role, 'PCA') && ~any(strcmp(scal(i).name, core)), scal(i).role = 'EXTENDED'; end
        end
    case 'full'
        i = find(strcmp({desc.name}, 'field_kg_LoC'));       % duplicates the split k1/k2 fields
        desc(i).role = 'EXTENDED';
    otherwise
        error('tpms_analyze_level:descriptorSet', 'cfg.descriptorSet must be ''core'' or ''full''.');
end

R.t = t; R.S = S; R.G = G; R.opt = opt; R.desc = desc; R.scalars = scal; R.lambda = lam;
R.descriptorSet = ds;
end

% -------------------------------------------------------------------------
function k = famKey(fam)
k = strrep(strrep(fam, '+', 'P'), '-', 'M');
end
function n = vnorm(A)
n = sqrt(sum(A.^2, 2));
end
