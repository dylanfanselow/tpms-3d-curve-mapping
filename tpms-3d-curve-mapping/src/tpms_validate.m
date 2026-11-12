function V = tpms_validate(fld, cfg)
%TPMS_VALIDATE  Self-checks for the Method-1 gyroid curve pipeline.
%
%   V = tpms_validate(fld, cfg)   prints a report and returns pass/fail + values.
%
%   1  Symbolic vs hand-derived derivatives (F, grad, Hess, 3rd derivatives)
%   2  Eigen-solve K, H vs closed-form level-set formulas
%   3  Hellmann-Feynman curvature derivatives vs finite differences on the surface
%   4  Analytic kappa_g (Codazzi / Liouville) vs finite differences of the traced
%      polyline, at two step sizes -> observed convergence order (expect ~2).
%      Both hmax and maxTurnDeg are halved: on lines of curvature the turn
%      limit, not hmax, usually sets the step.
%   5  Identities that make some curve quantities REDUNDANT with k1, k2:
%        lines of curvature: kappa_n = k_i, tau_g = 0
%        asymptotic curves : kappa_n = 0,  |tau_g| = sqrt(-K)  (Beltrami-Enneper)
%   6  Curves stay on the surface (|F-t|/|gradF|)
%   7  Surface area per cell converges with the sampling grid
%   8  Cubic symmetry: area-weighted orientation tensor of d1 = I/3
%   9  Inversion symmetry x -> -x maps F=t to F=-t and swaps k1<->k2,
%      asym+<->asym-: the field statistics must match across the swap
%  10  Electroplating numbers: Wa and Da^-1 match the written formulas;
%      Wa = WaTarget at i_avg = i_max; Da^-1 = target at v = v_min;
%      geometric similarity (cell x2, a_v /2) halves Wa_geom, keeps Da_geom
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.

tol = cfg.validationTol;
V = struct(); rep = {};
t = cfg.validation_t;

% 1 ------------------------------------------------------------------------
X = 2*pi*rand(2000,3);
if strcmp(fld.method, 'symbolic')
    fm = tpms_build_fields('gyroid', 'manual');
    e = max([max(abs(fld.F0(X)-fm.F0(X))), max(max(abs(fld.grad(X)-fm.grad(X)))), ...
             max(max(abs(fld.hess(X)-fm.hess(X)))), max(max(abs(fld.d3(X)-fm.d3(X))))]);
    rep = addRow(rep, '1 symbolic vs manual derivatives (max abs diff)', e, e < 1e-12);
    V.sym_vs_manual = e;
else
    rep = addRow(rep, '1 symbolic vs manual derivatives', NaN, true, 'skipped (manual mode)');
end

% 2 ------------------------------------------------------------------------
S = tpms_surface_sample(fld, t, cfg.nGrid);
G = tpms_geometry(fld, S.v, true);
eK = max(abs(G.K - G.K_closed)); eH = max(abs(G.H - G.H_closed));
rep = addRow(rep, '2 K eigen vs closed form (max abs)', eK, eK < 1e-12);
rep = addRow(rep, '2 H eigen vs closed form (max abs)', eH, eH < 1e-12);

% 3 ------------------------------------------------------------------------
ii = 1:max(1,floor(size(S.v,1)/200)):size(S.v,1);
x  = S.v(ii,:); Gi = tpms_geometry(fld, x, true); ep = 1e-5;
dirs = {Gi.d1, Gi.d2}; err = 0;
ok = Gi.gap > 0.05*median(G.gap);
for j = 1:2
    xp = tpms_project(fld, x + ep*dirs{j}, t);
    xm = tpms_project(fld, x - ep*dirs{j}, t);
    Gp = tpms_geometry(fld, xp); Gm = tpms_geometry(fld, xm);
    err = max(err, max(abs((Gp.k1(ok)-Gm.k1(ok))/(2*ep) - Gi.dk(ok,j))));
    err = max(err, max(abs((Gp.k2(ok)-Gm.k2(ok))/(2*ep) - Gi.dk(ok,2+j))));
end
rel = err / max(abs(Gi.dk(:)));
rep = addRow(rep, '3 dk (Hellmann-Feynman) vs FD (max rel)', rel, rel < 1e-6);

% 4, 5, 6 --------------------------------------------------------------------
opt  = tpms_default_trace_opts(cfg, G);
allowed = G.gap > 0.2*median(G.gap) & G.K < -0.2*median(abs(G.K));
seeds = S.v(tpms_seed_points(S.v, 10, allowed), :);
fams = {'k1','k2','asym+','asym-'};
hs = [cfg.hmax, cfg.hmax/2];
for f = 1:4
    e2 = zeros(1,2);
    for q = 1:2
        o = opt; o.hmax = hs(q); o.h0 = hs(q)/2; o.Lmax = 2*pi;
        o.maxTurnDeg = cfg.maxTurnDeg * hs(q)/hs(1);   % refine BOTH step limits
        C = tpms_trace_curves(fld, t, seeds, fams{f}, o);
        D = tpms_curve_descriptors(fld, t, C, o);
        k = D.check;
        e2(q) = median(abs(k.kg_ana - k.kg_num)) / median(abs(k.kg_ana));
        if q == 1
            isLoC = f <= 2;
            if isLoC, nm = {'kappa_n - k_i', '|tau_g|'}; else, nm = {'|kappa_n|', '| |tau_g| - sqrt(-K) |'}; end
            r1 = max(k.kn_resid); r2 = max(k.tg_resid); dr = max(k.drift);
            rep = addRow(rep, sprintf('5 %-5s %s (max)', fams{f}, nm{1}), r1, r1 < 1e-10);
            rep = addRow(rep, sprintf('5 %-5s %s (max)', fams{f}, nm{2}), r2, r2 < 1e-10);
            rep = addRow(rep, sprintf('6 %-5s surface drift (max)', fams{f}), dr, dr < 1e-10);
        end
    end
    ord = log2(e2(1)/e2(2));
    rep = addRow(rep, sprintf('4 %-5s kappa_g ana vs FD rel err @h=%.3g / %.3g', fams{f}, hs(1), hs(2)), ...
                 e2, e2(2) < 1e-2 && ord > 1.5, sprintf('order %.2f', ord));
    V.kg_order.(strrep(strrep(fams{f}, '+', 'P'), '-', 'M')) = ord;   % valid MATLAB field names
end

% 7 ------------------------------------------------------------------------
ng = round(cfg.nGrid*[0.5 1 1.5]);
A = zeros(size(ng));
for q = 1:numel(ng), Sq = tpms_surface_sample(fld, t, ng(q)); A(q) = Sq.A/(2*pi)^2; end
dA = abs(A(end) - A(end-1))/A(end);
rep = addRow(rep, sprintf('7 area/a^2 at nGrid=%s', mat2str(ng)), A, dA < 5e-3);

% 8 ------------------------------------------------------------------------
w = S.w/sum(S.w);
Mt = zeros(3);
for i = 1:3, for j = 1:3, Mt(i,j) = sum(w.*G.d1(:,i).*G.d1(:,j)); end, end
dev = max(abs(Mt(:) - reshape(eye(3)/3, [], 1)));
rep = addRow(rep, '8 orientation tensor of d1 vs I/3 (max abs)', dev, dev < 2e-2);

% 9 ------------------------------------------------------------------------
if t ~= 0
    Sn = tpms_surface_sample(fld, -t, cfg.nGrid);
    Gn = tpms_geometry(fld, Sn.v, true);
    pairs = {'k1 vs -k2',        G.k1,             -Gn.k2;
             '|kg1| vs |kg2|',   abs(G.kg1),       abs(Gn.kg2);
             '|kgA+| vs |kgA-|', abs(G.kgAp),      abs(Gn.kgAm);
             '|gradS H|',        vnorm(G.gradS_H), vnorm(Gn.gradS_H)};
    for p = 1:size(pairs,1)
        a = tpms_wstats(pairs{p,2}, S.w); b = tpms_wstats(pairs{p,3}, Sn.w);
        d = abs(a.median - b.median) / max(abs(a.median), 1e-12);
        rep = addRow(rep, sprintf('9 inversion t<->-t  %s (median rel diff)', pairs{p,1}), d, d < 2e-2);
    end
end

% 10 -----------------------------------------------------------------------
% tests the formulas with the paper's placeholder values, independent of cfg.plating
pc = struct('N_field',4,'N_flow',8,'T',298,'alpha_c',0.485,'kappa',6.4,'i_avg',3.5, ...
            'n',2,'C_o',10,'v',0.05,'WaTarget',1,'DaInvTarget',[440 2280]);
e0 = 0.75; av0 = 1200; L0 = 0.005;
P  = tpms_plating_numbers(e0, av0, L0, pc);
VT = 8.314*pc.T/(pc.alpha_c*96485);
WaRef = VT*pc.kappa*e0/(pc.i_avg*av0*(pc.N_field*L0)^2);
DaRef = pc.n*96485*pc.v*pc.C_o/(pc.i_avg*av0*pc.N_flow*L0);
q = pc; q.i_avg = P.i_max;   Pq = tpms_plating_numbers(e0, av0, L0, q);
r = pc; r.v = P.v_min(1);    Pr = tpms_plating_numbers(e0, av0, L0, r);
Ps = tpms_plating_numbers(e0, av0/2, 2*L0, pc);
eP = max(abs([P.Wa/WaRef - 1, P.DaInv/DaRef - 1, Pq.Wa/pc.WaTarget - 1, ...
              Pr.DaInv/pc.DaInvTarget(1) - 1, Ps.Wa_geom/P.Wa_geom - 0.5, Ps.Da_geom/P.Da_geom - 1]));
rep = addRow(rep, '10 electroplating numbers self-consistency (max rel)', eP, eP < 1e-12);

% report ---------------------------------------------------------------------
fprintf('\n==== Validation (gyroid, t = %.3g) ====\n', t);
allPass = true;
for r = 1:numel(rep)
    R = rep{r};
    fprintf('  [%s] %-58s %s  %s\n', ternary(R.pass,'PASS','FAIL'), R.name, fmtv(R.value), R.note);
    allPass = allPass && R.pass;
end
fprintf('  Overall: %s\n\n', ternary(allPass, 'ALL PASS', 'SOME CHECKS FAILED'));
V.rows = rep; V.allPass = allPass; V.t = t; V.area_per_a2 = A;
if ~allPass && tol.stopOnFail
    error('tpms_validate:fail', 'Validation failed; see report above.');
end
end

% -------------------------------------------------------------------------
function rep = addRow(rep, name, value, pass, note)
if nargin < 5, note = ''; end
rep{end+1} = struct('name', name, 'value', value, 'pass', logical(pass), 'note', note);
end
function s = fmtv(v)
s = strjoin(arrayfun(@(q) sprintf('%.3g', q), v, 'UniformOutput', false), ' / ');
end
function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
function n = vnorm(A)
n = sqrt(sum(A.^2, 2));
end
