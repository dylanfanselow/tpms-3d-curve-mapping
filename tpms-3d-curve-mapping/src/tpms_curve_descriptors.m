function [D, perCurve] = tpms_curve_descriptors(fld, t, C, opt)
%TPMS_CURVE_DESCRIPTORS  Space-curve descriptors along traced surface curves.
%
%   [D, perCurve] = tpms_curve_descriptors(fld, t, C, opt)
%
%   For each node of each curve, in the Darboux frame (T, u = n x T, n):
%     kn     normal curvature     = <T', n>     analytic: -T'*HessF*T/|gradF|
%     kg     geodesic curvature   = <T', u>     analytic: Codazzi (lines of curvature)
%                                               or Liouville (asymptotic curves)
%     tg     geodesic torsion     = <u', n>     analytic: -(n x T)'*HessF*T/|gradF|
%     kappa  Frenet curvature     = sqrt(kn^2 + kg^2)
%     tau    Frenet torsion       = tg - dphi/ds,   phi = atan2(kg, kn)
%   plus finite-difference versions (kg_num, kn_num, tau_num) of the traced
%   polyline, used only as validation of the analytic values.
%
%   Straightness over a window W of arc length:  chord/W = |x(s+W) - x(s)| / W
%   (0 = curve returns to where it started, 1 = straight line), evaluated for
%   every window in opt.windows. Tortuosity = 1/straightness at opt.tortWindow.
%   Straightness is dimensionless and cell-size independent when W is a
%   fixed fraction of the cell. Closed curves are unrolled periodically first.
%
%   Direction memory (persistence length): arc length over which the tangent
%   autocorrelation <T(s).T(s+lag)> falls to 1/e. Per-curve sums are kept so
%   confidence intervals can be bootstrapped over curves.
%
%   Nodes too close to a singularity of the field (umbilic for lines of
%   curvature, parabolic line for asymptotic curves) are flagged invalid for
%   the statistics, because kg diverges there by construction; the fraction
%   excluded is reported.
%
%   OUTPUT
%     D.nodes    struct of pooled node arrays (valid nodes only) + arc weights w
%     D.windows  pooled windowed tortuosities at opt.tortWindow
%     D.straight cell, pooled straightness for each window in opt.windows
%     D.acLags, D.acSum, D.acCnt  tangent autocorrelation (per-curve rows)
%     D.persistence  pooled persistence length (rad units; NaN if > max lag)
%     D.check    validation metrics (analytic vs numeric, redundancy identities)
%     perCurve   struct array: length, status, closed, translation, mean |kg|, ...
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.

isLoC = any(strcmp(C(1).family, {'k1','k2'}));
L     = fld.period;
fam   = C(1).family;

pool = struct('kg',[],'kn',[],'tg',[],'kappa',[],'tau',[],'w',[], ...
              'bz',[],'absK',[]);
chk  = struct('kn_resid',[],'tg_resid',[],'kg_ana',[],'kg_num',[], ...
              'tau_ana',[],'tau_num',[],'drift',[],'nNodes',0,'nExcluded',0);
Wl   = opt.windows(:).';
iDef = find(abs(Wl - opt.tortWindow) < 1e-9, 1);
if isempty(iDef), Wl = [Wl opt.tortWindow]; iDef = numel(Wl); end
straight = cell(1, numel(Wl));
acLags = 0:opt.acStep:opt.acMaxLag;
acSum = zeros(0, numel(acLags)); acCnt = zeros(0, numel(acLags));
perCurve = struct('length',{},'statusStart',{},'statusEnd',{},'closed',{}, ...
                  'translation',{},'closeGap',{},'mean_abs_kg',{},'mean_kappa',{}, ...
                  'mean_tau',{},'nNodes',{},'median_abs_kg',{},'median_abs_tau',{}, ...
                  'tau_wsum',{},'tau_wabs',{},'tau_wsign',{},'wsum',{},'straight1',{});

for m = 1:numel(C)
    x = C(m).x; s = C(m).s;
    if size(x,1) < 5, continue; end
    G = tpms_geometry(fld, x, true);
    switch fam
        case 'k1',    Df = G.d1; kgF = G.kg1;
        case 'k2',    Df = G.d2; kgF = G.kg2;
        case 'asym+', Df = G.ap; kgF = G.kgAp;
        case 'asym-', Df = G.am; kgF = G.kgAm;
    end
    chord = nudiff(x, s);
    sig = sign(sum(Df.*chord, 2)); sig(sig == 0) = 1;
    T   = sig .* Df;
    n   = G.n;
    u   = cross(n, T, 2);

    kn  = -sum(T.*hmul(G.Hs, T), 2) ./ G.gn;
    tg  = -sum(u.*hmul(G.Hs, T), 2) ./ G.gn;
    kg  = sig .* kgF;
    kap = sqrt(kn.^2 + kg.^2);
    phi = unwrap(2*atan2(kg, kn))/2;          % mod pi: N and -N give the same torsion
    tau = tg - nudiff(phi, s);

    Tp      = nudiff(T, s);
    Tpp     = nudiff(Tp, s);
    kg_num  = sum(u.*Tp, 2);
    kn_num  = sum(n.*Tp, 2);
    tau_num = sum(cross(T, Tp, 2).*Tpp, 2) ./ sum(Tp.^2, 2);

    w = zeros(size(s));
    ds = diff(s);
    w(1:end-1) = w(1:end-1) + ds/2;
    w(2:end)   = w(2:end)   + ds/2;

    if isLoC
        valid = G.gap >= opt.umbExcl;
        chk.kn_resid = [chk.kn_resid; abs(kn - ifelse(strcmp(fam,'k1'), G.k1, G.k2))];
        chk.tg_resid = [chk.tg_resid; abs(tg)];
    else
        valid = G.K <= -opt.parExcl;
        chk.kn_resid = [chk.kn_resid; abs(kn)];
        chk.tg_resid = [chk.tg_resid; abs(abs(tg) - sqrt(max(-G.K,0)))];
    end
    % interior nodes only for the finite-difference comparison
    int = valid; int([1 2 end-1 end]) = false;
    chk.kg_ana  = [chk.kg_ana;  kg(int)];
    chk.kg_num  = [chk.kg_num;  kg_num(int)];
    chk.tau_ana = [chk.tau_ana; tau(int & kap > 1e-3)];
    chk.tau_num = [chk.tau_num; tau_num(int & kap > 1e-3)];
    chk.drift   = [chk.drift; abs(fld.F0(x) - t) ./ G.gn];
    chk.nNodes    = chk.nNodes + numel(s);
    chk.nExcluded = chk.nExcluded + sum(~valid);

    pool.kg    = [pool.kg;    kg(valid)];
    pool.kn    = [pool.kn;    kn(valid)];
    pool.tg    = [pool.tg;    tg(valid)];
    pool.kappa = [pool.kappa; kap(valid)];
    pool.tau   = [pool.tau;   tau(valid)];
    pool.w     = [pool.w;     w(valid)];
    if isfield(opt,'buildDir'), bdir = opt.buildDir; else, bdir = [0 0 1]; end
    pool.bz    = [pool.bz;    abs(T(valid,:)*bdir(:))];
    pool.absK  = [pool.absK;  abs(G.K(valid))];

    str = windowStraightness(x, s, C(m), Wl, L);
    for q = 1:numel(Wl), straight{q} = [straight{q}; str{q}]; end
    [as, ac] = tangentAutocorr(T, s, C(m).closed, acLags, opt.acStep);
    acSum(end+1,:) = as; acCnt(end+1,:) = ac; %#ok<AGROW>

    pc.length = C(m).length; pc.statusStart = C(m).statusStart;
    pc.statusEnd = C(m).statusEnd; pc.closed = C(m).closed;
    pc.translation = C(m).translation; pc.closeGap = C(m).closeGap;
    pc.mean_abs_kg = wmean(abs(kg(valid)), w(valid));
    pc.mean_kappa  = wmean(kap(valid), w(valid));
    pc.mean_tau    = wmean(tau(valid), w(valid));
    pc.nNodes      = numel(s);
    skg = tpms_wstats(abs(kg(valid)), w(valid));  pc.median_abs_kg  = skg.median;
    stu = tpms_wstats(abs(tau(valid)), w(valid)); pc.median_abs_tau = stu.median;
    pc.tau_wsum    = sum(w(valid).*tau(valid));
    pc.tau_wabs    = sum(w(valid).*abs(tau(valid)));
    pc.tau_wsign   = sum(w(valid).*sign(tau(valid)));   % length twisting + minus length twisting -
    pc.wsum        = sum(w(valid));
    if isempty(str{iDef}), pc.straight1 = NaN; else, pc.straight1 = mean(str{iDef}); end
    perCurve(end+1) = pc; %#ok<AGROW>
end

D.family  = fam;
D.nodes   = pool;
D.windows  = 1 ./ max(straight{iDef}, eps);
D.straight = straight;
D.windowsUsed = Wl;
D.acLags = acLags; D.acSum = acSum; D.acCnt = acCnt;
D.persistence = tpms_persistence(acLags, sum(acSum,1), sum(acCnt,1));
D.check   = chk;
end

% -------------------------------------------------------------------------
function out = windowStraightness(x, s, c, Wl, L)
% chord/W for every start node and every window W (cell array, one per W)
out = cell(1, numel(Wl));
Wmax = max(Wl);
if c.closed
    Lc = s(end);
    if Lc <= 0, return; end
    nrep = ceil(Wmax/Lc) + 1;
    xx = x; ss = s;
    for k = 1:nrep
        xx = [xx; x(2:end,:) + k*L*c.translation]; %#ok<AGROW>
        ss = [ss; s(2:end) + k*Lc];                %#ok<AGROW>
    end
    s0all = s(1:end-1);
else
    xx = x; ss = s; s0all = s;
end
for q = 1:numel(Wl)
    W = Wl(q);
    if c.closed, s0 = s0all; else, s0 = s0all(s0all <= ss(end) - W); end
    if isempty(s0), out{q} = zeros(0,1); continue; end
    xe = interp1(ss, xx, s0 + W, 'linear');
    xs = interp1(ss, xx, s0,     'linear');
    out{q} = sqrt(sum((xe - xs).^2, 2)) / W;
end
end

function [aSum, aCnt] = tangentAutocorr(T, s, closed, lags, ds)
% sums of T(s).T(s+lag) on a uniform arc-length grid
aSum = zeros(1, numel(lags)); aCnt = zeros(1, numel(lags));
if s(end) < 2*ds, return; end
su = (0:ds:s(end)).';
Tu = interp1(s, T, su, 'linear');
Tu = Tu ./ sqrt(sum(Tu.^2, 2));
nl = numel(lags);
if closed
    n0 = numel(su) - 1;                       % one period, then unroll
    nrep = ceil(nl/max(n0,1)) + 1;
    Tu = repmat(Tu(1:n0,:), nrep + 1, 1);
    for k = 1:nl
        d = sum(Tu(1:n0,:) .* Tu(k:k+n0-1,:), 2);
        aSum(k) = sum(d); aCnt(k) = n0;
    end
else
    n = size(Tu,1);
    for k = 1:min(nl, n)
        d = sum(Tu(1:n-k+1,:) .* Tu(k:n,:), 2);
        aSum(k) = sum(d); aCnt(k) = numel(d);
    end
end
end

function d = nudiff(f, s)
% second-order finite difference on a non-uniform grid, column-wise
n = numel(s);
d = zeros(size(f));
if n < 3, return; end
h1 = s(2:n-1) - s(1:n-2);
h2 = s(3:n)   - s(2:n-1);
d(2:n-1,:) = (-h2./(h1.*(h1+h2))).*f(1:n-2,:) + ((h2-h1)./(h1.*h2)).*f(2:n-1,:) ...
           + (h1./(h2.*(h1+h2))).*f(3:n,:);
a = s(2)-s(1); b = s(3)-s(2);
d(1,:) = -(2*a+b)/(a*(a+b))*f(1,:) + (a+b)/(a*b)*f(2,:) - a/(b*(a+b))*f(3,:);
a = s(n-1)-s(n-2); b = s(n)-s(n-1);
d(n,:) = b/(a*(a+b))*f(n-2,:) - (a+b)/(a*b)*f(n-1,:) + (2*b+a)/(b*(a+b))*f(n,:);
end

function Hu = hmul(Hs, u)
Hu = [Hs(:,1).*u(:,1) + Hs(:,4).*u(:,2) + Hs(:,6).*u(:,3), ...
      Hs(:,4).*u(:,1) + Hs(:,2).*u(:,2) + Hs(:,5).*u(:,3), ...
      Hs(:,6).*u(:,1) + Hs(:,5).*u(:,2) + Hs(:,3).*u(:,3)];
end

function m = wmean(v, w)
if isempty(v) || sum(w) == 0, m = NaN; else, m = sum(v.*w)/sum(w); end
end

function out = ifelse(c, a, b)
if c, out = a; else, out = b; end
end
