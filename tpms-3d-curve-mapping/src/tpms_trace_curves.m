function C = tpms_trace_curves(fld, t, seeds, family, opt)
%TPMS_TRACE_CURVES  Trace integral curves of a surface line field on F0 = t.
%
%   C = tpms_trace_curves(fld, t, seeds, family, opt)
%
%   family : 'k1'    lines of curvature along d1 (max principal curvature)
%            'k2'    lines of curvature along d2
%            'asym+' asymptotic curves, direction  cos(thA) d1 + sin(thA) d2
%            'asym-' asymptotic curves, direction  cos(thA) d1 - sin(thA) d2
%   seeds  : M-by-3 points on the surface.
%
%   All seeds are traced together (vectorised). Each seed is integrated in
%   both directions, because the fields are LINE fields (no preferred sign):
%   at every RK4 stage the field vector is sign-aligned with the current
%   tangent. After every step the point is Newton-projected back onto F0 = t.
%
%   Step control: a step is rejected (and h halved) if the tangent turns by
%   more than opt.maxTurnDeg, the field is undefined, or projection fails.
%   After easy steps h grows back toward opt.hmax.
%
%   Stop reasons (statusStart / statusEnd):
%     'umbilic'       k1 - k2 < opt.umbTol          (lines of curvature)
%     'parabolic'     K > -opt.parTol               (asymptotic curves)
%     'closed'        passed within closeTol of the seed (modulo the lattice;
%                     segment distance, step-size independent), same tangent
%     'max_length'    arc length reached opt.Lmax
%     'max_steps'     storage limit reached
%     'step_underflow' h fell below opt.hmin (direction field singular)
%     'seed_invalid'  field undefined / singular at the seed
%
%   OUTPUT  C(m) struct: family, seed, x (n-by-3, unwrapped coordinates),
%           s (arc length from start), statusStart, statusEnd, closed,
%           translation (lattice vector of a closed curve, in cells),
%           closeGap (distance to the seed at closure, modulo the lattice), length.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.

isLoC = any(strcmp(family, {'k1','k2'}));
L     = tpms_period3(fld.period);   % [Px Py Pz]
M     = size(seeds,1);
nH    = 2*M;
X0    = [seeds; seeds];
sg    = [ones(M,1); -ones(M,1)];

[D0, st0] = fieldDir(fld, X0, family);
T0 = sg .* D0;

P   = nan(nH, opt.maxSteps, 3);
S   = nan(nH, opt.maxSteps);
cnt = ones(nH,1);
P(:,1,:) = reshape(X0, nH, 1, 3);
S(:,1)   = 0;

x = X0; T = T0; s = zeros(nH,1); h = opt.h0*ones(nH,1);
status = repmat({''}, nH, 1);
trans  = zeros(nH,3);
cgap   = nan(nH,1);
act    = true(nH,1);

bad0 = ~all(isfinite(D0),2) | isSingular(st0, isLoC, opt);
status(bad0) = {'seed_invalid'};
act(bad0) = false;

stride = nH*opt.maxSteps;
while any(act)
    ia = find(act);
    xa = x(ia,:); Ta = T(ia,:); ha = h(ia);

    q1 = alignedDir(fld, xa,                 family, Ta);
    q2 = alignedDir(fld, xa + 0.5*ha.*q1,    family, Ta);
    q3 = alignedDir(fld, xa + 0.5*ha.*q2,    family, Ta);
    q4 = alignedDir(fld, xa + ha.*q3,        family, Ta);
    xn = xa + (ha/6).*(q1 + 2*q2 + 2*q3 + q4);
    fin = all(isfinite(xn),2);
    xn(~fin,:) = xa(~fin,:);
    [xn, res] = tpms_project(fld, xn, t, 1e-13, 10);
    [Dn, stn] = fieldDir(fld, xn, family);
    Tn   = alignTo(Dn, Ta);
    cosT = max(-1, min(1, sum(Ta.*Tn, 2)));
    turn = acos(cosT)*180/pi;

    rej = ~fin | ~all(isfinite(Tn),2) | ~(turn <= opt.maxTurnDeg) | ~(res < opt.projTol);

    % ---- rejected: halve the step ------------------------------------------
    ir = ia(rej);
    h(ir) = h(ir)/2;
    under = ir(h(ir) < opt.hmin);
    status(under) = {'step_underflow'};
    act(under) = false;

    % ---- accepted: store node ----------------------------------------------
    lc = find(~rej);
    ic = ia(lc);
    if isempty(ic), continue; end
    x(ic,:) = xn(lc,:);
    T(ic,:) = Tn(lc,:);
    s(ic)   = s(ic) + ha(lc);
    cnt(ic) = cnt(ic) + 1;
    lin = sub2ind([nH opt.maxSteps], ic, cnt(ic));
    for d = 1:3
        P(lin + (d-1)*stride) = xn(lc,d);
    end
    S(lin) = s(ic);
    grow = turn(lc) < opt.maxTurnDeg/4;
    h(ic(grow)) = min(2*h(ic(grow)), opt.hmax);

    % ---- stop checks -------------------------------------------------------
    stc.gap = stn.gap(lc); stc.K = stn.K(lc);
    sing = isSingular(stc, isLoC, opt);
    if isLoC, lbl = 'umbilic'; else, lbl = 'parabolic'; end
    status(ic(sing)) = {lbl};

    % closure: distance from the seed (modulo the lattice) to the whole step
    % segment, so detection does not depend on where the nodes happen to land
    dx    = x(ic,:) - X0(ic,:);
    shift = L.*round(dx./L);
    dprev = xa(lc,:) - X0(ic,:) - shift;
    seg   = x(ic,:) - xa(lc,:);
    tseg  = max(0, min(1, -sum(dprev.*seg, 2) ./ max(sum(seg.^2, 2), eps)));
    dw    = sqrt(sum((dprev + tseg.*seg).^2, 2));
    near  = dw < opt.closeTol & s(ic) > opt.minLoop ...
            & sum(T(ic,:).*T0(ic,:), 2) > 0.98 & ~sing;
    status(ic(near)) = {'closed'};
    trans(ic(near),:) = round(dx(near,:)./L);
    cgap(ic(near)) = dw(near);

    longE = s(ic) >= opt.Lmax & ~sing & ~near;
    status(ic(longE)) = {'max_length'};
    full = cnt(ic) >= opt.maxSteps & ~sing & ~near & ~longE;
    status(ic(full)) = {'max_steps'};

    act(ic(sing | near | longE | full)) = false;
end

% ---- merge the two half-curves of every seed -------------------------------
C = struct('family', {}, 'seed', {}, 'x', {}, 's', {}, 'statusStart', {}, ...
           'statusEnd', {}, 'closed', {}, 'translation', {}, 'closeGap', {}, 'length', {});
for m = 1:M
    f = m; b = m + M;
    xf = squeeze(P(f, 1:cnt(f), :)); sf = S(f, 1:cnt(f)).';
    xb = squeeze(P(b, 1:cnt(b), :)); sb = S(b, 1:cnt(b)).';
    if cnt(f) == 1, xf = reshape(xf, 1, 3); end
    if cnt(b) == 1, xb = reshape(xb, 1, 3); end
    c.family = family; c.seed = seeds(m,:);
    if strcmp(status{f}, 'closed') || strcmp(status{b}, 'closed')
        if strcmp(status{f}, 'closed'), k = f; xx = xf; ss = sf;
        else,                           k = b; xx = xb; ss = sb; end
        c.x = xx; c.s = ss;
        c.statusStart = 'closed'; c.statusEnd = 'closed';
        c.closed = true; c.translation = trans(k,:); c.closeGap = cgap(k);
    else
        c.x = [flipud(xb(2:end,:)); xf];
        c.s = [-flipud(sb(2:end)); sf];
        c.s = c.s - c.s(1);
        c.statusStart = status{b}; c.statusEnd = status{f};
        c.closed = false; c.translation = [NaN NaN NaN]; c.closeGap = NaN;
    end
    c.length = c.s(end);
    C(end+1) = c; %#ok<AGROW>
end
end

% -------------------------------------------------------------------------
function [D, st] = fieldDir(fld, X, family)
G = tpms_geometry(fld, X, false);
switch family
    case 'k1',    D = G.d1;
    case 'k2',    D = G.d2;
    case 'asym+', D = G.ap;
    case 'asym-', D = G.am;
    otherwise, error('tpms_trace_curves:family', 'Unknown family ''%s''.', family);
end
st.gap = G.gap; st.K = G.K;
end

function D = alignedDir(fld, X, family, Tref)
D = alignTo(fieldDir(fld, X, family), Tref);
end

function D = alignTo(D, Tref)
sgn = sign(sum(D.*Tref, 2));
sgn(sgn == 0) = 1;
D = sgn .* D;
end

function tf = isSingular(st, isLoC, opt)
if isLoC
    tf = ~(st.gap >= opt.umbTol);
else
    tf = ~(st.K <= -opt.parTol);
end
end
