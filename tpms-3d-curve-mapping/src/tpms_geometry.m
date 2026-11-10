function G = tpms_geometry(fld, X, wantDerivs)
%TPMS_GEOMETRY  Analytic differential geometry of the level set at points X.
%
%   G = tpms_geometry(fld, X)          curvatures and direction fields
%   G = tpms_geometry(fld, X, true)    also curvature derivatives (uses d3)
%
%   X (N-by-3) should lie on the surface (see tpms_project).
%
%   CONVENTIONS
%     n  = grad F / |grad F|                (unit normal, fixed global orientation)
%     S  = -P*HessF*P / |grad F|            (shape operator, S = -dn)
%     kappa_n(u) = <gamma'', n> = -u'*HessF*u/|grad F|   for a unit tangent u
%     k1 >= k2 are the eigenvalues of S on the tangent plane
%     d1 = principal direction of k1,  d2 = n x d1   (so (d1, d2, n) is right-handed)
%     H = (k1+k2)/2,  K = k1*k2
%   With this convention a sphere F = |x|^2 - R^2 has k1 = k2 = -1/R.
%   Signs of k1, k2, H flip with the normal orientation; K does not.
%
%   OUTPUT (all N-by-1 unless noted)
%     g (N-by-3), gn = |grad F|, n (N-by-3), Hs (N-by-6)
%     k1, k2, H, K, gap = k1 - k2
%     d1, d2 (N-by-3)
%     K_closed, H_closed      closed-form checks (independent of eigen solve)
%     thA                     asymptotic angle from d1, where K < 0 (NaN elsewhere)
%     ap, am (N-by-3)         asymptotic directions  cos(thA)*d1 +/- sin(thA)*d2
%   with wantDerivs:
%     dk (N-by-4)   [e1(k1) e2(k1) e1(k2) e2(k2)]   directional derivatives
%     kg1, kg2      geodesic curvature of the k1- and k2-line fields (Codazzi)
%     kgAp, kgAm    geodesic curvature of the asymptotic fields (Liouville)
%     gradS_H, gradS_K (N-by-3)  surface gradients of H and K
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.

if nargin < 3, wantDerivs = false; end
N  = size(X,1);
g  = fld.grad(X);
Hs = fld.hess(X);
gn = sqrt(sum(g.^2, 2));
n  = g ./ gn;

% --- tangent basis and 2x2 shape operator ---------------------------------
[t1, t2] = tangentBasis(n);
a = -quadH(Hs, t1, t1) ./ gn;
b = -quadH(Hs, t1, t2) ./ gn;
c = -quadH(Hs, t2, t2) ./ gn;
m = (a + c)/2;
r = sqrt(((a - c)/2).^2 + b.^2);
k1 = m + r;
k2 = m - r;
th = 0.5*atan2(2*b, a - c);
d1 = cos(th).*t1 + sin(th).*t2;
d2 = cross(n, d1, 2);

G.X = X; G.g = g; G.gn = gn; G.n = n; G.Hs = Hs;
G.k1 = k1; G.k2 = k2; G.H = m; G.K = k1.*k2; G.gap = k1 - k2;
G.d1 = d1; G.d2 = d2;

% --- closed-form checks -----------------------------------------------------
xx = Hs(:,1); yy = Hs(:,2); zz = Hs(:,3); xy = Hs(:,4); yz = Hs(:,5); xz = Hs(:,6);
A11 = yy.*zz - yz.^2;  A22 = xx.*zz - xz.^2;  A33 = xx.*yy - xy.^2;
A12 = xz.*yz - xy.*zz; A13 = xy.*yz - xz.*yy; A23 = xy.*xz - xx.*yz;
gx = g(:,1); gy = g(:,2); gz = g(:,3);
gAg = A11.*gx.^2 + A22.*gy.^2 + A33.*gz.^2 + 2*(A12.*gx.*gy + A13.*gx.*gz + A23.*gy.*gz);
G.K_closed = gAg ./ gn.^4;
G.H_closed = (quadH(Hs, g, g) - gn.^2 .* (xx + yy + zz)) ./ (2*gn.^3);

% --- asymptotic directions (exist only where K < 0) -------------------------
hyp = (k1 > 0) & (k2 < 0);
thA = nan(N,1);
thA(hyp) = atan2(sqrt(k1(hyp)), sqrt(-k2(hyp)));
G.thA = thA;
G.ap  = cos(thA).*d1 + sin(thA).*d2;
G.am  = cos(thA).*d1 - sin(thA).*d2;

if ~wantDerivs, return; end

% --- directional derivatives of k1, k2 (Hellmann-Feynman on S) --------------
% e_v(k_i) = -(n'Hv) k_i/|g| + 2 (d_i'Hv)(n'H d_i)/|g|^2 - T3(d_i,d_i,v)/|g|
D3 = fld.d3(X);
Hd1 = hmul(Hs, d1); Hd2 = hmul(Hs, d2);
nHd1 = sum(n.*Hd1, 2); nHd2 = sum(n.*Hd2, 2);
ek = @(k, d, Hd, nHd, v, nHv) -nHv.*k./gn + 2*sum(Hd.*v,2).*nHd./gn.^2 ...
                              - t3contract(D3, d, d, v)./gn;
dk = [ek(k1, d1, Hd1, nHd1, d1, nHd1), ...   % e1(k1)
      ek(k1, d1, Hd1, nHd1, d2, nHd2), ...   % e2(k1)
      ek(k2, d2, Hd2, nHd2, d1, nHd1), ...   % e1(k2)
      ek(k2, d2, Hd2, nHd2, d2, nHd2)];      % e2(k2)
G.dk  = dk;

% Codazzi: geodesic curvature of lines of curvature (tangent +d1, resp. +d2)
G.kg1 = dk(:,2) ./ G.gap;
G.kg2 = dk(:,3) ./ G.gap;

gradS_k1 = dk(:,1).*d1 + dk(:,2).*d2;
gradS_k2 = dk(:,3).*d1 + dk(:,4).*d2;
G.gradS_H = 0.5*(gradS_k1 + gradS_k2);
G.gradS_K = k2.*gradS_k1 + k1.*gradS_k2;

% Liouville: kappa_g = dtheta/ds + cos(theta)*kg1 + sin(theta)*kg2
% theta = +thA (family a+), -thA (family a-);  tan(thA) = sqrt(k1/(-k2))
cA = cos(thA); sA = sin(thA);
p  = sqrt(-k2./k1); q = sqrt(k1./(-k2));        % NaN/Inf outside K<0: masked later
dthA = @(ek1, ek2) (p.*ek1 + q.*ek2) ./ (2*G.gap);
ekp1 = cA.*dk(:,1) + sA.*dk(:,2);  ekp2 = cA.*dk(:,3) + sA.*dk(:,4);   % along a+
ekm1 = cA.*dk(:,1) - sA.*dk(:,2);  ekm2 = cA.*dk(:,3) - sA.*dk(:,4);   % along a-
G.kgAp =  dthA(ekp1, ekp2) + cA.*G.kg1 + sA.*G.kg2;
G.kgAm = -dthA(ekm1, ekm2) + cA.*G.kg1 - sA.*G.kg2;
G.kgAp(~hyp) = NaN; G.kgAm(~hyp) = NaN;
end

% -------------------------------------------------------------------------
function [t1, t2] = tangentBasis(n)
N = size(n,1);
[~, imin] = min(abs(n), [], 2);
E = zeros(N,3);
E(sub2ind([N 3], (1:N)', imin)) = 1;
t1 = cross(n, E, 2);
t1 = t1 ./ sqrt(sum(t1.^2, 2));
t2 = cross(n, t1, 2);
end

function Hu = hmul(Hs, u)
% Hs = [xx yy zz xy yz xz]
Hu = [Hs(:,1).*u(:,1) + Hs(:,4).*u(:,2) + Hs(:,6).*u(:,3), ...
      Hs(:,4).*u(:,1) + Hs(:,2).*u(:,2) + Hs(:,5).*u(:,3), ...
      Hs(:,6).*u(:,1) + Hs(:,5).*u(:,2) + Hs(:,3).*u(:,3)];
end

function q = quadH(Hs, u, v)
q = sum(u .* hmul(Hs, v), 2);
end

function s = t3contract(D3, u, v, w)
% sum_ijk T_ijk u_i v_j w_k for the symmetric 3rd-derivative tensor stored as
% [xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz]
persistent M
if isempty(M)
    M = zeros(3,3,3);
    lst = [1 1 1; 2 2 2; 3 3 3; 1 1 2; 1 1 3; 1 2 2; 2 2 3; 1 3 3; 2 3 3; 1 2 3];
    for i = 1:3, for j = 1:3, for k = 1:3
        key = sort([i j k]);
        M(i,j,k) = find(all(lst == key, 2));
    end, end, end
end
s = zeros(size(u,1), 1);
for i = 1:3, for j = 1:3, for k = 1:3
    s = s + D3(:, M(i,j,k)) .* u(:,i) .* v(:,j) .* w(:,k);
end, end, end
end
