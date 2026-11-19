function S = tpms_surface_sample(fld, t, nGrid)
%TPMS_SURFACE_SAMPLE  Quadrature points on one unit cell of F0 = t.
%
%   S = tpms_surface_sample(fld, t, nGrid)
%
%   A marching-cubes mesh (isosurface) is used ONLY to obtain well-spread
%   sample points and area weights. Every vertex is then Newton-projected
%   onto the exact level set, and all geometry is evaluated analytically at
%   the projected vertices -- the mesh never supplies curvature.
%
%   S.v (N-by-3) vertices on the surface, S.f faces, S.w (N-by-1) vertex area
%   weights (1/3 of incident triangle area), S.A total area in one cell,
%   S.projRes projection residuals.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.
if nargin < 3 || isempty(nGrid), nGrid = 64; end
P  = tpms_period3(fld.period);                 % [2pi 2pi 2pi] for a cubic cell
ng = round(nGrid * P / (2*pi));               % = nGrid on every axis for a cube
[X, Y, Z] = meshgrid(linspace(0, P(1), ng(1)), linspace(0, P(2), ng(2)), linspace(0, P(3), ng(3)));
V = reshape(fld.F0([X(:) Y(:) Z(:)]), size(X));
[f, v] = isosurface(X, Y, Z, V, t);
if isempty(f)
    error('tpms_surface_sample:empty', 'No surface for t = %g.', t);
end
[v, res] = tpms_project(fld, v, t, 1e-13, 20);
e1 = v(f(:,2),:) - v(f(:,1),:);
e2 = v(f(:,3),:) - v(f(:,1),:);
fa = 0.5*sqrt(sum(cross(e1, e2, 2).^2, 2));
w  = accumarray(f(:), repmat(fa/3, 3, 1), [size(v,1) 1]);
keep = w > 0;                              % drop unreferenced vertices
map = zeros(size(v,1),1); map(keep) = 1:nnz(keep);
S.v = v(keep,:); S.w = w(keep); S.f = map(f);
S.A = sum(fa); S.projRes = res(keep); S.t = t; S.nGrid = nGrid;
end
