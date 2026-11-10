function [X, res] = tpms_project(fld, X, t, tol, maxIt)
%TPMS_PROJECT  Newton projection of points onto the level set F0 = t.
%   [X, res] = tpms_project(fld, X, t)   res = |F0 - t| / |grad F| (first-order distance)
%
%   Update:  X <- X - (F0(X) - t) * grad F / |grad F|^2   (along the normal).
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.
if nargin < 4 || isempty(tol),   tol = 1e-12; end
if nargin < 5 || isempty(maxIt), maxIt = 20;  end
for it = 1:maxIt
    f  = fld.F0(X) - t;
    g  = fld.grad(X);
    g2 = sum(g.^2, 2);
    res = abs(f) ./ sqrt(g2);
    if all(res < tol), return; end
    X = X - (f ./ g2) .* g;
end
f   = fld.F0(X) - t;
res = abs(f) ./ sqrt(sum(fld.grad(X).^2, 2));
end
