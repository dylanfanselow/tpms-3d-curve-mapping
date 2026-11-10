function fld = tpms_build_fields(tpmsType, method)
%TPMS_BUILD_FIELDS  Analytic derivative fields of a TPMS level-set function.
%
%   fld = tpms_build_fields('gyroid')             % Symbolic Toolbox + matlabFunction
%   fld = tpms_build_fields('gyroid', 'manual')   % hand-derived (independent cross-check)
%
%   The surface is the level set  F0(x,y,z) = t  (t is supplied later, so one
%   field struct serves every level-set value). Coordinates are in radians;
%   one unit cell is [0, 2*pi]^3.
%
%   Every handle takes X (N-by-3) and returns an N-by-m array:
%     fld.F0(X)    N-by-1
%     fld.grad(X)  N-by-3    [Fx Fy Fz]
%     fld.hess(X)  N-by-6    [Fxx Fyy Fzz Fxy Fyz Fxz]
%     fld.d3(X)    N-by-10   [Fxxx Fyyy Fzzz Fxxy Fxxz Fxyy Fyyz Fxzz Fyzz Fxyz]
%
%   Third derivatives are needed because the geodesic curvature of lines of
%   curvature depends on surface derivatives of k1, k2 (Codazzi equations).
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.

if nargin < 2 || isempty(method), method = 'symbolic'; end
tpmsType = lower(tpmsType);
if ~strcmp(tpmsType, 'gyroid')
    error('tpms_build_fields:type', ...
        'Only ''gyroid'' is implemented in this iteration (got ''%s'').', tpmsType);
end

fld.type   = tpmsType;
fld.method = lower(method);
fld.period = 2*pi;
fld.expr   = 'sin(x)*cos(y) + sin(y)*cos(z) + sin(z)*cos(x)';

switch fld.method
    case 'symbolic'
        syms x y z real
        F = sin(x)*cos(y) + sin(y)*cos(z) + sin(z)*cos(x);
        v = [x y z];
        g = [diff(F,x), diff(F,y), diff(F,z)];
        H = [diff(g(1),x), diff(g(2),y), diff(g(3),z), ...
             diff(g(1),y), diff(g(2),z), diff(g(1),z)];
        idx3 = [1 1 1; 2 2 2; 3 3 3; 1 1 2; 1 1 3; 1 2 2; 2 2 3; 1 3 3; 2 3 3; 1 2 3];
        T3 = sym(zeros(1,10));
        for k = 1:10
            T3(k) = diff(diff(diff(F, v(idx3(k,1))), v(idx3(k,2))), v(idx3(k,3)));
        end
        fld.F0   = makeEvaluator(F,  x, y, z);
        fld.grad = makeEvaluator(g,  x, y, z);
        fld.hess = makeEvaluator(H,  x, y, z);
        fld.d3   = makeEvaluator(T3, x, y, z);

    case 'manual'
        fld.F0   = @gyrF;
        fld.grad = @gyrGrad;
        fld.hess = @gyrHess;
        fld.d3   = @gyrD3;

    otherwise
        error('tpms_build_fields:method', 'Unknown method ''%s''.', method);
end
end

% -------------------------------------------------------------------------
function ev = makeEvaluator(exprVec, x, y, z)
% Build one numeric handle per component. Components that simplify to a
% constant (e.g. Fxyz = 0 for the gyroid) return a scalar from
% matlabFunction, so each is expanded to N-by-1 at evaluation time.
m  = numel(exprVec);
hc = cell(1, m);
for k = 1:m
    hc{k} = matlabFunction(exprVec(k), 'Vars', {x, y, z});
end
ev = @(X) evalHandles(hc, X);
end

function out = evalHandles(hc, X)
N   = size(X,1);
out = zeros(N, numel(hc));
for k = 1:numel(hc)
    out(:,k) = hc{k}(X(:,1), X(:,2), X(:,3)) + zeros(N,1);
end
end

% ---- hand-derived gyroid derivatives (independent of the Symbolic path) ----
function f = gyrF(X)
sx = sin(X(:,1)); cx = cos(X(:,1)); sy = sin(X(:,2)); cy = cos(X(:,2));
sz = sin(X(:,3)); cz = cos(X(:,3));
f = sx.*cy + sy.*cz + sz.*cx;
end

function g = gyrGrad(X)
sx = sin(X(:,1)); cx = cos(X(:,1)); sy = sin(X(:,2)); cy = cos(X(:,2));
sz = sin(X(:,3)); cz = cos(X(:,3));
g = [cx.*cy - sx.*sz, cy.*cz - sx.*sy, cz.*cx - sy.*sz];
end

function H = gyrHess(X)
sx = sin(X(:,1)); cx = cos(X(:,1)); sy = sin(X(:,2)); cy = cos(X(:,2));
sz = sin(X(:,3)); cz = cos(X(:,3));
H = [-sx.*cy - sz.*cx, -sx.*cy - sy.*cz, -sy.*cz - sz.*cx, ...
     -cx.*sy,          -cy.*sz,          -sx.*cz];
end

function T = gyrD3(X)
sx = sin(X(:,1)); cx = cos(X(:,1)); sy = sin(X(:,2)); cy = cos(X(:,2));
sz = sin(X(:,3)); cz = cos(X(:,3));
T = [ sx.*sz - cx.*cy, ...   % Fxxx
      sx.*sy - cy.*cz, ...   % Fyyy
      sy.*sz - cz.*cx, ...   % Fzzz
      sx.*sy, ...            % Fxxy
     -cz.*cx, ...            % Fxxz
     -cx.*cy, ...            % Fxyy
      sy.*sz, ...            % Fyyz
      sx.*sz, ...            % Fxzz
     -cy.*cz, ...            % Fyzz
      zeros(size(sx))];      % Fxyz
end
