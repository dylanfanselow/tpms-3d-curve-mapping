function idx = tpms_seed_points(v, M, allowed, firstIdx, period)
%TPMS_SEED_POINTS  Farthest-point sampling of M seeds among allowed vertices.
%   Periodic (minimum-image) distances on the 2*pi cell are used so seeds
%   do not crowd at the cell faces. Deterministic for a given firstIdx.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.
if nargin < 5 || isempty(period), L = 2*pi; else, L = tpms_period3(period); end
cand = find(allowed);
if isempty(cand), idx = []; return; end
M = min(M, numel(cand));
if nargin < 4 || isempty(firstIdx), firstIdx = cand(1); end
idx = zeros(M,1); idx(1) = firstIdx;
dmin = inf(numel(cand),1);
for k = 2:M
    d = v(cand,:) - v(idx(k-1),:);
    d = d - L.*round(d./L);
    dmin = min(dmin, sum(d.^2,2));
    [~, j] = max(dmin);
    idx(k) = cand(j);
end
end
