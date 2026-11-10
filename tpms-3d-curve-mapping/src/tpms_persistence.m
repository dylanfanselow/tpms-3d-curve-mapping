function ell = tpms_persistence(lags, aSum, aCnt)
%TPMS_PERSISTENCE  Arc length at which the mean tangent autocorrelation
%   <T(s).T(s+lag)> first drops below 1/e ("direction memory" of a curve
%   family). Returns NaN if it never drops within the largest lag.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.
C = aSum ./ max(aCnt, 1);
C(aCnt == 0) = NaN;
k = find(C < exp(-1), 1);
if isempty(k) || k == 1, ell = NaN; return; end
ell = lags(k-1) + (C(k-1) - exp(-1)) * (lags(k) - lags(k-1)) / (C(k-1) - C(k));
end
