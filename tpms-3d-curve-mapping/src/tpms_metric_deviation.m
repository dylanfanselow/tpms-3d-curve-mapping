function [dev, du, ok, hw] = tpms_metric_deviation(m, m0, study)
%TPMS_METRIC_DEVIATION  Deviation of a summary metric from a reference value.
%   [dev, unit, ok] = tpms_metric_deviation(m, m0, study)
%     bounded metrics (0-1, -1..+1): absolute difference, tolerance study.tolAbs
%     percentages: difference in points, tolerance study.tolPct
%     everything else: relative %, tolerance study.tolRel
%   hw = half-width of m's 95%% CI in the same units as dev (0 if no CI)
%   d = tpms_metric_deviation(row, [], study)  -> |dev| / tolerance for a
%   sampling-study row (fields dev, devUnit), floored at 1e-3 for log plots.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.
if isempty(m0)
    switch m.devUnit
        case 'abs',    dev = abs(m.dev)/study.tolAbs;
        case 'points', dev = abs(m.dev)/study.tolPct;
        otherwise,     dev = abs(m.dev)/study.tolRel;
    end
    if ~isfinite(dev), dev = 1e-3; end
    dev = max(dev, 1e-3); du = ''; ok = true; hw = 0;
    return
end
hwr = (m.hi - m.lo)/2;
if ~strcmp(m.rangeType, 'CI95') || ~isfinite(hwr), hwr = 0; end
if any(strcmp(m.unit, {'0-1', '-1..+1'}))
    dev = m.value - m0.value; du = 'abs'; ok = abs(dev) <= study.tolAbs; hw = hwr;
elseif strcmp(m.unit, '%')
    dev = m.value - m0.value; du = 'points'; ok = abs(dev) <= study.tolPct; hw = hwr;
else
    dev = 100*(m.value - m0.value)/abs(m0.value); du = '%'; ok = abs(dev) <= study.tolRel;
    hw = 100*hwr/abs(m0.value);
end
cm = strcmp(m.rangeType, 'LOWER'); c0 = strcmp(m0.rangeType, 'LOWER');
if cm || c0                       % censored (longer than max lag): agree only if both are
    ok = (cm == c0); if ~ok, dev = NaN; end
end
if ~isfinite(dev) && (isfinite(m.value) == isfinite(m0.value)) && ~(cm ~= c0), ok = true; end
end
