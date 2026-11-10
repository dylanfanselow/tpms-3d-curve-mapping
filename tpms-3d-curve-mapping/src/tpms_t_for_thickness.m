function [t, info] = tpms_t_for_thickness(targetThickness, cellSize, nGrid)
%TPMS_T_FOR_THICKNESS  Level-set value t that gives a target sheet wall thickness.
%
%   t = tpms_t_for_thickness(0.4, 5)        % 0.4 mm walls in a 5 mm cell
%   [t, info] = tpms_t_for_thickness(...)    % info: volume fraction, S/V, ...
%
%   Uses the exact (along-normal) mean thickness of tpms_wall_thickness,
%   interpolated on a fine t grid. Thickness and cellSize in the same unit.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.
if nargin < 3 || isempty(nGrid), nGrid = 48; end
fld = tpms_build_fields('gyroid', 'manual');
cfg.cellSize = cellSize; cfg.nGrid = nGrid; cfg.vfGrid = 96;
tg = 0:0.02:1.3;
TT = tpms_wall_thickness(fld, tg, cfg);
th = [TT.thick_mean];
if targetThickness < 0 || targetThickness > max(th)
    error('tpms_t_for_thickness:range', ...
        'Target %.3g is outside the achievable range 0 .. %.3g for cell size %.3g.', ...
        targetThickness, max(th), cellSize);
end
t = interp1(th, tg, targetThickness, 'pchip');
info.t = t;
info.thickness_over_cell = targetThickness/cellSize;
info.vf_sheet = interp1(tg, [TT.vf_sheet], t, 'pchip');
info.SV_sheet = interp1(tg(2:end), [TT(2:end).SV_sheet], t, 'pchip');
info.table = TT;
end
