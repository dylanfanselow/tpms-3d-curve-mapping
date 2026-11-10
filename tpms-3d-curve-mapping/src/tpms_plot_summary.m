function tpms_plot_summary(Ms, tList, thickFine, cfg)
%TPMS_PLOT_SUMMARY  Trend figure (each summary metric vs t, with 95% CI)
%   and wall-thickness figure. Light style regardless of MATLAB theme.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.
keys = {Ms{1}.key}; n = numel(keys);
nc = 4; nr = ceil(n/nc);
fg = newFig(cfg, [50 50 1700 380*nr]);
for k = 1:n
    v = nan(1,numel(tList)); lo = v; hi = v; m0 = [];
    for j = 1:numel(tList)
        i = find(strcmp({Ms{j}.key}, keys{k}), 1); if isempty(i), continue; end
        m = Ms{j}(i); if isempty(m0), m0 = m; end
        v(j) = m.value; lo(j) = m.lo; hi(j) = m.hi;
    end
    ax = subplot(nr, nc, k, 'Parent', fg); hold(ax, 'on');
    if any(isfinite(lo)) && strcmp(m0.rangeType, 'CI95')
        e = errorbar(ax, tList, v, v - lo, hi - v, 'o-');
    elseif any(isfinite(lo))
        e = errorbar(ax, tList, v, v - lo, hi - v, 's-');
    else
        e = plot(ax, tList, v, 'o-');
    end
    set(e, 'Color', [0.05 0.30 0.85], 'LineWidth', 1.8, 'MarkerFaceColor', [0.05 0.30 0.85]);
    title(ax, wrapLabel(m0.label)); xlabel(ax, 't'); ylabel(ax, m0.unit);
    styleAxes(ax);
end
print(fg, fullfile(cfg.outDir, 'summary_trends.png'), '-dpng', '-r130');
if strcmp(cfg.figVisible, 'off'), close(fg); end

% wall thickness lookup
if isempty(cfg.cellSize), uL = 'cell lengths'; f = 1/(2*pi); else, uL = cfg.lengthUnit; f = 1; end
tg = [thickFine.t]; rel = [thickFine.thick_fail_frac] < 0.01;
fg = newFig(cfg, [50 50 1300 520]);
ax = subplot(1,2,1, 'Parent', fg); hold(ax, 'on');
fill_ = [tg(rel), fliplr(tg(rel))];
patch(ax, fill_, f*[[thickFine(rel).thick_p5], fliplr([thickFine(rel).thick_p95])], [0.75 0.82 0.95], 'EdgeColor', 'none');
h1 = plot(ax, tg(rel), f*[thickFine(rel).thick_mean], '-', 'Color', [0.05 0.30 0.85], 'LineWidth', 2.2);
h2 = plot(ax, tg(rel), f*[thickFine(rel).thick_firstorder_mean], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
xlabel(ax, 't'); ylabel(ax, ['wall thickness (' uL ')']);
if isempty(cfg.cellSize), ttl = 'Sheet wall thickness vs t';
else, ttl = sprintf('Sheet wall thickness vs t (cell %g %s)', cfg.cellSize, cfg.lengthUnit); end
title(ax, ttl); styleAxes(ax);
styleLegend(legend(ax, [h1 h2], {'exact (band = 5-95%)', 'first-order 2t/|\nablaF|'}, 'Location', 'northwest'));
ax = subplot(1,2,2, 'Parent', fg); hold(ax, 'on');
h1 = plot(ax, tg, 100*[thickFine.vf_sheet], '-', 'Color', [0.85 0.10 0.10], 'LineWidth', 2.2);
h2 = plot(ax, tg, 100*[thickFine.vf_network], '-', 'Color', [0.00 0.55 0.15], 'LineWidth', 2.2);
xlabel(ax, 't'); ylabel(ax, 'solid volume fraction (%)'); title(ax, 'Solid fraction vs t'); styleAxes(ax);
styleLegend(legend(ax, [h1 h2], {'sheet', 'network'}, 'Location', 'northwest'));
print(fg, fullfile(cfg.outDir, 'wall_thickness_vs_t.png'), '-dpng', '-r130');
if strcmp(cfg.figVisible, 'off'), close(fg); end
end

% -------------------------------------------------------------------------
function s = wrapLabel(s)
i = strfind(s, ': ');
if ~isempty(i), s = {s(1:i(1)-1), s(i(1)+2:end)}; end
end
function fg = newFig(cfg, pos)
fg = figure('Visible', cfg.figVisible, 'Color', 'w', 'Position', pos, 'InvertHardcopy', 'off');
try, theme(fg, 'light'); catch, end
end
function styleAxes(ax)
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', ...
        'FontSize', 13, 'LineWidth', 1, 'Box', 'on');
set(get(ax, 'Title'),  'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold');
set(get(ax, 'XLabel'), 'Color', 'k', 'FontSize', 14);
set(get(ax, 'YLabel'), 'Color', 'k', 'FontSize', 14);
end
function styleLegend(lg)
set(lg, 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k', 'FontSize', 12);
end
