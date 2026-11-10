function tpms_plot_results(R, cfg)
%TPMS_PLOT_RESULTS  Figures for one level set: curve networks and distributions.
%   Saved as PNG in cfg.outDir. Curves are wrapped into the [0, 2*pi]^3 cell.
%   All figures are forced to a light style (white background, black text),
%   so they read the same whether MATLAB is in light or dark theme.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.
tag = sprintf('t%+.2f', R.t);
S = R.S; G = R.G;
fams = fieldnames(R.curves);
col = struct('k1',[0.85 0.10 0.10], 'k2',[0.05 0.30 0.85], ...
             'asymP',[0.00 0.55 0.15], 'asymM',[0.55 0.10 0.65]);
lbl = struct('k1','k1-lines', 'k2','k2-lines', 'asymP','asym+', 'asymM','asym-');

% 1-2: curve networks on the surface ----------------------------------------------
groups = {{'k1','k2'}, 'lines_of_curvature'; {'asymP','asymM'}, 'asymptotic_curves'};
for g = 1:size(groups,1)
    fg = newFig(cfg, [100 100 950 850]);
    ax = axes('Parent', fg);
    patch(ax, 'Faces', S.f, 'Vertices', S.v, 'FaceColor', [0.80 0.80 0.82], ...
          'EdgeColor', 'none', 'FaceAlpha', 0.35);
    hold(ax, 'on')
    hl = {}; nm = {};
    for k = 1:numel(groups{g,1})
        key = groups{g,1}{k};
        if ~isfield(R.curves, key), continue; end
        C = R.curves.(key);
        if isempty(C), continue; end
        for m = 1:numel(C)
            xw = wrapCurve(C(m).x);
            h = plot3(ax, xw(:,1), xw(:,2), xw(:,3), '-', 'Color', col.(key), 'LineWidth', 1.2);
        end
        hl{end+1} = h; nm{end+1} = lbl.(key); %#ok<AGROW>
    end
    axis(ax, 'equal'); axis(ax, [0 2*pi 0 2*pi 0 2*pi]); view(ax, 35, 25);
    camlight(ax, 'headlight'); lighting(ax, 'gouraud');
    xlabel(ax, 'x'); ylabel(ax, 'y'); zlabel(ax, 'z');
    title(ax, sprintf('Gyroid t = %.2f: %s', R.t, strrep(groups{g,2}, '_', ' ')));
    styleAxes(ax);
    if ~isempty(hl), styleLegend(legend(ax, [hl{:}], nm, 'Location', 'northeast')); end
    saveFig(fg, cfg, sprintf('%s_%s.png', groups{g,2}, tag));
end

% 3: surface coloured by |kappa_g| of the k1-line field ------------------------------
fg = newFig(cfg, [100 100 950 850]);
ax = axes('Parent', fg);
c = abs(G.kg1); c(G.gap < R.opt.umbExcl) = NaN;
patch(ax, 'Faces', S.f, 'Vertices', S.v, 'FaceVertexCData', c, 'FaceColor', 'interp', 'EdgeColor', 'none');
axis(ax, 'equal'); axis(ax, [0 2*pi 0 2*pi 0 2*pi]); view(ax, 35, 25);
caxis(ax, [0 prctile_(c, 95)]);
xlabel(ax, 'x'); ylabel(ax, 'y'); zlabel(ax, 'z');
title(ax, {sprintf('Gyroid t = %.2f', R.t), '|geodesic curvature| of k1-lines'});
styleAxes(ax);
cb = colorbar(ax);
set(cb, 'Color', 'k', 'FontSize', 14, 'LineWidth', 1);
ylabel(cb, '|\kappa_g|', 'Color', 'k', 'FontSize', 15);
saveFig(fg, cfg, sprintf('field_kg_k1lines_%s.png', tag));

% 4: distributions ------------------------------------------------------------------
fg = newFig(cfg, [100 100 1650 560]);
ax = subplot(1,3,1, 'Parent', fg); hold(ax, 'on')
[hl, nm] = deal({}, {});
for k = 1:numel(fams)
    d = findDesc(R, ['curve_' fams{k} '_kg']);
    if ~isempty(d), hl{end+1} = histLine(ax, d.values, col.(fams{k})); nm{end+1} = lbl.(fams{k}); end %#ok<AGROW>
end
title(ax, '|\kappa_g| along curves'); xlabel(ax, '|\kappa_g|'); ylabel(ax, 'fraction of nodes');
styleAxes(ax); if ~isempty(hl), styleLegend(legend(ax, [hl{:}], nm, 'Location', 'northeast')); end

ax = subplot(1,3,2, 'Parent', fg); hold(ax, 'on')
[hl, nm] = deal({}, {});
for k = 1:numel(fams)
    d = findDesc(R, ['curve_' fams{k} '_tort']);
    if ~isempty(d), hl{end+1} = histLine(ax, d.values, col.(fams{k})); nm{end+1} = lbl.(fams{k}); end %#ok<AGROW>
end
title(ax, sprintf('Windowed tortuosity (W = %.2f)', R.opt.tortWindow));
xlabel(ax, 'W / chord'); ylabel(ax, 'fraction of windows');
styleAxes(ax); if ~isempty(hl), styleLegend(legend(ax, [hl{:}], nm, 'Location', 'north')); end

ax = subplot(1,3,3, 'Parent', fg); hold(ax, 'on')
[hl, nm] = deal({}, {});
for k = {'k1','k2'}
    d = findDesc(R, ['curve_' k{1} '_tau']);
    if ~isempty(d), hl{end+1} = histLine(ax, d.values, col.(k{1})); nm{end+1} = lbl.(k{1}); end %#ok<AGROW>
end
title(ax, 'Frenet torsion, lines of curvature'); xlabel(ax, '\tau'); ylabel(ax, 'fraction of nodes');
styleAxes(ax); if ~isempty(hl), styleLegend(legend(ax, [hl{:}], nm, 'Location', 'northeast')); end
saveFig(fg, cfg, sprintf('distributions_%s.png', tag));
end

% -------------------------------------------------------------------------
function fg = newFig(cfg, pos)
fg = figure('Visible', cfg.figVisible, 'Color', 'w', 'Position', pos, 'InvertHardcopy', 'off');
try, theme(fg, 'light'); catch, end          % R2025a+ dark theme -> light (no-op on older releases)
end

function styleAxes(ax)
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', ...
        'FontSize', 14, 'LineWidth', 1, 'Box', 'on', 'GridColor', [0.2 0.2 0.2]);
set(get(ax, 'Title'),  'Color', 'k', 'FontSize', 17, 'FontWeight', 'bold');
set(get(ax, 'XLabel'), 'Color', 'k', 'FontSize', 15);
set(get(ax, 'YLabel'), 'Color', 'k', 'FontSize', 15);
set(get(ax, 'ZLabel'), 'Color', 'k', 'FontSize', 15);
end

function styleLegend(lg)
set(lg, 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k', 'FontSize', 13, 'Box', 'on');
end

function xw = wrapCurve(x)
xw = mod(x, 2*pi);
jump = any(abs(diff(xw)) > pi, 2);
if any(jump)
    idx = find(jump);
    out = zeros(size(xw,1) + numel(idx), 3); r = 0; last = 0;
    for j = idx.'
        seg = xw(last+1:j, :); out(r+1:r+size(seg,1), :) = seg; r = r + size(seg,1);
        out(r+1, :) = NaN; r = r + 1; last = j;
    end
    seg = xw(last+1:end, :); out(r+1:r+size(seg,1), :) = seg;
    xw = out;
end
end

function d = findDesc(R, name)
i = find(strcmp({R.desc.name}, name), 1);
if isempty(i), d = []; else, d = R.desc(i); end
end

function h = histLine(ax, v, c)
v = v(isfinite(v));
if isempty(v), v = 0; end
hi = prctile_(v, 99); lo = min(v);
if hi <= lo, hi = lo + 1; end
edges = linspace(lo, hi, 60);
cnt = histc(v, edges); cnt = cnt(:)'/max(sum(cnt),1);
h = plot(ax, edges, cnt, '-', 'Color', c, 'LineWidth', 2);
end

function p = prctile_(v, q)
v = sort(v(isfinite(v))); if isempty(v), p = 1; return; end
p = v(max(1, round(q/100*numel(v))));
if p <= 0, p = max(v) + eps; end
end

function saveFig(fg, cfg, name)
print(fg, fullfile(cfg.outDir, name), '-dpng', '-r150');
if strcmp(cfg.figVisible, 'off'), close(fg); end
end
