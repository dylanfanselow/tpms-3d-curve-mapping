function opt = tpms_default_trace_opts(cfg, G)
%TPMS_DEFAULT_TRACE_OPTS  Tracer options from cfg, with singularity thresholds
%   set relative to the curvature scale of the surface being analysed (G).
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.
medGap = median(G.gap);
medK   = median(abs(G.K));
opt.h0         = cfg.h0;
opt.hmin       = cfg.hmin;
opt.hmax       = cfg.hmax;
opt.maxTurnDeg = cfg.maxTurnDeg;
opt.Lmax       = cfg.Lmax;
opt.maxSteps   = cfg.maxSteps;
opt.projTol    = 1e-10;
opt.closeTol   = cfg.closeTol;
opt.minLoop    = 0.5;
opt.umbTol     = cfg.umbStopRel * medGap;   % stop lines of curvature here
opt.umbExcl    = cfg.umbExclRel * medGap;   % exclude from statistics below this
opt.parTol     = cfg.parStopRel * medK;     % stop asymptotic curves here
opt.parExcl    = cfg.parExclRel * medK;
opt.tortWindow = cfg.tortWindow;
if isfield(cfg, 'windowsCells'), opt.windows = 2*pi*cfg.windowsCells; else, opt.windows = cfg.tortWindow; end
if isfield(cfg, 'acStep'),   opt.acStep = cfg.acStep;     else, opt.acStep = 0.05;      end
if isfield(cfg, 'acMaxLag'), opt.acMaxLag = cfg.acMaxLag; else, opt.acMaxLag = 4*pi;   end
opt.medGap     = medGap;
opt.medK       = medK;
end
