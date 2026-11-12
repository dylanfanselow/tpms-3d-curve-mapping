function PL = tpms_electroplating(thick, cfg, outFile)
%TPMS_ELECTROPLATING  Lattice Wagner / inverse Damkohler numbers versus t.
%
%   PL = tpms_electroplating(thick, cfg)            struct array, one entry per t
%   PL = tpms_electroplating(thick, cfg, outFile)   also writes a CSV
%
%   thick is the output of tpms_wall_thickness (lengths in cfg.lengthUnit).
%   Geometry per lattice type:
%     sheet    eps = 1 - vf_sheet,   wetted area = BOTH faces (2A per cell);
%              both gyroid channel systems are open to the electrolyte
%     network  eps = 1 - vf_network, wetted area = the single surface (A per cell)
%   cfg.plating holds the process inputs (see tpms_plating_numbers).
%   Returns [] if cfg.cellSize is empty (physical size needed) or if
%   cfg.plating is missing.
%
%   For RANKING designs only: equations reconstructed from Clemens et al.,
%   ACS Appl. Eng. Mater. 2025, 3(3), 600-612 (absolute values may be off by ~1.5x).
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.

PL = [];
if isempty(cfg.cellSize) || ~isfield(cfg, 'plating') || isempty(cfg.plating), return; end
toM = lengthToMeters(cfg.lengthUnit);
L   = cfg.cellSize * toM;
pc  = cfg.plating;

for k = 1:numel(thick)
    r = thick(k);
    if strcmp(cfg.latticeType, 'sheet')
        por = 1 - r.vf_sheet;   ad = r.area_density_sheet;
    else
        por = 1 - r.vf_network; ad = r.area_density_network;
    end
    a_v = ad / toM;                               % 1/lengthUnit -> 1/m
    P = tpms_plating_numbers(por, a_v, L, pc);
    e.t = r.t; e.porosity = por; e.a_v = a_v; e.avL = P.avL;
    e.Wa_geom = P.Wa_geom; e.Da_geom = P.Da_geom;
    e.Wa = P.Wa; e.DaInv = P.DaInv; e.i_max = P.i_max;
    e.v_min = P.v_min; e.DaInvTarget = P.DaInvTarget;
    if isempty(PL), PL = e; else, PL(end+1) = e; end %#ok<AGROW>
end

if nargin >= 3 && ~isempty(outFile) && ~isempty(PL)
    fid = fopen(outFile, 'w');
    fprintf(fid, ['# Lattice electroplating numbers (Clemens et al., ACS Appl. Eng. Mater. 2025, 3(3), 600-612). ' ...
        'Reconstructed equations: RANK designs only; absolute values may be off by ~1.5x.\n']);
    fprintf(fid, ['# Assumptions: Tafel kinetics; uniform electrolyte conductivity; perfectly conductive substrate; ' ...
        'convection-dominated ion supply. Wa predicts uniformity only when DaInv is high.\n']);
    fprintf(fid, '# lattice=%s cell=%g %s N_field=%g N_flow=%g T=%g K alpha_c=%g kappa=%g S/m i_avg=%g A/m2 n=%g C_o=%g mol/m3 v=%g m/s\n', ...
        cfg.latticeType, cfg.cellSize, cfg.lengthUnit, pg(pc,'N_field'), pg(pc,'N_flow'), pg(pc,'T'), ...
        pg(pc,'alpha_c'), pg(pc,'kappa'), pg(pc,'i_avg'), pg(pc,'n'), pg(pc,'C_o'), pg(pc,'v'));
    fprintf(fid, 't,porosity,a_v_per_m,a_v_times_L,Wa_geom_per_m,Da_geom,Wa,DaInv,i_max_mA_per_cm2_at_Wa_%g', pg(pc,'WaTarget'));
    tg = PL(1).DaInvTarget;
    for j = 1:numel(tg), fprintf(fid, ',v_min_m_per_s_for_DaInv_%g', tg(j)); end
    fprintf(fid, '\n');
    for e = PL
        fprintf(fid, '%.4g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g,%.5g', e.t, e.porosity, e.a_v, e.avL, ...
            e.Wa_geom, e.Da_geom, e.Wa, e.DaInv, e.i_max/10);        % A/m^2 -> mA/cm^2
        fprintf(fid, ',%.5g', e.v_min);
        fprintf(fid, '\n');
    end
    fclose(fid);
end
end

% -------------------------------------------------------------------------
function f = lengthToMeters(u)
switch lower(strtrim(u))
    case 'm',                f = 1;
    case 'cm',               f = 1e-2;
    case 'mm',               f = 1e-3;
    case {'um', 'micron'},   f = 1e-6;
    case {'in', 'inch'},     f = 0.0254;
    otherwise
        error('tpms_electroplating:unit', ...
            'Unknown cfg.lengthUnit ''%s'' (use m, cm, mm, um or in).', u);
end
end

function v = pg(s, f)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = NaN; end
end
