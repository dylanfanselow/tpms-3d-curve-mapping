function P = tpms_library_plating(fld, cfg)
%TPMS_LIBRARY_PLATING  Electroplating numbers for one lattice at a matched solid fraction.
%
%   P = tpms_library_plating(fld, cfg)
%
%   At t = 0 a sheet lattice has zero wall thickness (porosity 1), so the
%   Wagner and inverse Damkohler numbers would only reflect surface area.
%   Lattices are therefore compared at the same solid fraction,
%   cfg.platingSolidFraction (default 0.16 = porosity 0.84, as in the reference
%   strut lattices of Clemens et al. 2025):
%     sheet    solid = |F| < t,  t = that quantile of |F| over the cell
%     network  solid =  F  < t,  t = that quantile of F
%   tpms_wall_thickness then gives the porosity, wetted area (both faces for a
%   sheet) and mean wall thickness at that t, and tpms_plating_numbers the
%   Wagner / inverse Damkohler numbers with the process inputs in cfg.plating.
%
%   NOTE: at equal porosity both Wa and Da^-1 scale as 1/a_v, so they rank
%   lattices identically; their absolute values answer different questions
%   (coating uniformity vs. flow needed).

if isfield(cfg, 'platingSolidFraction'), sf = cfg.platingSolidFraction; else, sf = 0.16; end
if isfield(cfg, 'vfGrid'), nv = cfg.vfGrid; else, nv = 96; end
Pv = tpms_period3(fld.period);                 % [2pi 2pi 2pi] for a cubic cell
nvv = round(nv * Pv / (2*pi));
[X, Y, Z] = ndgrid(((1:nvv(1)) - 0.5) * Pv(1) / nvv(1), ((1:nvv(2)) - 0.5) * Pv(2) / nvv(2), ...
                   ((1:nvv(3)) - 0.5) * Pv(3) / nvv(3));
Fv = fld.F0([X(:) Y(:) Z(:)]);
if strcmp(cfg.latticeType, 'sheet'), q = sort(abs(Fv)); else, q = sort(Fv); end
t = q(max(1, round(sf*numel(q))));

T = tpms_wall_thickness(fld, t, cfg);
toM = lengthToMeters(cfg.lengthUnit);
if strcmp(cfg.latticeType, 'sheet')
    por = 1 - T.vf_sheet;   ad = T.area_density_sheet;   th = T.thick_mean;
else
    por = 1 - T.vf_network; ad = T.area_density_network; th = NaN;
end
a_v = ad / toM;
N = tpms_plating_numbers(por, a_v, cfg.cellSize*toM, cfg.plating);

P.solidFractionTarget = sf;
P.t = t;
P.porosity = por;
P.a_v = a_v;
P.thickness = th;                          % in cfg.lengthUnit
P.thicknessFailFrac = T.thick_fail_frac;
P.Wa_geom = N.Wa_geom; P.Wa = N.Wa; P.i_max = N.i_max;
P.Da_geom = N.Da_geom; P.DaInv = N.DaInv;
P.v_min = N.v_min; P.DaInvTarget = N.DaInvTarget;
end

function f = lengthToMeters(u)
switch lower(strtrim(u))
    case 'm', f = 1; case 'cm', f = 1e-2; case 'mm', f = 1e-3;
    case {'um', 'micron'}, f = 1e-6; case {'in', 'inch'}, f = 0.0254;
    otherwise, error('tpms_library_plating:unit', 'Unknown length unit ''%s''.', u);
end
end
