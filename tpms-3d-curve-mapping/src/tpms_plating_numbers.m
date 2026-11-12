function P = tpms_plating_numbers(por, a_v, L, pc)
%TPMS_PLATING_NUMBERS  Lattice Wagner and inverse Damkohler numbers for electroplating.
%
%   P = tpms_plating_numbers(por, a_v, L, pc)
%
%   Lattice versions of two dimensionless numbers after Clemens et al.,
%   ACS Appl. Eng. Mater. 2025, 3(3), 600-612, doi:10.1021/acsaenm.4c00561.
%   The equations are RECONSTRUCTED from the paper's variable definitions:
%   the scaling is right but absolute values may be off by a constant
%   (about 1.5x when checked). Use them to RANK designs, not as absolute
%   predictions.
%
%   This function depends only on bulk geometry, so it can be used for any
%   lattice (TPMS or strut, analytic or STL-based):
%     por   porosity eps = void volume / bounding volume      (0..1)
%     a_v   wetted surface area / bounding volume              [1/m]
%     L     unit cell size                                     [m]
%   por, a_v and L may be vectors of the same size (one entry per design).
%
%   pc (process config, any field may be empty):
%     N_field      cells along the current direction (toward the anode);
%                  if plated from several sides, the depth from the nearest exposed face
%     N_flow       cells along the electrolyte flow direction
%     T [K], alpha_c (cathodic charge-transfer coefficient), kappa [S/m],
%     i_avg [A/m^2], n (electrons per ion), C_o [mol/m^3], v [m/s]
%     WaTarget     target Wagner number for the "speed limit" i_max
%     DaInvTarget  target inverse Damkohler number(s) for the flow requirement v_min
%
%   OUTPUT (NaN where the needed inputs are missing)
%     Wa_geom   = eps / (a_v (N_field L)^2)                         [1/m]
%     Da_geom   = 1 / (a_v N_flow L)                                [-]
%     avL       = a_v L, dimensionless surface density               [-]
%     Wa        = Wa_geom * (R T/(alpha_c F)) * kappa / i_avg      [-]   coating uniformity
%     DaInv     = Da_geom * n F v C_o / i_avg                      [-]   ion supply vs consumption
%     i_max     = (R T/(alpha_c F)) * kappa * Wa_geom / WaTarget   [A/m^2]
%     v_min     = DaInvTarget * i_avg / (Da_geom n F C_o)          [m/s], one column per target
%
%   ASSUMPTIONS: Tafel kinetics, uniform electrolyte conductivity (no
%   tortuosity correction beyond eps), perfectly conductive substrate,
%   convection-dominated ion supply. Wa predicts uniformity only when DaInv
%   is high (enough flow); otherwise depletion effects dominate.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.

Rg = 8.314;  Fc = 96485;
g  = @(f) getf(pc, f);

P.Wa_geom = por ./ (a_v .* (g('N_field') .* L).^2);
P.Da_geom = 1 ./ (a_v .* g('N_flow') .* L);
P.avL     = a_v .* L;

VT = Rg * g('T') / (g('alpha_c') * Fc);           % Tafel voltage scale [V]
P.Wa    = P.Wa_geom .* VT .* g('kappa') ./ g('i_avg');
P.DaInv = P.Da_geom .* g('n') .* Fc .* g('v') .* g('C_o') ./ g('i_avg');
P.i_max = VT .* g('kappa') .* P.Wa_geom ./ g('WaTarget');

tg = getf(pc, 'DaInvTarget');
P.DaInvTarget = tg(:).';
P.v_min = nan(numel(por), numel(tg));
for k = 1:numel(tg)
    P.v_min(:,k) = tg(k) .* g('i_avg') ./ (P.Da_geom(:) .* g('n') .* Fc .* g('C_o'));
end
end

% -------------------------------------------------------------------------
function v = getf(s, f)
% missing or empty field -> NaN, so dependent outputs become NaN
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = NaN; end
end
