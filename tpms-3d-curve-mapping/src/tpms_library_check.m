function D = tpms_library_check(fld, S, G)
%TPMS_LIBRARY_CHECK  Sanity diagnostics for one level-set surface.
%
%   D = tpms_library_check(fld, S, G)
%     S  tpms_surface_sample output, G  tpms_geometry output at S.v
%
%   D.genus          genus in one 2*pi box from Gauss-Bonnet:
%                    integral of K dA = 2*pi*chi,  genus = 1 - chi/2
%   D.genusRaw       unrounded value (should be close to an integer)
%   D.nTrans         translations t (multiples of pi/4, within the box) with
%                    F(x+t) = F(x): surface symmetries that keep its two sides
%   D.nSwap          translations with F(x+t) = -F(x): map the surface onto
%                    itself but swap its sides (the two labyrinths are congruent)
%   D.genusPrim      genus per primitive cell = 1 + (genus - 1)/nTrans, the value
%                    usually quoted in the literature (side-swapping translations
%                    are excluded: they give a non-orientable quotient)
%   D.minimality     median |H| / median curvedness (0 for an exact minimal
%                    surface; small for a good nodal approximation)
%   D.smoothness     min |gradF| / median |gradF| on the surface
%                    (near 0 = singular point, e.g. surfaces touching)
%   D.cubicDev       max deviation of the area-weighted orientation tensor of
%                    the max-bending direction from I/3 (0 for cubic symmetry)
%   D.areaPerCell    surface area in one box divided by (2*pi)^2, i.e. per
%                    unit cell of side 1
%   D.halfPeriodAxes axes along which F repeats after pi (box holds 2 cells there)
%   D.curvednessChangeLength  median curvedness / median |surface gradient of
%                    curvedness| (radian units). Unlike the mean-curvature
%                    version it stays meaningful for minimal surfaces, where H = 0
%                    and its gradient only reflects the nodal-approximation error.

w = S.w;
D.genusRaw   = 1 - sum(G.K .* w)/(2*pi)/2;
D.genus      = round(D.genusRaw);
C            = sqrt((G.k1.^2 + G.k2.^2)/2);
sH = tpms_wstats(abs(G.H), w); sC = tpms_wstats(C, w);
D.minimality = sH.median / sC.median;
D.smoothness = min(G.gn) / median(G.gn);
wn = w/sum(w); Mt = zeros(3);
for i = 1:3, for j = 1:3, Mt(i,j) = sum(wn .* G.d1(:,i) .* G.d1(:,j)); end, end
D.cubicDev   = max(abs(Mt(:) - reshape(eye(3)/3, [], 1)));
Pv = tpms_period3(fld.period);                 % [2pi 2pi 2pi] for a cubic cell
D.areaPerCell = S.A / (Pv(1)*Pv(2));           % per a x b face; (2pi)^2 for a cube
X = Pv .* rand(400, 3);
f0 = fld.F0(X); ax = '';
names = 'xyz';
for i = 1:3
    Xs = X; Xs(:,i) = Xs(:,i) + Pv(i)/2;
    if max(abs(fld.F0(Xs) - f0)) < 1e-9, ax(end+1) = names(i); end %#ok<AGROW>
end
D.halfPeriodAxes = ax;
[i1, i2, i3] = ndgrid(0:7, 0:7, 0:7);
T = [i1(:) i2(:) i3(:)] .* (Pv/8);              % multiples of pi/4 for a cube
D.nTrans = 0; D.nSwap = 0;
for q = 1:size(T,1)
    fq = fld.F0(X + T(q,:));
    if max(abs(fq - f0)) < 1e-9, D.nTrans = D.nTrans + 1;
    elseif max(abs(fq + f0)) < 1e-9, D.nSwap = D.nSwap + 1; end
end
D.genusPrim = 1 + (D.genus - 1)/D.nTrans;
gk1 = G.dk(:,1).*G.d1 + G.dk(:,2).*G.d2;               % surface gradient of k1
gk2 = G.dk(:,3).*G.d1 + G.dk(:,4).*G.d2;               % surface gradient of k2
gC  = (G.k1.*gk1 + G.k2.*gk2) ./ (2*max(C, eps));      % gradient of curvedness
sG  = tpms_wstats(sqrt(sum(gC.^2, 2)), w);
D.curvednessChangeLength = sC.median / sG.median;
end
