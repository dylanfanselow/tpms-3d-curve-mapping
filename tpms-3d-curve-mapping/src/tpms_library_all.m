function L = tpms_library_all(includeNonCubic)
%TPMS_LIBRARY_ALL  Cubic library plus (optionally) the non-cubic add-on.
%   L = tpms_library_all(true)   every entry gets the fields period ([] = cubic
%   2*pi box) and cell ('cubic' or the non-cubic cell type). The cubic entries
%   are exactly those of tpms_library, in the same order.
if nargin < 1, includeNonCubic = true; end
L = tpms_library();
[L.period] = deal([]);
[L.cell] = deal('cubic');
if includeNonCubic
    N = tpms_library_noncubic();
    L = [L, N];
end
end
