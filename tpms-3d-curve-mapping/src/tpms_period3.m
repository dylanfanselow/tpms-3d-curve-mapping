function P = tpms_period3(p)
%TPMS_PERIOD3  Period of the repeating box as a 1-by-3 vector [Px Py Pz].
%   P = tpms_period3(fld.period)  accepts a scalar (cubic cell, the default
%   2*pi) or a 3-vector (non-cubic cell). A cubic period gives [2*pi 2*pi 2*pi],
%   so every per-axis calculation reduces to the original cubic one.
if isscalar(p), P = [p p p]; else, P = p(:).'; end
end
