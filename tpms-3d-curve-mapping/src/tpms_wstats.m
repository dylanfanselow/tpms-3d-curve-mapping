function st = tpms_wstats(x, w)
%TPMS_WSTATS  Weighted distribution statistics (mean, std, CV, p5, median, p95, RMS).
%   Same statistic set as the STL-based tool, so the columns line up for PCA.
%   Weighted quantiles use the midpoint rule on the cumulative weight.
%
%   Copyright (C) 2026 YOUR NAME
%   SPDX-License-Identifier: GPL-3.0-or-later
%   See LICENSE in the repository root for the full license text.
x = x(:); if nargin < 2 || isempty(w), w = ones(size(x)); end
w = w(:);
ok = isfinite(x) & isfinite(w) & w > 0;
x = x(ok); w = w(ok);
st = struct('mean',NaN,'std',NaN,'CV',NaN,'p5',NaN,'median',NaN,'p95',NaN,'RMS',NaN,'n',numel(x));
if isempty(x), return; end
W  = sum(w);
mu = sum(w.*x)/W;
sd = sqrt(sum(w.*(x-mu).^2)/W);
[xs, o] = sort(x); ws = w(o);
cw = (cumsum(ws) - ws/2)/W;
q  = @(p) quantileMid(xs, cw, p);
st.mean = mu; st.std = sd; st.CV = sd/abs(mu);
st.p5 = q(0.05); st.median = q(0.5); st.p95 = q(0.95);
st.RMS = sqrt(sum(w.*x.^2)/W);
end

function v = quantileMid(xs, cw, p)
if numel(xs) == 1, v = xs; return; end
if p <= cw(1), v = xs(1); return; end
if p >= cw(end), v = xs(end); return; end
[cwu, iu] = unique(cw);
v = interp1(cwu, xs(iu), p, 'linear');
end
