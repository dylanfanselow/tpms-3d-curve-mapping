function P = tpms_library_duplicates(lib)
%TPMS_LIBRARY_DUPLICATES  Find library equations that describe the same lattice.
%
%   P = tpms_library_duplicates(lib)   struct array of duplicate pairs
%
%   Two equations describe the same lattice if F_b(freq*(Q*x + t)) = c*F_a(x)
%   for a cubic rotation/reflection Q (48), a translation t (multiples of pi/4,
%   512), a scale c (any non-zero constant, sign included) and freq = 1 or 2
%   (frequency 2 = the same lattice at half the cell size).
%   Stage 1 screens pairs with an invariant: the distribution of F/std(F) at
%   200,000 random points, which no such transformation changes. Stage 2 tests
%   all transformations exactly for the pairs that pass (residual < 1e-8).
%
%   Only lattices with the same repeating box are compared (entries without a
%   'period' field are cubic, 2*pi); axis swaps are limited to those that keep
%   the box shape, and points and translations scale with the box.
%
%   P(k).a, P(k).b   keys;  P(k).how  the transformation found

n = numel(lib);
f = cell(n,1);
Pd = zeros(n,3);
for i = 1:n
    if isfield(lib, 'period') && ~isempty(lib(i).period), Pd(i,:) = tpms_period3(lib(i).period);
    else, Pd(i,:) = [2*pi 2*pi 2*pi]; end
end
for i = 1:n
    e = regexprep(lib(i).expr, '(?<!\.)([\*/\^])', '.$1');   % vectorise
    f{i} = str2func(['@(x,y,z) ' e ' + 0*x']);
end
st = rng; rng(7);
R0 = 2*pi*rand(200000, 3);
Pt = 2*pi*rand(150, 3);
rng(st);
qi = round(linspace(0.02, 0.98, 49) * size(R0,1));
FP = zeros(n, 49); FN = FP;
for i = 1:n
    Ri = R0 .* (Pd(i,:)/(2*pi));                      % unchanged for a cubic box
    v = f{i}(Ri(:,1), Ri(:,2), Ri(:,3)); v = v - mean(v);
    if std(v) > 0, v = v/std(v); end
    s1 = sort(v); s2 = sort(-v);
    FP(i,:) = s1(qi); FN(i,:) = s2(qi);
end
pm = perms(1:3);
[s1, s2, s3] = ndgrid([1 -1], [1 -1], [1 -1]); sg = [s1(:) s2(:) s3(:)];
[t1, t2, t3] = ndgrid(0:7, 0:7, 0:7); T = [t1(:) t2(:) t3(:)] * pi/4;
P = struct('a', {}, 'b', {}, 'how', {});
for i = 1:n
    for j = i+1:n
        if any(Pd(i,:) ~= Pd(j,:)), continue; end      % different boxes: different lattices
        if min(max(abs(FP(i,:) - FP(j,:))), max(abs(FN(i,:) - FP(j,:)))) > 0.03, continue; end
        sc = Pd(i,:)/(2*pi);
        Pi = Pd(i,:);
        keep = all(Pi(pm) == Pi, 2);                     % shape-preserving axis swaps
        how = testPair(f{i}, f{j}, Pt .* sc, pm(keep,:), sg, T .* sc);
        if isempty(how), how = testPair(f{j}, f{i}, Pt .* sc, pm(keep,:), sg, T .* sc); end
        if ~isempty(how)
            P(end+1) = struct('a', lib(i).key, 'b', lib(j).key, 'how', how); %#ok<AGROW>
        end
    end
end
end

% -------------------------------------------------------------------------
function how = testPair(fa_h, fb_h, Pt, pm, sg, T)
how = '';
fa = fa_h(Pt(:,1), Pt(:,2), Pt(:,3)).';                   % 1 x 150
nt = size(T,1); np = size(Pt,1);
for freq = 1:2
    for p = 1:size(pm,1)
        for s = 1:size(sg,1)
            Q = Pt(:, pm(p,:)) .* sg(s,:);
            X = freq*(repmat(Q(:,1).', nt, 1) + T(:,1));
            Y = freq*(repmat(Q(:,2).', nt, 1) + T(:,2));
            Z = freq*(repmat(Q(:,3).', nt, 1) + T(:,3));
            fb = reshape(fb_h(X(:), Y(:), Z(:)), nt, np);
            den = sum(fb.^2, 2); ok = den > 1e-12;
            c = zeros(nt,1); c(ok) = (fb(ok,:)*fa.') ./ den(ok);
            r = max(abs(fa - c.*fb), [], 2) / max(max(abs(fa)), 1e-12);
            r(~ok) = Inf;
            [rmin, k] = min(r);
            if rmin < 1e-8
                how = sprintf('frequency x%d, axes [%d %d %d], signs [%d %d %d], shift [%g %g %g]*pi, scale %.3g', ...
                    freq, pm(p,:), sg(s,:), T(k,:)/pi, c(k));
                return
            end
        end
    end
end
end
