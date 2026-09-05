% step2_global_corr.m
% -------------------------------------------------------------------------
% GLOBAL-level partial correlations, ALL combos (open/closed x dwpli/aec).
%
% FDR scheme (2026-08-12 decision, sole criterion at this level):
%   family = one global metric x one band;
%   pool   = 6 AQ variables (bands are NOT pooled across).
%   16 families per combo (4 metrics x 4 bands).
%
% Output (per combo): results_global_<cond>_<conn>.csv
% -------------------------------------------------------------------------

cfg = pipeline_config; outDir = cfg.corrDir;
alphaFDR = 0.05;
minN     = 10;

conds = {'open','closed'};
conns = {'dwpli','aec'};

for ci1 = 1:numel(conds)
    for ci2 = 1:numel(conns)
        cond = conds{ci1}; conn = conns{ci2};
        S = load(fullfile(outDir, sprintf('design_data_%s_%s.mat', cond, conn))); D = S.D;
        aqNames = D.aqNames; bands = D.bands; globalQ = D.globalQ;
        Z = [D.age, D.gender];
        nAQ = numel(aqNames); nB = numel(bands);

        nR = numel(globalQ) * nB * nAQ;
        T = table('Size',[nR 11], 'VariableTypes', ...
            [repmat({'string'},1,3), repmat({'double'},1,6), repmat({'string'},1,2)], ...
            'VariableNames', {'metric','band','aqvar','r','ci_lo','ci_hi','p','n','q4','level','cond_conn'});
        T.level = repmat("global", nR, 1);
        T.cond_conn = repmat(string(sprintf('%s_%s', cond, conn)), nR, 1);

        row = 0;
        for q = 1:numel(globalQ)
            for b = 1:nB
                fn = sprintf('%s_%s_%s', conn, bands{b}, globalQ{q});
                x = D.global.(fn);
                for a = 1:nAQ
                    row = row + 1;
                    [r,p,n] = partialCorr(x, D.(aqNames{a}), Z, minN);
                    [lo,hi] = ciFisher(r, n, size(Z,2));
                    T.metric(row)=string(fn); T.band(row)=string(bands{b});
                    T.aqvar(row)=string(aqNames{a});
                    T.r(row)=r; T.ci_lo(row)=lo; T.ci_hi(row)=hi; T.p(row)=p; T.n(row)=n;
                end
            end
        end

        % FDR: family = metric x band, pool = 6 AQ
        T.q4 = nan(nR,1);
        for q = 1:numel(globalQ)
            for b = 1:nB
                fn = sprintf('%s_%s_%s', conn, bands{b}, globalQ{q});
                idx = find(T.metric == string(fn));
                T.q4(idx) = fdrBH(T.p(idx));
            end
        end

        outFile = sprintf('results_global_%s_%s.csv', cond, conn);
        writetable(T, fullfile(outDir, outFile));

        sig = sortrows(T(T.q4 < alphaFDR, :), 'q4');
        fprintf('\n== GLOBAL %s %s: q4<0.05 hits = %d ==\n', cond, conn, height(sig));
        if height(sig)>0, disp(sig(:, {'metric','aqvar','r','p','q4'})); end
    end
end
fprintf('Done (global).\n');

%% ===== helpers =====
function [r, p, n] = partialCorr(x, y, Z, minN)
ok = ~isnan(x) & ~isnan(y) & all(~isnan(Z), 2);
n  = sum(ok);
if n < minN + size(Z,2) + 2, r = NaN; p = NaN; return; end
x = x(ok); y = y(ok); Z = Z(ok,:);
X = [ones(n,1), Z];
rx = x - X * (X \ x);
ry = y - X * (X \ y);
r = sum(rx.*ry) / sqrt(sum(rx.^2) * sum(ry.^2));
df = n - 2 - size(Z,2);
t  = r * sqrt(df / max(1 - r^2, eps));
p  = 2 * tcdf(-abs(t), df);
end

function [lo, hi] = ciFisher(r, n, k)
if isnan(r) || n <= k + 3, lo = NaN; hi = NaN; return; end
z  = atanh(max(min(r, 0.999999), -0.999999));
se = 1 / sqrt(n - k - 3);
lo = tanh(z - 1.96*se);
hi = tanh(z + 1.96*se);
end

function q = fdrBH(p)
q = NaN(size(p));
valid = ~isnan(p);
pv = p(valid);
m = numel(pv);
if m == 0, return; end
[sv, ord] = sort(pv);
qv = sv .* m ./ (1:m)';
qv = cummin(qv(end:-1:1)); qv = qv(end:-1:1);
qv = min(qv, 1);
tmp = NaN(m,1); tmp(ord) = qv;
q(valid) = tmp;
end
