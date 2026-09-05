% step3_network_corr.m
% -------------------------------------------------------------------------
% NETWORK-level partial correlations, ALL combos (open/closed x dwpli/aec).
% Network value = mean of member nodes' AUC metric (whole-graph nodal
% measures grouped by Schaefer-7 label; NOT within-network subgraphs).
%
% FDR scheme (2026-08-12 decision, sole criterion at this level):
%   family = one nodal quantity x one band;
%   pool   = 7 networks x 6 AQ variables = 42 (bands NOT pooled across).
%   8 families per combo (2 quantities x 4 bands).
%
% Output (per combo): results_network_<cond>_<conn>.csv
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
        aqNames = D.aqNames; bands = D.bands; nodalQ = D.nodalQ;
        netIdx = D.netIdx; NETWORKS = D.NETWORKS;
        Z = [D.age, D.gender];
        nAQ = numel(aqNames); nB = numel(bands); nNet = numel(NETWORKS);

        nR = numel(nodalQ) * nB * nNet * nAQ;
        T = table('Size',[nR 12], 'VariableTypes', ...
            [repmat({'string'},1,4), repmat({'double'},1,6), repmat({'string'},1,2)], ...
            'VariableNames', {'quantity','band','network','aqvar','r','ci_lo','ci_hi','p','n','q4','level','cond_conn'});
        T.level = repmat("network", nR, 1);
        T.cond_conn = repmat(string(sprintf('%s_%s', cond, conn)), nR, 1);

        row = 0;
        for q = 1:numel(nodalQ)
            for b = 1:nB
                M = D.nodal.(nodalQ{q}).(sprintf('%s_%s', conn, bands{b}));   % N x 100
                for k = 1:nNet
                    v = mean(M(:, netIdx == k), 2, 'omitnan');
                    for a = 1:nAQ
                        row = row + 1;
                        [r,p,n] = partialCorr(v, D.(aqNames{a}), Z, minN);
                        [lo,hi] = ciFisher(r, n, size(Z,2));
                        T.quantity(row)=string(nodalQ{q}); T.band(row)=string(bands{b});
                        T.network(row)=string(NETWORKS{k}); T.aqvar(row)=string(aqNames{a});
                        T.r(row)=r; T.ci_lo(row)=lo; T.ci_hi(row)=hi; T.p(row)=p; T.n(row)=n;
                    end
                end
            end
        end

        % FDR: family = quantity x band, pool = 7 networks x 6 AQ = 42
        T.q4 = nan(nR,1);
        for q = 1:numel(nodalQ)
            for b = 1:nB
                idx = find(T.quantity == string(nodalQ{q}) & T.band == string(bands{b}));
                T.q4(idx) = fdrBH(T.p(idx));
            end
        end

        outFile = sprintf('results_network_%s_%s.csv', cond, conn);
        writetable(T, fullfile(outDir, outFile));

        sig = sortrows(T(T.q4 < alphaFDR, :), 'q4');
        fprintf('\n== NETWORK %s %s: q4<0.05 hits = %d ==\n', cond, conn, height(sig));
        if height(sig)>0, disp(sig(:, {'quantity','band','network','aqvar','r','p','q4'})); end
    end
end
fprintf('Done (network).\n');

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
