% step4_node_corr.m
% -------------------------------------------------------------------------
% NODE-level partial correlations, ALL combos (open/closed x dwpli/aec).
%
% FDR scheme (2026-08-12 decision, sole criterion at this level):
%   family = one nodal quantity x one Schaefer network x one band;
%   pool   = nodes within that network x 6 AQ variables (bands NOT pooled).
%   56 families per combo (2 quantities x 7 networks x 4 bands);
%   pool sizes = n_nodes(net) x 6 (Vis 102 / SomMot 84 / DorsAttn 90 /
%   SalVentAttn 72 / Limbic 30 / Cont 78 / Default 144).
%
% Output (per combo): results_node_<cond>_<conn>.csv
%                     results_node_uncorrected_sig_<cond>_<conn>.csv
% -------------------------------------------------------------------------

cfg = pipeline_config; outDir = cfg.corrDir;
alphaFDR = 0.05;
minN     = 10;

conds = {'open','closed'};
conns = {'dwpli','aec'};

atlasCsv = cfg.atlasCsv;
A = readtable(atlasCsv, 'PreserveVariableNames', true);
roiNames = string(A.("ROI Name"));

for ci1 = 1:numel(conds)
    for ci2 = 1:numel(conns)
        cond = conds{ci1}; conn = conns{ci2};
        S = load(fullfile(outDir, sprintf('design_data_%s_%s.mat', cond, conn))); D = S.D;
        aqNames = D.aqNames; bands = D.bands; nodalQ = D.nodalQ;
        netIdx = D.netIdx; NETWORKS = D.NETWORKS;
        Z = [D.age, D.gender];
        nAQ = numel(aqNames); nB = numel(bands); nNet = numel(NETWORKS);

        nR = numel(nodalQ) * nB * 100 * nAQ;
        T = table('Size',[nR 12], 'VariableTypes', ...
            [repmat({'string'},1,3), {'double'}, {'string'}, repmat({'double'},1,6), {'string'}], ...
            'VariableNames', {'quantity','band','aqvar','node','roi','r','ci_lo','ci_hi','p','n','q4','level'});
        T.level = repmat("node", nR, 1);

        row = 0;
        for q = 1:numel(nodalQ)
            for b = 1:nB
                M = D.nodal.(nodalQ{q}).(sprintf('%s_%s', conn, bands{b}));   % N x 100
                for node = 1:100
                    for a = 1:nAQ
                        row = row + 1;
                        [r,p,n] = partialCorr(M(:,node), D.(aqNames{a}), Z, minN);
                        [lo,hi] = ciFisher(r, n, size(Z,2));
                        T.quantity(row)=string(nodalQ{q}); T.band(row)=string(bands{b});
                        T.aqvar(row)=string(aqNames{a}); T.node(row)=node;
                        T.roi(row)=roiNames(node);
                        T.r(row)=r; T.ci_lo(row)=lo; T.ci_hi(row)=hi; T.p(row)=p; T.n(row)=n;
                    end
                end
            end
        end

        % FDR: family = quantity x network x band, pool = nodes-in-net x 6 AQ
        T.q4 = nan(nR,1);
        for q = 1:numel(nodalQ)
            for k = 1:nNet
                nodesK = find(netIdx == k);
                for b = 1:nB
                    idx = find(T.quantity == string(nodalQ{q}) & T.band == string(bands{b}) ...
                        & ismember(T.node, nodesK));
                    T.q4(idx) = fdrBH(T.p(idx));
                end
            end
        end

        outFile = sprintf('results_node_%s_%s.csv', cond, conn);
        writetable(T, fullfile(outDir, outFile));

        Tu = sortrows(T(T.p < 0.05, :), 'p');
        writetable(Tu, fullfile(outDir, sprintf('results_node_uncorrected_sig_%s_%s.csv', cond, conn)));

        sig = sortrows(T(T.q4 < alphaFDR, :), 'q4');
        fprintf('\n== NODE %s %s: q4<0.05 hits = %d (uncorrected p<0.05: %d) ==\n', ...
            cond, conn, height(sig), height(Tu));
        if height(sig)>0
            disp(sig(:, {'quantity','band','node','roi','aqvar','r','p','q4'}));
        end
        [~,ord] = sort(T.p);
        fprintf('--- strongest 3 uncorrected (context) ---\n');
        disp(T(ord(1:3), {'quantity','band','node','roi','aqvar','r','p','q4'}));
    end
end
fprintf('Done (node).\n');

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
