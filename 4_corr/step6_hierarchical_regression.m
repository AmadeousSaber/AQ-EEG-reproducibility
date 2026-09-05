% step6_hierarchical_regression.m
% -------------------------------------------------------------------------
% Story-driven hierarchical regression, epoch2s_analysis_open(0.15_0.3).
% Conditions (open / closed) analyzed SEPARATELY. Narrative units rebuilt
% from the density-fixed, band-separated-FDR hit set (2026-08-12).
%
% Upgrades vs the parent version (epoch2s_analysis/4_corr/step6):
%   1. NESTED leave-one-out CV: per fold, the unit's FDR-family screen is
%      re-run on the training data, selected predictors are PCA-fused
%      (z-scored, sign set within training), and the held-out subject is
%      predicted. Gives an honest predictive q2 (vs age+gender baseline),
%      bounding the circularity of "screen-then-regress on the same data".
%   2. Freedman-Lane residual permutation for dR2 (main and interaction)
%      instead of permuting the brain rows unconditionally.
%   3. Predictors z-scored before PCA (old version used raw scales).
%   4. Influence diagnostics: max Cook's D; dR2 and p after excluding the
%      most influential subject.
%   5. BH across confirmatory tests for BOTH main (q_unit) and interaction
%      (q_inter). The open theta-imagination trend is reported as a
%      separate exploratory row, outside the confirmatory family.
%
% Models per unit x condition:
%   M1: aq ~ age + gender
%   M2: + BRAIN (single predictor as-is, or PCA-PC1 of z-scored predictors)
%   M3: + gender x BRAIN
%
% Outputs: results_regression_summary.csv ; plot_data_regression.mat
% -------------------------------------------------------------------------

rng(20260812);

cfg      = pipeline_config;
corrDir  = cfg.corrDir;
atlasCsv = cfg.atlasCsv;
B     = 2000;      % Freedman-Lane permutations
Bboot = 2000;      % bootstrap resamples for dR2 CI
alpha = 0.05;

NETWORKS = {'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
A = readtable(atlasCsv,'PreserveVariableNames',true);
netIdx = zeros(100,1);
for n = 1:100
    tok = strsplit(A.("ROI Name"){n},'_'); netIdx(n) = find(strcmp(NETWORKS,tok{3}));
end
NET = @(nm) find(strcmp(NETWORKS,nm));

%% ===== units: fixed predictors (full-sample hits) + scopes (nested screen) ==
% scope = {level, conn, band, qty, net} — defines the FDR family re-run per CV fold.
units = {};
units{end+1} = mkUnit('U1 detail x AEC visual-beta degree + gamma global','detail','confirm', ...
    struct('open',  {{ mkP('network','aec','beta','auc_degree',NET('Vis'),'net:Vis:beta:deg(AEC)'), ...
                       mkP('node','aec','beta','auc_degree',55,'node55(RH_Vis_5):beta:deg'), ...
                       mkP('node','aec','beta','auc_degree',54,'node54(RH_Vis_4):beta:deg'), ...
                       mkP('node','aec','beta','auc_degree',58,'node58(RH_Vis_8):beta:deg'), ...
                       mkP('node','aec','beta','auc_degree',5 ,'node5(LH_Vis_5):beta:deg') }}, ...
           'closed',{{ mkP('global','aec','gamma','auc_gcc',[],'g:gamma:gcc(AEC)'), ...
                       mkP('global','aec','gamma','auc_sw' ,[],'g:gamma:sw(AEC)'), ...
                       mkP('node','aec','beta','auc_degree',65,'node65(RH_SomMot_7):beta:deg') }}), ...
    struct('open',  {{ {'network','aec','beta','auc_degree',NET('Vis')}, {'node','aec','beta','auc_degree',NET('Vis')} }}, ...
           'closed',{{ {'global','aec','gamma','auc_gcc',[]}, {'global','aec','gamma','auc_sw',[]}, {'node','aec','beta','auc_degree',NET('SomMot')} }}));
units{end+1} = mkUnit('U2 switching x attention/control node degree','switching','confirm', ...
    struct('open',  {{ mkP('node','aec','beta','auc_degree',82,'node82(RH_Cont_Par_2):beta:deg') }}, ...
           'closed',{{ mkP('node','dwpli','alpha','auc_degree',89,'node89(RH_Cont_pCun_1):alpha:deg'), ...
                       mkP('node','dwpli','alpha','auc_cc',19,'node19(LH_DorsAttn_Post_4):alpha:cc') }}), ...
    struct('open',  {{ {'node','aec','beta','auc_degree',NET('Cont')} }}, ...
           'closed',{{ {'node','dwpli','alpha','auc_degree',NET('Cont')}, {'node','dwpli','alpha','auc_cc',NET('DorsAttn')} }}));
units{end+1} = mkUnit('U3 communication x node42 / Cont network','communication','confirm', ...
    struct('open',  {{ mkP('node','dwpli','alpha','auc_cc',42,'node42(LH_Default_PFC_1):alpha:cc') }}, ...
           'closed',{{ mkP('network','aec','beta','auc_degree',NET('Cont'),'net:Cont:beta:deg(AEC)') }}), ...
    struct('open',  {{ {'node','dwpli','alpha','auc_cc',NET('Default')} }}, ...
           'closed',{{ {'network','aec','beta','auc_degree',NET('Cont')} }}));
units{end+1} = mkUnit('U4 imagination x DorsAttn FEF gamma cc','imagination','confirm', ...
    struct('open',  {{}}, ...
           'closed',{{ mkP('node','dwpli','gamma','auc_cc',23,'node23(LH_DorsAttn_FEF_1):gamma:cc') }}), ...
    struct('open',  {{}}, ...
           'closed',{{ {'node','dwpli','gamma','auc_cc',NET('DorsAttn')} }}));
units{end+1} = mkUnit('U0 imagination x theta segregation [EXPLORATORY TREND]','imagination','explor', ...
    struct('open',  {{ mkP('global','dwpli','theta','auc_gcc',[],'g:theta:gcc'), ...
                       mkP('global','dwpli','theta','auc_sw' ,[],'g:theta:sw'), ...
                       mkP('network','dwpli','theta','auc_cc',NET('Default'),'net:Default:theta:cc'), ...
                       mkP('network','dwpli','theta','auc_cc',NET('SalVentAttn'),'net:SalVentAttn:theta:cc') }}, ...
           'closed',{{}}), ...
    struct('open',  {{ {'global','dwpli','theta','auc_gcc',[]}, {'global','dwpli','theta','auc_sw',[]}, ...
                       {'network','dwpli','theta','auc_cc',NET('Default')}, {'network','dwpli','theta','auc_cc',NET('SalVentAttn')} }}, ...
           'closed',{{}}));

%% ===== run =====
rows = {}; plotS = struct(); pi = 0;
conds = {'open','closed'};

for c = 1:numel(conds)
    cond = conds{c};
    S = load(fullfile(corrDir, sprintf('design_data_%s_dwpli.mat', cond))); Dw = S.D;
    S = load(fullfile(corrDir, sprintf('design_data_%s_aec.mat', cond)));   Da = S.D;
    assert(isequal(Dw.subj_ids, Da.subj_ids), 'subject mismatch %s', cond);
    Dc = mergeDesign(Dw, Da);
    age = Dc.age(:); gender01 = Dc.gender(:) - 1;
    fprintf('\n==================== %s  (N=%d) ====================\n', upper(cond), numel(age));

    for u = 1:numel(units)
        U = units{u};
        preds  = U.preds.(cond);
        scopes = U.scopes.(cond);
        if isempty(preds)
            fprintf('  [%s] skipped in %s\n', U.name, cond);
            continue;
        end
        y = Dc.(U.aqvar)(:);
        ok = ~isnan(y) & ~isnan(age) & ~isnan(gender01);
        Pmat = zeros(numel(y), numel(preds));
        for k = 1:numel(preds), Pmat(:,k) = getCol(Dc, preds{k}, netIdx); end
        ok = ok & all(~isnan(Pmat),2);
        n  = sum(ok);
        yk = y(ok); ak = age(ok); gk = gender01(ok);
        Pk = Pmat(ok,:);
        Pk = (Pk - mean(Pk,1)) ./ max(std(Pk,0,1), eps);   % z-score predictors

        % ---- fused BRAIN component ----
        if size(Pk,2) > 1
            [coef, score] = pca(Pk);
            brain = score(:,1);
            if corr(brain, yk) < 0, brain = -brain; coef(:,1) = -coef(:,1); end
            comp_label = 'PC1[';
            for k = 1:numel(preds), comp_label = [comp_label sprintf('%s=%.2f ', preds{k}.tag, coef(k,1))]; end %#ok<AGROW>
            comp_label = [comp_label ']'];
        else
            brain = Pk(:,1);
            if corr(brain, yk) < 0, brain = -brain; end
            comp_label = preds{1}.tag;
        end

        X0 = [ones(n,1) ak gk];
        X1 = [X0 brain];
        X2 = [X1 brain.*gk];
        o0 = ols(yk, X0); o1 = ols(yk, X1); o2 = ols(yk, X2);
        dR2m = o1.r2 - o0.r2;
        dR2i = o2.r2 - o1.r2;
        p_main_asy = o1.p(end);
        p_int_asy  = o2.p(end);
        std_beta   = o1.b(end) * std(brain) / std(yk);

        % ---- Freedman-Lane permutation on dR2 ----
        b0fl = X0 \ yk;  fit0 = X0*b0fl;  e0 = yk - fit0;
        b1fl = X1 \ yk;  fit1 = X1*b1fl;  e1 = yk - fit1;
        null_m = zeros(B,1); null_i = zeros(B,1);
        for b = 1:B
            pm = randperm(n);
            null_m(b) = ols(fit0 + e0(pm), X1).r2 - o0.r2;
            pi2 = randperm(n);
            null_i(b) = ols(fit1 + e1(pi2), X2).r2 - o1.r2;
        end
        p_main_perm = (sum(null_m >= dR2m) + 1)/(B+1);
        p_int_perm  = (sum(null_i >= dR2i) + 1)/(B+1);

        % ---- bootstrap percentile CI for dR2m ----
        boots = zeros(Bboot,1);
        for b = 1:Bboot
            bi = randi(n, n, 1);
            boots(b) = ols(yk(bi), X1(bi,:)).r2 - ols(yk(bi), X0(bi,:)).r2;
        end
        ci_lo = prctile(boots, 2.5); ci_hi = prctile(boots, 97.5);

        % ---- influence diagnostics ----
        [maxD, dR2_excl, p_excl] = cookCheck(yk, X0, brain);

        % ---- nested LOO predictive q2 ----
        [q2_base, q2_unit] = nestedLOO(Dc, scopes, U.aqvar, age, gender01, netIdx);

        % ---- simple slopes if interaction notable ----
        sm = NaN; pm = NaN; sf = NaN; pf = NaN;
        if min(p_int_asy, p_int_perm) < 0.10
            for s = 0:1
                idx = (gk == s);
                os = ols(yk(idx), [ones(sum(idx),1) ak(idx) brain(idx)]);
                if s == 0, sm = os.b(end); pm = os.p(end); else, sf = os.b(end); pf = os.p(end); end
            end
        end

        pi = pi+1;
        plotS(pi).unit = U.name; plotS(pi).cond = cond; plotS(pi).aqvar = U.aqvar;
        plotS(pi).comp_label = comp_label;
        plotS(pi).resid_brain = residv(brain, X0);
        plotS(pi).resid_aq    = residv(yk, X0);
        plotS(pi).gender01 = gk;
        plotS(pi).dR2m = dR2m; plotS(pi).p_main_perm = p_main_perm;
        plotS(pi).q2_unit = q2_unit; plotS(pi).q2_base = q2_base;

        rows = [rows; {cond, U.name, U.family, U.aqvar, n, comp_label, ...
            o0.r2, o1.r2, dR2m, std_beta, p_main_asy, p_main_perm, ci_lo, ci_hi, ...
            o2.r2, dR2i, p_int_asy, p_int_perm, sm, pm, sf, pf, ...
            maxD, dR2_excl, p_excl, q2_base, q2_unit, q2_unit - q2_base}]; %#ok<AGROW>

        fprintf(['  [%s] %s | N=%d | dR2=%.3f (p_asy=%.3f p_FL=%.3f, bootCI [%.2f %.2f])' ...
                 ' | cookD=%.2f -> excl dR2=%.3f (p=%.3f) | q2: base=%.3f unit=%.3f (dq2=%.3f)'], ...
            U.name, U.aqvar, n, dR2m, p_main_asy, p_main_perm, ci_lo, ci_hi, ...
            maxD, dR2_excl, p_excl, q2_base, q2_unit, q2_unit - q2_base);
        if min(p_int_asy, p_int_perm) < 0.10
            fprintf(' | INTER dR2=%.3f (p_asy=%.3f p_FL=%.3f) M=%.2f(p%.3f) F=%.2f(p%.3f)\n', ...
                dR2i, p_int_asy, p_int_perm, sm, pm, sf, pf);
        else
            fprintf(' | inter n.s.\n');
        end
    end
end

varn = {'cond','unit','family','aqvar','N','brain_component','R2_base','R2_main','dR2_main', ...
        'std_beta','dR2_main_p_asy','dR2_main_p_permFL','dR2_boot_lo','dR2_boot_hi', ...
        'R2_inter','dR2_inter','inter_p_asy','inter_p_permFL', ...
        'slope_male','slope_male_p','slope_female','slope_female_p', ...
        'maxCookD','dR2_excl_maxD','dR2_excl_maxD_p','q2_base','q2_unit','dq2'};
T = cell2table(rows,'VariableNames',varn);

% BH across CONFIRMATORY tests only (exploratory rows reported uncorrected)
isConf = strcmp(T.family, 'confirm');
q1 = fdrBH(T.dR2_main_p_permFL(isConf));
q2i = fdrBH(T.inter_p_permFL(isConf));
T.q_unit = nan(height(T),1); T.q_inter = nan(height(T),1);
T.q_unit(isConf) = q1; T.q_inter(isConf) = q2i;

writetable(T, fullfile(corrDir,'results_regression_summary.csv'));
save(fullfile(corrDir,'plot_data_regression.mat'),'plotS','-v7.3');
fprintf('\nSaved results_regression_summary.csv and plot_data_regression.mat\n');

%% ============================ helpers ============================
function u = mkUnit(name, aqvar, family, preds, scopes)
u = struct('name',{name},'aqvar',{aqvar},'family',{family},'preds',{preds},'scopes',{scopes});
end

function p = mkP(level,conn,band,qty,net,tag)
p = struct('level',level,'conn',conn,'band',band,'qty',qty,'net',net,'tag',tag);
end

function D = mergeDesign(Dw, Da)
D = Dw;
fF = fieldnames(Da.global);
for i = 1:numel(fF), D.global.(fF{i}) = Da.global.(fF{i}); end
for q = fieldnames(Da.nodal).'
    qty = q{1};
    bF = fieldnames(Da.nodal.(qty));
    for i = 1:numel(bF), D.nodal.(qty).(bF{i}) = Da.nodal.(qty).(bF{i}); end
end
end

function v = getCol(Dc, p, netIdx)
cb = sprintf('%s_%s', p.conn, p.band);
switch p.level
    case 'global'
        v = Dc.global.(sprintf('%s_%s_%s', p.conn, p.band, p.qty));
    case 'network'
        M = Dc.nodal.(p.qty).(cb);
        v = mean(M(:, netIdx == p.net), 2, 'omitnan');
    case 'node'
        v = Dc.nodal.(p.qty).(cb)(:, p.net);
end
v = v(:);
end

function o = ols(y, X)
n = size(X,1); p = size(X,2);
b = X \ y; r = y - X*b;
ssres = r'*r; sstot = sum((y-mean(y)).^2);
r2 = 1 - ssres/max(sstot,eps);
df = n - p; s2 = ssres/max(df,1);
try covb = s2*inv(X'*X); se = sqrt(max(diag(covb),0)'); catch, se = nan(p,1); end
t = b ./ max(se,eps); pv = 2*tcdf(-abs(t), max(df,1));
o = struct('b',b,'r2',r2,'se',se,'t',t,'p',pv,'df',df);
end

function r = residv(y, X)
r = y - X*(X\y);
end

function q = fdrBH(p)
p = p(:);
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

function [r, p] = parcorrRows(x, y, Z)
% partial corr of column vector x with y, controlling Z (rows already valid)
X = [ones(numel(x),1) Z];
rx = x - X*(X\x); ry = y - X*(X\y);
r = sum(rx.*ry)/sqrt(sum(rx.^2)*sum(ry.^2));
df = numel(x) - 2 - size(Z,2);
t = r*sqrt(df/max(1-r^2,eps));
p = 2*tcdf(-abs(t), df);
end

function cols = screenScope(Dc, scope, aqNames, aqTarget, rows, netIdx, NETWORKS)
% Re-run ONE scope's FDR family on `rows` (training subjects); return cell
% array of full-length predictor columns whose (candidate x aqTarget) test
% passes q<0.05. Mirrors the family definitions of step2/3/4.
[level, conn, band, qty, net] = deal(scope{1}, scope{2}, scope{3}, scope{4}, scope{5});
Z = [Dc.age(rows), Dc.gender(rows)];
cb = sprintf('%s_%s', conn, band);
cols = {};
switch level
    case 'global'
        x = Dc.global.(sprintf('%s_%s_%s', conn, band, qty));
        ps = nan(numel(aqNames),1);
        for a = 1:numel(aqNames)
            [~, ps(a)] = parcorrRows(x(rows), Dc.(aqNames{a})(rows), Z);
        end
        q = fdrBH(ps);   % family: metric x band, pool = 6 AQ
        ia = find(strcmp(aqNames, aqTarget));
        if q(ia) < 0.05, cols{end+1} = x; end

    case 'network'
        M = Dc.nodal.(qty).(cb);
        ps = nan(numel(NETWORKS), numel(aqNames));
        for k = 1:numel(NETWORKS)
            v = mean(M(:, netIdx == k), 2, 'omitnan');
            for a = 1:numel(aqNames)
                [~, ps(k,a)] = parcorrRows(v(rows), Dc.(aqNames{a})(rows), Z);
            end
        end
        q = fdrBH(ps(:));   % family: qty x band, pool = 7 nets x 6 AQ
        q = reshape(q, size(ps));
        ia = find(strcmp(aqNames, aqTarget));
        if q(net, ia) < 0.05
            cols{end+1} = mean(M(:, netIdx == net), 2, 'omitnan');
        end

    case 'node'
        M = Dc.nodal.(qty).(cb);
        nodesK = find(netIdx == net);
        ps = nan(numel(nodesK), numel(aqNames));
        for j = 1:numel(nodesK)
            for a = 1:numel(aqNames)
                [~, ps(j,a)] = parcorrRows(M(rows, nodesK(j)), Dc.(aqNames{a})(rows), Z);
            end
        end
        q = fdrBH(ps(:));   % family: qty x net x band, pool = nodes x 6 AQ
        q = reshape(q, size(ps));
        ia = find(strcmp(aqNames, aqTarget));
        hit = find(q(:, ia) < 0.05);
        for j = hit.'
            cols{end+1} = M(:, nodesK(j));
        end
end
end

function [q2_base, q2_unit] = nestedLOO(Dc, scopes, aqTarget, age, gender01, netIdx)
% Honest predictive q2: per fold, re-screen scopes on the training subjects,
% PCA-fuse selected predictors (z-scored in training), fit M1/M2 on training,
% predict the held-out subject. Folds with no selected predictor fall back
% to the baseline model for that subject.
aqNames = Dc.aqNames; NETWORKS = Dc.NETWORKS;
y = Dc.(aqTarget)(:);
ok = ~isnan(y) & ~isnan(age) & ~isnan(gender01);
idxAll = find(ok);
yh = y(idxAll); ah = age(idxAll); gh = gender01(idxAll);
n = numel(yh);
pred_base = nan(n,1); pred_unit = nan(n,1);
for i = 1:n
    tr = setdiff((1:n)', i);
    X0tr = [ones(numel(tr),1) ah(tr) gh(tr)];
    b0 = X0tr \ yh(tr);
    pred_base(i) = [1 ah(i) gh(i)] * b0;

    cols = {};
    for s = 1:numel(scopes)
        cols = [cols screenScope(Dc, scopes{s}, aqNames, aqTarget, idxAll(tr), netIdx, NETWORKS)]; %#ok<AGROW>
    end
    if isempty(cols)
        pred_unit(i) = pred_base(i); continue;
    end
    Ptr = nan(numel(tr), numel(cols)); Pte = nan(1, numel(cols));
    for k = 1:numel(cols)
        v = cols{k}(idxAll);
        mu = mean(v(tr),'omitnan'); sd = std(v(tr),'omitnan'); if sd < eps, sd = 1; end
        Ptr(:,k) = (v(tr)-mu)/sd; Pte(k) = (v(i)-mu)/sd;
    end
    good = ~any(isnan(Ptr),1) & ~isnan(Pte);
    Ptr = Ptr(:, good); Pte = Pte(:, good);
    if isempty(Ptr)
        pred_unit(i) = pred_base(i); continue;
    end
    if size(Ptr,2) > 1
        [coef, score] = pca(Ptr);
        br_tr = score(:,1);
        br_te = (Pte - mean(Ptr,1)) * coef(:,1);
        if corr(br_tr, yh(tr)) < 0, br_tr = -br_tr; br_te = -br_te; end
    else
        br_tr = Ptr(:,1); br_te = Pte(1);
        if corr(br_tr, yh(tr)) < 0, br_tr = -br_tr; br_te = -br_te; end
    end
    X1tr = [X0tr br_tr];
    b1 = X1tr \ yh(tr);
    pred_unit(i) = [1 ah(i) gh(i) br_te] * b1;
end
SST = sum((yh - mean(yh)).^2);
q2_base = 1 - sum((yh - pred_base).^2)/SST;
q2_unit = 1 - sum((yh - pred_unit).^2)/SST;
end

function [maxD, dR2_excl, p_excl] = cookCheck(y, X0, brain)
X1 = [X0 brain];
n = numel(y); p = size(X1,2);
b = X1\y; r = y - X1*b; mse = (r'*r)/(n-p);
H = X1*((X1'*X1)\X1'); h = diag(H);
D = (r.^2/(p*mse)).*(h./max(1-h,eps).^2);
[maxD, imax] = max(D);
keep = true(n,1); keep(imax) = false;
o0 = ols(y(keep), X0(keep,:)); o1 = ols(y(keep), X1(keep,:));
dR2_excl = o1.r2 - o0.r2; p_excl = o1.p(end);
end
