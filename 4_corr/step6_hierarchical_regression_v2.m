function step6_hierarchical_regression_v2(overwrite)
% Seven exploratory, selection-conditioned follow-up models (2026-09-10).
%
% Usage from the project root:
%   addpath('4_corr'); step6_hierarchical_regression_v2
% Explicitly refresh existing v2 outputs:
%   step6_hierarchical_regression_v2(true)
%
% Requires Statistics and Machine Learning Toolbox, pipeline_config.m, the
% atlas CSV, and four authorized design_data_{open,closed}_{aec,dwpli}.mat
% files in cfg.corrDir. Participant-level inputs are not publicly distributed.
% This file does not depend on the work/ audit directory or original step6.
%
% Primary residual-permutation statistic: full-minus-reduced R2, with BOTH
% models refitted to the SAME pseudo-outcome (10,000 permutations).
% Partial F is a companion statistic using exactly the same permutations.
% BH is applied separately across the seven main effects and interactions.
% Bootstrap uses 5,000 paired row resamples and re-estimates predictor
% standardization/PCA within each draw, retaining the original fixed features.
% A percentile interval for nonnegative delta R2 is descriptive uncertainty,
% NOT by itself a zero-effect test, and does not account for feature selection.
% Original conditional LOO scopes, folds, screening and q2 denominator remain
% unchanged; this is not a fully nested validation of the discovery process.
%
% Only filenames containing regression_v2 are written. The historical script
% and outputs remain unchanged. plot_data_regression_v2.mat contains individual
% residuals and must remain private; only aggregate CSV outputs may be shared.
if nargin<1, overwrite=false; end
assert(islogical(overwrite) && isscalar(overwrite),'overwrite must be logical.');
cfg=pipeline_config; corrDir=cfg.corrDir;
outputNames={'results_regression_v2_summary.csv','plot_data_regression_v2.mat', ...
    'results_regression_v2_permutations.csv','results_regression_v2_bootstrap.csv'};
for i=1:numel(outputNames)
    assert(overwrite || ~isfile(fullfile(corrDir,outputNames{i})), ...
        'V2 output exists. Use step6_hierarchical_regression_v2(true) to refresh v2 only.');
end
B=10000; Bboot=5000;
NETWORKS={'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
A=readtable(cfg.atlasCsv,'PreserveVariableNames',true);
netIdx=zeros(100,1);
for i=1:100
    tok=strsplit(A.("ROI Name"){i},'_');
    netIdx(i)=find(strcmp(NETWORKS,tok{3}));
end
NET=@(nm) find(strcmp(NETWORKS,nm));

units = {};
units{end+1} = mkUnit('U1 detail x AEC visual-beta degree + gamma global','detail','exploratory_fixed_units', ...
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
units{end+1} = mkUnit('U2 switching x attention/control node degree','switching','exploratory_fixed_units', ...
    struct('open',  {{ mkP('node','aec','beta','auc_degree',82,'node82(RH_Cont_Par_2):beta:deg') }}, ...
           'closed',{{ mkP('node','dwpli','alpha','auc_degree',89,'node89(RH_Cont_pCun_1):alpha:deg'), ...
                       mkP('node','dwpli','alpha','auc_cc',19,'node19(LH_DorsAttn_Post_4):alpha:cc') }}), ...
    struct('open',  {{ {'node','aec','beta','auc_degree',NET('Cont')} }}, ...
           'closed',{{ {'node','dwpli','alpha','auc_degree',NET('Cont')}, {'node','dwpli','alpha','auc_cc',NET('DorsAttn')} }}));
units{end+1} = mkUnit('U3 communication x node42 / Cont network','communication','exploratory_fixed_units', ...
    struct('open',  {{ mkP('node','dwpli','alpha','auc_cc',42,'node42(LH_Default_PFC_1):alpha:cc') }}, ...
           'closed',{{ mkP('network','aec','beta','auc_degree',NET('Cont'),'net:Cont:beta:deg(AEC)') }}), ...
    struct('open',  {{ {'node','dwpli','alpha','auc_cc',NET('Default')} }}, ...
           'closed',{{ {'network','aec','beta','auc_degree',NET('Cont')} }}));
units{end+1} = mkUnit('U4 imagination x DorsAttn FEF gamma cc','imagination','exploratory_fixed_units', ...
    struct('open',  {{}}, ...
           'closed',{{ mkP('node','dwpli','gamma','auc_cc',23,'node23(LH_DorsAttn_FEF_1):gamma:cc') }}), ...
    struct('open',  {{}}, ...
           'closed',{{ {'node','dwpli','gamma','auc_cc',NET('DorsAttn')} }}));

% Build exactly the original seven fixed units, in the audited stream order.
models={};
for cond={'open','closed'}
    condition=cond{1};
    S=load(fullfile(corrDir,sprintf('design_data_%s_dwpli.mat',condition))); Dw=S.D;
    S=load(fullfile(corrDir,sprintf('design_data_%s_aec.mat',condition))); Da=S.D;
    assert(isequal(Dw.subj_ids,Da.subj_ids),'Subject mismatch.');
    Dc=mergeDesign(Dw,Da); age=Dc.age(:); sex=Dc.gender(:)-1;
    for i=1:numel(units)
        U=units{i}; preds=U.preds.(condition); scopes=U.scopes.(condition);
        if isempty(preds), continue; end
        y=Dc.(U.aqvar)(:);
        P=zeros(numel(y),numel(preds));
        for k=1:numel(preds), P(:,k)=getCol(Dc,preds{k},netIdx); end
        ok=~isnan(y)&~isnan(age)&~isnan(sex)&all(~isnan(P),2);
        y=y(ok); P=P(ok,:); X0=[ones(sum(ok),1),age(ok),sex(ok)];
        assert(all(isfinite([y,P,X0]),'all'),'Nonfinite fixed-unit values.');
        assert(all(ismember(X0(:,3),[0,1])),'Unexpected sex encoding.');
        [brain,coef]=brainComponent(P,y);
        X1=[X0,brain]; X2=[X1,brain.*X0(:,3)];
        label=preds{1}.tag;
        if size(P,2)>1
            label='PC1[';
            for k=1:numel(preds)
                label=[label,sprintf('%s=%.2f ',preds{k}.tag,coef(k))]; %#ok<AGROW>
            end
            label=[label,']'];
        end
        [qb,qu]=conditionalLOO(Dc,scopes,U.aqvar,age,sex,netIdx);
        [maxD,exclDelta,exclP]=cookCheck(y,X0,brain);
        models{end+1}=struct('id',[condition,'_',U.aqvar], ...
            'cond',condition,'unit',U.name,'family',U.family,'aqvar',U.aqvar, ...
            'y',y,'Praw',P,'brain',brain,'label',label,'X0',X0,'X1',X1,'X2',X2, ...
            'o0',ols(y,X0),'o1',ols(y,X1),'o2',ols(y,X2), ...
            'q2base',qb,'q2unit',qu,'maxCookD',maxD, ...
            'deltaExcl',exclDelta,'pExcl',exclP); %#ok<AGROW>
    end
end
assert(numel(models)==7,'Expected exactly seven exploratory fixed units.');

% Audit-identical order: all seven permutation blocks before any bootstrap.
% For each unit, replicate b draws main randperm first, interaction second.
rng(20260910,'twister');
permRows={};
for u=1:7
    M=models{u}; y=M.y; n=numel(y);
    Q={basis(M.X0),basis(M.X1),basis(M.X2)};
    fitted={Q{1}*(Q{1}'*y),Q{2}*(Q{2}'*y)};
    residual={y-fitted{1},y-fitted{2}};
    pseudo={zeros(n,B),zeros(n,B)};
    for b=1:B
        pm=randperm(n); pi=randperm(n);
        pseudo{1}(:,b)=fitted{1}+residual{1}(pm);
        pseudo{2}(:,b)=fitted{2}+residual{2}(pi);
    end
    for effect=1:2
        labels={'main','interaction'}; methods={'paired_delta_R2','partial_F'};
        [observed,null]=pairedStatistics(y,pseudo{effect},Q{effect},Q{effect+1});
        for method=1:2
            count=sum(null(method,:)>=observed(method));
            [lo,hi]=wilson(count,B);
            permRows(end+1,:)={M.id,M.cond,M.aqvar,labels{effect},methods{method}, ...
                n,B,observed(method),count,(count+1)/(B+1),lo,hi}; %#ok<AGROW>
        end
    end
    fprintf('V2 paired permutations complete: %s.\n',M.id);
end
permutations=cell2table(permRows,'VariableNames',{'id','cond','aqvar','effect', ...
    'method','n','B','observed_statistic','exceedances','p_plus_one', ...
    'tail_probability_wilson95_lo','tail_probability_wilson95_hi'});
permutations.q_BH_seven=nan(height(permutations),1);
for effect={'main','interaction'}
    for method={'paired_delta_R2','partial_F'}
        mask=strcmp(permutations.effect,effect{1})&strcmp(permutations.method,method{1});
        assert(sum(mask)==7);
        permutations.q_BH_seven(mask)=fdrBH(permutations.p_plus_one(mask));
    end
end

% Fixed-component comparison is retained as an aggregate reproducibility check;
% the summary reports the re-estimated-component interval.
rng(20260911,'twister');
bootRows={};
for u=1:7
    M=models{u}; n=numel(M.y); draws=nan(Bboot,2);
    baselineBad=0; fixedBad=0; refitBad=0; zeroScale=0; pcTies=0;
    for b=1:Bboot
        ix=randi(n,n,1); y=M.y(ix); X0=M.X0(ix,:);
        if rank(X0)<size(X0,2) || sum((y-mean(y)).^2)<=eps
            baselineBad=baselineBad+1; continue;
        end
        Xfixed=M.X1(ix,:);
        if rank(Xfixed)<size(Xfixed,2)
            fixedBad=fixedBad+1;
        else
            draws(b,1)=deltaR2(y,X0,Xfixed);
        end
        P=M.Praw(ix,:); zeroScale=zeroScale+any(std(P,0,1)<=eps);
        [brain,~,tie]=brainComponent(P,y); pcTies=pcTies+tie;
        Xrefit=[X0,brain];
        if rank(Xrefit)<size(Xrefit,2)
            refitBad=refitBad+1;
        else
            draws(b,2)=deltaR2(y,X0,Xrefit);
        end
    end
    paired=all(isfinite(draws),2); assert(any(paired),'No valid bootstrap draws.');
    fixedCI=prctile(draws(paired,1),[2.5,97.5]);
    refitCI=prctile(draws(paired,2),[2.5,97.5]);
    dif=draws(paired,2)-draws(paired,1);
    bootRows(end+1,:)={M.id,M.cond,M.aqvar,n,size(M.Praw,2),Bboot,sum(paired), ...
        baselineBad,fixedBad,refitBad,zeroScale,pcTies,fixedCI(1),fixedCI(2), ...
        refitCI(1),refitCI(2),mean(dif),max(abs(dif)), ...
        sum(draws(paired,1)<-1e-12),sum(draws(paired,2)<-1e-12)}; %#ok<AGROW>
    fprintf('V2 bootstrap complete: %s (%d/%d paired valid).\n',M.id,sum(paired),Bboot);
end
bootstraps=cell2table(bootRows,'VariableNames',{'id','cond','aqvar','n', ...
    'feature_count','Bboot','paired_valid_draws','baseline_degenerate_draws', ...
    'fixed_component_rank_deficient_draws','refit_component_rank_deficient_draws', ...
    'draws_with_zero_scale_feature','draws_with_PC1_PC2_tie','fixed_component_ci_lo', ...
    'fixed_component_ci_hi','refit_component_ci_lo','refit_component_ci_hi', ...
    'mean_paired_delta_difference','max_abs_paired_delta_difference', ...
    'fixed_negative_draws_beyond_1e_12','refit_negative_draws_beyond_1e_12'});

rows={}; plotS=struct([]);
for u=1:7
    M=models{u}; y=M.y; n=numel(y);
    pm=selectPermutation(permutations,M.id,'main','paired_delta_R2');
    pi=selectPermutation(permutations,M.id,'interaction','paired_delta_R2');
    fm=selectPermutation(permutations,M.id,'main','partial_F');
    fi=selectPermutation(permutations,M.id,'interaction','partial_F');
    boot=bootstraps(strcmp(bootstraps.id,M.id),:); assert(height(boot)==1);
    deltaMain=M.o1.r2-M.o0.r2; deltaInter=M.o2.r2-M.o1.r2;
    sm=NaN; spm=NaN; sf=NaN; spf=NaN;
    % Compatibility descriptive slopes only. This gate is not an additional
    % confirmatory family or a substitute for the corrected interaction test.
    if min(M.o2.p(end),pi.p_plus_one)<0.10
        for sex=0:1
            ix=M.X0(:,3)==sex;
            os=ols(y(ix),[ones(sum(ix),1),M.X0(ix,2),M.brain(ix)]);
            if sex==0, sm=os.b(end); spm=os.p(end);
            else, sf=os.b(end); spf=os.p(end); end
        end
    end
    rows(end+1,:)={M.cond,M.unit,M.family,M.aqvar,n,M.label, ...
        M.o0.r2,M.o1.r2,deltaMain,M.o1.b(end)*std(M.brain)/std(y), ...
        M.o1.p(end),pm.p_plus_one,boot.refit_component_ci_lo,boot.refit_component_ci_hi, ...
        M.o2.r2,deltaInter,M.o2.p(end),pi.p_plus_one,sm,spm,sf,spf, ...
        M.maxCookD,M.deltaExcl,M.pExcl,M.q2base,M.q2unit,M.q2unit-M.q2base, ...
        pm.q_BH_seven,pi.q_BH_seven,fm.p_plus_one,fm.q_BH_seven, ...
        fi.p_plus_one,fi.q_BH_seven,B,Bboot,boot.paired_valid_draws}; %#ok<AGROW>
    plotS(u).unit=M.unit; plotS(u).cond=M.cond; plotS(u).aqvar=M.aqvar;
    plotS(u).comp_label=M.label;
    plotS(u).resid_brain=residv(M.brain,M.X0);
    plotS(u).resid_aq=residv(y,M.X0); plotS(u).gender01=M.X0(:,3);
    plotS(u).dR2m=deltaMain; plotS(u).p_main_perm=pm.p_plus_one;
    plotS(u).q2_unit=M.q2unit; plotS(u).q2_base=M.q2base;
    plotS(u).q_unit=pm.q_BH_seven; plotS(u).q_inter=pi.q_BH_seven;
    plotS(u).dR2_boot_lo=boot.refit_component_ci_lo;
    plotS(u).dR2_boot_hi=boot.refit_component_ci_hi;
end
varn={'cond','unit','family','aqvar','N','brain_component','R2_base','R2_main', ...
    'dR2_main','std_beta','dR2_main_p_asy','dR2_main_p_permFL','dR2_boot_lo','dR2_boot_hi', ...
    'R2_inter','dR2_inter','inter_p_asy','inter_p_permFL','slope_male','slope_male_p', ...
    'slope_female','slope_female_p','maxCookD','dR2_excl_maxD','dR2_excl_maxD_p', ...
    'q2_base','q2_unit','dq2','q_unit','q_inter','dR2_main_p_permF','q_unit_F', ...
    'inter_p_permF','q_inter_F','B_permutations','B_bootstrap','bootstrap_valid_draws'};
T=cell2table(rows,'VariableNames',varn);
metadata=struct('analysisVersion','step6_v2_20260910','matlabVersion',version, ...
    'permutationSeed',20260910,'bootstrapSeed',20260911,'B',B,'Bboot',Bboot, ...
    'primaryStatistic','paired full-minus-reduced R2 on each pseudo-outcome', ...
    'companionStatistic','partial F on identical pseudo-outcomes', ...
    'bootstrapScope','fixed features; standardization and PC1 re-estimated within each draw', ...
    'ciInterpretation','descriptive percentile interval; not a zero-effect or selection-adjusted test', ...
    'validationScope','conditional LOO within full-sample-defined scopes and AQ targets', ...
    'privacy','plotS contains individual residuals; do not publicly distribute');
writetable(T,fullfile(corrDir,outputNames{1}));
save(fullfile(corrDir,outputNames{2}),'plotS','metadata','-v7.3');
writetable(permutations,fullfile(corrDir,outputNames{3}));
writetable(bootstraps,fullfile(corrDir,outputNames{4}));
fprintf('Saved four v2 outputs; historical script and outputs were not changed.\n');
end

function T=selectPermutation(table,id,effect,method)
T=table(strcmp(table.id,id)&strcmp(table.effect,effect)&strcmp(table.method,method),:);
assert(height(T)==1,'Ambiguous permutation result.');
end

function [brain,coef,tie]=brainComponent(P,y)
P=(P-mean(P,1))./max(std(P,0,1),eps); tie=false;
if size(P,2)>1
    [coeff,score,latent]=pca(P);
    if isempty(score), brain=zeros(size(P,1),1); coef=zeros(size(P,2),1);
    else, brain=score(:,1); coef=coeff(:,1); end
    if numel(latent)>1
        tie=abs(latent(1)-latent(2))<=1e-10*max(latent(1),eps);
    end
else
    brain=P(:,1); coef=1;
end
if corr(brain,y)<0, brain=-brain; coef=-coef; end
end

function Q=basis(X)
assert(rank(X)==size(X,2),'Rank-deficient fixed design.');
[Q,~]=qr(X,0);
end

function value=sse(Y,Q)
r=Y-Q*(Q'*Y); value=sum(r.^2,1);
end

function value=deltaR2(y,X0,X1)
value=(sse(y,basis(X0))-sse(y,basis(X1)))/sum((y-mean(y)).^2);
end

function [observed,null]=pairedStatistics(y,Y,Q0,Q1)
n=numel(y); df1=size(Q1,2)-size(Q0,2); df2=n-size(Q1,2);
sst=sum((y-mean(y)).^2); ss0=sse(y,Q0); ss1=sse(y,Q1);
observed=[(ss0-ss1)/sst;((ss0-ss1)/df1)/(ss1/df2)];
sstStar=sum((Y-mean(Y,1)).^2,1);
s0=sse(Y,Q0); s1=sse(Y,Q1);
assert(all(sstStar>eps)&all(s1>eps),'Degenerate pseudo-outcome.');
null=[(s0-s1)./sstStar;((s0-s1)/df1)./(s1/df2)];
end

function [lo,hi]=wilson(k,n)
% Monte Carlo interval for underlying exceedance probability k/B.
z=1.959963984540054; phat=k/n; den=1+z^2/n;
center=(phat+z^2/(2*n))/den;
half=z*sqrt(phat*(1-phat)/n+z^2/(4*n^2))/den;
lo=max(0,center-half); hi=min(1,center+half);
end

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

function o = ols(y,X)
% Column-shaped standard errors; QR solves avoid an inverse of X'*X.
n=size(X,1); p=size(X,2);
assert(rank(X)==p && n>p,'OLS requires a full-rank design with residual df.');
b=X\y; r=y-X*b; ssres=r'*r; sstot=sum((y-mean(y)).^2);
df=n-p; s2=ssres/df;
[~,R]=qr(X,0); Rsolve=R\eye(p);
se=sqrt(max(s2*sum(Rsolve.^2,2),0));
t=b./max(se,eps); pv=2*tcdf(-abs(t),df);
o=struct('b',b,'r2',1-ssres/max(sstot,eps),'se',se,'t',t,'p',pv,'df',df);
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

function [q2_base, q2_unit] = conditionalLOO(Dc, scopes, aqTarget, age, gender01, netIdx)
% Conditional internal q2: full-sample-defined scopes and target stay fixed.
% This is NOT validation of the complete discovery/selection process.
% Per fold, re-screen those scopes on the training subjects,
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
[Q,~] = qr(X1,0); h = sum(Q.^2,2);
D = (r.^2/(p*mse)).*(h./max(1-h,eps).^2);
[maxD, imax] = max(D);
keep = true(n,1); keep(imax) = false;
o0 = ols(y(keep), X0(keep,:)); o1 = ols(y(keep), X1(keep,:));
dR2_excl = o1.r2 - o0.r2; p_excl = o1.p(end);
end


