function step7_interaction_stability()
% Reproducible targeted post hoc EO communication-by-sex stability audit.
% Requires private design_data_open_dwpli.mat from step1; public aggregates
% cannot substitute for participant-level design inputs. Nothing is uploaded.
% Fixed feature: node 42 (LH_Default_PFC_1), dwPLI alpha clustering.
% Locked 2026-09-10: all 39 deletions, 10000 matched FL permutations/fit,
% seed 20260913; 5000 sex-stratified bootstrap draws, seed 20260914;
% fixed-design interaction-null Gaussian simulation 1000 x 999, seed 20260915.
% Never outcome-flip the reported brain direction. Interaction beta is the
% female-minus-male slope difference per brain SD; AQ remains in score units.
% HC3 uses an approximate residual-df t reference. These are fixed-feature
% sensitivity diagnostics, not selection-adjusted or causal inference.
% No BH across the 39 deletion diagnostics. Retain every case in final models.
% All outputs, including private delete-one traces, remain local. Publication
% requires separate authorization and selection of aggregate-only files.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root,'-begin');
cfg = pipeline_config();
assert(strcmpi(cfg.studentRoot,root), 'Wrong pipeline_config resolved.');
assert(strcmpi(cfg.corrDir,fullfile(root,'4_corr')), ...
    'Expected project-local 4_corr path.');
auditDir = fullfile(cfg.corrDir,'interaction_stability_v1');
if ~isfolder(auditDir), mkdir(auditDir); end
outputNames = {'observed_interaction.csv','full_sample_permutation.csv', ...
    'private_delete_one_trace.csv','delete_one_summary.csv', ...
    'maximum_cook_deletion.csv','bootstrap_summary.csv', ...
    'bootstrap_draws.csv','interaction_null_summary.csv', ...
    'interaction_null_pvalues.csv','interaction_null_disagreements.csv', ...
    'run_metadata.csv'};
for j = 1:numel(outputNames)
    assert(~isfile(fullfile(auditDir,outputNames{j})), ...
        'Refusing to overwrite existing output %s.',outputNames{j});
end
inputPath = fullfile(cfg.corrDir,'design_data_open_dwpli.mat');
assert(isfile(inputPath), ...
    'Private design_data_open_dwpli.mat is required; aggregate results are insufficient.');
source = load(inputPath,'D');
D = source.D;
p = double(D.nodal.auc_cc.dwpli_alpha(:,42));
valid = all(isfinite([D.communication(:),D.age(:),D.gender(:),p]),2);
y = double(D.communication(valid)); y=y(:);
rawBrain = p(valid);
age = double(D.age(valid)); age=age(:);
sex = double(D.gender(valid))-1; sex=sex(:);
n = numel(y);
assert(n==39 && size(rawBrain,2)==1);
assert(all(isfinite([y,rawBrain,age,sex]),'all'));
assert(all(ismember(sex,[0,1])) && all([sum(sex==0),sum(sex==1)]>1));
% The producing step1 source explicitly maps 1=male, 2=female. Preserve that
% mapping, then reconstruct the original outcome-oriented single component
% only to verify fitted-model invariance; all reported coefficients use raw z.
mappingSource = fileread(fullfile(cfg.corrDir,'step1_extract_design.m'));
assert(contains(mappingSource,'1 = male, 2 = female'), ...
    'Re-check source gender mapping before labeling groups.');
originalBrain = (rawBrain-mean(rawBrain))/max(std(rawBrain),eps);
if corr(originalBrain,y)<0, originalBrain=-originalBrain; end
U = struct('brain',originalBrain,'X2', ...
    [ones(n,1),age,sex,originalBrain,sex.*originalBrain]);
selfTest();
started = tic;
[Xr,Xf,z] = designs(age,sex,rawBrain);
observed = fitDiagnostics(y,Xr,Xf);
orientation = sign(z' * U.brain(:));
assert(ismember(orientation,[-1,1]));
assert(norm(U.brain(:)-orientation*z)<1e-10);
qOriginal = basis(U.X2);
assert(abs(sum((y-qOriginal*(qOriginal'*y)).^2)-observed.rssFull)<1e-9, ...
    'Original component orientation must not change fitted model.');
[maxCook,maxCookIndex] = max(observed.cook);
observedTable = table(n,sum(sex==0),sum(sex==1),orientation, ...
    observed.r2Reduced,observed.r2Full,observed.delta,observed.F, ...
    observed.beta(end),observed.beta(end)/std(y), ...
    observed.se(end),observed.t(end),observed.p(end), ...
    observed.classicCI(1),observed.classicCI(2), ...
    observed.hc3SE,observed.hc3t,observed.hc3p, ...
    observed.hc3CI(1),observed.hc3CI(2),observed.df,maxCook, ...
    sum(observed.cook>4/n),max(observed.leverage), ...
    'VariableNames',{'n','male_n','female_n','original_brain_orientation_vs_raw', ...
    'R2_reduced','R2_full','delta_R2_interaction','partial_F', ...
    'beta_interaction_AQ_points_per_brain_SD','beta_interaction_divided_by_AQ_SD', ...
    'conventional_SE','conventional_t','conventional_p','conventional_CI_lo', ...
    'conventional_CI_hi','HC3_SE','HC3_t','HC3_p_approx','HC3_CI_lo', ...
    'HC3_CI_hi','full_residual_df','max_Cook_D','count_Cook_D_over_4_over_n', ...
    'max_leverage'});

% Observed full sample first, then each deletion in private original row order.
B = 10000;
rng(20260913,'twister');
fullPerm = permutationP(y,Xr,Xf,B);
fullPermutationTable = table(["paired_delta_R2";"partial_F"], ...
    fullPerm.observed(:),fullPerm.exceedances(:),fullPerm.p(:), ...
    fullPerm.lo(:),fullPerm.hi(:),repmat(B,2,1), ...
    'VariableNames',{'method','observed_statistic','exceedances','p_plus_one', ...
    'tail_probability_Wilson95_lo','tail_probability_Wilson95_hi','permutations'});
deleteRows = cell(n,17);
for d = 1:n
    keep = true(n,1); keep(d)=false;
    [Xrd,Xfd] = designs(age(keep),sex(keep),rawBrain(keep));
    fd = fitDiagnostics(y(keep),Xrd,Xfd);
    pd = permutationP(y(keep),Xrd,Xfd,B);
    deleteRows(d,:) = {d,d==maxCookIndex,observed.cook(d), ...
        sum(keep),fd.beta(end),fd.beta(end)/std(y(keep)),fd.delta,fd.F, ...
        pd.exceedances(1),pd.p(1),pd.lo(1),pd.hi(1), ...
        pd.exceedances(2),pd.p(2),pd.lo(2),pd.hi(2), ...
        sign(fd.beta(end))==sign(observed.beta(end))};
end
deletions = cell2table(deleteRows,'VariableNames',{'private_deletion_order', ...
    'is_maximum_Cook_deletion','original_Cook_D','retained_n', ...
    'beta_interaction_AQ_points_per_brain_SD','beta_interaction_divided_by_AQ_SD', ...
    'delta_R2_interaction','partial_F','paired_delta_exceedances','paired_delta_p', ...
    'paired_delta_tail_Wilson95_lo','paired_delta_tail_Wilson95_hi', ...
    'partial_F_exceedances','partial_F_p','partial_F_tail_Wilson95_lo', ...
    'partial_F_tail_Wilson95_hi','interaction_sign_retained'});
maximumCookDeletion = removevars(deletions(deletions.is_maximum_Cook_deletion,:), ...
    'private_deletion_order');
deletionSummary = table(n,min(deletions.beta_interaction_AQ_points_per_brain_SD), ...
    median(deletions.beta_interaction_AQ_points_per_brain_SD), ...
    max(deletions.beta_interaction_AQ_points_per_brain_SD), ...
    sum(deletions.interaction_sign_retained), ...
    min(deletions.delta_R2_interaction),median(deletions.delta_R2_interaction), ...
    max(deletions.delta_R2_interaction),min(deletions.paired_delta_p), ...
    median(deletions.paired_delta_p),max(deletions.paired_delta_p), ...
    sum(deletions.paired_delta_p<=.05),min(deletions.partial_F_p), ...
    median(deletions.partial_F_p),max(deletions.partial_F_p), ...
    sum(deletions.partial_F_p<=.05), ...
    'VariableNames',{'deletions','beta_min','beta_median','beta_max', ...
    'same_direction_count','delta_R2_min','delta_R2_median','delta_R2_max', ...
    'paired_delta_p_min','paired_delta_p_median','paired_delta_p_max', ...
    'paired_delta_p_le_05_count','partial_F_p_min','partial_F_p_median', ...
    'partial_F_p_max','partial_F_p_le_05_count'});
fprintf('All %d deletion diagnostics completed; elapsed %.1f s.\n',n,toc(started));

% Fixed-sex-count resampling; never outcome-align the brain sign.
Bboot = 5000;
rng(20260914,'twister');
male = find(sex==0); female = find(sex==1);
boot = nan(Bboot,3);
zeroScale = 0; rankDeficient = 0; degenerateOutcome = 0;
for b = 1:Bboot
    take = [male(randi(numel(male),numel(male),1)); ...
        female(randi(numel(female),numel(female),1))];
    by = y(take); bp = rawBrain(take);
    if std(bp)<=eps, zeroScale=zeroScale+1; continue; end
    if std(by)<=eps, degenerateOutcome=degenerateOutcome+1; continue; end
    [Br,Bf] = designs(age(take),sex(take),bp);
    if rank(Br)<size(Br,2) || rank(Bf)<size(Bf,2)
        rankDeficient=rankDeficient+1; continue;
    end
    coefficients = Bf\by;
    boot(b,:) = [coefficients(end),coefficients(end)/std(by), ...
        incrementalDelta(by,Br,Bf)];
end
good = all(isfinite(boot),2);
assert(sum(good)>0,'No valid bootstrap draws.');
ci = prctile(boot(good,1),[2.5 97.5]);
ciFully = prctile(boot(good,2),[2.5 97.5]);
bootstrapSummary = table(Bboot,sum(good),zeroScale,rankDeficient,degenerateOutcome, ...
    observed.beta(end),ci(1),ci(2), ...
    mean(sign(boot(good,1))==sign(observed.beta(end))), ...
    mean(boot(good,1)>0),mean(boot(good,1)<0), ...
    ciFully(1),ciFully(2), ...
    'VariableNames',{'bootstrap_draws','valid_draws','zero_brain_scale_draws', ...
    'rank_deficient_draws','degenerate_outcome_draws','observed_beta_brain_standardized', ...
    'percentile95_beta_lo','percentile95_beta_hi','descriptive_sign_retention', ...
    'fraction_beta_positive','fraction_beta_negative', ...
    'percentile95_beta_divided_by_AQ_SD_lo','percentile95_beta_divided_by_AQ_SD_hi'});
bootstrapDraws = table((1:Bboot)',good,boot(:,1),boot(:,2),boot(:,3), ...
    'VariableNames',{'replicate','valid','beta_brain_standardized', ...
    'beta_divided_by_AQ_SD','delta_R2_interaction'});
fprintf('Bootstrap %d/%d valid; elapsed %.1f s.\n',sum(good),Bboot,toc(started));

% Focal interaction-null calibration retains age, sex AND brain main effect.
S = 1000; Bnull = 999;
rng(20260915,'twister');
Qr = basis(Xr); Qf = basis(Xf);
nullMean = Xr*(Xr\y);
sigma = sqrt(sum((y-nullMean).^2)/(n-size(Xr,2)));
nullP = nan(S,3);
for s = 1:S
    simY = nullMean+sigma*randn(n,1);
    [obs,perms] = commonPermutationStatistics(simY,Qr,Qf,Bnull);
    nullP(s,:) = (1+sum(perms>=obs,1))/(Bnull+1);
    if mod(s,250)==0
        fprintf('Interaction-null calibration %d/%d; elapsed %.1f s.\n',s,S,toc(started));
    end
end
method = ["original_constant_reduced_R2";"paired_delta_R2";"partial_F"];
reject = nullP<=.05;
counts = sum(reject,1)';
[lo,hi] = wilson95(counts,S);
nullSummary = table(method,repmat(S,3,1),repmat(Bnull,3,1), ...
    counts,counts/S,lo,hi,repmat(sigma,3,1), ...
    'VariableNames',{'method','simulations','inner_permutations','rejections', ...
    'rejection_rate','Wilson95_lo','Wilson95_hi','generating_residual_SD'});
nullPvalues = table((1:S)',nullP(:,1),nullP(:,2),nullP(:,3), ...
    'VariableNames',{'simulation','p_original','p_paired_delta_R2','p_partial_F'});
pairs = [1 2;1 3;2 3]; pairRows=cell(3,7);
for k=1:3
    a=pairs(k,1); b=pairs(k,2); ra=reject(:,a); rb=reject(:,b);
    pairRows(k,:)={method(a),method(b),sum(ra&rb),sum(~ra&~rb), ...
        sum(ra&~rb),sum(~ra&rb),sum(ra~=rb)};
end
nullDisagreements=cell2table(pairRows,'VariableNames',{'method_a','method_b', ...
    'both_reject','neither_reject','a_only_reject','b_only_reject','disagreements'});
metadata = table(string(version),20260913,B,20260914,Bboot,20260915,S,Bnull, ...
    "0=male, 1=female; original toGender code 1=male,2=female minus 1", ...
    "raw node42 dwPLI alpha cc positive orientation; brain re-z-scored in retained/resampled sample", ...
    "female-minus-male slope difference; AQ points per sample brain SD; age slope shared", ...
    "targeted post hoc fixed-feature audit, not discovery-adjusted or causal", ...
    "HC3 t with n-rank(Xfull) df is approximate diagnostic inference", ...
    "delete-one p values are unadjusted sensitivity diagnostics, not 39 discoveries", ...
    "interaction null mean retains X1 age+sex+brain; iid Gaussian errors only", ...
    toc(started), 'VariableNames',{'matlab_version','deletion_seed','permutations', ...
    'bootstrap_seed','bootstrap_draws','interaction_null_seed','null_simulations', ...
    'null_inner_permutations','sex_mapping','brain_direction_and_scaling', ...
    'coefficient_interpretation','scope','HC3_reference', ...
    'deletion_inference_scope','null_scope','elapsed_seconds'});
tables = {observedTable,fullPermutationTable,deletions,deletionSummary, ...
    maximumCookDeletion,bootstrapSummary,bootstrapDraws,nullSummary,nullPvalues, ...
    nullDisagreements,metadata};
for j=1:numel(tables)
    writetable(tables{j},fullfile(auditDir,outputNames{j}));
end
fprintf('Focal private interaction audit finished. Originals unchanged.\n');
end

function [Xr,Xf,z] = designs(age,sex,brain)
scale=std(brain);
assert(scale>eps,'Degenerate brain scale.');
z=(brain-mean(brain))/scale;
Xr=[ones(numel(brain),1),age,sex,z];
Xf=[Xr,sex.*z];
end

function Q = basis(X)
assert(rank(X)==size(X,2),'Rank-deficient design.');
[Q,~]=qr(X,0);
end

function result = fitDiagnostics(y,Xr,Xf)
Qr=basis(Xr); Qf=basis(Xf); n=numel(y); p=size(Xf,2);
beta=Xf\y;
res=y-Xf*beta;
rssFull=sum(res.^2); rssReduced=sum((y-Qr*(Qr'*y)).^2);
sst=sum((y-mean(y)).^2); df=n-p; mse=rssFull/df;
% Stable inverse via QR, with independent sandwich test below.
[~,R]=qr(Xf,0); inverseR=R\eye(p); bread=inverseR*inverseR';
covariance=mse*bread;
se=sqrt(diag(covariance)); t=beta./se;
leverage=sum(Qf.^2,2);
assert(all(leverage<1),'Unit leverage invalidates HC3.');
weights=res.^2./(1-leverage).^2;
HC3=bread*(Xf'*(Xf.*weights))*bread;
hc3SE=sqrt(max(HC3(end,end),0));
hc3t=beta(end)/hc3SE;
crit=tinv(.975,df);
result=struct('beta',beta,'se',se,'t',t,'p',2*tcdf(-abs(t),df), ...
    'classicCI',beta(end)+[-1,1]*crit*se(end), ...
    'hc3SE',hc3SE,'hc3t',hc3t,'hc3p',2*tcdf(-abs(hc3t),df), ...
    'hc3CI',beta(end)+[-1,1]*crit*hc3SE,'df',df, ...
    'rssFull',rssFull,'r2Reduced',1-rssReduced/sst,'r2Full',1-rssFull/sst, ...
    'delta',(rssReduced-rssFull)/sst,'F',(rssReduced-rssFull)/(rssFull/df), ...
    'leverage',leverage,'cook',(res.^2/(p*mse)).*leverage./(1-leverage).^2);
end

function value = incrementalDelta(y,Xr,Xf)
Qr=basis(Xr); Qf=basis(Xf);
value=(sum((y-Qr*(Qr'*y)).^2)-sum((y-Qf*(Qf'*y)).^2))/sum((y-mean(y)).^2);
end

function result = permutationP(y,Xr,Xf,B)
[obs,perms]=commonPermutationStatistics(y,basis(Xr),basis(Xf),B);
obs=obs(2:3); perms=perms(:,2:3);
counts=sum(perms>=obs,1);
[lo,hi]=wilson95(counts,B);
result=struct('observed',obs,'exceedances',counts,'p',(counts+1)/(B+1), ...
    'lo',lo,'hi',hi);
end

function [observed,null] = commonPermutationStatistics(y,Qr,Qf,B)
n=numel(y); df=n-size(Qf,2);
assert(size(Qf,2)-size(Qr,2)==1 && df>0);
fit=Qr*(Qr'*y); residual=y-fit;
Y=zeros(n,B);
for b=1:B, Y(:,b)=fit+residual(randperm(n)); end
[r0,r1,F]=statistics(y,Qr,Qf,df);
[r0p,r1p,Fp]=statistics(Y,Qr,Qf,df);
observed=[r1-r0,r1-r0,F];
null=[(r1p-r0)',(r1p-r0p)',Fp'];
end

function [r0,r1,F] = statistics(Y,Qr,Qf,df)
ss0=sum((Y-Qr*(Qr'*Y)).^2,1);
ss1=sum((Y-Qf*(Qf'*Y)).^2,1);
sst=sum((Y-mean(Y,1)).^2,1);
assert(all(sst>eps) && all(ss1>eps));
r0=1-ss0./sst; r1=1-ss1./sst;
F=(ss0-ss1)./(ss1/df);
end

function [lo,hi] = wilson95(k,n)
z=1.959963984540054; p=k/n; den=1+z^2/n;
center=(p+z^2/(2*n))/den;
half=z*sqrt(p.*(1-p)/n+z^2/(4*n^2))/den;
lo=max(0,center-half); hi=min(1,center+half);
end

function selfTest()
t=(1:21)'; sex=mod(t,2); brain=sin(t); age=19+t/10;
[Xr,Xf]=designs(age,sex,brain);
y=.2*t+cos(t/2)+sin(t).*sex;
fit=fitDiagnostics(y,Xr,Xf);
beta=Xf\y; residual=y-Xf*beta;
% A separate direct pseudoinverse sandwich checks the QR-based HC3 formula.
pinvX=pinv(Xf); leverage=diag(Xf*pinvX);
covariance=pinvX*diag((residual./(1-leverage)).^2)*pinvX';
assert(abs(fit.hc3SE-sqrt(covariance(end,end)))<1e-10);
Qr=basis(Xr); Qf=basis(Xf); Y=[y,y+cos(t)];
[r0,r1,F]=statistics(Y,Qr,Qf,numel(y)-size(Xf,2));
ss0=sum((Y-Xr*(Xr\Y)).^2,1); ss1=sum((Y-Xf*(Xf\Y)).^2,1);
sst=sum((Y-mean(Y,1)).^2,1);
assert(max(abs(r0-(1-ss0./sst)))<1e-12);
assert(max(abs(r1-(1-ss1./sst)))<1e-12);
assert(max(abs(F-(ss0-ss1)./(ss1/16)))<1e-10);
assert(abs(fit.delta-incrementalDelta(y,Xr,Xf))<1e-12);
end


