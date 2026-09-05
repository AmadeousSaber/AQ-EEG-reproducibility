function check_setup(requireGraphToolboxes)
%check_setup Validate frozen data, portable paths, and optional toolboxes.

if nargin < 1
    requireGraphToolboxes = false;
end
assert(islogical(requireGraphToolboxes) && isscalar(requireGraphToolboxes), ...
    'requireGraphToolboxes must be a logical scalar.');

cfg = pipeline_config('canonical');
nFail = 0;
fprintf('===== student_pipeline_3.0 setup check =====\n');

% Core portable inputs.
nFail = nFail + report(isfile(cfg.aqXlsx), 'de-identified AQ workbook');
nFail = nFail + report(isfile(cfg.atlasCsv), 'local Schaefer-100 atlas');
nFail = nFail + report(isfolder(cfg.canonicalGraphRoot), ...
    'canonical graph root');
nFail = nFail + report(isfolder(cfg.frozenDerivativeRoot), ...
    'frozen derivative root');
nFail = nFail + report(isfolder(cfg.referenceOutputDir), ...
    'reference output root');

T = readtable(cfg.aqXlsx, 'PreserveVariableNames', true);
expectedAqHeaders = {'subject_id','age','gender','total','social skill', ...
    'attention switching','attention to detail','communication','imagination'};
aqHeadersOk = all(ismember(expectedAqHeaders, T.Properties.VariableNames));
aqContentOk = false;
if height(T) == 41 && aqHeadersOk
    aqIDs = double(T.('subject_id'));
    aqNumeric = T{:, {'age','total','social skill','attention switching', ...
        'attention to detail','communication','imagination'}};
    gender = string(T.('gender'));
    aqContentOk = isequal(aqIDs(:), (1:41)') && ...
        isnumeric(aqNumeric) && all(isfinite(aqNumeric), 'all') && ...
        all(ismember(gender, ["男","女","1","2"]));
end
nFail = nFail + report(height(T) == 41 && aqHeadersOk && aqContentOk, ...
    'AQ workbook: IDs 1-41 and expected de-identified values');

for c = 1:numel(cfg.conds)
    cond = cfg.conds{c};
    frozen = fullfile(cfg.frozenDerivativeRoot, cond);
    postIcaMat = dir(fullfile(frozen, 'sub-*', 'eeg', '*_eeg.mat'));
    postIcaSet = dir(fullfile(frozen, 'sub-*', 'eeg', '*_eeg.set'));
    postIcaFdt = dir(fullfile(frozen, 'sub-*', 'eeg', '*_eeg.fdt'));
    matStems = sort(erase({postIcaMat.name}, '.mat'));
    setStems = sort(erase({postIcaSet.name}, '.set'));
    fdtStems = sort(erase({postIcaFdt.name}, '.fdt'));
    sourceFiles = dir(fullfile(frozen, 'EEG_features', 'source', ...
        '*_source_*.mat'));
    connFiles = dir(fullfile(frozen, 'EEG_features', 'connectivity', '*.mat'));
    graphFiles = dir(fullfile(cfg.canonicalGraphRoot, cond, ...
        'EEG_features', 'graph_measures', '*.mat'));

    nFail = nFail + report(numel(postIcaMat) == 41 && ...
        numel(postIcaSet) == 41 && numel(postIcaFdt) == 41 && ...
        isequal(matStems, setStems, fdtStems), ...
        sprintf('%s frozen post-ICA EEG: 41 MAT/SET/FDT triplets', cond));
    nFail = nFail + report(numel(sourceFiles) == 164, ...
        sprintf('%s frozen source matrices: 164', cond));
    nFail = nFail + report(numel(connFiles) == ...
        cfg.expectedConnectivityFiles.(cond), ...
        sprintf('%s frozen connectivity files: %d', cond, ...
        cfg.expectedConnectivityFiles.(cond)));
    nFail = nFail + report(numel(graphFiles) == ...
        cfg.expectedGraphFiles.(cond), ...
        sprintf('%s canonical graph files: %d', cond, ...
        cfg.expectedGraphFiles.(cond)));
end

sampleGraph = fullfile(cfg.canonicalGraphRoot, 'open', 'EEG_features', ...
    'graph_measures', 'sub-01_task-rest_graph_dwpli_theta.mat');
if isfile(sampleGraph)
    S = load(sampleGraph, 'graph_measures');
    fields = {'auc_gcc','auc_geff','auc_cpl','auc_sw','auc_degree','auc_cc'};
    graphOk = isfield(S, 'graph_measures') && ...
        all(isfield(S.graph_measures, fields)) && ...
        isequal(double(S.graph_measures.sparsity_range(:))', ...
        cfg.sparsityRange);
else
    graphOk = false;
end
nFail = nFail + report(graphOk, ...
    'canonical graph schema and true 15-30% threshold labels');

referenceCsv = dir(fullfile(cfg.referenceOutputDir, '*.csv'));
nFail = nFail + report(numel(referenceCsv) == 18, ...
    'reference statistical outputs: 18 CSV files');
nFail = nFail + report(exist('tcdf', 'file') == 2 && ...
    exist('pca', 'file') == 2 && exist('prctile', 'file') == 2, ...
    'MATLAB Statistics and Machine Learning Toolbox functions');

if requireGraphToolboxes
    nFail = nFail + report(isfolder(cfg.toolkitDir), ...
        sprintf('external DISCOVER-EEG: %s', cfg.toolkitDir));
    nFail = nFail + report(isfolder(cfg.bctDir), ...
        sprintf('external BCT: %s', cfg.bctDir));
    if isfolder(cfg.toolkitDir) && isfolder(cfg.bctDir)
        addpath(genpath(fullfile(cfg.toolkitDir, 'custom_functions')));
        addpath(cfg.bctDir);
        rng(1, 'twister');
        W = rand(100); W = (W + W') / 2; W(1:101:end) = 0;
        [A, ~] = graph_threshold_network(W, 0.15);
        realizedDensity = mean(sum(A)) / 99;
        nFail = nFail + report(abs(realizedDensity - 0.15) < 0.005, ...
            sprintf('patched density function: requested .15, got %.4f', ...
            realizedDensity));
    end
end

if nFail > 0
    error('student_pipeline_3.0 setup check failed: %d issue(s).', nFail);
end
fprintf('Setup check passed. ICA will not be rerun.\n');
end

function n = report(ok, label)
if ok
    fprintf('[OK]   %s\n', label);
    n = 0;
else
    fprintf('[FAIL] %s\n', label);
    n = 1;
end
end
