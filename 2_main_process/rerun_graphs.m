function rerun_graphs(forceRecompute)
%rerun_graphs Regenerate graph measures from frozen connectivity matrices.
%
% The input connectivity matrices are copied into this pipeline and are
% unaffected by the historical density bug. Outputs are written only to
% work/graph_results, leaving canonical graph files immutable.
%
% rerun_graphs()       skips existing work files.
% rerun_graphs(true)   recomputes and overwrites work files.

if nargin < 1
    forceRecompute = false;
end
assert(islogical(forceRecompute) && isscalar(forceRecompute), ...
    'forceRecompute must be a logical scalar.');

addpath(fileparts(fileparts(mfilename('fullpath'))));
cfg = pipeline_config('work');
assert(isfolder(cfg.toolkitDir), 'DISCOVER-EEG not found: %s', cfg.toolkitDir);
assert(isfolder(cfg.bctDir), 'BCT not found: %s', cfg.bctDir);

addpath(genpath(fullfile(cfg.toolkitDir, 'custom_functions')));
addpath(cfg.bctDir);

% Functional guard against the historical ~2x density bug.
rng(1, 'twister');
W = rand(100); W = (W + W') / 2; W(1:101:end) = 0;
[A, ~] = graph_threshold_network(W, 0.15);
realizedDensity = mean(sum(A)) / 99;
assert(abs(realizedDensity - 0.15) < 0.005, ...
    ['graph_threshold_network failed the true-density check: requested ' ...
     '0.15, realized %.4f. Use the patched toolkit.'], realizedDensity);

seedFile = fullfile(cfg.manifestDir, 'graph_seed_manifest.csv');
assert(isfile(seedFile), 'Missing graph seed manifest: %s', seedFile);
S = readtable(seedFile, 'TextType', 'string');
required = {'condition','subject_id','connectivity','band','seed','filename'};
required = [required, {'reproduction_mode','historical_task_index', ...
    'rng_algorithm'}];
assert(all(ismember(required, S.Properties.VariableNames)), ...
    'graph_seed_manifest.csv has an unexpected schema.');

tasks = cell(0, 6);
copiedCanonical = 0;
for i = 1:height(S)
    cond = char(S.condition(i));
    bidsID = char(S.subject_id(i));
    connMeasure = char(S.connectivity(i));
    freqBand = char(S.band(i));
    inputDir = fullfile(cfg.frozenDerivativeRoot, cond, ...
        'EEG_features', 'connectivity');
    outputDir = fullfile(cfg.workGraphRoot, cond, ...
        'EEG_features', 'graph_measures');
    if ~isfolder(outputDir), mkdir(outputDir); end

    inputFile = fullfile(inputDir, sprintf('%s_%s_%s.mat', ...
        bidsID, connMeasure, freqBand));
    outputFile = fullfile(outputDir, char(S.filename(i)));
    assert(isfile(inputFile), 'Missing frozen connectivity file: %s', inputFile);
    if S.reproduction_mode(i) == "copy_canonical"
        canonicalFile = fullfile(cfg.canonicalGraphRoot, cond, ...
            'EEG_features', 'graph_measures', char(S.filename(i)));
        assert(isfile(canonicalFile), ...
            'Missing frozen stochastic reference: %s', canonicalFile);
        if forceRecompute || ~isfile(outputFile)
            copyfile(canonicalFile, outputFile, 'f');
            copiedCanonical = copiedCanonical + 1;
        end
        continue;
    end
    assert(S.reproduction_mode(i) == "recompute", ...
        'Unknown reproduction mode in row %d.', i);
    assert(~isnan(S.seed(i)), 'Missing seed in manifest row %d.', i);
    if isfile(outputFile) && ~forceRecompute
        continue;
    end
    tasks(end+1, :) = {inputDir, outputDir, bidsID, freqBand, ...
        connMeasure, S.seed(i)}; %#ok<AGROW>
end

nTasks = size(tasks, 1);
fprintf('Graph source: frozen connectivity inside student_pipeline_3.0\n');
fprintf('Output root: %s\n', cfg.workGraphRoot);
fprintf('Canonical stochastic files copied: %d\n', copiedCanonical);
fprintf('Tasks to run: %d of %d\n', nTasks, height(S));
if nTasks == 0
    fprintf('All work graph files already exist. Nothing to do.\n');
    return;
end

pool = gcp('nocreate');
if isempty(pool)
    parpool('local', min(12, feature('numcores')));
end

tAll = tic;
parfor t = 1:nTasks
    addpath(genpath(fullfile(cfg.toolkitDir, 'custom_functions')));
    addpath(cfg.bctDir);
    [connPath, graphPath, bidsID, freqBand, connMeasure, seed] = ...
        tasks{t, :};

    params = struct();
    params.ConnectivityPath = connPath;
    params.GraphPath = graphPath;
    params.SparsityRange = cfg.sparsityRange;
    params.ConnMatrixThreshold = 0.3;
    params.NRandNetworks = cfg.nRandNetworks;
    params.UseWeightedGraph = true;
    params.UseMST = false;

    % Match the historical parfor worker generator. Do not force Twister:
    % rng(seed) resets the worker's default parallel generator and is what
    % the authoritative corrected-graph script used.
    rng(seed);
    t0 = tic;
    compute_graph_measures(params, bidsID, freqBand, connMeasure);
    fprintf('[%d/%d] %s %s %s seed=%d (%.1fs)\n', ...
        t, nTasks, bidsID, connMeasure, freqBand, seed, toc(t0));
end
fprintf('Completed %d graph files in %.1f min.\n', nTasks, toc(tAll)/60);
end
