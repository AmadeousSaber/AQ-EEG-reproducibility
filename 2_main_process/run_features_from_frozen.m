function run_features_from_frozen(cond, forceRecompute)
%run_features_from_frozen Recompute source/connectivity without rerunning ICA.
%
% This optional audit path starts from the frozen post-ICA FieldTrip files
% under data/frozen_derivatives/<condition>/sub-*/eeg/*_eeg.mat. It writes
% only to work/recomputed_features and never modifies frozen data.
%
% The canonical analysis does not require this slow step. Use it only to
% audit the post-ICA -> source -> connectivity lineage.

if nargin < 1 || isempty(cond), cond = 'open'; end
if nargin < 2, forceRecompute = false; end
assert(ismember(cond, {'open','closed'}), ...
    'cond must be ''open'' or ''closed''.');
assert(islogical(forceRecompute) && isscalar(forceRecompute), ...
    'forceRecompute must be a logical scalar.');

addpath(fileparts(fileparts(mfilename('fullpath'))));
cfg = pipeline_config('canonical');
assert(isfolder(cfg.toolkitDir), 'DISCOVER-EEG not found: %s', cfg.toolkitDir);
assert(isfolder(cfg.fieldtripPath), 'FieldTrip not found: %s', cfg.fieldtripPath);

addpath(cfg.fieldtripPath);
ft_defaults;
addpath(genpath(fullfile(cfg.toolkitDir, 'custom_functions')));
addpath(cfg.bctDir);

inputRoot = fullfile(cfg.frozenDerivativeRoot, cond);
paramsFile = fullfile(inputRoot, 'params.json');
assert(isfile(paramsFile), 'Missing frozen parameter file: %s', paramsFile);
params = jsondecode(fileread(paramsFile));

outputRoot = fullfile(cfg.workRoot, 'recomputed_features', cond);
params.PreprocessedDataPath = inputRoot;
params.PowerPath = fullfile(outputRoot, 'EEG_features', 'power');
params.SourcePath = fullfile(outputRoot, 'EEG_features', 'source');
params.ConnectivityPath = fullfile(outputRoot, 'EEG_features', 'connectivity');
params.GraphPath = fullfile(outputRoot, 'EEG_features', 'graph_measures');
params.HeadModelPath = resolveDependencyPath(cfg.toolkitDir, params.HeadModelPath);
params.SurfaceModelPath = resolveDependencyPath(cfg.toolkitDir, params.SurfaceModelPath);
params.AtlasPath = resolveDependencyPath(cfg.toolkitDir, params.AtlasPath);
for pathCell = {params.PowerPath, params.SourcePath, ...
        params.ConnectivityPath, params.GraphPath}
    if ~isfolder(pathCell{1}), mkdir(pathCell{1}); end
end

referenceConn = fullfile(inputRoot, 'EEG_features', 'connectivity');
files = dir(fullfile(referenceConn, '*_dwpli_theta.mat'));
assert(~isempty(files), 'No frozen connectivity files found in %s', referenceConn);
bidsIDs = sort(erase({files.name}, '_dwpli_theta.mat'));
freqs = cfg.bands;

fprintf('Recomputing from frozen post-ICA data: %s (%d subjects)\n', ...
    cond, numel(bidsIDs));
fprintf('Output root: %s\n', outputRoot);
for i = 1:numel(bidsIDs)
    bidsID = bidsIDs{i};
    frozenMat = fullfile(inputRoot, bidsID(1:6), 'eeg', ...
        [bidsID '_eeg.mat']);
    assert(isfile(frozenMat), 'Missing frozen post-ICA data: %s', frozenMat);
    for f = 1:numel(freqs)
        band = freqs{f};
        sourceFile = fullfile(params.SourcePath, ...
            sprintf('%s_source_%s.mat', bidsID, band));
        dwpliFile = fullfile(params.ConnectivityPath, ...
            sprintf('%s_dwpli_%s.mat', bidsID, band));
        aecFile = fullfile(params.ConnectivityPath, ...
            sprintf('%s_aec_%s.mat', bidsID, band));

        if forceRecompute || ~isfile(sourceFile)
            compute_spatial_filter(params, bidsID, band);
        end
        if forceRecompute || ~isfile(dwpliFile)
            compute_dwpli(params, bidsID, band);
        end
        if forceRecompute || ~isfile(aecFile)
            compute_aec(params, bidsID, band);
        end
        fprintf('[%d/%d] %s %s complete\n', ...
            (i-1)*numel(freqs)+f, numel(bidsIDs)*numel(freqs), bidsID, band);
    end
end
fprintf('Frozen post-ICA feature audit complete.\n');
end

function pathOut = resolveDependencyPath(toolkitRoot, pathIn)
pathIn = char(pathIn);
if isfile(pathIn)
    pathOut = pathIn;
elseif isfile(fullfile(toolkitRoot, pathIn))
    pathOut = fullfile(toolkitRoot, pathIn);
else
    pathOut = which(pathIn);
    assert(~isempty(pathOut), 'Dependency file not found: %s', pathIn);
end
end
