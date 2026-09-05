function cfg = pipeline_config(graphSource)
%pipeline_config Self-contained data paths for student_pipeline_3.0.
%
% Data and analysis inputs are resolved relative to this folder. Only the
% external software/toolbox paths may need editing or environment variables.
%
% graphSource:
%   'canonical' (default) - use frozen, density-corrected graph files.
%   'work'                - use graphs regenerated into work/graph_results.

if nargin < 1 || isempty(graphSource)
    graphSource = getenv('AQ_PIPELINE_GRAPH_SOURCE');
end
if isempty(graphSource)
    graphSource = 'canonical';
end
graphSource = lower(char(graphSource));
assert(ismember(graphSource, {'canonical','work'}), ...
    'graphSource must be ''canonical'' or ''work''.');

cfg = struct();
cfg.pipelineVersion = 'student_pipeline_3.0';
cfg.freezeDate = '2026-08-18';
cfg.studentRoot = fileparts(mfilename('fullpath'));

% ---- Frozen project data (portable, relative paths) --------------------
cfg.dataRoot = fullfile(cfg.studentRoot, 'data');
cfg.frozenDerivativeRoot = fullfile(cfg.dataRoot, 'frozen_derivatives');
cfg.connRoot = cfg.frozenDerivativeRoot;
cfg.canonicalGraphRoot = fullfile(cfg.dataRoot, 'canonical_graphs');
cfg.workRoot = fullfile(cfg.studentRoot, 'work');
cfg.workGraphRoot = fullfile(cfg.workRoot, 'graph_results');
cfg.corrDir = fullfile(cfg.studentRoot, '4_corr');
cfg.referenceOutputDir = fullfile(cfg.studentRoot, 'reference_outputs', '4_corr');

cfg.aqXlsx = fullfile(cfg.studentRoot, 'analysis_inputs', ...
    'aq_analysis_deidentified.xlsx');
cfg.atlasCsv = fullfile(cfg.studentRoot, 'analysis_inputs', ...
    'Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv');
cfg.manifestDir = fullfile(cfg.studentRoot, 'manifests');

if strcmp(graphSource, 'canonical')
    cfg.graphRoot = cfg.canonicalGraphRoot;
else
    cfg.graphRoot = cfg.workGraphRoot;
end
cfg.graphSource = graphSource;

% ---- External software only --------------------------------------------
% Environment variables take precedence, which avoids editing this file on
% shared machines. Machine-specific paths are not distributed in this release.
cfg.toolkitDir = envOrDefault('DISCOVER_EEG_HOME', ...
    '');
cfg.bctDir = envOrDefault('BCT_HOME', ...
    '');
cfg.eeglabPath = envOrDefault('EEGLAB_HOME', ...
    '');
cfg.fieldtripPath = envOrDefault('FIELDTRIP_HOME', ...
    '');

% ---- Fixed analysis constants -----------------------------------------
cfg.conds = {'open','closed'};
cfg.conns = {'dwpli','aec'};
cfg.bands = {'theta','alpha','beta','gamma'};
cfg.sparsityRange = [0.15 0.20 0.25 0.30];
cfg.nRandNetworks = 100;
cfg.expectedSubjects = struct('open', 39, 'closed', 41);
cfg.expectedGraphFiles = struct('open', 312, 'closed', 328);
cfg.expectedConnectivityFiles = struct('open', 312, 'closed', 328);
end

function value = envOrDefault(name, defaultValue)
value = getenv(name);
if isempty(value)
    value = defaultValue;
end
end
