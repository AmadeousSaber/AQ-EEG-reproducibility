function run_pipeline(mode)
%run_pipeline One-entry driver for student_pipeline_3.0.
%
% run_pipeline('stats')  Fast, exact reproduction from canonical graphs.
% run_pipeline('graphs') Recompute graphs from frozen connectivity, compare
%                        with canonical graphs, then rerun statistics.
%
% ICA is never rerun by either mode.

if nargin < 1 || isempty(mode)
    mode = 'stats';
end
mode = lower(char(mode));
assert(ismember(mode, {'stats','graphs'}), ...
    'mode must be ''stats'' or ''graphs''.');

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir, fullfile(rootDir, '2_main_process'), ...
    fullfile(rootDir, '4_corr'), fullfile(rootDir, 'tools'));

check_setup(strcmp(mode, 'graphs'));
switch mode
    case 'stats'
        run_all('canonical');
        verify_results;
    case 'graphs'
        rerun_graphs(false);
        compare_graphs;
        run_all('work');
        verify_results;
end
end
