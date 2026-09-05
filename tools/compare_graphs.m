function compare_graphs(tolerance)
%compare_graphs Compare regenerated work graphs with canonical graph inputs.

if nargin < 1, tolerance = 1e-12; end
cfg = pipeline_config('work');
nFiles = 0;
maxAbsDiff = 0;
failures = strings(0,1);

for c = 1:numel(cfg.conds)
    cond = cfg.conds{c};
    canonicalDir = fullfile(cfg.canonicalGraphRoot, cond, ...
        'EEG_features', 'graph_measures');
    workDir = fullfile(cfg.workGraphRoot, cond, ...
        'EEG_features', 'graph_measures');
    files = dir(fullfile(canonicalDir, '*.mat'));
    for i = 1:numel(files)
        canonicalFile = fullfile(canonicalDir, files(i).name);
        workFile = fullfile(workDir, files(i).name);
        nFiles = nFiles + 1;
        if ~isfile(workFile)
            failures(end+1) = "missing: " + string(workFile); %#ok<AGROW>
            continue;
        end
        A = load(canonicalFile, 'graph_measures');
        B = load(workFile, 'graph_measures');
        fieldsA = sort(fieldnames(A.graph_measures));
        fieldsB = sort(fieldnames(B.graph_measures));
        if ~isequal(fieldsA, fieldsB)
            failures(end+1) = string(files(i).name) + ...
                " graph_measures schema differs"; %#ok<AGROW>
            continue;
        end
        fields = fieldsA;
        for f = 1:numel(fields)
            field = fields{f};
            x = A.graph_measures.(field);
            y = B.graph_measures.(field);
            if ~isequal(size(x), size(y))
                failures(end+1) = string(files(i).name) + ...
                    " size mismatch " + field; %#ok<AGROW>
                continue;
            end
            if ~isnumeric(x) && ~islogical(x)
                if ~isequaln(x, y)
                    failures(end+1) = string(files(i).name) + ...
                        " differs in " + field; %#ok<AGROW>
                end
                continue;
            end
            x = double(x);
            y = double(y);
            d = abs(x - y);
            sameNaN = isnan(x) & isnan(y);
            sameInf = isinf(x) & isinf(y) & sign(x) == sign(y);
            d(sameNaN | sameInf) = 0;
            if any(xor(isnan(x), isnan(y)), 'all') || ...
                    any(xor(isinf(x), isinf(y)), 'all') || ...
                    any(d > tolerance, 'all')
                failures(end+1) = string(files(i).name) + ...
                    " differs in " + field; %#ok<AGROW>
            end
            if ~isempty(d)
                maxAbsDiff = max(maxAbsDiff, max(d, [], 'all', 'omitnan'));
            end
        end
    end
end

fprintf('Compared %d graph files; max absolute numeric difference = %.3g\n', ...
    nFiles, maxAbsDiff);
if ~isempty(failures)
    disp(failures(1:min(20,end)));
    error('Graph comparison failed for %d field/file checks.', ...
        numel(failures));
end
fprintf('All regenerated graph_measures fields match canonical files.\n');
end
