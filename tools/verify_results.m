function verify_results
%verify_results Compare generated design matrices and CSVs with references.

cfg = pipeline_config;
nFail = 0;

expectedDesign = cell(0, 1);
expectedCsv = {'plot_data_all.csv'; 'results_regression_summary.csv'};
for cond = cfg.conds
    for conn = cfg.conns
        expectedDesign{end+1, 1} = sprintf('design_data_%s_%s.mat', cond{1}, conn{1}); %#ok<AGROW>
        for level = {'global','network','node','node_uncorrected_sig'}
            expectedCsv{end+1, 1} = sprintf('results_%s_%s_%s.csv', ...
                level{1}, cond{1}, conn{1}); %#ok<AGROW>
        end
    end
end

designFiles = dir(fullfile(cfg.referenceOutputDir, 'design_data_*.mat'));
assert(isequal(sort({designFiles.name}'), sort(expectedDesign)), ...
    ['Full pipeline verification needs four private reference design matrices. ' ...
     'For the public aggregate release, use scripts/verify_public_results.py.']);
for i = 1:numel(designFiles)
    referenceFile = fullfile(cfg.referenceOutputDir, designFiles(i).name);
    generatedFile = fullfile(cfg.corrDir, designFiles(i).name);
    ok = isfile(generatedFile);
    if ok
        R = load(referenceFile, 'D');
        G = load(generatedFile, 'D');
        ok = isequaln(R.D, G.D);
    end
    nFail = nFail + report(ok, ['design ' designFiles(i).name]);
end

csvFiles = dir(fullfile(cfg.referenceOutputDir, '*.csv'));
assert(isequal(sort({csvFiles.name}'), sort(expectedCsv)), ...
    'Full pipeline verification requires the complete private reference set.');
nByteIdentical = 0;
for i = 1:numel(csvFiles)
    referenceFile = fullfile(cfg.referenceOutputDir, csvFiles(i).name);
    generatedFile = fullfile(cfg.corrDir, csvFiles(i).name);
    ok = isfile(generatedFile);
    byteIdentical = false;
    if ok
        byteIdentical = strcmp(fileread(referenceFile), fileread(generatedFile));
        if byteIdentical
            ok = true;
        else
            R = readtable(referenceFile, 'PreserveVariableNames', true, ...
                'TextType', 'string');
            G = readtable(generatedFile, 'PreserveVariableNames', true, ...
                'TextType', 'string');
            ok = isequaln(R, G);
        end
    end
    nByteIdentical = nByteIdentical + double(byteIdentical);
    nFail = nFail + report(ok, ['CSV ' csvFiles(i).name]);
end

fprintf('%d/%d CSV files are byte-identical to the reference set.\n', ...
    nByteIdentical, numel(csvFiles));
if nFail > 0
    error('Result verification failed: %d artifact(s) differ.', nFail);
end
fprintf('All design matrices and statistical CSVs match their references.\n');
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
