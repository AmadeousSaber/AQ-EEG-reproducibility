% step1_extract_design.m
% -------------------------------------------------------------------------
% student_pipeline_3.0: design matrices for ALL four combinations
%   condition    : open + closed
%   connectivity : dwpli + aec
%
% Source: graph root selected by pipeline_config ('canonical' or 'work').
% Canonical files are frozen at verified true 15-30% densities.
%
% Output (per combo): design_data_<cond>_<conn>.mat
%   D.global.<conn_band_qty>                     (N x 1)
%   D.nodal.<auc_degree|auc_cc>.<conn_band>      (N x 100)
%   D.age / D.gender / D.<aqvar> / D.subj_ids
%   D.netIdx (100x1 Schaefer-7 membership), D.NETWORKS, D.aqNames, D.bands
% -------------------------------------------------------------------------

cfg      = pipeline_config;     % 学生管线：路径集中在 pipeline_config.m
outDir   = cfg.corrDir;
atlasCsv = cfg.atlasCsv;
aqXlsx   = cfg.aqXlsx;

conds   = cfg.conds;
conns   = cfg.conns;
bands   = cfg.bands;
globalQ = {'auc_gcc','auc_geff','auc_cpl','auc_sw'};
nodalQ  = {'auc_degree','auc_cc'};
aqNames = {'total','social','switching','detail','communication','imagination'};

% --- AQ table (shared across combos) ---
T = readtable(aqXlsx, 'PreserveVariableNames', true);
hdr = T.Properties.VariableNames;
COL = struct('id','subject_id','age','age','gender','gender','total','total', ...
    'social','social skill','switching','attention switching', ...
    'detail','attention to detail','communication','communication', ...
    'imagination','imagination');
ci = struct();
for f = fieldnames(COL).'
    k = f{1};
    ci.(k) = find(strcmp(hdr, COL.(k)), 1);
    assert(~isempty(ci.(k)), 'AQ column not found: %s', COL.(k));
end
aqId     = toNum(T{:, ci.id});
aqAge    = toNum(T{:, ci.age});
aqGender = toGender(T{:, ci.gender});
AQ = struct();
for k = aqNames, AQ.(k{1}) = toNum(T{:, ci.(k{1})}); end

% --- Schaefer-7 network membership (shared) ---
NETWORKS = {'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
A = readtable(atlasCsv, 'PreserveVariableNames', true);
roiNames = A.("ROI Name");
netIdx = zeros(100,1);
for n = 1:100
    tok = strsplit(roiNames{n}, '_');
    netIdx(n) = find(strcmp(NETWORKS, tok{3}));
end

for ci1 = 1:numel(conds)
    for ci2 = 1:numel(conns)
        cond = conds{ci1}; conn = conns{ci2};
        graphDir = fullfile(cfg.graphRoot, cond, ...
            'EEG_features', 'graph_measures');
        pat = sprintf('*_graph_%s_*.mat', conn);
        assert(isfolder(graphDir) && ~isempty(dir(fullfile(graphDir, pat))), ...
            'Missing regenerated graph files: %s', fullfile(graphDir, pat));

        % subjects with all 4 band files for this conn
        files = dir(fullfile(graphDir, pat));
        count = containers.Map();
        for i = 1:numel(files)
            tok = regexp(files(i).name, '^(.+)_graph_', 'tokens');
            if isempty(tok), continue; end
            sid = tok{1}{1};
            count(sid) = getOr(count, sid, 0) + 1;
        end
        allSids = keys(count);
        sids = {};
        for i = 1:numel(allSids)
            if count(allSids{i}) == numel(bands)
                sids{end+1} = allSids{i}; %#ok<AGROW>
            end
        end
        sids = sort(sids);
        N = numel(sids);

        subNum = cellfun(@(s) str2double(regexp(s, '\d+', 'match', 'once')), sids);
        [tf, iaq] = ismember(subNum, aqId);
        assert(all(tf), 'subjects missing in AQ sheet: %s', strjoin(sids(~tf), ', '));

        D = struct();
        D.subj_ids = sids;
        G = struct(); Nd = struct();
        for q = globalQ
            for b = 1:numel(bands)
                G.(sprintf('%s_%s_%s', conn, bands{b}, q{1})) = NaN(N,1);
            end
        end
        for q = nodalQ
            for b = 1:numel(bands)
                Nd.(q{1}).(sprintf('%s_%s', conn, bands{b})) = NaN(N,100);
            end
        end
        for k = 1:N
            for b = 1:numel(bands)
                Sg = load(fullfile(graphDir, sprintf('%s_graph_%s_%s.mat', sids{k}, conn, bands{b})));
                gm = Sg.graph_measures;
                for q = globalQ
                    G.(sprintf('%s_%s_%s', conn, bands{b}, q{1}))(k) = gm.(q{1});
                end
                for q = nodalQ
                    Nd.(q{1}).(sprintf('%s_%s', conn, bands{b}))(k,:) = gm.(q{1});
                end
            end
        end
        D.global = G; D.nodal = Nd;
        D.age    = aqAge(iaq);
        D.gender = aqGender(iaq);
        for k = aqNames, D.(k{1}) = AQ.(k{1})(iaq); end

        D.netIdx = netIdx; D.NETWORKS = NETWORKS;
        D.aqNames = aqNames; D.bands = bands;
        D.globalQ = globalQ; D.nodalQ = nodalQ;
        D.connMeas = conn; D.cond = cond;

        outFile = fullfile(outDir, sprintf('design_data_%s_%s.mat', cond, conn));
        save(outFile, 'D');
        fprintf('Saved %s (N=%d, %s, %s)\n', outFile, N, cond, conn);
    end
end

%% ===== helpers =====
function v = toNum(x)
if isnumeric(x), v = double(x); return; end
v = str2double(string(x));
end

function g = toGender(x)
% Spreadsheet encodes 男/女 (or 1/2). Output: 1 = male, 2 = female.
if isnumeric(x)
    g = double(x); return;
end
s = string(x);
g = nan(size(s));
g(s == "男" | s == "1") = 1;
g(s == "女" | s == "2") = 2;
end

function v = getOr(map, key, default)
if isKey(map, key), v = map(key); else, v = default; end
end
