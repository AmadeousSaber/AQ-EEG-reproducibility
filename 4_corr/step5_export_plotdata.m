% step5_export_plotdata.m
% -------------------------------------------------------------------------
% Export partial-residual scatter data (brain & AQ residualised on age+gender)
% for every q4<0.05 hit across ALL combos (open/closed x dwpli/aec) and all
% three spatial levels. Output: plot_data_all.csv (tidy; one row per subject
% per hit; columns include cond + conn + level).
% -------------------------------------------------------------------------

cfg = pipeline_config; outDir = cfg.corrDir;
alphaFDR = 0.05;

conds = {'open','closed'};
conns = {'dwpli','aec'};

rows = {}; hid = 0;

for ci1 = 1:numel(conds)
    for ci2 = 1:numel(conns)
        cond = conds{ci1}; conn = conns{ci2};
        S = load(fullfile(outDir, sprintf('design_data_%s_%s.mat', cond, conn))); D = S.D;
        Z = [D.age, D.gender];
        netIdx = D.netIdx; NETWORKS = D.NETWORKS;

        % ---- global hits ----
        Tg = readtable(fullfile(outDir, sprintf('results_global_%s_%s.csv', cond, conn)));
        sg = Tg(Tg.q4 < alphaFDR, :);
        for i = 1:height(sg)
            hid = hid + 1;
            fn = toStr(sg.metric(i));
            x = D.global.(fn);
            y = D.(toStr(sg.aqvar(i)));
            label = sprintf('%s x %s (global)', fn, toStr(sg.aqvar(i)));
            rows = appendRows(rows, hid, cond, conn, 'global', label, ...
                toStr(sg.aqvar(i)), sg.r(i), sg.p(i), sg.q4(i), x, y, Z, D.gender);
        end

        % ---- network hits ----
        Tn = readtable(fullfile(outDir, sprintf('results_network_%s_%s.csv', cond, conn)));
        sn = Tn(Tn.q4 < alphaFDR, :);
        for i = 1:height(sn)
            hid = hid + 1;
            fn = sprintf('%s_%s', conn, toStr(sn.band(i)));
            M  = D.nodal.(toStr(sn.quantity(i))).(fn);
            v  = mean(M(:, netIdx == find(strcmp(NETWORKS, toStr(sn.network(i))))), 2, 'omitnan');
            y = D.(toStr(sn.aqvar(i)));
            label = sprintf('%s %s %s x %s (network)', fn, toStr(sn.quantity(i)), toStr(sn.network(i)), toStr(sn.aqvar(i)));
            rows = appendRows(rows, hid, cond, conn, 'network', label, ...
                toStr(sn.aqvar(i)), sn.r(i), sn.p(i), sn.q4(i), v, y, Z, D.gender);
        end

        % ---- node hits ----
        Td = readtable(fullfile(outDir, sprintf('results_node_%s_%s.csv', cond, conn)));
        sd = Td(Td.q4 < alphaFDR, :);
        for i = 1:height(sd)
            hid = hid + 1;
            fn = sprintf('%s_%s', conn, toStr(sd.band(i)));
            M  = D.nodal.(toStr(sd.quantity(i))).(fn);
            v  = M(:, sd.node(i));
            y = D.(toStr(sd.aqvar(i)));
            label = sprintf('%s %s node%d x %s (node)', fn, toStr(sd.quantity(i)), sd.node(i), toStr(sd.aqvar(i)));
            rows = appendRows(rows, hid, cond, conn, 'node', label, ...
                toStr(sd.aqvar(i)), sd.r(i), sd.p(i), sd.q4(i), v, y, Z, D.gender);
        end
    end
end

T = cell2table(rows, 'VariableNames', {'hit_id','cond','conn','level','label','aqvar','r','p','q4','resid_brain','resid_aq','gender'});
writetable(T, fullfile(outDir, 'plot_data_all.csv'));
fprintf('Exported %d hits -> plot_data_all.csv\n', hid);

%% ===== helpers =====
function s = toStr(v)
% robust scalar text extraction from table cell (string or cellstr)
if iscell(v), v = v{1}; end
s = char(string(v));
end

function rows = appendRows(rows, hid, cond, conn, level, label, aqvar, r, p, q4, x, y, Z, gender)
ok = ~isnan(x) & ~isnan(y) & all(~isnan(Z),2);
x = x(ok); y = y(ok); Zk = Z(ok,:); g = gender(ok);
X = [ones(numel(x),1), Zk];
rx = x - X*(X\x);
ry = y - X*(X\y);
for k = 1:numel(x)
    rows(end+1,:) = {hid, cond, conn, level, label, aqvar, r, p, q4, rx(k), ry(k), g(k)}; %#ok<AGROW>
end
end
