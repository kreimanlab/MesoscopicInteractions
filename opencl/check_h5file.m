close all;
clear;

[~,hname] = system('hostname');
if strcmp(strip(hname),'ubuntu_1604')
    resultsDir = '/nas_share/RawData/data/results';
else
    resultsDir = '/mnt/cuenap2/data/results';
end

o = getOutputFilenames(resultsDir, true);

for i = 1:o.n_f
    try
        graph_f = [resultsDir,'/',o.SubjectsAllGraph{i}];
        %fprintf('%s\n',graph_f)
        h5read(graph_f,'/R');
    catch
        sid = strsplit(o.SubjectsAllGraph{i},'_');
        sid = sid{1};
        metric = strsplit(graph_f,'.h5');
        metric = strsplit(metric{1},'-');
        metric = metric{end};
        dev = mod(i,2)+1;
        fprintf('\tpython3 graph.py %s %s %i\n',sid,metric,dev)
    end

    try
        perm_f = [resultsDir,'/',o.SubjectsAllPerm{i}];
        %fprintf('%s\n',perm_f)
        h5read(perm_f,'/R');
    catch
        sid = strsplit(o.SubjectsAllPerm{i},'_');
        sid = sid{1};
        metric = strsplit(perm_f,'.h5');
        metric = strsplit(metric{1},'-');
        metric = metric{end};
        dev = mod(i,2)+1;
        fprintf('\tpython3 perm.py %s %s 10000 %i\n',sid,metric,dev)
    end

    try
        dist_f = [resultsDir,'/',o.SubjectsAllDists{i}];
        %fprintf('%s\n',dist_f)
        load(dist_f);
    catch
        fprintf('dist: %s\n',dist_f)
    end
end
