close all; clear;

%Subs = {'m00006','m00019','m00023','m00024','m00026','m00030','m00037','m00038','m00043','m00060','m00068','m00083'};
%Subs = {'m00006'};

% Use all folders for SubjectsAll
[~,host] = system('hostname');
if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix    
    resultsDir = '/mnt/cuenap2/data/results';
    h5Dir = '/mnt/cuenap2/scripts/synth/out';
end

if contains(host,'ubuntu_1604')
    resultsDir = '/nas_share/RawData/scripts/opencl/results_res5hz';
    h5Dir = '/nas_share/RawData/data/h5_notch20';
end

extra_check = false;

o = getOutputFilenames(resultsDir,false);
n_f = o.n_f;
SubjectsAllPerm = o.SubjectsAllPerm;
SubjectsAllGraph = o.SubjectsAllGraph;
SubjectsAllDists = o.SubjectsAllDists;
% fprintf('[*] Looking for results in: %s\n',resultsDir);
% fprintf('[*] Looking for h5eeg in: %s\n',h5Dir);
% D = dir(resultsDir);
% Disdir = [D.isdir];
% Dname = {D.name};
% Dname = Dname(~Disdir);
% SubjectsAllPerm = Dname(startsWith(Dname,'m') & endsWith(Dname,'.h5') & contains(Dname,'perm'));
% SubjectsAllGraph = Dname(startsWith(Dname,'m') & endsWith(Dname,'.h5') & contains(Dname,'graph'));
% SubjectsAllDists = Dname(startsWith(Dname,'m') & endsWith(Dname,'.mat') & contains(Dname,'dists'));
% n_f = length(SubjectsAllPerm);

if (extra_check)
% Trim graph files with no perm partner
SubjectsAllGraphN = {};
count = 1;
for i = 1:length(SubjectsAllGraph)
    sag = SubjectsAllGraph{i};
    sag = strsplit(sag,'.h5');
    sad = replace(sag{1},'graph','dists');
    sag = replace(sag{1},'graph','perm');
    if ((sum(startsWith(SubjectsAllPerm,sag)) > 0) && (sum(startsWith(SubjectsAllDists,sad)) == 0))
        SubjectsAllGraphN{count} = SubjectsAllGraph{i};
        count = count + 1;
    end
end
SubjectsAllGraph = SubjectsAllGraphN;

% Trim perm files with no graph partner
SubjectsAllPermN = {};
count = 1;
for i = 1:length(SubjectsAllPerm)
    sap = SubjectsAllPerm{i};
    sap = strsplit(sap,'.h5');
    sad = replace(sap{1},'perm','dists');
    sap = replace(sap{1},'perm','graph');
    sap = strsplit(sap,'-');
    sap = strjoin(sap(1:(end-1)),'-');
    if ((sum(startsWith(SubjectsAllGraph,sap)) > 0) && (sum(startsWith(SubjectsAllDists,sad)) == 0))
        SubjectsAllPermN{count} = SubjectsAllPerm{i};
        count = count + 1;
    end
end
SubjectsAllPerm = SubjectsAllPermN;

end

n_f = length(SubjectsAllPerm);
for i = 1:n_f
    fprintf('%s : %s\n',SubjectsAllPerm{i},SubjectsAllGraph{i})
end
%return

parfor iSub = 1:n_f
    ss = strsplit(SubjectsAllPerm{iSub},'_');
    ss2 = strsplit(ss{2},'-');
    sid = ss{1};
    metric = ss2{2};
    suffix = strjoin({ss2{2},ss2{end}(1:(end-3))},'-');
    permf = sprintf('%s/%s',resultsDir,SubjectsAllPerm{iSub});
    graphf = sprintf('%s/%s_graph-%s.h5',resultsDir,sid,metric);
    fprintf('[#] Processing: %s\n',permf)
    
    % Check if files exist
    h5fname = sprintf('%s/%s.h5',h5Dir,sid);
    if (~exist(graphf,'file'))
        fprintf('[!] skip: graph file could not be read: %s\n',graphf);
    elseif (~exist(permf,'file'))
        fprintf('[!] skip: perm file could not be read: %s\n',permf);
    elseif (~exist(h5fname,'file'))
        fprintf('[!] skip: h5 file could not be read: %s\n',h5fname);
    else
        % Read h5eeg
        try
            h5fname = sprintf('%s/%s.h5',h5Dir,sid);
            % h5fname = sprintf('/Volumes/RawData/scripts/synth/out/%s.h5',sid);
            fs = h5readatt(h5fname,'/h5eeg/eeg','rate');
            n_chan = double(h5readatt(h5fname,'/h5eeg/eeg','n_chan'));
            n_samples = double(h5readatt(h5fname,'/h5eeg/eeg','n_samples'));
            chan_labels = h5readatt(h5fname,'/h5eeg/eeg','labels');
            h5_n_samples = double(h5readatt(h5fname,'/h5eeg','n_samples'));
            bip = h5readatt(h5fname,'/h5eeg/eeg','bip');
            arts = h5read(h5fname,'/h5eeg/artifacts');
            width = double(h5readatt(h5fname,'/h5eeg/artifacts','width'));

            % Read results
            R = double(h5read(permf,'/R'));
            %if ( startsWith(metric,'pc') || startsWith(metric,'sc') )
            %    R = sqrt(R);
            %end
            starts = h5read(permf,'/starts');
            chan1 = h5read(permf,'/chan1');
            chan2 = h5read(permf,'/chan2');
            r_rows = h5read(permf,'/r_rows');
            r_cols = h5read(permf,'/r_cols');
            w = double(h5read(permf,'/w'));
            alpha = double(h5read(permf,'/alpha'));

            DIST = 'tLocationScale';
            [n_comb,n_perm] = size(R);
            Dists = cell(1,n_comb);
            for iDist = 1:n_comb
                tic;
                Rclean = R(iDist,:)';
                cleanIdx = ~isinf(Rclean);
                Rclean = Rclean(cleanIdx);
                Dists{iDist} = fitdist(Rclean,DIST);
                fprintf('[%i/%i] %s-%s (%i/%i) n:%i, ETA:%.2f hrs\n',...
                    iSub,n_f,sid,metric,iDist,n_comb,sum(cleanIdx),toc*((n_f-iSub)*n_comb+(n_comb-iDist))/3600)
            end
            %save(sprintf('%s/%s_%s_dists',resultsDir,sid,suffix),'Dists')
            saveFitPerm(sprintf('%s/%s_dists-%s',resultsDir,sid,suffix),Dists);
        catch e
            %throw(e);
            z23 = 1 + 1;
        end
    end
    
end
fprintf('[*] Done.\n')


function [] = saveFitPerm( fn, d )
    save(fn,'d','-v7.3');
end
