close all; clear;

P_VALUE = 0.05/60;
openclDir = '.';
N_plot = 5;
DIST_THRESH_MM = 15;

% Use all folders for SubjectsAll
if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix    
    resultsDir = '/mnt/cuenap2/data/results';
    h5Dir = '/mnt/cuenap2/scripts/synth/out';
end

[~,host] = system('hostname');
if contains(host,'ubuntu_1604')
    resultsDir = '/nas_share/RawData/data/results';
    h5Dir = '/nas_share/RawData/scripts/synth/out';
    %coregDir = '/mnt/cuenap_ssd/coregistration';
    %outDir = '/nas_share/RawData/data/results/figures';
end

o = getOutputFilenames(resultsDir,true);
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
% n_f = length(SubjectsAllPerm);

for iSub = 1:n_f
    
    ss = strsplit(SubjectsAllPerm{iSub},'_');
    ss2 = strsplit(ss{2},'-');
    sid = ss{1};
    metric = ss2{2};
    suffix = strjoin({ss2{2},ss2{end}(1:(end-3))},'-');
    permf = sprintf('%s/%s',resultsDir,SubjectsAllPerm{iSub});
    graphf = sprintf('%s/%s_graph-%s.h5',resultsDir,sid,metric);
    distsf = sprintf('%s/%s_dists-%s.mat',resultsDir,sid,suffix);
    fprintf('[#] Processing: %s\n',permf)
    
    % make sure output directory exists
    system(sprintf('mkdir %s/figures',openclDir));
    system(sprintf('mkdir %s/figures/%s_%s',openclDir,sid,suffix));
    
    % Read h5eeg
    h5fname = sprintf('%s/%s.h5',h5Dir,sid);
    fs = h5readatt(h5fname,'/h5eeg/eeg','rate');
    n_chan = double(h5readatt(h5fname,'/h5eeg/eeg','n_chan'));
    n_samples = double(h5readatt(h5fname,'/h5eeg/eeg','n_samples'));
    chan_labels = h5readatt(h5fname,'/h5eeg/eeg','labels');
    h5_n_samples = double(h5readatt(h5fname,'/h5eeg','n_samples'));
    bip = h5readatt(h5fname,'/h5eeg/eeg','bip');
    arts = h5read(h5fname,'/h5eeg/artifacts');
    width = double(h5readatt(h5fname,'/h5eeg/artifacts','width'));
    bipChan = h5readatt(h5fname,'/h5eeg/eeg','bip');
    [n_bchan,~] = size(bipChan);
    
    % Read graph results
    R = h5read(graphf,'/R');
    if ( startsWith(metric,'pc') || startsWith(metric,'sc') )
        R = sqrt(R);
    end
    chan1 = h5read(graphf,'/chan1');
    chan2 = h5read(graphf,'/chan2');
    chan1 = double(chan1)+1;
    chan2 = double(chan2)+1;
    r_cols = h5read(graphf,'/r_cols');
    r_rows = h5read(graphf,'/r_rows');
    w = double(h5read(graphf,'/w'));
    %return
    
    % Load null distributions
    % ==== PLACEHOLDER ====
    fprintf('[*] Loading distributions..\n')
    df = load(distsf,'d');
    Dists = df.d;
%     load('dist.mat','dist')
%     Dists = cell(1,length(chan1));
%     for i = 1:length(chan1)
%         Dists{i} = dist;
%     end

    % Make distance matrix
    Dmat = zeros(n_bchan,n_bchan);
    for i2 = 1:length(Dmat)
        for j2 = 1:length(Dmat)
            a = reshape(bipChan(i2,4:end),2,[]);
            b = reshape(bipChan(j2,4:end),2,[]);
            dist2 = zeros(2,2);
            for a1 = 1:2
                for b1 = 1:2
                    dist2(a1,b1) = sqrt(sum((a(a1,:) - b(b1,:)).^2));
                end
            end
            Dmat(i2,j2) = min(dist2(:));
        end
    end
    
    % Threshold R by distance
    [~,R_n] = size(R);
    co = 0;
    for i3 = 1:(n_bchan-1)
        for j3 = (i3+1):n_bchan
            co = co + 1;
            distance = Dmat(i3,j3);
            if (distance < DIST_THRESH_MM)
                R(co,:) = nan(size(R(co,:)));
            end
        end
    end
    
    % Randomly choose pair to plot
    rIdx = false(size(chan1));
    rIdx(1) = true;
    rIdx = rIdx(randperm(length(rIdx)));
    rIdxs = cell(1,1);
    rIdxs{1} = rIdx;
    
    % Plot top,mid,low
    Rmag = nanmean(abs(R),2);
    [~,sortIdx] = sort(Rmag);
    sortIdx = sortIdx(1:(find(isnan(Rmag(sortIdx)),1)-1));
    %sortIdx = sortIdx(~isnan(Rmag(sortIdx)));
    counter = 1;
    rIdxs = {};
    for i = 1:N_plot
        % Lowest correlations
        rIdx = false(size(chan1));
        rIdx(sortIdx(i)) = true;
        rIdxs{counter} = rIdx;
        counter = counter + 1;
        
        % Medium correlations
        rIdx = false(size(chan1));
        rIdx(sortIdx(round(length(Rmag)/2)+i)) = true;
        rIdxs{counter} = rIdx;
        counter = counter + 1;
        
        % Highest correlations
        rIdx = false(size(chan1));
        rIdx(sortIdx(end-i-1)) = true;
        rIdxs{counter} = rIdx;
        counter = counter + 1;
    end
    %return
    
    for plotIdx = 1:length(rIdxs)
        rIdx = rIdxs{plotIdx};

        % Plot correlation vs time
        h = figure('Position',[0 0 2000 800]);
        YMIN = -1;
        YMAX = 1;
        Adj = R(rIdx,:);
        T = ((1:length(Adj))-1)*(w)/(3600*24);

        % Calculate artifact locations
        arts_index = arts(2,(arts(1,:)~=0))/(3600*24*fs);
        offset = 0;
        AdjRej = nan(size(Adj));
        AdjPlot = Adj;
        for iArt = 1:length(arts_index)
            rejectIdx = find(arts_index(iArt) < T,1,'first');
            AdjRej(rejectIdx) = Adj(rejectIdx);
            AdjPlot(rejectIdx) = NaN;
        end

        % Significance
        %P_VALUE = P_VALUE_FINAL;%/length(chan1);
        p = cdf(Dists{rIdx},AdjPlot);
        sigIdx = ((p<(P_VALUE/2))|(p>(1-P_VALUE/2)));

        % Plot file markers
        subplot(2,1,1)
        xoffset = 0.02;
        for iSamp = 1:(length(h5_n_samples)-1)
            xpos = (sum(h5_n_samples(1:iSamp))/(3600*24*fs));
            text(xpos - xoffset,YMIN,'//'); hold on;
            plot(xpos*ones(1,2),[YMIN YMAX],'--','Color',0.5*ones(1,3)); hold on;
        end
        plot(T,AdjRej,'.','color',0.5*ones(1,3),'MarkerSize',4); hold on;
        plot(T,AdjPlot,'black.'); hold on;
        plot(T(sigIdx),AdjPlot(sigIdx),'red.'); hold on;
        axis([min(T),max(T),YMIN,YMAX])

        % Process title text
        Dist = Dmat(chan1(rIdx),chan2(rIdx));
%         Dist = zeros(1,4);
%         Dist(1) = sqrt(sum(( bip(chan1(rIdx),4:6) - bip(chan2(rIdx),4:6) ).^2));
%         Dist(2) = sqrt(sum(( bip(chan1(rIdx),4:6) - bip(chan2(rIdx),7:9) ).^2));
%         Dist(3) = sqrt(sum(( bip(chan1(rIdx),7:9) - bip(chan2(rIdx),4:6) ).^2));
%         Dist(4) = sqrt(sum(( bip(chan1(rIdx),7:9) - bip(chan2(rIdx),7:9) ).^2));

        e1Name = [chan_labels{bip(chan1(rIdx),1)},'-',chan_labels{bip(chan1(rIdx),2)}];
        e2Name = [chan_labels{bip(chan2(rIdx),1)},'-',chan_labels{bip(chan2(rIdx),2)}];
        title(sprintf('%s: %s - %s (%.0f mm)',sid,e1Name,e2Name,min(Dist)))
        xlabel('DAYS')
        ylabel('SPEARMAN BIP CORR COEF')

        % Plot sampled drom null
        subplot(2,1,2)
        
        for iSamp = 1:(length(h5_n_samples)-1)
            xpos = (sum(h5_n_samples(1:iSamp))/(3600*24*fs));
            plot(xpos*ones(1,2),[YMIN YMAX],'--','Color',0.5*ones(1,3)); hold on;
        end
        
        AdjNeg = random(Dists{rIdx},1,length(AdjPlot));
        pNeg = cdf(Dists{rIdx},AdjNeg);
        sigIdxNeg = ((pNeg<(P_VALUE/2))|(pNeg>(1-P_VALUE/2)));
        plot(T,AdjNeg,'black.'); hold on;
        plot(T(sigIdxNeg),AdjNeg(sigIdxNeg),'red.'); hold on;
        axis([min(T),max(T),YMIN,YMAX])
        xlabel('DAYS')
        ylabel('SPEARMAN BIP CORR COEF')
        title(sprintf('%s: %s - %s (%.0f mm) CHANCE',sid,e1Name,e2Name,min(Dist)))

        % Save figure
        print(h,sprintf('%s/figures/%s_%s/%s_%s_%s_%i_%i.eps',...
            openclDir,sid,suffix,sid,e1Name,e2Name,round(min(Dist)),round(1000*nanmean(AdjPlot))),'-depsc')
        close(h);
        
        
        % Plot timeseries
        
    end
    fprintf('[!] Finished patient: %s\n',sid);
end
fprintf('[!] Done.\n')
