close all; clear;

%Subs = {'m00006','m00019','m00023','m00024','m00026','m00030','m00037','m00038','m00043','m00060','m00068','m00083'};
%Subs = {'m00006','m00019','m00023','m00024','m00026','m00030','m00037'};
Subs = {'m00006','m00019','m00023'};
P_VALUE_FINAL = 0.05;
DIST_THRESH_MM = 15;
h5eegDir = '.';
openclDir = '.';

mmpLabelsAll = {};
for iSub = 1:length(Subs)
    h5fname = sprintf('%s/%s.h5',h5eegDir,Subs{iSub});
    labels_mmp = h5readatt(h5fname,'/h5eeg/eeg','labels_mmp');
    mmpLabelsT = unique(labels_mmp);
    mmpLabels = cell(size(mmpLabelsT));
    for i = 1:length(mmpLabelsT)
        mmplabel = strsplit(mmpLabelsT{i},'_');
        if (length(mmplabel) == 1)
            mmpLabels{i} = mmplabel{1};
        else
            mmpLabels{i} = mmplabel{2};
        end
    end
    mmpLabelsAll = unique({mmpLabelsAll{:},mmpLabels{:}});
end

n_atlas = length(mmpLabelsAll);
Adj_cov = zeros(n_atlas,n_atlas);
Adj_sub = cell(n_atlas,n_atlas);
Adj_score = cell(n_atlas,n_atlas);
Adj_scoreAvg = zeros(n_atlas,n_atlas);
Adj_scoreNeg = cell(n_atlas,n_atlas);
Adj_scoreNegAvg = zeros(n_atlas,n_atlas);
for iSub = 1:length(Subs)
    
    system(sprintf('mkdir %s/figures',openclDir));
    
    % Read h5eeg
    h5fname = sprintf('%s/%s.h5',h5eegDir,Subs{iSub});
    fs = h5readatt(h5fname,'/h5eeg/eeg','rate');
    n_chan = double(h5readatt(h5fname,'/h5eeg/eeg','n_chan'));
    n_samples = double(h5readatt(h5fname,'/h5eeg/eeg','n_samples'));
    chan_labels = h5readatt(h5fname,'/h5eeg/eeg','labels');
    h5_n_samples = double(h5readatt(h5fname,'/h5eeg','n_samples'));
    bip = h5readatt(h5fname,'/h5eeg/eeg','bip');
    labels_mmp = h5readatt(h5fname,'/h5eeg/eeg','labels_mmp');
    labels_m132 = h5readatt(h5fname,'/h5eeg/eeg','labels_m132');
    arts = h5read(h5fname,'/h5eeg/artifacts');
    width = double(h5readatt(h5fname,'/h5eeg/artifacts','width'));
    
    % Read graph results
    R = h5read(sprintf('%s/out/%s_graph.h5',openclDir,Subs{iSub}),'/R');
    chan1 = h5read(sprintf('%s/out/%s_graph.h5',openclDir,Subs{iSub}),'/chan1');
    chan2 = h5read(sprintf('%s/out/%s_graph.h5',openclDir,Subs{iSub}),'/chan2');
    chan1 = double(chan1)+1;
    chan2 = double(chan2)+1;
    r_cols = h5read(sprintf('%s/out/%s_graph.h5',openclDir,Subs{iSub}),'/r_cols');
    r_rows = h5read(sprintf('%s/out/%s_graph.h5',openclDir,Subs{iSub}),'/r_rows');
    w = double(h5read(sprintf('%s/out/%s_graph.h5',openclDir,Subs{iSub}),'/w'));
        
    % Load null distributions
    % ==== PLACEHOLDER ====
    fprintf('[*] Loading distributions..\n')
    load(sprintf('%s/out/%s_dists.mat',openclDir,Subs{iSub}),'Dists');
    AtlasLabels = labels_mmp;
    
    % Calculate artifact locations
    
    P_VALUE = P_VALUE_FINAL;
    Rmag = mean(abs(R),2);
    CT = zeros(size(Rmag));
    RmagNeg = zeros(size(Rmag));
    CTneg = zeros(size(Rmag));
    for dIdx = 1:length(Dists)
        Adj = R(dIdx,:);
        T = ((1:length(Adj))-1)*(w)/(3600*24);
        arts_index = arts(2,(arts(1,:)~=0))/(3600*24*fs);
        offset = 0;
        AdjRej = nan(size(Adj));
        AdjPlot = Adj;
        AdjNeg = random(Dists{dIdx},1,length(AdjPlot));
        RmagNeg(dIdx) = mean(abs(AdjNeg));
%         for iArt = 1:length(arts_index)
%             rejectIdx = find(arts_index(iArt) < T,1,'first');
%             AdjRej(rejectIdx) = Adj(rejectIdx);
%             AdjPlot(rejectIdx) = NaN;
%         end
        p = cdf(Dists{dIdx},AdjPlot);
        sigIdx = ((p<(P_VALUE/2))|(p>(1-P_VALUE/2)));
        CT(dIdx) = sum(sigIdx)/length(sigIdx);
        pNeg = cdf(Dists{dIdx},AdjNeg);
        sigIdxNeg = ((pNeg<(P_VALUE/2))|(pNeg>(1-P_VALUE/2)));
        CTneg(dIdx) = sum(sigIdxNeg)/length(sigIdxNeg);
    end

    
    n_bchan = double(r_cols);
    bipA = zeros(n_bchan,2);
    for j = 1:n_bchan
        e = bip(j,1:2);
        ss1 = strsplit(AtlasLabels{e(1)},'_');
        ss2 = strsplit(AtlasLabels{e(2)},'_');

        % unknown electrodes
        if (length(ss1) == 1)
            atIdx1 = 0;
        else
            atIdx1 = find(strcmp(mmpLabelsAll,ss1{2}));
        end
        if (length(ss2) == 1)
            atIdx2 = 0;
        else
            atIdx2 = find(strcmp(mmpLabelsAll,ss2{2}));
        end
        bipA(j,:) = [atIdx1,atIdx2];
    end
    

    % Map to coverage matrix
    for m = 1:n_bchan
        for n = 1:n_bchan
            for a1 = 1:2
                for a2 = 1:2
                    % Distance thresholding
                    cBip = [bip(m,4:end),bip(n,4:end)];
                    cBip = reshape(cBip,4,3);
                    Adist = zeros(2,2);
                    for d1 = 1:2
                        for d2 = 1:2
                            Adist(d1,d2) = sqrt(sum((cBip(d1,:) - cBip(d2+2,:)).^2));
                        end
                    end
%                     Adist = zeros(4,4);
%                     for d1 = 1:4
%                         for d2 = 1:4
%                             Adist(d1,d2) = sqrt(sum((cBip(d1,:) - cBip(d2,:)).^2));
%                         end
%                     end
                    skip = (min(Adist(~(Adist==0))) > DIST_THRESH_MM);
                    if (~skip)
                        try
                            
                            m2 = m;
                            n2 = n;
                            if (m > n)
                                mtemp = m;
                                m2 = n;
                                n2 = mtemp;
                            end
                            combIdx = (chan1==m2) & (chan2==n2);
                            Adj_score{bipA(m,a1),bipA(n,a2)} = [Adj_score{bipA(m,a1),bipA(n,a2)},CT(combIdx)+Rmag(combIdx)];
                            Adj_scoreAvg(bipA(m,a1),bipA(n,a2)) = (Adj_scoreAvg(bipA(m,a1),bipA(n,a2)) + CT(combIdx)+Rmag(combIdx));
                            Adj_scoreNeg{bipA(m,a1),bipA(n,a2)} = [Adj_scoreNeg{bipA(m,a1),bipA(n,a2)},CTneg(combIdx)+RmagNeg(combIdx)];
                            Adj_scoreNegAvg(bipA(m,a1),bipA(n,a2)) = (Adj_scoreNegAvg(bipA(m,a1),bipA(n,a2)) + CTneg(combIdx)+RmagNeg(combIdx));
                            
                            Adj_cov(bipA(m,a1),bipA(n,a2)) = Adj_cov(bipA(m,a1),bipA(n,a2)) + 1;
                            Adj_sub{bipA(m,a1),bipA(n,a2)} = [Adj_sub{bipA(m,a1),bipA(n,a2)},iSub];
                            
                        catch
                            try
                                m2 = m;
                                n2 = n;
                                if (m < n)
                                    mtemp = m;
                                    m2 = n;
                                    n2 = mtemp;
                                end
                                combIdx = (chan1==m2) & (chan2==n2);
                                Adj_score{bipA(m,a1),bipA(n,a2)} = [Adj_score{bipA(m,a1),bipA(n,a2)},CT(combIdx)+Rmag(combIdx)];
                                Adj_scoreAvg(bipA(m,a1),bipA(n,a2)) = (Adj_scoreAvg(bipA(m,a1),bipA(n,a2)) + CT(combIdx)+Rmag(combIdx));
                                Adj_scoreNeg{bipA(m,a1),bipA(n,a2)} = [Adj_scoreNeg{bipA(m,a1),bipA(n,a2)},CTneg(combIdx)+RmagNeg(combIdx)];
                                Adj_scoreNegAvg(bipA(m,a1),bipA(n,a2)) = (Adj_scoreNegAvg(bipA(m,a1),bipA(n,a2)) + CTneg(combIdx)+RmagNeg(combIdx));
                            catch
                                % placeholder
                            end
                            %fprintf('Skip. \n')
                        end
                    end
                end
            end
        end
    end

    fprintf('[!] Finished patient: %s\n',Subs{iSub});
end

% --- CONFIG ---
fsize = 5;
fsize_title = 8;

h = figure;
pctCov = 100*(sum(sum(Adj_cov ~= 0))/(180*180));
rIdx = ~(sum((Adj_cov == 0)) == n_atlas);
imagesc(Adj_cov(rIdx,rIdx))
cc = parula(10*64);
cc = [[1 1 1];cc];
colormap(cc);
set(gca,'xtick',[1:length(mmpLabelsAll(rIdx))],'Ticklength',[2e-3 1e-5],...
    'xticklabel',mmpLabelsAll(rIdx),'fontsize',fsize,'TickDir','out');
xtickangle(90);
set(gca,'ytick',[1:length(mmpLabelsAll(rIdx))],'Ticklength',[2e-3 1e-5],...
    'yticklabel',mmpLabelsAll(rIdx),'fontsize',fsize,'TickDir','out');
title(sprintf('HCP-MMP parcelation coverage (%.2f %%, n = %i)',pctCov,length(Subs)),'fontsize',8)
c = colorbar; c.FontSize = fsize_title;
print(h,sprintf('figures/mmp_Adj_bipcov-%i',length(Subs)),'-depsc')

% Number of subjects per ROI pair    
Adj_nsub = zeros(size(Adj_sub));
for m = 1:length(Adj_sub)
    for n = 1:length(Adj_sub)
        Adj_nsub(m,n) = length(unique(Adj_sub{m,n}));
    end
end

h = figure;
rIdx = ~(sum((Adj_nsub == 0)) == n_atlas);
imagesc(Adj_nsub(rIdx,rIdx))
cc = parula(10*64);
cc = [[1 1 1];cc];
colormap(cc);
%fsize = 3;
set(gca,'xtick',[1:length(mmpLabelsAll(rIdx))],'Ticklength',[2e-3 1e-5],...
    'xticklabel',mmpLabelsAll(rIdx),'fontsize',fsize,'TickDir','out');
xtickangle(90);
set(gca,'ytick',[1:length(mmpLabelsAll(rIdx))],'Ticklength',[2e-3 1e-5],...
    'yticklabel',mmpLabelsAll(rIdx),'fontsize',fsize,'TickDir','out');
title(sprintf('HCP-MMP parcelation number of patients (%.2f %% coverage, n = %i)',pctCov,length(Subs)),'fontsize',8)
c = colorbar; c.FontSize = fsize_title;
print(h,sprintf('figures/mmp_Adj_bipnsub-%i',length(Subs)),'-depsc')

h = figure;
AdjPlot = (Adj_scoreAvg./Adj_cov);
imagesc(AdjPlot(rIdx,rIdx));
colormap(cc);
%caxis([min(AdjPlot(:)) max(AdjPlot(:))]);
caxis([0 2])
colorbar;
set(gca,'xtick',[1:length(mmpLabelsAll(rIdx))],'Ticklength',[2e-3 1e-5],...
    'xticklabel',mmpLabelsAll(rIdx),'fontsize',fsize,'TickDir','out');
xtickangle(90);
set(gca,'ytick',[1:length(mmpLabelsAll(rIdx))],'Ticklength',[2e-3 1e-5],...
    'yticklabel',mmpLabelsAll(rIdx),'fontsize',fsize,'TickDir','out');
title('Composite Interactivity')
print(h,sprintf('figures/mmp_Adj_scoreAvgnsub-%i',length(Subs)),'-depsc')

h = figure;
AdjPlotNeg = (Adj_scoreNegAvg./Adj_cov);
imagesc(AdjPlotNeg(rIdx,rIdx));
%cc = parula(10*64);
colormap(cc);
%caxis([min(AdjPlot(:)) max(AdjPlot(:))]);
caxis([0 2])
colorbar;
set(gca,'xtick',[1:length(mmpLabelsAll(rIdx))],'Ticklength',[2e-3 1e-5],...
    'xticklabel',mmpLabelsAll(rIdx),'fontsize',fsize,'TickDir','out');
xtickangle(90);
set(gca,'ytick',[1:length(mmpLabelsAll(rIdx))],'Ticklength',[2e-3 1e-5],...
    'yticklabel',mmpLabelsAll(rIdx),'fontsize',fsize,'TickDir','out');
title('Composite Interactivity CHANCE')
print(h,sprintf('figures/mmp_Adj_scoreAvgNegnsub-%i',length(Subs)),'-depsc')

fprintf('[!] Done.\n')