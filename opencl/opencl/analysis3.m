close all;

% 23 subject fit: coverage = 1.675 * n + 0.2905  

P_VALUE_FINAL = 0.01;
P_VALUE = 0.05/60;%P_VALUE_FINAL;%/n_w;
DIST_THRESH_MM = 15;
MIN_UNIQUE_SUB = 2; % number of unique subjects
FRAC_SIG_THRESH = 0.75; % fraction of roi pairs that are significant has to be at least this
%coregDir = './coreg';
%outDir = './figures';
N_ATLAS = 20;


% % Use all folders for SubjectsAll
% D = dir(coregDir);
% Disdir = [D.isdir];
% Dname = {D.name};
% Dname = Dname(Disdir);
% SubjectsAll = Dname(startsWith(Dname,'m'));
% %SubjectsAll = {'sub21','sub33','sub40'};

% Use all folders for SubjectsAll
if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
    coregDir = '/Volumes/cuenap_ssd/coregistration';
    outDir = '/Volumes/RawData/data/results/figures';
elseif isunix    
    resultsDir = '/mnt/cuenap2/data/results';
    h5Dir = '/mnt/cuenap2/scripts/synth/out';
    coregDir = '/mnt/cuenap_ssd/coregistration';
    outDir = '/mnt/cuenap2/data/results/figures';
end

% special cases
[~,host] = system('hostname');
if contains(host,'ubuntu_1604')
    resultsDir = '/nas_share/RawData/data/results';
    h5Dir = '/nas_share/RawData/scripts/synth/out';
    coregDir = '/mnt/cuenap_ssd/coregistration';
    outDir = '/nas_share/RawData/data/results/figures';
end

o = getOutputFilenames(resultsDir,true);
n_f = o.n_f;
SubjectsAllPermAll = o.SubjectsAllPerm;
SubjectsAllGraphAll = o.SubjectsAllGraph;
SubjectsAllDists = o.SubjectsAllDists;


% exclude mSu
popi1 = true(size(SubjectsAllPermAll));
popi2 = true(size(SubjectsAllGraphAll));
popi3 = true(size(SubjectsAllDists));
for i = 1:length(SubjectsAllPermAll)
    
    if startsWith(SubjectsAllPermAll{i},'mSu')
        popi1(i) = false;
    end

    if startsWith(SubjectsAllGraphAll{i},'mSu')
        popi2(i) = false;
    end
    
    if startsWith(SubjectsAllDists{i},'mSu')
        popi3(i) = false;
    end
end
SubjectsAllPermAll = SubjectsAllPermAll(popi1);
SubjectsAllGraphAll = SubjectsAllGraphAll(popi2);
SubjectsAllDists = SubjectsAllDists(popi3);
n_f = length(SubjectsAllPermAll);


% Separate by metric
Metric = cell(1,n_f);
SubjectsAllAll = cell(1,n_f);
for i = 1:n_f
    ss = strsplit(SubjectsAllPermAll{i},'-');
    ss2 = strsplit(ss{1},'_');
    Metric{i} = ss{2};
    SubjectsAllAll{i} = ss2{1};
end
Ms = unique(Metric);
Ms = Ms(randperm(length(Ms)));

for iM = 1:length(Ms)
    
    % Index into metric
    metric = Ms{iM};
    mIdx = strcmp(Metric,Ms{iM});
    SubjectsAllPerm = SubjectsAllPermAll(mIdx);
    SubjectsAllGraph = SubjectsAllGraphAll(mIdx);
    SubjectsAll = SubjectsAllAll(mIdx);

    P = cell(size(SubjectsAll));

    AdjNcov = cell(length(P),N_ATLAS);
    AdjCTaAll = cell(length(P),N_ATLAS);
    AdjMagaAll = cell(length(P),N_ATLAS);
    AdjXsub = cell(N_ATLAS,2);
    
    for i = 1:length(P)
        
        % construct filenames
        ss = strsplit(SubjectsAllPerm{i},'_');
        ss2 = strsplit(ss{2},'-');
        sid = SubjectsAll{i};
        suffix = strjoin({ss2{2},ss2{end}(1:(end-3))},'-');
        permf = sprintf('%s/%s',resultsDir,SubjectsAllPerm{i});
        graphf = sprintf('%s/%s_graph-%s.h5',resultsDir,sid,metric);
        distsf = sprintf('%s/%s_dists-%s.mat',resultsDir,sid,suffix);
        h5fname = sprintf('%s/%s.h5',h5Dir,sid);
        fprintf('[%i/%i] Processing: %s\n',i,length(P),permf)
        
        % Check all files are present
        allp = (exist(permf,'file') & exist(graphf,'file') & exist(distsf,'file') & exist(h5fname,'file'));
        if (~ allp)
            % can't compute this metric without all files
            break
        end
        
        
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
        
        % load null distributions
        fprintf('[*] Loading distributions..\n')
        df = load(distsf,'d');
        Dists = df.d;
        
        % load parcellation
        fprintf('[*] Loading parcellations..\n')
        P{i} = load(sprintf('%s/%s/label/all_parcellation.mat',coregDir,SubjectsAll{i}));
        n_atlas = length(P{i}.AtlLabels);

        % load bip montage
        %bipChan = bip(SubjectsAll{i});
        fprintf('[*] Loading bipolar coordinates..\n')
        bipChan = h5readatt(h5fname,'/h5eeg/eeg','bip');
        chan_labels = h5readatt(h5fname,'/h5eeg/eeg','labels');
        [n_bchan,~] = size(bipChan);
        % bip names
        bip_labels = cell(1,n_bchan);
        for iBn = 1:n_bchan
            bip_labels{iBn} = sprintf('%s-%s',chan_labels{iBn},chan_labels{iBn});
        end

        % Significance
        fprintf('[*] Thresholding for statistical significance..\n')
        [n_comb,n_w] = size(R);
        Rct = zeros(1,n_comb);
        Rmag = zeros(1,n_comb);
        RNegct = zeros(1,n_comb);
        RNegmag = zeros(1,n_comb);
        
        AdjCT = nan(n_bchan,n_bchan);
        AdjMag = nan(n_bchan,n_bchan);
        AdjNegCT = nan(n_bchan,n_bchan);
        AdjNegMag = nan(n_bchan,n_bchan);
        
        parfor rIdx = 1:n_comb
            Adj = R(rIdx,:);
            AdjNeg = random(Dists{rIdx},size(Adj));
            p = cdf(Dists{rIdx},Adj);
            sigIdx = ((p<(P_VALUE/2))|(p>(1-P_VALUE/2)));
            pNeg = cdf(Dists{rIdx},AdjNeg);
            sigIdxNeg = ((pNeg<(P_VALUE/2))|(pNeg>(1-P_VALUE/2)));
            
            % --- interaction strength metrics ---
            Rct(rIdx) = sum(sigIdx)/n_w;
            Rmag(rIdx) = abs(nanmean(Adj));
            RNegct(rIdx) = sum(sigIdxNeg)/n_w;
            RNegmag(rIdx) = abs(nanmean(AdjNeg));
        end
        
        % Channel adjacency matrix
        fprintf('[*] Building electrode adjacency matrix..\n')
        RNegctS = sort(RNegct);
        ct_thresh = RNegctS(end-round(P_VALUE_FINAL*(n_comb)-1)); %P_VALUE;%binopdf(5,nchoosek(n_bchan,2),P_VALUE_FINAL); sqrt((P_VALUE*(1-P_VALUE)));
        for rIdx = 1:n_comb
            % Apply significance
            if (Rct(rIdx) > ct_thresh)
                AdjCT(chan1(rIdx),chan2(rIdx)) = Rct(rIdx);
                AdjMag(chan1(rIdx),chan2(rIdx)) = Rmag(rIdx);
                AdjCT(chan2(rIdx),chan1(rIdx)) = Rct(rIdx);
                AdjMag(chan2(rIdx),chan1(rIdx)) = Rmag(rIdx);
            else
                AdjCT(chan1(rIdx),chan2(rIdx)) = 0;
                AdjMag(chan1(rIdx),chan2(rIdx)) = 0;
                AdjCT(chan2(rIdx),chan1(rIdx)) = 0;
                AdjMag(chan2(rIdx),chan1(rIdx)) = 0;
            end
            if (RNegct(rIdx) > ct_thresh)
                AdjNegCT(chan1(rIdx),chan2(rIdx)) = RNegct(rIdx);
                AdjNegMag(chan1(rIdx),chan2(rIdx)) = RNegmag(rIdx);
                AdjNegCT(chan2(rIdx),chan1(rIdx)) = RNegct(rIdx);
                AdjNegMag(chan2(rIdx),chan1(rIdx)) = RNegmag(rIdx);
            else
                AdjNegCT(chan1(rIdx),chan2(rIdx)) = 0;
                AdjNegMag(chan1(rIdx),chan2(rIdx)) = 0;
                AdjNegCT(chan2(rIdx),chan1(rIdx)) = 0;
                AdjNegMag(chan2(rIdx),chan1(rIdx)) = 0;
            end
            
        end
        
        % Map onto atlas
        fprintf('[*] Building atlas adjacency matrix..\n')
        
        for j = 1:n_atlas
            % strip extras for mmp
            if strcmp(P{i}.AtlNames{j},'HCP-MMP1')
                P{i}.AtlLabels{j} = replace(P{i}.AtlLabels{j},'_ROI','');
                P{i}.AtlLabels{j} = replace(P{i}.AtlLabels{j},'L_','');
                P{i}.AtlLabels{j} = replace(P{i}.AtlLabels{j},'R_','');
                P{i}.AtlROIs{j}.LH.struct_names = replace(P{i}.AtlROIs{j}.LH.struct_names,'_ROI','');
                P{i}.AtlROIs{j}.LH.struct_names = replace(P{i}.AtlROIs{j}.LH.struct_names,'L_','');
                P{i}.AtlROIs{j}.LH.struct_names = replace(P{i}.AtlROIs{j}.LH.struct_names,'R_','');
            end

            % rename
            Labels = P{i}.AtlLabels{j};
            Rois = P{i}.AtlROIs{j}.LH.struct_names;
            n_rois = length(Rois);
            n_chan = length(Labels);

            % build adjacency matrix for atlas
            AdjN = zeros(n_rois,n_rois);
            
            AdjCTa = cell(n_rois,n_rois);
            AdjMaga = cell(n_rois,n_rois);
            AdjNegCTa = cell(n_rois,n_rois);
            AdjNegMaga = cell(n_rois,n_rois);
            
            for i1 = 1:(n_bchan-1)
                for i2 = (i1+1):n_bchan
                    % Check if electrodes are too close
                    a = reshape(bipChan(i1,4:end),2,[]);
                    b = reshape(bipChan(i2,4:end),2,[]);
                    dist2 = zeros(2,2);
                    for a1 = 1:2
                        for b1 = 1:2
                            dist2(a1,b1) = sqrt(sum((a(a1,:) - b(b1,:)).^2));
                        end
                    end
                    if (min(dist2(:)) > DIST_THRESH_MM)
                        % electrodes are far enough apart
                        e1s = bipChan(i1,1:2);
                        e2s = bipChan(i2,1:2);
                        for a2 = 1:2
                            for b2 = 1:2
                                % electrode indices
                                e1 = e1s(a2);
                                e2 = e2s(b2);

                                % find ROI mapping for electrodes
                                roiIdx1 = contains(Rois,Labels{e1});
                                roiIdx2 = contains(Rois,Labels{e2});
                                if ((sum(roiIdx1) == 1) && (sum(roiIdx2) == 1))
                                    AdjN(roiIdx1,roiIdx2) = AdjN(roiIdx1,roiIdx2) + 1;
                                    AdjN(roiIdx2,roiIdx1) = AdjN(roiIdx2,roiIdx1) + 1;
                                    
                                    AdjCTa{roiIdx1,roiIdx2} = [AdjCTa{roiIdx1,roiIdx2},AdjCT(i1,i2)];
                                    AdjCTa{roiIdx2,roiIdx1} = [AdjCTa{roiIdx2,roiIdx1},AdjCT(i1,i2)];
                                    AdjMaga{roiIdx1,roiIdx2} = [AdjMaga{roiIdx1,roiIdx2},AdjMag(i1,i2)];
                                    AdjMaga{roiIdx2,roiIdx1} = [AdjMaga{roiIdx2,roiIdx1},AdjMag(i1,i2)];
                                    
                                    AdjNegCTa{roiIdx1,roiIdx2} = [AdjNegCTa{roiIdx1,roiIdx2},AdjNegCT(i1,i2)];
                                    AdjNegCTa{roiIdx2,roiIdx1} = [AdjNegCTa{roiIdx2,roiIdx1},AdjNegCT(i1,i2)];
                                    AdjNegMaga{roiIdx1,roiIdx2} = [AdjNegMaga{roiIdx1,roiIdx2},AdjNegMag(i1,i2)];
                                    AdjNegMaga{roiIdx2,roiIdx1} = [AdjNegMaga{roiIdx2,roiIdx1},AdjNegMag(i1,i2)];
                                end
                            end
                        end
                    end
                end
            end
            
%             % Render per subject atlas adj matrix
%             AdjCTaAvg = cell(n_rois,n_rois);
%             AdjMagaAvg = cell(n_rois,n_rois);
%             AdjNegCTaAvg = cell(n_rois,n_rois);
%             AdjNegMagaAvg = cell(n_rois,n_rois);
%             for ii1 = 1:n_rois
%                 for ii2 = 1:n_rois
%                     AdjCTaAvg(ii1,ii2) = AdjCTa{ii1,ii2};
%                     AdjMagaAvg(ii1,ii2) = AdjMaga{ii1,ii2};
%                 end
%             end
%             AdjCTaAvg(isnan(AdjCTaAvg)) = 0;
%             AdjMagaAvg(isnan(AdjMagaAvg)) = 0;

            AdjCTaAll{i,j} = AdjCTa;
            AdjMagaAll{i,j} = AdjMaga;
            AdjNcov{i,j} = AdjN;

            %fprintf('\t%i of %i subjects, %i of %i atlases\n',i,length(P),j,n_atlas)
            %figure; imagesc(AdjN); colorbar; title(sprintf('%s: %s',SubjectsAll{i},P{i}.AtlNames{j}))
        end
    end
    
    
    
    %return
 
    for j = 1:n_atlas
        Rois = P{1}.AtlROIs{j}.LH.struct_names;
        n_rois = length(Rois);
        
        Usub = zeros(n_rois,n_rois);
        FracSig = zeros(n_rois,n_rois);
        AdjCTAvg = nan(n_rois,n_rois);
        AdjCTAvgAll = cell(n_rois,n_rois);
        AdjMagAvg = nan(n_rois,n_rois);
        AdjMagAvgAll = cell(n_rois,n_rois);
        for iii1 = 1:n_rois
            for iii2 = 1:n_rois
                for i = 1:length(P)
                    ct = AdjCTaAll{i,j}{iii1,iii2};
                    FracSig(iii1,iii2) = sum(ct ~= 0)/numel(ct);
                    AdjCTAvgAll{iii1,iii2} = [AdjCTAvgAll{iii1,iii2}, ct];
                    AdjMagAvgAll{iii1,iii2} = [AdjMagAvgAll{iii1,iii2}, AdjMagaAll{i,j}{iii1,iii2}];
                    Usub(iii1,iii2) = Usub(iii1,iii2) + (AdjNcov{i,j}(iii1,iii2) ~= 0);
                end
            end
        end
        %FracSig(isnan(FracSig)) = 0;
        condSig = (FracSig >= FRAC_SIG_THRESH);
        cond = (Usub >= MIN_UNIQUE_SUB);
        condNeg = (Usub < MIN_UNIQUE_SUB) & (Usub > 0);
        
        for iii1 = 1:n_rois
            for iii2 = 1:n_rois
                if (cond(iii1,iii2))
                    AdjCTAvg(iii1,iii2) = nanmean(AdjCTAvgAll{iii1,iii2});
                    magT = AdjMagAvgAll{iii1,iii2};
                    AdjMagAvg(iii1,iii2) = nanmean(magT(~isinf(magT)));
                end
            end
        end
        
        % Set nonsignificant areas to -1
        AdjCTAvg(~condSig & ~isnan(FracSig)) = -1;
        AdjMagAvg(~condSig & ~isnan(FracSig)) = -1;
        
        % Set 0 coverage areas to zero
        AdjCTAvg(condNeg) = 0;
        AdjMagAvg(condNeg) = 0;
        
        % Save
        AdjXsub{j,1} = AdjCTAvg;
        AdjXsub{j,2} = AdjMagAvg;

        % --- Plot CT ------------------------------------
        if (j == 6) % mmp
            roi_font_size = 3;
        else
            roi_font_size = 8;
        end
        colorbar_font_size = 8;
        h = figure;
        set(h,'Position',[0 0 1080 900])
        v = AdjCTAvg;
        map = corrcmap(100);
        minv = min(v(:));
        maxv = max(v(:));
        ncol = size(map,1);
        s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
        Im = ind2rgb(s,map);
        for i3 = 1:n_rois
            for j3 = 1:n_rois
                % Areas not significantly interacting
                if (v(i3,j3) == -1)
                    Im(i3,j3,:) = [1 1 1];
                end
                % Areas with not enough subject coverage
                if (v(i3,j3) == 0)
                    Im(i3,j3,:) = [0.65 0.65 0.65];
                end
                % Areas not covered by any electrodes
                if isnan(v(i3,j3))
                    Im(i3,j3,:) = [0.4 0.4 0.4];
                end
            end
        end
        % Diagonal
        for i3 = 1:n_rois
            Im(i3,i3,:) = [0 0 0];
        end
        imagesc(Im);
        colormap(map);
        caxis([minv maxv]);
        colorbar;
        fsize = roi_font_size;
        
        set(gca,'xtick',1:length(Rois),'Ticklength',[2e-3 1e-5],...
            'xticklabel',replace(Rois,'_','-'),'fontsize',fsize,'TickDir','out');
        xtickangle(90);
        set(gca,'ytick',1:length(Rois),'Ticklength',[2e-3 1e-5],...
            'yticklabel',replace(Rois,'_','-'),'fontsize',fsize,'TickDir','out');
        title(sprintf('Cross-subject Time Consistency: %s (%i), n = %i',P{1}.AtlNames{j},length(Rois),length(P)),'fontsize',8)
        c = colorbar; c.FontSize = colorbar_font_size;
        system(sprintf('mkdir %s/adjacency_matrices',outDir));
        system(sprintf('mkdir %s/adjacency_matrices/xsub',outDir));
        system(sprintf('mkdir %s/adjacency_matrices/xsub/%s',outDir,suffix));
        print(h,sprintf('%s/adjacency_matrices/xsub/%s/AdjCTxsub_%s_%s_n-%i',outDir,suffix,suffix,P{i}.AtlNames{j},length(P)),'-depsc')
        close(h)
        
        % --- Plot Mag ------------------------------------
        if (j == 6) % mmp
            roi_font_size = 3;
        else
            roi_font_size = 8;
        end
        colorbar_font_size = 8;
        h = figure;
        set(h,'Position',[0 0 1080 900])
        v = AdjMagAvg;
        map = corrcmap(100);
        minv = min(v(:));
        maxv = max(v(:));
        ncol = size(map,1);
        s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
        Im = ind2rgb(s,map);
        for i3 = 1:n_rois
            for j3 = 1:n_rois
                % Areas not significantly interacting
                if (v(i3,j3) == -1)
                    Im(i3,j3,:) = [1 1 1];
                end
                % Areas with not enough subject coverage
                if (v(i3,j3) == 0)
                    Im(i3,j3,:) = [0.65 0.65 0.65];
                end
                % Areas not covered by any electrodes
                if isnan(v(i3,j3))
                    Im(i3,j3,:) = [0.4 0.4 0.4];
                end
            end
        end
        % Diagonal
        for i3 = 1:n_rois
            Im(i3,i3,:) = [0 0 0];
        end
        imagesc(Im);
        colormap(map);
        caxis([minv maxv]);
        colorbar;
        fsize = roi_font_size;
        
        set(gca,'xtick',1:length(Rois),'Ticklength',[2e-3 1e-5],...
            'xticklabel',replace(Rois,'_','-'),'fontsize',fsize,'TickDir','out');
        xtickangle(90);
        set(gca,'ytick',1:length(Rois),'Ticklength',[2e-3 1e-5],...
            'yticklabel',replace(Rois,'_','-'),'fontsize',fsize,'TickDir','out');
        title(sprintf('Cross-subject Average Magnitude: %s (%i), n = %i',P{1}.AtlNames{j},length(Rois),length(P)),'fontsize',8)
        c = colorbar; c.FontSize = colorbar_font_size;
        print(h,sprintf('%s/adjacency_matrices/xsub/%s/AdjMagxsub_%s_%s_n-%i',outDir,suffix,suffix,P{i}.AtlNames{j},length(P)),'-depsc')
        close(h)
    end
    
    %return

%     %%
%     % Build atlas coverage matrices
%     for j = 1:n_atlas
%         Labels = P{1}.AtlLabels{j};
%         Rois = P{1}.AtlROIs{j}.LH.struct_names;
%         n_rois = length(Rois);
% 
%         AdjNa = zeros(n_rois,n_rois);
%         AdjNa2 = zeros(n_rois,n_rois);
%         AdjAtCT = cell(n_rois,n_rois);
%         AdjAtMag = cell(n_rois,n_rois);
%         for i = 1:length(P)
%             AdjNa = AdjNa + AdjNcov{i,j};
%             AdjNa2 = AdjNa2 + (AdjNcov{i,j} ~= 0);
%             for iii1 = 1:n_rois
%                 for iii2 = 1:n_rois
%                     el = AdjCTaAvgs{i,j}(iii1,iii2);
%                     if (el ~= 0)
%                         AdjAtCT{iii1,iii2} = [AdjAtCT{iii1,iii2},el];
%                     end
%                     el2 = AdjMagaAvgs{i,j}(iii1,iii2);
%                     if (el2 ~= 0)
%                         AdjAtMag{iii1,iii2} = [AdjAtMag{iii1,iii2},el2];
%                     end
%                     
%                 end
%             end
%         end
% 
%         % plot massive consistency across time matrix
%         h = figure;
%         set(gcf, 'PaperUnits', 'inches');
%         set(gcf, 'PaperPosition', [0 0 36 36]);
%         bins = linspace(0,1,11);
%         for jj1 = 1:n_rois
%             for jj2 = 1:n_rois
%                 if (~isempty(AdjAtCT{jj1,jj2}))
%                     subplot(n_rois,n_rois,(jj1-1)*n_rois+jj2);
%                     [f,x] = hist(AdjAtCT{jj1,jj2},bins);
%                     plot(x,f,'black-');
%                     title(sprintf('%s/%s',replace(Rois{jj1},'_','-'),replace(Rois{jj2},'_','-')))
%                     set(gca,'fontsize',4)
%                     xlabel('CorrCoeff')
%                     ylabel('Number of subjects')
%                     axis tight;
%                 end
%             end
%         end
%         print(h,sprintf('%s/AdjCT_xsub-%s-%i_%s',outDir,P{i}.AtlNames{j},length(P),suffix),'-depsc')
%         close(h)
%         
%         
%         h = figure;
%         set(gcf, 'PaperUnits', 'inches');
%         set(gcf, 'PaperPosition', [0 0 36 36]);
%         bins = linspace(0,1,11);
%         for jj1 = 1:n_rois
%             for jj2 = 1:n_rois
%                 if (~isempty(AdjAtMag{jj1,jj2}))
%                     subplot(n_rois,n_rois,(jj1-1)*n_rois+jj2);
%                     [f,x] = hist(AdjAtMag{jj1,jj2},bins);
%                     plot(x,f,'black-');
%                     title(sprintf('%s/%s',replace(Rois{jj1},'_','-'),replace(Rois{jj2},'_','-')))
%                     set(gca,'fontsize',4)
%                     xlabel('CorrCoeff')
%                     ylabel('Number of subjects')
%                     axis tight;
%                 end
%             end
%         end
%         print(h,sprintf('%s/AdjMag_xsub-%s-%i_%s',outDir,P{i}.AtlNames{j},length(P),suffix),'-depsc')
%         close(h)
        
        
        
%         % percent coverage
%         pctCov1 = 100*sum(sum(AdjNa2 > 0))/numel(AdjNa2);
% 
%         % percent coverage with at least 2 subjects
%         pctCov2 = 100*sum(sum(AdjNa2 > 1))/numel(AdjNa2);
        
% 
%         h = figure;
%         rIdx = ~(sum((AdjNa == 0)) == n_rois);
%         imagesc(AdjNa(rIdx,rIdx))
%         cc = parula(10*64);
%         cc = [[1 1 1];cc];
%         colormap(cc);
%         fsize = 3;
%         set(gca,'xtick',1:length(Rois(rIdx)),'Ticklength',[2e-3 1e-5],...
%             'xticklabel',replace(Rois(rIdx),'_','-'),'fontsize',fsize,'TickDir','out');
%         xtickangle(90);
%         set(gca,'ytick',1:length(Rois(rIdx)),'Ticklength',[2e-3 1e-5],...
%             'yticklabel',replace(Rois(rIdx),'_','-'),'fontsize',fsize,'TickDir','out');
%         title(sprintf('%s (%i) coverage for %i patients, %.1f%% coverage, %.1f%% double coverage',...
%             P{i}.AtlNames{j},length(Rois),length(P),pctCov1,pctCov2),'fontsize',6)
%         c = colorbar; c.FontSize = 8;
%         print(h,sprintf('%s/coverage-%s-%i_%s',outDir,P{i}.AtlNames{j},length(P),suffix),'-depsc')
%         close(h)
%         fprintf('[*] Finished %s: %s\n',SubjectsAll{i},P{i}.AtlNames{j})
% 
%         h = figure;
%         rIdx = ~(sum((AdjNa2 == 0)) == n_rois);
%         imagesc(AdjNa2(rIdx,rIdx))
%         cc = parula(10*64);
%         cc = [[1 1 1];cc];
%         colormap(cc);
%         fsize = 3;
%         set(gca,'xtick',1:length(Rois(rIdx)),'Ticklength',[2e-3 1e-5],...
%             'xticklabel',replace(Rois(rIdx),'_','-'),'fontsize',fsize,'TickDir','out');
%         xtickangle(90);
%         set(gca,'ytick',1:length(Rois(rIdx)),'Ticklength',[2e-3 1e-5],...
%             'yticklabel',replace(Rois(rIdx),'_','-'),'fontsize',fsize,'TickDir','out');
%         title(sprintf('%s (%i) patient-patient overlap for %i patients, %.1f%% coverage, %.1f%% double coverage',...
%             P{i}.AtlNames{j},length(Rois),length(P),pctCov1,pctCov2),'fontsize',6)
%         c = colorbar; c.FontSize = 8;
%         print(h,sprintf('%s/coverage2-%s-%i_%s',outDir,P{i}.AtlNames{j},length(P),suffix),'-depsc')
%         close(h)
%         fprintf('[*] Finished %s: %s\n',SubjectsAll{i},P{i}.AtlNames{j})
%     end

end

fprintf('[!] All Done.\n')
