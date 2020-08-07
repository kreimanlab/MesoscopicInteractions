classdef Analysis
    %Analysis
    
    properties
        col
        const
        r
        h5eeg
    end
    
    methods
        %
        %   Analysis class constructor
        %
        function self = Analysis(results_dir, h5_dir)
            % initialize colors
            self.col.ADJ_N_SIG = 1 * [1 1 1]; %0.667 * [0.97 1 0.97];
            self.col.ADJ_D_THRESH = 0.333 * [1 1 1];
            self.col.ADJ_NAN = 2*0.333 * [1 1 1];
            self.col.ADJ_DIAG = 0.333 * [1 1 1];
            self.col.COLORMAP = corrcmap(100);
            self.col.SIG_ISSIG = self.col.COLORMAP(end,:);%0.3*[1 0 0];
            self.col.SIG_NOTSIG = 0.667*[1 1 1]; % 0.6
            self.col.SIG_ART = 0.9*[1 1 1];
            
            % initialize constants
            % 0.05 false positives per 60 minutes
            self.const.P_VALUE = 0.05;%0.05/(60);
            % at most 1 connection marked false positive per 51 subjects
            self.const.P_VALUE_CT = 0.05;%(0.05/nchoosek(176,2))/51;
            % at most 1 false positive per 20 atlases
            self.const.P_VALUE_CS = (0.05/nchoosek(180,2))/20;
            % Maximum metric score allowed (absolute value)
            self.const.R_MAX = 0.8;
            % Minimum distance between electrodes to threshold
            self.const.DIST_THRESH_MM = 20;
            self.const.FIG_FILETYPE = 'pdf';
            self.const.FIG_DPI = '-r400';
            
            % initialize frequency bands
            self.const.DEL_S = 0.5;
            self.const.DEL_E = 3;
            self.const.THE_S = 3;
            self.const.THE_E = 8;
            self.const.ALP_S = 8;
            self.const.ALP_E = 12;
            self.const.BET_S = 12;
            self.const.BET_E = 25;
            self.const.GAM_S = 25;
            self.const.GAM_E = 100;
            self.const.HG_S = 70;
            self.const.HG_E = 200;

            % get result files
            o = getOutputFilenames(results_dir,true);         
            subjects = cell(1,o.n_f);
            metrics = cell(1,o.n_f);
            for iSub = 1:(o.n_f)
                ss = strsplit(o.SubjectsAllPerm{iSub},'_');
                ss2 = strsplit(ss{2},'-');
                metrics{iSub} = ss2{2};
                sid = ss{1};
                subjects{iSub} = sid;
            end
            self.r.n_f = o.n_f;
            self.r.perm_f = o.SubjectsAllPerm;
            self.r.graph_f = o.SubjectsAllGraph;
            self.r.dist_f = o.SubjectsAllDists;
            self.r.subjects = subjects;
            self.r.metrics = metrics;
            self.r.results_dir = results_dir;
            

            % get subject files
            subjectsU = unique(subjects);
            n_sub = length(subjectsU);
            self.h5eeg.n_sub = n_sub;
            self.h5eeg.subjects = subjectsU;
            
            % init variables
            self.h5eeg.filenames = cell(1,n_sub);
            self.h5eeg.fs = cell(1,n_sub);
            for iSub = 1:n_sub
                % Read h5eeg
                h5fname = sprintf('%s/%s.h5',h5_dir,subjectsU{iSub});
                fs = h5readatt(h5fname,'/h5eeg/eeg','rate');
                n_chan = double(h5readatt(h5fname,'/h5eeg/eeg','n_chan'));
                n_samples = double(h5readatt(h5fname,'/h5eeg/eeg','n_samples'));
                chan_labels = h5readatt(h5fname,'/h5eeg/eeg','labels');
                h5_n_samples = double(h5readatt(h5fname,'/h5eeg','n_samples'));
                bip = h5readatt(h5fname,'/h5eeg/eeg','bip');
                %arts = h5read(h5fname,'/h5eeg/artifacts');
                %width = double(h5readatt(h5fname,'/h5eeg/artifacts','width'));
                %bipChan = h5readatt(h5fname,'/h5eeg/eeg','bip');
                [n_bchan,~] = size(bip);
                
                bip_labels = cell(1,n_bchan);
                for iBn = 1:n_bchan
                    bip_labels{iBn} = sprintf('%s-%s',chan_labels{bip(iBn,1)},chan_labels{bip(iBn,2)});
                end
                
                % Save variables
                self.h5eeg.filenames{iSub} = h5fname;
                self.h5eeg.fs{iSub} = fs;
                self.h5eeg.n_chan{iSub} = n_chan;
                self.h5eeg.n_bchan{iSub} = n_bchan;
                self.h5eeg.n_comb{iSub} = nchoosek(n_bchan,2);
                self.h5eeg.n_samples{iSub} = n_samples;
                self.h5eeg.chan_labels{iSub} = chan_labels;
                self.h5eeg.h5_n_samples{iSub} = h5_n_samples;
                self.h5eeg.bip{iSub} = bip;
                self.h5eeg.bip_labels{iSub} = bip_labels;
                
                %fprintf('read: %i of %i\n',iSub,length(subjects))
            end
            self.h5eeg.h5_dir = h5_dir;
            
        end
        
        %
        %   getCT - calculate consistency across time
        %       
        %   inputs:
        %       sIdx -      subject index
        %       metric -    metric name
        %       is_null -   where or not to sample from null distribution
        %
        %   outputs:
        %       AdjCT -     the consistency across time matrix taking
        %                   values from 0 to 1
        %       AdjMag -    the average magnitude of the given metric from
        %                   0 to 1
        %
        function [CT] = getCT(self, sIdx, metric, is_null)
            % Find file matrix
            fIdx = NaN;
            for i = 1:self.r.n_f
                pcond = contains(self.r.dist_f{i},[self.h5eeg.subjects{sIdx},'_']) ...
                    & contains(self.r.dist_f{i},['dists-',metric,'-']);
                if (pcond)
                    fIdx = i;
                end
            end
            if (isnan(fIdx))
                fprintf(2,'W: Analysis.getCT(%i,''%s'') could not find metric and subject combo.\n',sIdx,metric)
                CT = [];
                return;
            end
            
            % Read distribution matrix
            D_filename = sprintf('%s/%s',self.r.results_dir,self.r.dist_f{fIdx});
            R_filename = sprintf('%s/%s',self.r.results_dir,self.r.graph_f{fIdx});
            fprintf('Reading dist file: %s.\n',D_filename)
            fprintf('Reading graph file: %s.\n',R_filename)
            D = load(D_filename);
            R = h5read(R_filename,'/R');
            w = double(h5read(R_filename,'/w'));
            
            % Read artifacts
            if (~ (is_null) )
                h5fname = self.h5eeg.filenames{sIdx};
                arts = h5read(h5fname,'/h5eeg/artifacts');
                arts_width = double(h5readatt(h5fname,'/h5eeg/artifacts','width'));
            end
            
            % Adjust for magnitude-square coherence
            if (startsWith(metric,'pc') || startsWith(metric,'sc'))
                R = sqrt(R);
            end
            
            % Precalculations
            [n_comb,n_w] = size(R);
            n_bchan = self.h5eeg.n_bchan{sIdx};
            P_VALUE = self.const.P_VALUE;
            AdjCT = nan(n_bchan,n_bchan);
            AdjMag = nan(n_bchan,n_bchan);
            AdjAll = cell(n_bchan,n_bchan);
            SigIdx = cell(n_bchan,n_bchan);
            ArtIdx = cell(n_bchan,n_bchan);

            % Parallel cdf calculation
            p_all = cell(1,n_comb);
            if (is_null)
                Adj_all = cell(1,n_comb);
            end
            parfor ip = 1:n_comb
                if (is_null)
                    Adj = random(D.d{ip},size(R(ip,:)));
                    Adj_all{ip} = Adj;
                else
                    Adj = R(ip,:);
                end
                p_all{ip} = cdf(D.d{ip},Adj);
            end

            rIdx = 1;
            % for optimizing artifact insignificance index
            if (~(is_null))
                arts_n_zero = (arts(1,:)~=0);
                arts_index = arts(2,arts_n_zero);
            end
            for i = 1:(n_bchan-1)
                for j = (i+1):n_bchan
                    if (is_null)
                        %Adj = random(D.d{rIdx},size(R(rIdx,:)));
                        Adj = Adj_all{rIdx};
                    else
                        Adj = R(rIdx,:);
                    end
                    %p = cdf(D.d{rIdx},Adj);
                    p = p_all{rIdx};
                    sigIdx = ((p<(P_VALUE/2))|(p>(1-P_VALUE/2)));
                    
                    % Threshold max metric and set as insignificant
                    sigIdx(abs(Adj) > self.const.R_MAX) = false;
                    
                    if (~ (is_null) )
                        %Set artifact areas as insignificant
                        %arts_index = arts(2,(arts(1,:)~=0));
                        %w_samp = round(w*self.h5eeg.fs{sIdx});
                        %r_index = unique(ceil(arts_index/w_samp));
                        
                        % optimized for speed
                        w_samp = round(w*self.h5eeg.fs{sIdx});
                        r_index = ceil(arts_index/w_samp);
                        sigIdx(r_index) = false;
                    end
                    
                    % --- interaction strength metrics ---
                    adjct_ij = sum(sigIdx)/n_w;
                    
                    % --- Magnitude ---------------------------------------
                    adjmag_ij = abs(nanmean(Adj(sigIdx)));
                    % -----------------------------------------------------
                    
                    AdjCT(i,j) = adjct_ij;
                    AdjCT(j,i) = adjct_ij;
                    AdjMag(i,j) = adjmag_ij;
                    AdjMag(j,i) = adjmag_ij;
                    AdjAll{i,j} = Adj;
                    SigIdx{i,j} = sigIdx;
                    if (~ (is_null))
                        ArtIdx{i,j} = r_index;
                    else
                        ArtIdx{i,j} = [];
                    end
                    
                    rIdx = rIdx + 1;
                end
            end
            
            CT.Adj = AdjAll;
            CT.AdjCT = AdjCT;
            CT.AdjMag = AdjMag;
            CT.SigIdx = SigIdx;
            CT.ArtIdx = ArtIdx;
            CT.n_w = n_w;
            CT.w = w;
            CT.metric = metric;
            CT.is_null = is_null;
            
        end
        
        %
        %   getAT - makes atlas-space adjacency matrices
        %
        function [AT] = getAT(self, sIdx, CT, dmat)
            % Get coregistration directory
            [~,subjects_dir] = system('echo $SUBJECTS_DIR');
            coregDir = strip(subjects_dir);
            
            % Load parcellations
            P = load(sprintf('%s/%s/label/all_parcellation.mat',...
                coregDir,self.h5eeg.subjects{sIdx}));
            n_atlas = length(P.AtlLabels);
            n_bchan = self.h5eeg.n_bchan{sIdx};
            
            % init
            AdjNcov = cell(1,n_atlas);
            AdjCTaAll = cell(1,n_atlas);
            AdjMagaAll = cell(1,n_atlas);
            AdjBchanAll = cell(1,n_atlas);
            AdjDistAll = cell(1,n_atlas);
            AdjSubAll = cell(1,n_atlas);
            N_roisAll = zeros(1,n_atlas);
            
            for j = 1:n_atlas
                % strip extras for mmp
                if strcmp(P.AtlNames{j},'HCP-MMP1')
                    P.AtlLabels{j} = replace(P.AtlLabels{j},'_ROI','');
                    P.AtlLabels{j} = replace(P.AtlLabels{j},'L_','');
                    P.AtlLabels{j} = replace(P.AtlLabels{j},'R_','');
                    P.AtlROIs{j}.LH.struct_names = replace(P.AtlROIs{j}.LH.struct_names,'_ROI','');
                    P.AtlROIs{j}.LH.struct_names = replace(P.AtlROIs{j}.LH.struct_names,'L_','');
                    P.AtlROIs{j}.LH.struct_names = replace(P.AtlROIs{j}.LH.struct_names,'R_','');
                end

                % rename
                Labels = P.AtlLabels{j};
                Rois = P.AtlROIs{j}.LH.struct_names;
                n_rois = length(Rois);
                N_roisAll(j) = n_rois;

                % build adjacency matrix for atlas
                AdjN = zeros(n_rois,n_rois);
                AdjCTa = cell(n_rois,n_rois);
                AdjMaga = cell(n_rois,n_rois);
                AdjBchan = cell(n_rois,n_rois);
                AdjDist = cell(n_rois,n_rois);
                AdjSub = cell(n_rois,n_rois);
                
                % Loop through bipolar pairs
                for i1 = 1:(n_bchan-1)
                    for i2 = (i1+1):n_bchan
                        
                        % if electrodes are far enough apart
                        if (dmat(i1,i2) > self.const.DIST_THRESH_MM)
                            e1s = self.h5eeg.bip{sIdx}(i1,1:2);
                            e2s = self.h5eeg.bip{sIdx}(i2,1:2);
                            
                            % loop through combinations of real electrodes
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

                                        AdjCTa{roiIdx1,roiIdx2} = [AdjCTa{roiIdx1,roiIdx2}, CT.AdjCT(i1,i2)];
                                        AdjCTa{roiIdx2,roiIdx1} = [AdjCTa{roiIdx2,roiIdx1}, CT.AdjCT(i1,i2)];
                                        AdjMaga{roiIdx1,roiIdx2} = [AdjMaga{roiIdx1,roiIdx2}, CT.AdjMag(i1,i2)];
                                        AdjMaga{roiIdx2,roiIdx1} = [AdjMaga{roiIdx2,roiIdx1}, CT.AdjMag(i1,i2)];
                                        AdjBchan{roiIdx1,roiIdx2} = [AdjBchan{roiIdx1,roiIdx2}, [i1;i2]];
                                        AdjBchan{roiIdx2,roiIdx1} = [AdjBchan{roiIdx2,roiIdx1}, [i1;i2]];
                                        %dmat(i1,i2)
                                        AdjDist{roiIdx1,roiIdx2} = [AdjDist{roiIdx1,roiIdx2}, dmat(i1,i2)];
                                        AdjDist{roiIdx2,roiIdx1} = [AdjDist{roiIdx2,roiIdx1}, dmat(i1,i2)];
                                        % subjects
                                        AdjSub{roiIdx1,roiIdx2} = [AdjSub{roiIdx1,roiIdx2}, self.h5eeg.subjects{sIdx}];
                                        AdjSub{roiIdx2,roiIdx1} = [AdjSub{roiIdx2,roiIdx1}, self.h5eeg.subjects{sIdx}];
                                        
%                                     else
%                                         fprintf(2,'Could not map: %s - %s\n',Labels{e1},Labels{e2})
                                    end
                                end
                            end
                        end
                    end
                end

                AdjCTaAll{1,j} = AdjCTa;
                AdjMagaAll{1,j} = AdjMaga;
                AdjNcov{1,j} = AdjN;
                AdjBchanAll{1,j} = AdjBchan;
                AdjDistAll{1,j} = AdjDist;
                AdjSubAll{1,j} = AdjSub;

                %fprintf('\t%i of %i subjects, %i of %i atlases\n',i,length(P),j,n_atlas)
                %figure; imagesc(AdjN); colorbar; title(sprintf('%s: %s',SubjectsAll{i},P{i}.AtlNames{j}))
            end
            
            AT.P = P;
            AT.AdjCT = AdjCTaAll;
            AT.AdjMag = AdjMagaAll;
            AT.Bchan = AdjBchanAll;
            AT.Dist = AdjDistAll;
            AT.Sub = AdjSubAll;
            AT.Ncov = AdjNcov;
            AT.n_atlas = n_atlas;
            AT.n_rois = N_roisAll;
        end
        
        %
        %   getDistanceMatrix - makes the inter-electrode distance matrix
        %       
        %   inputs:
        %       sIdx -  subject index
        %
        %   outputs:
        %       Dmat -  the distance matrix where the ij-th element is the
        %               shortest distance between bipolar channels i, j
        %
        function Dmat = getDistanceMatrix(self, sIdx)
            % Make distance matrix
            Dmat = zeros(self.h5eeg.n_bchan{sIdx},self.h5eeg.n_bchan{sIdx});
            for i2 = 1:length(Dmat)
                for j2 = 1:length(Dmat)
                    a = reshape(self.h5eeg.bip{sIdx}(i2,4:end),2,[]);
                    b = reshape(self.h5eeg.bip{sIdx}(j2,4:end),2,[]);
                    dist2 = zeros(2,2);
                    for a1 = 1:2
                        for b1 = 1:2
                            dist2(a1,b1) = sqrt(sum((a(a1,:) - b(b1,:)).^2));
                        end
                    end
                    Dmat(i2,j2) = min(dist2(:));
                end
            end
        end
        
        %
        %   saveFig - saves figure
        %       
        %   inputs:
        %       h -         figure handle
        %       filename -  name of output figure to save
        %
        function saveFig(self, h, filename)
            print(h, filename, ['-d',self.const.FIG_FILETYPE], self.const.FIG_DPI);
        end
        
        %
        %   plotAdjElectrode - plots adjacency matrix in electrode space
        %       
        %   inputs:
        %       sIdx -  subject index
        %       Adj -   Adjacency matrix to plot
        %       ttl -   text to print as title
        %
        %   outputs:
        %       h -     figure handle
        %
        function [ax2] = plotAdjElectrode(self, ax, sIdx, Adj, ttl)
            %h = figure;
            
            % generate image from matrix
            v = Adj;
            map = self.col.COLORMAP;
            minv = min(v(:));
            maxv = max(v(:));
            ncol = size(map,1);
            s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
            Im = ind2rgb(s,map);
            
            % plot image
            imagesc(Im, 'Parent', ax);
            
            % label axes
            bip_labels = self.h5eeg.bip_labels{sIdx};
            if (length(bip_labels) <= 25)
                fsize = 8;
            elseif (length(bip_labels) <= 50)
                fsize = 8;
            elseif (length(bip_labels) <= 75)
                fsize = 7;
            elseif (length(bip_labels) <= 100)
                fsize = 6;
            else
                fsize = 5;
            end
            
            stagger_labels = true;
            % Axis label staggering
            xtick_vec = 1:length(bip_labels);
            ytick_vec = 1:length(bip_labels);
            ax1_xind = true(1,length(xtick_vec));
            ax1_yind = true(1,length(ytick_vec));
            ax2_xind = true(1,length(xtick_vec));
            ax2_yind = true(1,length(ytick_vec));
            if (stagger_labels)
                odd_i = downsample(1:length(ax1_xind),2,1);
                even_i = downsample(1:length(ax1_xind),2,0);
                ax1_xind(even_i) = false;
                ax2_xind(odd_i) = false;
                ax1_yind(odd_i) = false;
                ax2_yind(even_i) = false;
            end
            
            % Show axis labels
            %ax1 = ax;
            set(ax,'xtick',xtick_vec(ax1_xind),'Ticklength',[2e-3 1e-5],...
                'xticklabel',bip_labels(ax1_xind),'fontsize',fsize,'TickDir','out');
            xtickangle(90);
            set(ax,'ytick',ytick_vec(ax1_yind),'Ticklength',[2e-3 1e-5],...
                'yticklabel',bip_labels(ax1_yind),'fontsize',fsize,'TickDir','out');

            ax2 = copyobj(ax,ax.Parent);
            set(ax2,'YAxisLocation','right');
            set(ax2,'XAxisLocation','top');
            set(ax2,'xtick',xtick_vec(ax2_xind),'Ticklength',[2e-3 1e-5],...
                'xticklabel',bip_labels(ax2_xind),'fontsize',fsize,'TickDir','out');
            set(ax2,'ytick',ytick_vec(ax2_yind),'Ticklength',[2e-3 1e-5],...
                'yticklabel',bip_labels(ax2_yind),'fontsize',fsize,'TickDir','out');

            % colormap
            colormap(ax,map)
            caxis(ax,[minv maxv]);
            colormap(ax2,map)
            caxis(ax2,[minv maxv]);
            
            % Title
            text(length(Adj)/2,-length(Adj)*0.22,ttl,'HorizontalAlignment','center')
            
        end
        
        %
        %   similar to plotAdjElectrode, but thresholded by distance
        %
        %   inputs Dmat and Adj must be same size
        %
        function [ax2] = plotAdjElectrodeDistThresh(self, ax, sIdx, Adj, ttl, Dmat, Adj4Sig, AdjSig, AT, aIdx)
            %h = figure;
            if (nargin > 6)
                do_sig_thresh = true;
            else
                do_sig_thresh = false;
            end
            
            if (nargin > 8)
                do_atlas = true;
            else
                do_atlas = false;
            end
            
            % generate image from matrix
            Adj(Dmat <= self.const.DIST_THRESH_MM) = mean(Adj(:));
            v = Adj;
            [v_n, v_m] = size(v);
            map = self.col.COLORMAP;
            minv = min(v(:));
            %maxv = self.const.R_MAX;
            maxv = max(v(:));
            ncol = size(map,1);
            s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
            Im = ind2rgb(s,map);
            n_thresh = 0;
            if (do_sig_thresh)
                n_sig_thresh = 0;
            end
            if (do_atlas)
                n_nomap = 0;
            end
            for i3 = 1:v_n
                for j3 = 1:v_m
                    % distance threshold
                    if (Dmat(i3,j3) <= self.const.DIST_THRESH_MM)
                        Im(i3,j3,:) = self.col.ADJ_D_THRESH;
                        n_thresh = n_thresh + 1;
                    else
                        if (do_sig_thresh)
                            if (Adj4Sig(i3,j3) <= AdjSig)
                                Im(i3,j3,:) = self.col.ADJ_N_SIG;
                                n_sig_thresh = n_sig_thresh + 1;
                            end
                        end
                        if (do_atlas)
                            if (isnan(Adj(i3,j3)))
                                Im(i3,j3,:) = self.col.ADJ_NAN;
                                n_nomap = n_nomap + 1;
                            end
                        end
                    end
                end
            end
            
            % plot image
            imagesc(Im, 'Parent', ax);
            
            % label axes
            if (length(Adj) == self.h5eeg.n_bchan{sIdx})
                bip_labels = self.h5eeg.bip_labels{sIdx};
            else
                bip_labels = AT.P.AtlROIs{aIdx}.RH.struct_names;
                for i2 = 1:length(bip_labels)
                    %ts = lower(replace(bip_labels{i2},'_','-'));
                    %ts(1) = upper(ts(1));
                    bip_labels{i2} = replace(bip_labels{i2},'_','-');
                    if (isempty(bip_labels{i2}))
                        bip_labels{i2} = 'UNKNOWN';
                    end
                end
                %fprintf(2,'Input adjacency matrix is not in bip space.\n')
            end
            if (length(bip_labels) <= 25)
                fsize = 9;
            elseif (length(bip_labels) <= 50)
                fsize = 9;
            elseif (length(bip_labels) <= 75)
                fsize = 7;
            elseif (length(bip_labels) <= 100)
                fsize = 6;
            else
                fsize = 5;
            end

            stagger_labels = true;
            
            % Axis label staggering
            xtick_vec = 1:length(bip_labels);
            ytick_vec = 1:length(bip_labels);
            ax1_xind = true(1,length(xtick_vec));
            ax1_yind = true(1,length(ytick_vec));
            ax2_xind = true(1,length(xtick_vec));
            ax2_yind = true(1,length(ytick_vec));
            if (stagger_labels)
                odd_i = downsample(1:length(ax1_xind),2,1);
                even_i = downsample(1:length(ax1_xind),2,0);
                ax1_xind(even_i) = false;
                ax2_xind(odd_i) = false;
                ax1_yind(odd_i) = false;
                ax2_yind(even_i) = false;
            end
            
            % Show axis labels
            %ax1 = ax;
            set(ax,'xtick',xtick_vec(ax1_xind),'Ticklength',[2e-3 1e-5],...
                'xticklabel',bip_labels(ax1_xind),'fontsize',fsize,'TickDir','out');
            xtickangle(90);
            set(ax,'ytick',ytick_vec(ax1_yind),'Ticklength',[2e-3 1e-5],...
                'yticklabel',bip_labels(ax1_yind),'fontsize',fsize,'TickDir','out');

            ax2 = copyobj(ax,ax.Parent);
            set(ax2,'YAxisLocation','right');
            set(ax2,'XAxisLocation','top');
            set(ax2,'xtick',xtick_vec(ax2_xind),'Ticklength',[2e-3 1e-5],...
                'xticklabel',bip_labels(ax2_xind),'fontsize',fsize,'TickDir','out');
            set(ax2,'ytick',ytick_vec(ax2_yind),'Ticklength',[2e-3 1e-5],...
                'yticklabel',bip_labels(ax2_yind),'fontsize',fsize,'TickDir','out');
            
            % colormap
            colormap(ax,map)
            colormap(ax2,map)
            if (minv == maxv)
                caxis(ax,[(minv - 1e-6), (maxv + 1e-6)]); 
                caxis(ax2,[(minv - 1e-6), (maxv + 1e-6)]); 
            else
                caxis(ax,[minv maxv]);
                caxis(ax2,[minv maxv]);
            end
            
            % Title
            if (do_sig_thresh)
                if (do_atlas)
                    pct_sig = 100*(1-n_sig_thresh/(v_n*v_m-n_thresh-n_nomap));
                else
                    pct_sig = 100*(1-n_sig_thresh/(v_n*v_m-n_thresh));
                end
                ttltxt = sprintf('%s - Threshold: %.2f mm (%.2f%%)(%.2f%% significant)',ttl,...
                    self.const.DIST_THRESH_MM, 100*n_thresh/(v_n*v_m), pct_sig );
            else
                ttltxt = sprintf('%s - Threshold: %.2f mm (%.2f%%)',ttl,...
                self.const.DIST_THRESH_MM, 100*n_thresh/(v_n*v_m));
            end
            
            text_shift = 0.22;
            if (do_atlas)
                text_shift = 0.31;
            end
            text(length(Adj)/2,-length(Adj)*text_shift,ttltxt,...
                'HorizontalAlignment','center')
        end
        
        function [] = plotLine(self, ax1, v1, sIdx)
            fs = self.h5eeg.fs{sIdx};
            T = 1:length(v1);
            T = (T - 1)/(fs);
            v1 = double(v1 - mean(v1));
            %h = figure;
            plot(T,v1,'black','Parent',ax1,'LineWidth',0.05);
            scale_font_size = 4;

            t_scale_incr = 5*10^(floor(log10(max(T)))-1);% seconds
            t_scale_frac = 2/10;
            v_scale_incr = 50; % microvolts
            v_scale_frac = 1;
            
            t_scale = round((max(T)*t_scale_frac)/t_scale_incr)*t_scale_incr;
            v_scale = round((max(abs(v1))*v_scale_frac)/v_scale_incr)*v_scale_incr; 
            if (v_scale == 0)
                v_scale = v_scale_incr;
            end
            if (t_scale == 0)
                t_scale = t_scale_incr;
            end
            scale_col = [0 0 0];
            minv = min(v1)-v_scale;
            maxv = max(v1)+v_scale;
            minv_s = minv;
            
            % Controls time gap between scale marker and start of signal
            t_shift = 0.01*(max(T)-min(T));
            T = T - t_shift;
            line([min(T) min(T)],[minv_s minv_s+v_scale],'color',scale_col)
            line([min(T) min(T)+t_scale],[minv_s minv_s],'color',scale_col)
            text(min(T),(minv_s+0*v_scale),sprintf('%i \x03BCV',v_scale),...
                'fontsize',scale_font_size,'HorizontalAlignment','left',...
                'VerticalAlignment','bottom','Rotation',90)
            if (t_scale >= 1)
                t_scale_fmt = '%.0f sec';
            else
                t_scale_fmt = '%.1f sec';
            end
            text(min(T)+0*t_scale,minv_s,sprintf(t_scale_fmt,t_scale),...
                'fontsize',scale_font_size,'VerticalAlignment','top',...
                'HorizontalAlignment','left')
            axis([min(T) max(T) 1.05*minv 1.05*maxv])
            set(gca,'Visible','off')
        end
        
        %
        % Plot significance across time
        %
        function plotSig(self, sIdx, ax, CT, bchan1, bchan2, dmat)
            fs = self.h5eeg.fs{sIdx};
            R = CT.Adj{bchan1,bchan2};
            sigidx = CT.SigIdx{bchan1,bchan2};
            artidx = CT.ArtIdx{bchan1,bchan2};
            T = 1:length(R);
            T = (T - 1)*(round(fs * CT.w)/(round(fs*3600*24)));
            
            
            % Flip sign
            if (nanmean(R) < 0)
                R = (-1)*R;
            end
            % Plot
%             plot(T(artidx),R(artidx),'.','color',self.col.SIG_ART,...
%                 'Parent',ax,'MarkerSize',3);
%             hold on;
            R(artidx) = NaN;
            plot(T(~sigidx),R(~sigidx),'.','color',self.col.SIG_NOTSIG,...
                'Parent',ax,'MarkerSize',3);
            hold on;
            plot(T(sigidx),R(sigidx),'.','color',self.col.SIG_ISSIG,...
                'Parent',ax,'MarkerSize',3);
            
            % Set axes
            Rmax = self.const.R_MAX;
            if ((~isempty(T)) && (~isempty(Rmax)))
                %disp([min(T),max(T),-Rmax,Rmax]);
                axis([min(T),max(T),-Rmax,Rmax])
            else
                axis tight;
            end
            xlabel('Days')
            ylabel(self.metric_tofull(CT.metric))
            
            % Tick marks
            fsize = 8;
            tick_len = [5e-3 1e-3];
            x_tick = min(T):1:max(T);
            x_tick_lab = cell(1,length(x_tick));
            for i = 1:length(x_tick)
                x_tick_lab{i} = sprintf('%.0f',x_tick(i));
            end
            n_ytick = 5;
            y_tick = linspace(-Rmax,Rmax,n_ytick);
            y_tick_lab = cell(1,n_ytick);
            for i = 1:n_ytick
                y_tick_lab{i} = sprintf('%.2f',y_tick(i));
            end
            set(ax,'xtick',x_tick,'Ticklength',tick_len,...
                'xticklabel',x_tick_lab,'fontsize',fsize,'TickDir','out');
            set(ax,'ytick',y_tick,'Ticklength',tick_len,...
                'yticklabel',y_tick_lab,'fontsize',fsize,'TickDir','out');
            
            % Show electrode
            if (CT.is_null)
                titlef = '%s : %s (%.0f mm) Randomly Sampled';
            else
                titlef = '%s : %s (%.0f mm)';
            end
            text(mean(T),Rmax,sprintf(titlef,...
                self.h5eeg.bip_labels{sIdx}{bchan1},...
                self.h5eeg.bip_labels{sIdx}{bchan2},...
                dmat(bchan1,bchan2)),'HorizontalAlignment','center',...
                'VerticalAlignment','bottom');
            
            % ax formatting
            box(ax,'off')
            set(ax,'Layer','top')
        end
        
        % filter Delta
        function vf = del(self, v1, sIdx)
            fs = self.h5eeg.fs{sIdx};
            v1 = double(v1 - mean(v1));
            [b,a] = butter(3, [self.const.DEL_S self.const.DEL_E]/fs);
            vf = filter(b,a,v1);
        end
        % filter Theta
        function vf = the(self, v1, sIdx)
            fs = self.h5eeg.fs{sIdx};
            v1 = double(v1 - mean(v1));
            [b,a] = butter(3, [self.const.THE_S self.const.THE_E]/fs);
            vf = filter(b,a,v1);
        end
        % filter Alpha
        function vf = alp(self, v1, sIdx)
            fs = self.h5eeg.fs{sIdx};
            v1 = double(v1 - mean(v1));
            [b,a] = butter(3, [self.const.ALP_S self.const.ALP_E]/fs);
            vf = filter(b,a,v1);
        end
        % filter Beta
        function vf = bet(self, v1, sIdx)
            fs = self.h5eeg.fs{sIdx};
            v1 = double(v1 - mean(v1));
            [b,a] = butter(3, [self.const.BET_S self.const.BET_E]/fs);
            vf = filter(b,a,v1);
        end
        % filter Gamma
        function vf = gam(self, v1, sIdx)
            fs = self.h5eeg.fs{sIdx};
            v1 = double(v1 - mean(v1));
            [b,a] = butter(3, [self.const.GAM_S self.const.GAM_E]/fs);
            vf = filter(b,a,v1);
        end
        % Hilbert envelope
        function ve = env(self, v1)
            ve = abs(hilbert(v1));
        end
        % Metric short name to full name
        function mfull = metric_tofull(self,metric)
            switch(metric)
                case ('s')
                    mfull = 'Spearman Correlation';
                case ('p')
                    mfull = 'Pearson Correlation';
                case ('sP')
                    mfull = 'Spearman Partial Correlation';
                case ('pP')
                    mfull = 'Pearson Partial Correlation';
                case ('sd')
                    mfull = 'Delta Envelope Spearman Correlation';
                case ('st')
                    mfull = 'Theta Envelope Spearman Correlation';
                case ('sa')
                    mfull = 'Alpha Envelope Spearman Correlation';
                case ('sb')
                    mfull = 'Beta Envelope Spearman Correlation';
                case ('sg')
                    mfull = 'Gamma Envelope Spearman Correlation';
                case ('pc')
                    mfull = 'Coherence';
                case ('sc')
                    mfull = 'Spearman Coherence';
                case ('pcBroadband')
                    mfull = 'Coherence';
                case ('pcDelta')
                    mfull = 'Delta Coherence';
                case ('pcTheta')
                    mfull = 'Theta Coherence';
                case ('pcAlpha')
                    mfull = 'Alpha Coherence';
                case ('pcBeta')
                    mfull = 'Beta Coherence';
                case ('pcGamma')
                    mfull = 'Gamma Coherence';
                case ('scBroadband')
                    mfull = 'Spearman Coherence';
                case ('scDelta')
                    mfull = 'Delta Spearman Coherence';
                case ('scTheta')
                    mfull = 'Theta Spearman Coherence';
                case ('scAlpha')
                    mfull = 'Alpha Spearman Coherence';
                case ('scBeta')
                    mfull = 'Beta Spearman Coherence';
                case ('scGamma')
                    mfull = 'Gamma Spearman Coherence';
                otherwise
                    mfull = 'UNKNOWN';
            end
        end
    end
    
end

