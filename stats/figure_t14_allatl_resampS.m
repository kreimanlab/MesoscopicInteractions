% Figure W3


close all;
clear;
clc;
rng shuffle;

n_perm = 10000;
perm_alpha_sec = 20;
cp_thresh_override = 0.05; 
n_pairs_thresh = 10;% 10; % at least this many electrode pairs to be considered
n_subs_thresh = 2;% 2; % at least this many subjects to be considered
n_subs_ct_thresh = 0;% 2; % significant CTs in region pair must be from at least this many subjects
fig_size_scale = 0.75;
stride_axislabel = 1; % stride for x and y axis labels
%n_resample = 20; % number of times to resample with n_pairs_thresh pairs and n_subs_thresh subjects

% n_pairs_thresh = 10 * 4; % at least this many electrode pairs to be considered
% n_subs_thresh = 2; % at least this many subjects to be considered
% [*] Original ROI coverage: 78.0952 %
% [*] Thresholded ROI coverage: 61.4286 %

fig_fmt = '-dpng';
trig_eps = true;
trig_mag_no_cp = true;
trig_plot_mag_cov = true;
system('mkdir figures');
system('mkdir figures/T14_allatl_resampS');



[~,host] = system('hostname');

if contains(host,'ubuntu_1604')
    dir_artL = '/nas_share/RawData/data/h5_notch20/art';
    dir_resL = '/nas_share/RawData/data/results/coh_w10';
    dir_corL = '/nas_share/RawData/data/coreg';
    dir_cacheL = './cache';
    subjects_dirL = '/mnt/cuenap_ssd/coregistration';
    dir_h5L = '/nas_share/RawData/data/h5_notch20';
else
%     dir_artL = '/media/klab/internal/data/h5_notch20/art';
%     dir_resL = '/media/klab/internal/data/results/coh_w10';
%     dir_corL = '/media/klab/internal/data/coreg';
    dir_artL = '/media/klab/KLAB101/h5_notch20/art_nosz';
    dir_resL = '/media/klab/KLAB101/results/coh_w10';
    dir_corL = '/media/klab/internal/data/coreg';
    dir_cacheL = './cache';
    subjects_dirL = '/mnt/cuenap_ssd/coregistration';
    dir_h5L = '/media/klab/KLAB101/h5_notch20';
end

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
metrics_suffix = {'0.5-125 Hz','3-8 Hz','8-12 Hz','12-30 Hz','30-100 Hz'};
%metrics = {'pcBroadband'};

% Get atlas names
CaAtl = load(sprintf('./cache/xsub_out_all_%i',1));
AtlNames = CaAtl.C.AtlNames;
for atl = 2 %1:20

    system(sprintf('mkdir figures/T14_allatl_resampS/atl%i_%s',atl,AtlNames{atl}));
    
    for iM = [1 5] %1:length(metrics) % [1 5] %
        metric = metrics{iM};

        % Load human cache
        Ca_hum = load(sprintf('%s/xsub_out_all_%i_atl%i.mat',dir_cacheL,iM,atl));
        roi_labels = Ca_hum.C.AtlROIs{atl}.RH.struct_names;
        n_rois = Ca_hum.n_rois;
        %AdjAA = cell(n_rois,n_rois);
        
                
        % Vectorize matrix
        n_roi_pairs = nchoosek(n_rois,2);
        Vec_ij = zeros(n_roi_pairs,2);
        c = 1;
        for i1 = 1:(n_rois-1)
            for i2 = (i1+1):n_rois
                Vec_ij(c,1) = i1;
                Vec_ij(c,2) = i2;
                c = c + 1;
            end
        end
        
        
        % Unknown areas filter
        rois = Ca_hum.rois;
        cond_empty = false(length(rois),1);
        for ce = 1:length(cond_empty)
            elem = rois{ce};
            cond_empty(ce) = (all(isspace(elem)) | isempty(elem));
        end
        cond_contains = (contains(lower(rois),'unknown')) | (contains(lower(rois),'???')) | (contains(lower(rois),'wall')) | (cond_empty);
        known_idx = (~ cond_contains);
        
        
        % Permutation parameters <-- important
        %DropoutN = [0 8]; %0:length(Ca_hum.Subjects);
        DropoutN = 0;
        n_perm = 100000;
        
        
        
        % Drop out subjects
        n_sub = length(Ca_hum.Subjects);
        i_nSub = 1;
        AdjMag_Hi = nan(n_roi_pairs,length(DropoutN));
        AdjMag_Lo = nan(n_roi_pairs,length(DropoutN));
        for nSub = DropoutN
            
            % Set number of resampling permutations
            if (nSub == 0)
                n_perm_ac = 1;
            else
                n_perm_ac = n_perm;
            end
            
            fprintf('[n_perm: %i] atl: %i\tmetric: %i\t#drop: %i\tof %i\n',...
                n_perm_ac,atl,iM,nSub,max(DropoutN));
            
            AdjMag3 = nan(n_rois,n_rois,n_perm_ac);
            AdjMag3_1 = nan(n_rois,n_rois,n_perm_ac);
            AdjMag3_2 = nan(n_rois,n_rois,n_perm_ac);
            tic;
            
            for ip = 1:n_perm_ac
                
                % Pick nSub subjects to drop out by randperm
                sIdx = 1:n_sub;
                sIdx = sIdx(randperm(n_sub));
                sub_ex = sIdx(1:nSub);
                
                % Rebuild AdjMag excluding these subjects
                Adj = nan(n_rois,n_rois);
                AdjMag = nan(n_rois,n_rois);
                AdjMag_1 = nan(n_rois,n_rois);
                AdjMag_2 = nan(n_rois,n_rois);
                AdjMag4cl = nan(n_rois,n_rois);
                AdjCP = nan(n_rois,n_rois);
                AdjMagVar = nan(n_rois,n_rois);
                AdjMagReS = nan(n_rois,n_rois);
                AdjMagL = cell(n_rois,n_rois);
                N_bchan = nan(n_rois,n_rois);
                cp_thresh = cp_thresh_override;
                Dmat = Inf(n_rois,n_rois); % set to inf to avoid removing from every instance
                DistsAtl = [];
                
                elec_count = 0;
                for i1 = 1:n_rois
                    for i2 = (i1):n_rois
                        
                        cond_known = (known_idx(i1) && known_idx(i2));
                        
                        if (cond_known)
                        
                            AA_sub_all = Ca_hum.AdjAtl_sid{i1,i2};
                            AA_all = Ca_hum.AdjAtl{i1,i2};
                            AA_dist_all = Ca_hum.adjct_dist{i1,i2};


                            % Filter sub
                            i_trim = false(size(AA_sub_all));
                            for i3 = 1:length(i_trim)
                                if (mod(elec_count,2) == 0)
                                    i_trim(i3) = true;
                                else
                                    i_trim(i3) = false;
                                end
                                elec_count = elec_count + 1;
                                %i_trim = i_trim | (AA_sub == sub_ex(i3));
                            end

                            % apply filter
                            AA_sub = AA_sub_all(i_trim);
                            AA = AA_all(i_trim);
                            AA_dist = AA_dist_all(i_trim);

                            % apply complement filter
                            AA_sub_c = AA_sub_all(~i_trim);
                            AA_c = AA_all(~i_trim);
                            AA_dist_c = AA_dist_all(~i_trim);


                            n_pairs = length(AA_sub_all);
                            n_subs = length(unique(AA_sub_all));
                            n_subs_ct = length(unique(AA_sub_all(AA_all ~= 0)));
                            % ROI pair coverage condition (grey)
                            if ( (~ isempty(AA_all))  && (n_pairs >= n_pairs_thresh) && (n_subs >= n_subs_thresh))
                                N_bchan(i1,i2) = length(AA_all);
                                frac_cp = sum(AA_all ~= 0)/length(AA_all);
                                AdjCP(i1,i2) = frac_cp;
                                % ROI pair significance condition (white)
                                if ( (frac_cp > cp_thresh) && (n_subs_ct >= n_subs_ct_thresh) )
                                    %return;
                                    Adj(i1,i2) = 1;
                                    AdjMag(i1,i2) = mean(AA_all(AA_all ~= 0));

                                    % Separate calculations for pair resamp
                                    AdjMag_1(i1,i2) = mean(AA(AA ~= 0));
                                    AdjMag_2(i1,i2) = mean(AA_c(AA_c ~= 0));

                                    AdjMagVar(i1,i2) = var(AA_all(AA_all ~= 0));
                                    AdjMagL{i1,i2} = AA_all(AA_all~=0);
                                    DistsAtl = [DistsAtl; [mean(AA_dist_all(AA_all ~= 0)), mean(AA_all(AA_all ~= 0))]];
                                else
                                    Adj(i1,i2) = 0;
                                    AdjMag(i1,i2) = 0;
                                    AdjMag_1(i1,i2) = 0;
                                    AdjMag_2(i1,i2) = 0;
                                    AdjMagVar(i1,i2) = 0;
                                    AdjMagL{i1,i2} = [];
                                end
                            end
                            % AdjMag for clustering
                            AdjMag4cl(i1,i2) = mean(AA(AA~=0));
                        end
                    end
                end
                
                AdjMag3(:,:,ip) = AdjMag;
                AdjMag3_1(:,:,ip) = AdjMag_1;
                AdjMag3_2(:,:,ip) = AdjMag_2;
            end
            
            
            % 95%CI for AdjMag3
            CI_alpha = 0.05;
            AdjMag_hi = nan(n_rois,n_rois);
            AdjMag_lo = nan(n_rois,n_rois);
            for i1 = 1:n_rois
                for i2 = (i1):n_rois
                    mags = squeeze(AdjMag3(i1,i2,:));
                    s_mags = sort(mags,'descend');
                    sIdx = (~isnan(s_mags)) & (s_mags ~= 0);
                    s_mags = s_mags(sIdx);
                    
                    n_mags = length(s_mags);
                    idx_hi = round((CI_alpha/2) * n_mags);
                    idx_lo = round((1-(CI_alpha/2)) * n_mags);
                    
                    % Exception for no resampling (nSub = 0)
                    if (idx_hi == 0)
                        idx_hi = 1;
                    end
                    
                    if (isempty(s_mags))
                        hi = NaN;
                        lo = NaN;
                    else
                        hi = s_mags(idx_hi);
                        lo = s_mags(idx_lo);
                    end
                    
                    AdjMag_hi(i1,i2) = hi;
                    AdjMag_lo(i1,i2) = lo;
                    
%                     if (i1 == 4) && (i2 == 10)
%                         return
%                     end
                end
            end
            
            % Vectorize CI
            AdjMag_hi_v = nan(n_roi_pairs,1);
            AdjMag_lo_v = nan(n_roi_pairs,1);
            for c = 1:n_roi_pairs
                i1 = Vec_ij(c,1);
                i2 = Vec_ij(c,2);
                AdjMag_hi_v(c) = AdjMag_hi(i1,i2);
                AdjMag_lo_v(c) = AdjMag_lo(i1,i2);
            end
            
            % Store
            AdjMag_Hi(:,i_nSub) = AdjMag_hi_v;
            AdjMag_Lo(:,i_nSub) = AdjMag_lo_v;
            
            t_sing = toc;
            fprintf('per atl-metric ETA: %.1f mins\n',(t_sing * (length(DropoutN) - i_nSub))/(60))
            
%             if (i_nSub == 3)
%                 return
%             end
            
            i_nSub = i_nSub + 1;
        end
        
        save(sprintf('./cache/figure_t14_allatl_resampS_im-%i',iM));
        
        
        %%
        % Important variables:
        %
        % AdjMag_Hi
        % AdjMag_Lo
        % Vec_ij
        % DropoutN
        %
        %clear;
        %load('./cache/figure_t14_allatl_resamp');
        load(sprintf('./cache/figure_t14_allatl_resampS_im-%i',iM));
        
        
        
        
        
        % Filter
        AdjMag3_k = AdjMag3(known_idx,known_idx,1);
        AdjMag3_1_k = AdjMag3_1(known_idx,known_idx,1);
        AdjMag3_2_k = AdjMag3_2(known_idx,known_idx,1);
        
        % Vectorize
        n_known = sum(known_idx);
        Av_0 = zeros(nchoosek(n_known,2),1);
        Av_1 = zeros(nchoosek(n_known,2),1);
        Av_2 = zeros(nchoosek(n_known,2),1);

        c = 0;
        for i1 = 1:(n_known - 1)
            for i2 = (i1 + 1):n_known
                Av_0(c+1) = AdjMag3_k(i1,i2);
                Av_1(c+1) = AdjMag3_1_k(i1,i2);
                Av_2(c+1) = AdjMag3_2_k(i1,i2);
                c = c + 1;
            end
        end
        
        
        for iz = 1:length(AdjMag_Hi)
            i1 = Vec_ij(iz,1);
            i2 = Vec_ij(iz,2);
            
            if ( (~ known_idx(i1)) || (~ known_idx(i1)) )
                AdjMag_Hi(iz) = NaN;
                AdjMag_Lo(iz) = NaN;
            end
        end
        
        
        
        pct_hi = nanmean(100 * ((AdjMag_Hi(:,2:end) ./ AdjMag_Hi(:,1)) - 1));
        pct_lo = nanmean(100 * (1 - (AdjMag_Lo(:,2:end) ./ AdjMag_Lo(:,1))));
        pct_hi_s = nanstd(100 * ((AdjMag_Hi(:,2:end) ./ AdjMag_Hi(:,1)) - 1));
        pct_lo_s = nanstd(100 * (1 - (AdjMag_Lo(:,2:end) ./ AdjMag_Lo(:,1))));
        
        fprintf('[*] Number of permutations: %i\n',n_perm);
        for i = 1:length(pct_hi)
            fprintf('# subject drop: %i, 95%% CI Hi: %.2f%% (+- %.2f%%), Lo:%.2f%% (+- %.2f%%)\n',DropoutN(i+1),pct_hi(i),pct_hi_s(i),pct_lo(i),pct_lo_s(i));
            
            [p,H,stats] = ranksum(v_A,v_A1);
        end
        
        
        % Remove zeros and nans
        key = AdjMag_Hi(:,1);
        idx_sig = (~isnan(key)) & ((key) ~= 0);
        AdjMag_Hi = AdjMag_Hi(idx_sig,:);
        AdjMag_Lo = AdjMag_Lo(idx_sig,:);
        Vec_ij = Vec_ij(idx_sig,:);
        
        
        % Print information
        fprintf('[*] metric: %i, atl: %i\n',iM,atl);
        fprintf('\tBip-pair full n=%i\n',elec_count);
        fprintf('\tFull\tn_sig=%i\n',sum(Av_0>0));
        fprintf('\tSplit1\tn_sig=%i\n',sum(Av_1>0));
        fprintf('\tSplit2\tn_sig=%i\n',sum(Av_2>0));
        idx_full = (Av_0 > 0);
        diff = Av_1(idx_full) - Av_2(idx_full);
        [p,h,s] = ranksum(Av_1(idx_full),Av_2(idx_full));
        fprintf('\tRanksum\tp=%.4d\n',p);
        [h,p,s] = ttest(Av_1(idx_full),Av_2(idx_full),'tail','both');
        fprintf('\tT-test\tp=%.4d\n',p);
        fprintf('\tdiff mean: %.4f\n',nanmean(diff));
        fprintf('\tdiff std: %.4f\n',nanstd(diff));
        fprintf('\tdiff n: %i\n',sum(~isnan(diff)));
        

        close all;
        h = figure('Position',[0 0 600 400],'visible','off');
        histogram(diff,31,'Normalization','pdf','FaceColor',[1 1 1]*0.5);
        set(gca,'TickDir','Out');
        
        xlabel('Coherence Difference','FontSize',10);
        ylabel('pdf');
        title('Dataset split by bipolar pairs');
        
        axis tight;
        xlim(max(abs(diff)) * [-1 1]);
        grid off;
        box off;
        print(h,sprintf('figures/figure_t14_allatl_resampS_metric%i_atl%i_diff',iM,atl),'-depsc');
        close(h);
        
        %--------------------------------------------------------
        
        
        
        % Sort
        key = AdjMag_Hi(:,1);
        [~,idx_sort] = sort(key,'descend');
        AdjMag_Hi = AdjMag_Hi(idx_sort,:);
        AdjMag_Lo = AdjMag_Lo(idx_sort,:);
        Vec_ij = Vec_ij(idx_sort,:);
        
        fprintf('[*] Number of edges compared: %i\n',length(AdjMag_Hi(:,1)));
        
        % X-axis
        X = 1:length(AdjMag_Hi);
        
        close all;
        h = figure('Position',[0 0 800 200],'visible','off');
        hold all;
        
        % Colormap
        try
            cmap = inferno(length(DropoutN) - 1);
            cmap = [(0.5*[1 1 1]); cmap];
        catch
            cmap = [1 1 1]*0.5;
        end
        % Plot
        DropoutNp = DropoutN(1:end);
        for iR = 1:length(DropoutNp)
            i = length(DropoutNp) - iR + 1;
            if (i == 1)
                mk_size = 5;
            else
                mk_size = 2;
            end
            Y = AdjMag_Hi(:,i);
            plot(X(Y~=0),Y(Y~=0),'.','Color',cmap(i,:),'MarkerSize',mk_size);
            Y = AdjMag_Lo(:,i);
            plot(X(Y~=0),Y(Y~=0),'.','Color',cmap(i,:),'MarkerSize',mk_size);
        end
        
        % X labels
        xlab_txt = cell(length(X),1);
        for iL = 1:length(X)
            xlab_txt{iL} = sprintf('%2i,%2i',Vec_ij(iL,1),Vec_ij(iL,2));
        end
        xticks(X);
        xticklabels(xlab_txt);
        xtickangle(90);
        
        % Y labels
        coh_max = nanmax(AdjMag_Hi(:));
        coh_min = nanmin(AdjMag_Lo(:));
        yt_vals = linspace(coh_min,coh_max,5);
        yt_txt = cell(size(yt_vals));
        yticks(yt_vals);
        for iY = 1:length(yt_vals)
            yt_txt{iY} = sprintf('%6.2f',yt_vals(iY));
        end
        yticklabels(yt_txt);
        
        % Axis font sizes
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize',3);
        set(xAX,'FontName','FixedWidth'); %FixedWidth
        set(xAX,'TickLength',[0.001, 0.001]);
        yAX = get(gca,'YAxis');
        set(yAX,'TickLength',[0.005, 0.005]);
        
        % Plot format
        set(gca,'TickDir','Out');
        xlabel('Desikan-Killiany area pairs','FontSize',10);
        ylabel('Coherence')
        grid off;
        axis tight;
        
        print(h,sprintf('figures/figure_t14_allatl_resampS_metric%i_atl%i',iM,atl),'-depsc');
        close(h);

    end

end

% 18 fevr. 2020

% [n_perm: 1] atl: 2	metric: 1	#drop: 0	of 0
% per atl-metric ETA: 0.0 mins
% [*] Number of permutations: 100000
% [*] metric: 1, atl: 2
% 	Bip-pair full n=130820
% 	Full	n_sig=193
% 	Split1	n_sig=190
% 	Split2	n_sig=189
% 	Ranksum	p=7.7738e-01
% 	T-test	p=1.2880e-01
% 	diff mean: -0.0053
% 	diff std: 0.0476
% 	diff n: 186
% [*] Number of edges compared: 193
% [n_perm: 1] atl: 2	metric: 5	#drop: 0	of 0
% per atl-metric ETA: 0.0 mins
% [*] Number of permutations: 100000
% [*] metric: 5, atl: 2
% 	Bip-pair full n=130820
% 	Full	n_sig=183
% 	Split1	n_sig=178
% 	Split2	n_sig=179
% 	Ranksum	p=7.5480e-01
% 	T-test	p=8.5822e-01
% 	diff mean: 0.0008
% 	diff std: 0.0554
% 	diff n: 174
% [*] Number of edges compared: 183
