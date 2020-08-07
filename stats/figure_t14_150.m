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
stride_axislabel = 5; % stride for x and y axis labels
%n_resample = 20; % number of times to resample with n_pairs_thresh pairs and n_subs_thresh subjects

% n_pairs_thresh = 10 * 4; % at least this many electrode pairs to be considered
% n_subs_thresh = 2; % at least this many subjects to be considered
% [*] Original ROI coverage: 78.0952 %
% [*] Thresholded ROI coverage: 61.4286 %

fig_fmt = '-dpng';
trig_eps = true;
trig_mag_no_cp = true;
trig_plot_mag_cov = true;
trig_override_rois = false;
system('mkdir figures');
system('mkdir figures/T14_150');


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

for iM = 1 %1:length(metrics) % [1 5] %
    metric = metrics{iM};
    
    % Load human cache
    Ca_hum = load(sprintf('%s/xsub_out_all_%i_150.mat',dir_cacheL,iM));
    n_rois = Ca_hum.n_rois;
    
    % Calculate final functional interaction matrix
    Adj = nan(n_rois,n_rois);
    AdjMag = nan(n_rois,n_rois);
    AdjNpairs = nan(n_rois,n_rois);
    AdjNpairs_sig = nan(n_rois,n_rois);
    AdjNusubs = nan(n_rois,n_rois);
    AdjNusubs_sig = nan(n_rois,n_rois);
    AdjMag4cl = nan(n_rois,n_rois);
    AdjCP = nan(n_rois,n_rois);
    AdjMagVar = nan(n_rois,n_rois);
    AdjMagReS = nan(n_rois,n_rois);
    AdjMagL = cell(n_rois,n_rois);
    AdjMagEs2 = cell(n_rois,n_rois);
    AdjMagEs2B = cell(n_rois,n_rois);
    %AdjCPVar = nan(n_rois,n_rois);
    N_bchan = nan(n_rois,n_rois);
    cp_thresh = cp_thresh_override;
    Dmat = Inf(n_rois,n_rois); % set to inf to avoid removing from every instance
    DistsAtl = [];
    
    for i1 = 1:n_rois
        for i2 = 1:n_rois
            AA = Ca_hum.AdjAtl{i1,i2};
            AA_dist = Ca_hum.adjct_dist{i1,i2};
            AA_sub = Ca_hum.AdjAtl_sid{i1,i2};
            AA_chans = double(reshape(Ca_hum.AdjAtl_chans{i1,i2},2,[])' + 1);
            n_pairs = length(AA_sub);
            n_subs = length(unique(AA_sub));
            n_subs_ct = length(unique(AA_sub(AA ~= 0)));
            % ROI pair coverage condition (grey)
            if ( (~ isempty(AA))  && (n_pairs >= n_pairs_thresh) && (n_subs >= n_subs_thresh))
                N_bchan(i1,i2) = length(AA);
                frac_cp = sum(AA ~= 0)/length(AA);
                AdjCP(i1,i2) = frac_cp;
                % ROI pair significance condition (white)
                if ( (frac_cp > cp_thresh) && (n_subs_ct >= n_subs_ct_thresh) )
                    %return;
                    Adj(i1,i2) = 1;
                    AdjMag(i1,i2) = mean(AA(AA ~= 0));
                    AdjMagVar(i1,i2) = var(AA(AA ~= 0));
                    AdjMagL{i1,i2} = AA(AA~=0);
                    DistsAtl = [DistsAtl; [mean(AA_dist(AA ~= 0)), mean(AA(AA ~= 0))]];
                    AdjMagEs2{i1,i2} = AA_sub(AA~=0);
                    AdjMagEs2B{i1,i2} = AA_chans(AA~=0,:);
                    AdjNpairs(i1,i2) = n_pairs;
                    AdjNpairs_sig(i1,i2) = length(AA(AA ~= 0));
                    AdjNusubs(i1,i2) = n_subs;
                    AdjNusubs_sig(i1,i2) = n_subs_ct;
%                     % resample
%                     AAs = nan(n_resample,1);
%                     AAr = AA;
%                     for i3 = 1:n_resample
%                         AAr = AAr(randperm(length(AAr)));
%                         AA_t = AAr(1:n_pairs_thresh);
%                         AAs(i3) = mean(AA_t(AA_t ~= 0));
%                     end
%                     AdjMagReS(i1,i2) = nanmean(AAs);
                else
                    Adj(i1,i2) = 0;
                    AdjMag(i1,i2) = 0;
                    AdjMagVar(i1,i2) = 0;
                    AdjMagL{i1,i2} = [0];
                    AdjMagEs2{i1,i2} = [];
                    AdjMagEs2B{i1,i2} = [];
                end
            end
            
            AdjMag4cl(i1,i2) = mean(AA(AA~=0));
        end
    end
    
    % --- debug ---
%     adjmagl = nan(size(AdjMagL));
%     for i = 1:length(AdjMagL)
%         for j = 1:length(AdjMagL)
%             if (AdjMag(i,j) == 0)
%                 adjmagl(i,j) = 0;
%             else
%                 adjmagl(i,j) = mean(AdjMagL{i,j});
%             end
%         end
%     end
%     figure; imagesc(AdjMag);
%     figure; imagesc(adjmagl);
%     return
    % --- debug ---
    
    frac_xsub_msu = sum(Adj(~isnan(Adj)))/numel(Adj(~isnan(Adj)));
    fprintf('%s - human fraction of ROIs significant: %.4f\n',metric,frac_xsub_msu)
    
    % --- PLOT  --------------------------------------------------
    colordef_adj;
    color_nocov = COLOR_ADJ_NOCOV;
    color_distt = color_nocov;
    color_isnan = color_nocov;
    color_not_sig0 = COLOR_ADJ_NOSIG; %0.6*[1 1 1]; %[0 0.5 0]; %0.4*[0.5 1 0.5];
    fontsz = 10;
    fsize = fontsz;

    % Calculate rows of nans
%     Adj2 = Adj;
%     ind_isnan = false(length(Adj2),1);
%     for i = 1:length(ind_isnan)
%         ind_isnan(i) = (sum(isnan(Adj2(i,:))) == length(ind_isnan));
%     end
    
    % overall 36-area isnan storage
    ind_isnan_master = false(n_rois,1);
    for iii = 3%1:3

        %##################################################################
        % Make figure
        %##################################################################
        
        for iii2 = 1:2

            h = figure('visible','off');
            %set(h,'Position',round(fig_size_scale*[0 0 1*1080 1*1080]));
            set(h,'PaperUnits','Inches');
            set(h,'PaperPosition',[0 0 8.5 6.8]);

            if (iii == 1)
                if (trig_mag_no_cp)
                    Adj_plt = AdjMag;
                else
                    Adj_plt = AdjCP;
                end
            elseif (iii == 2)
                Adj_plt = Adj;
            elseif (iii == 3)
                if (trig_mag_no_cp)
                    if (iii2 == 1)
                        Adj_plt = AdjMag;
                    elseif (iii2 == 2)
                        Adj_plt = AdjMagVar;
                    end
                else
                    Adj_plt = AdjCP;
                end
                %Adj_plt(Adj == 0) = 0;
                color_not_sig = color_not_sig0;
                Adj_plt_var = AdjMagVar;
            end

            % ROI list is from:
            %   Ca_hum.rois <- xsub_out_all_1_150.mat <-
            %   all_parcellation_150.mat:LH.struct_names <- cluster_i <-
            %   cache/fig_cluster3_cluster_i.mat <- clustering old Adj 150
            Cai = load('cache/fig_cluster3_cluster_i.mat');
            
            rois = Ca_hum.rois;
            
            % --- DEBUG ---
%             fprintf('Check ROI names:\n')
%             for itmp = 1:length(rois)
%                 fprintf('\troi start: %i (int)\tend: %s (str)\n',Cai.cluster_i(itmp),rois{itmp});
%             end
            
            if (trig_override_rois)
                % False
                CaT14 = load(sprintf('./cache/figure_t14_%i_150',iM));
                
                % Apply clustering for deprecated method in fig_cluster3.m
                cluster_i_override = Cai.cluster_i;
                rois_sav = rois;
                rois = rois(cluster_i_override);
                fprintf('[*] ROI override (cluster_i ID):\n')
                for itmp = 1:length(rois)
                    fprintf('\t%s -> %s (%i)\n',rois_sav{itmp},rois{itmp},cluster_i_override(itmp));
                end
            end
            
            Dmat_plt = Dmat;
            rois_plt = rois;

            % checkpoint:
            % Adjacency matrix Adj_plt: 36 by 36
            % Adjacency rois_plt: 36 by 1
            % Distance matrix Dmat_plt: 

            Adj_plt2_cl = AdjMag4cl;
            Adj_plt_npairs = AdjNpairs;
            Adj_plt_npairs_sig = AdjNpairs_sig;
            Adj_plt_nusubs = AdjNusubs;
            Adj_plt_nusubs_sig = AdjNusubs_sig;
            
            % Filter out unknown from adjacency to plot
            known_idx = (~ strcmp(rois,'unknown'));
            ind_isnan_master(~known_idx) = true;
            Adj_plt = Adj_plt(known_idx,known_idx);
            Adj_plt_var = Adj_plt_var(known_idx,known_idx);
            Adj_plt_npairs = Adj_plt_npairs(known_idx,known_idx);
            Adj_plt_npairs_sig = Adj_plt_npairs_sig(known_idx,known_idx);
            Adj_plt_nusubs = Adj_plt_nusubs(known_idx,known_idx);
            Adj_plt_nusubs_sig = Adj_plt_nusubs_sig(known_idx,known_idx);
            rois_plt = rois_plt(known_idx);
            Dmat_plt = Dmat_plt(known_idx,known_idx);
            AdjMagEs2 = AdjMagEs2(known_idx,known_idx);
            AdjMagEs2B = AdjMagEs2B(known_idx,known_idx);
            AdjMagL = AdjMagL(known_idx,known_idx);
            Adj_plt2_cl= Adj_plt2_cl(known_idx,known_idx);

            for i = 1:length(rois_plt)
                rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '_','/');
                %rois_plt{i} = convertRoiDK(rois_plt{i});
            end



            % Distance threshold: add nans to Adj_plt
            %Adj_plt(Dmat_plt <= dist_thresh) = nanmean(Adj_plt(:));
            
            dist_thresh = Ca_hum.dist_thresh;
            Adj_plt(Dmat_plt <= dist_thresh) = nan;
            Adj_plt_var(Dmat_plt <= dist_thresh) = nan;
            Adj_plt_npairs(Dmat_plt <= dist_thresh) = nan;
            Adj_plt_npairs_sig(Dmat_plt <= dist_thresh) = nan;
            Adj_plt_nusubs(Dmat_plt <= dist_thresh) = nan;
            Adj_plt_nusubs_sig(Dmat_plt <= dist_thresh) = nan;
            Adj_plt2_cl(Dmat_plt <= dist_thresh) = nan;
            Adj_plt2 = Adj_plt;
            Adj_plt2_var = Adj_plt_var;
            Adj_plt2_npairs = Adj_plt_npairs;
            Adj_plt2_npairs_sig = Adj_plt_npairs_sig;
            Adj_plt2_nusubs = Adj_plt_nusubs;
            Adj_plt2_nusubs_sig = Adj_plt_nusubs_sig;
            
            % --- debug ---
%             adjmagl = nan(size(AdjMagL));
%             for i = 1:length(AdjMagL)
%                 for j = 1:length(AdjMagL)
%                     if (AdjMag(i,j) == 0)
%                         adjmagl(i,j) = 0;
%                     else
%                         adjmagl(i,j) = mean(AdjMagL{i,j});
%                     end
%                 end
%             end
%             figure; imagesc(AdjMag);
%             figure; imagesc(adjmagl);
%             return
            % --- debug ---
            
            %AdjMagEs2{Dmat_plt <= dist_thresh} = nan;
            %AdjMagL{Dmat_plt <= dist_thresh} = nan;

            %return
%             % Filter out nans nodes
%             cov_idx = ~ all(isnan(Adj_plt));
%             Adj_plt = Adj_plt(cov_idx,cov_idx);
%             rois_plt = rois_plt(cov_idx);
%             Dmat_plt = Dmat_plt(cov_idx,cov_idx);
%             ind_isnan_master(~cov_idx) = true;

            % Get number of nodes
            [n_roi_p,~] = size(Adj_plt);

            % Main setting for plot
            v = Adj_plt;
            [v_n, v_m] = size(v);
            %map = corrcmap(100);
            map = COLOR_ADJ_CMAP; %inferno(100);
            if (iii == 3)
                minv = min(v(v~=0));
                maxv = max(v(v~=0));
            else
                minv = min(v(:));
                maxv = max(v(:));
            end
            ncol = size(map,1);
            s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
            Im = ind2rgb(s,map);
            for i = 1:n_roi_p
                for j = 1:n_roi_p
                    if (isnan(Adj_plt(i,j)))
                        Im(i,j,:) = color_isnan;
                    end
                    if (Dmat_plt(i,j) <= dist_thresh)
                        Im(i,j,:) = color_distt;
                    end
                    if ((iii == 3) && (Adj_plt(i,j) == 0))
                        Im(i,j,:) = color_not_sig;
                    end
                end
            end

            %Im = Im(~ind_isnan,~ind_isnan,:); % remove rows of nans
            %rois_plt = rois_plt(~ind_isnan);
            %Adj_plt2 = Adj_plt(~ind_isnan,~ind_isnan);
            %rois_plt2 = rois_plt;
            %Adj_plt2 = Adj_plt;
            %rois_plt2 = rois_plt;




            % --- Clustering (output: cluster_i) --------------------------------------

            if ((iM == 1) && (iii2 == 1))
                %adjct_dist_cl = adjct_dist(ind_hum2mac,ind_hum2mac);

                % TODO: RERUN xsub_out_stats to enable adjct_dist clustering
                %adjct_dist = Ca_hum.adjct_dist;
                %adjct_dist_cl = adjct_dist;
                %adjct_dist_cl = adjct_dist_cl(~ind_isnan_master,~ind_isnan_master);
                %adjct_dist = adjct_dist(~ind_isnan)
                n_rois_cl = length(rois_plt);
                Y = zeros(1,nchoosek(n_rois_cl,2));
                yc = 1;
                for m1 = 1:(n_rois_cl-1)
                    for m2 = (m1+1):n_rois_cl

                        % ========= CLUSTER BY DISTANCE ===========================
                        %ad = adjct_dist_cl{m2,m1};
        %                 if (isempty(ad))
        %                     Y(yc) = Inf;%Inf; 600;
        %                 else
        %                     Y(yc) = mean(ad);
        %                 end
        %                 yc = yc + 1;

                        % ========= CLUSTER BY COHERENCE ==========================
                        ad = Adj_plt(m2,m1);
                        if (isnan(ad))
                            Y(yc) = 0;
                        else
                            Y(yc) = ad;
                        end
                        yc = yc + 1;


                    end
                end
                roi_dist = squareform(Y);
                % average   2778.832643576550 mm
                % centroid  2773.791496333409 mm
                % complete  2808.894068525726 mm
                % median    2767.246370358684 mm
                % single    2836.305543514634 mm
                % ward      2778.832643576550 mm
                % weighted  2778.918494966365 mm

%                 Z = linkage(Y,'complete'); %,'centroid'
                 % Imputed
                %Z = linkage(Y,'complete'); %,'centroid'
                
                % Marginalized
                M = Adj_plt2_cl;
                M = (M + M.')/2;

                Z = linkage(M,'ward',@nanDist4clustering);


                cluster_i = optimalleaforder(Z,Y); % ,'transformation','inverse'
                roi_dist = roi_dist(cluster_i,cluster_i);
                %roi_dist(isinf(roi_dist)) = ;
                clash = nansum(nansum(triu(roi_dist,1) - triu(roi_dist,2)));
                fprintf('cluster clash: %.12f mm\n',clash)
            %cluster_i = optimalleaforder(Z,Y);
            end

            % -------------------------------------------------------------------------

            
            % bypass clustering
            %cluster_i = 1:length(cluster_i);
            %Cai = load('cache/fig_cluster3_cluster_i.mat');
            %cluster_i = Cai.cluster_i;

            % Apply clustering
            %cluster_i = 1:length(cluster_i);
            Im = Im(cluster_i,cluster_i,:);
            Adj_plt2 = Adj_plt2(cluster_i,cluster_i);
            Adj_plt2_var = Adj_plt2_var(cluster_i,cluster_i);
            Adj_plt2_npairs = Adj_plt2_npairs(cluster_i,cluster_i);
            Adj_plt2_npairs_sig = Adj_plt2_npairs_sig(cluster_i,cluster_i);
            Adj_plt2_nusubs = Adj_plt2_nusubs(cluster_i,cluster_i);
            Adj_plt2_nusubs_sig = Adj_plt2_nusubs_sig(cluster_i,cluster_i);
            %rois_plt = rois_plt(1:length(cluster_i));
            %rois_plt_sav = rois_plt;
            rois_plt = rois_plt(cluster_i);
            AdjMagEs2 = AdjMagEs2(cluster_i,cluster_i);
            AdjMagEs2B = AdjMagEs2B(cluster_i,cluster_i);
            AdjMagL = AdjMagL(cluster_i,cluster_i);
            
            
            %rois_plt_sav = rois_plt;
            %rois_plt = rois_plt(cluster_i);
            
            %rois_plt_all_parcellation = rois_plt;
            %rois_plt = rois_plt_sav;

            % Manual adjustment
            n_wrap = 0; %round((0.2258) * n_rois_cl); % 9 - Translate the whole clustering and wrap around
            cluster_i_manual = [(n_rois_cl - (n_wrap-1)):n_rois_cl, 1:(n_rois_cl - (n_wrap))];
            Im = Im(cluster_i_manual,cluster_i_manual,:);
            rois_plt = rois_plt(cluster_i_manual);
            Adj_plt2 = Adj_plt2(cluster_i_manual,cluster_i_manual);
            Adj_plt2_var = Adj_plt2_var(cluster_i_manual,cluster_i_manual);
            Adj_plt2_npairs = Adj_plt2_npairs(cluster_i_manual,cluster_i_manual);
            Adj_plt2_npairs_sig = Adj_plt2_npairs_sig(cluster_i_manual,cluster_i_manual);
            Adj_plt2_nusubs = Adj_plt2_nusubs(cluster_i_manual,cluster_i_manual);
            Adj_plt2_nusubs_sig = Adj_plt2_nusubs_sig(cluster_i_manual,cluster_i_manual);
            AdjMagEs2 = AdjMagEs2(cluster_i_manual,cluster_i_manual);
            AdjMagEs2B = AdjMagEs2B(cluster_i_manual,cluster_i_manual);
            AdjMagL = AdjMagL(cluster_i_manual,cluster_i_manual);
            
            %return

            % ROIs renaming, keeping old name
            rois_plt_all_parcellation = rois_plt;
            rois_plt = cell(length(rois_plt),1);
            for itmp = 1:length(rois_plt)
                rois_plt{itmp} = sprintf('%i',itmp);
            end
            
            rois_plt_all_parcellation_cluster_i = zeros(1,length(rois_plt_all_parcellation));
            for irpa = 1:length(rois_plt_all_parcellation_cluster_i)
                rois_plt_all_parcellation_cluster_i(irpa) = find(strcmp(rois_plt_all_parcellation,sprintf('%i',irpa)));
            end
            
            % --- NOTE ----------------------------------------------------
            %
            % rois_plt_all_parcellation - roi names found in
            %   all_parcellation_150.mat files, which is the same as
            %   cluster_i found in:
            %       CaT14 = load(sprintf('./cache/figure_t14_%i_150',1));
            %
            % rois_plt - renamed roi names after clustering by function
            %
            %--------------------------------------------------------------
%             if (iii2 == 1)
%                 save(sprintf('./cache/figure_t14_%i_150',iM),'Adj_plt','Adj_plt2','cluster_i','rois_plt','rois_plt_all_parcellation');
%             end
            
            imagesc(Im);
            %imagesc(Im, 'Parent', ax);

            % 
            % % --- Axis label staggering---------------------
            % bip_labels = rois_plt;
            % stagger_labels = true;
            % xtick_vec = 1:length(bip_labels);
            % ytick_vec = 1:length(bip_labels);
            % ax1_xind = true(1,length(xtick_vec));
            % ax1_yind = true(1,length(ytick_vec));
            % ax2_xind = true(1,length(xtick_vec));
            % ax2_yind = true(1,length(ytick_vec));
            % if (stagger_labels)
            %     odd_i = downsample(1:length(ax1_xind),2,1);
            %     even_i = downsample(1:length(ax1_xind),2,0);
            %     ax1_xind(even_i) = false;
            %     ax2_xind(odd_i) = false;
            %     ax1_yind(odd_i) = false;
            %     ax2_yind(even_i) = false;
            % end
            % 
            % % Show axis labels
            % %ax1 = ax;
            % ax = gca;
            % set(ax,'xtick',xtick_vec(ax1_xind),'Ticklength',[0.001, 0.001],...
            %     'xticklabel',bip_labels(ax1_xind),'fontsize',fsize,'TickDir','out');
            % xtickangle(90);
            % set(ax,'ytick',ytick_vec(ax1_yind),'Ticklength',[0.001, 0.001],...
            %     'yticklabel',bip_labels(ax1_yind),'fontsize',fsize,'TickDir','out');
            % 
            % ax2 = copyobj(ax,ax.Parent);
            % set(ax2,'XAxisLocation','top');
            % set(ax2,'YAxisLocation','right');
            % set(ax2,'xtick',xtick_vec(ax2_xind),'xticklabel',bip_labels(ax2_xind));
            % set(ax2,'ytick',ytick_vec(ax2_yind),'yticklabel',bip_labels(ax2_yind));
            % return

            % --- END Axis label staggering---------------------

             % colormap
            % colormap(ax,map)
            % colormap(ax2,map)
            % if (minv == maxv)
            %     caxis(ax,[(minv - 1e-6), (maxv + 1e-6)]); 
            %     caxis(ax2,[(minv - 1e-6), (maxv + 1e-6)]); 
            % else
            %     caxis(ax,[minv maxv]);
            %     caxis(ax2,[minv maxv]);
            % end

            set(gca,'YTickMode','manual')
            set(gca,'XTickMode','manual')
            
            sub_idx = stride_axislabel:stride_axislabel:length(rois_plt);
            rty = 1:length(rois_plt);
            yticks(rty(sub_idx));
            yticklabels(rois_plt(sub_idx));
            rtx = 1:length(rois_plt);
            xticks(rtx(sub_idx));
            xticklabels(rois_plt(sub_idx));
            xtickangle(90);
            
%             yticks(1:length(rois_plt));
%             yticklabels(rois_plt);
%             xticks(1:length(rois_plt));
%             xticklabels(rois_plt);
%             xtickangle(90);
            
            set(gca,'tickdir','out');
            set(gca,'fontsize',fontsz); %fontsz*(31/n_rois_cl)^(0.7)
            set(gca,'TickLength',[0.001, 0.001])
            daspect([1 1 1]);
            set(gca,'FontName','Arial');
            
            % Move axis if needed
            ax = gca;
            tp = ax.Position;
            tp(1) = tp(1) - 0.05;
            ax.Position = tp;
    %         
    %         if (strcmp(metric,'pcBroadband'))
    %             title('Functional Interactions')
    %         else
    %             title(sprintf('Functional Interactions - %s',metric(3:end)))
    %         end

            %return

            if (iii == 1)
                if (trig_mag_no_cp)
                    ss = 'Mag';
                else
                    ss = 'CP';
                end
            elseif (iii == 2)
                ss = 'bin';
            elseif (iii == 3)
                if (trig_mag_no_cp)
                    ss = 'Mag-bin';
                    %ss = sprintf('%s_cov-',ss,round(1000*));
    %                 mag_min = min(AdjMag(AdjMag ~= 0 ));
    %                 mag_max = max(AdjMag(AdjMag ~= 0 ));
                else
                    ss = 'CP-bin';
                end
            end
            
            
            %set(gca,'FontName','Times New Roman');

    %         colormap(map);
    %         if (((~isnan(minv)) && (~isnan(maxv))) && (minv ~= maxv))
    %             ctick = linspace(minv,maxv,5);
    %             ctickT = cell(size(ctick));
    %             for ict = 1:length(ctickT)
    %                 ctickT{ict} = sprintf('%.2f',ctick(ict));
    %             end
    %             cb = colorbar('Ytick',ctick,'TickLabels',ctickT);
    %             set(cb,'TickLength',0);
    %             caxis([minv maxv]);
    %         else
    %             cb = colorbar;
    %             set(cb,'TickLength',0);
    %         end

            % New colormap handling
            %-----------------------------------------
            colormap(map);
            if (((~isnan(minv)) && (~isnan(maxv))) && (minv ~= maxv))
                cytick = linspace(minv,maxv,3);
                cytick_str = cell(1,length(cytick));
                for kk = 1:length(cytick)
                    cytick_str{kk} = sprintf('%.2f',cytick(kk));
                end
                cb = colorbar('ytick',cytick,'yticklabel',cytick_str,'FontSize',fontsz);
                caxis([minv maxv]);
            else
                cb = colorbar('ytick',minv);
            end
            set(cb,'TickDir','out');
            set(cb,'TickLength',0);
            
            %--- display frequency band range ----
            metricTxt = metric(3:end);
            if (strcmp(metricTxt,'Broadband'))
                if (iii2 == 1)
                    ylabel(cb,sprintf('Coherence'));
                elseif (iii2 == 2)
                    ylabel(cb,sprintf('Variance of Coherence'));
                end
            else
                if (iii2 == 1)
                    ylabel(cb,sprintf('%s Coherence (%s)',metric(3:end),metrics_suffix{iM}));
                elseif (iii2 == 2)
                    ylabel(cb,sprintf('Variance of %s Coherence (%s)',metric(3:end),metrics_suffix{iM}));
                end
            end
            %-------------------------------------
            
            colormap(map);
            cb.Location = 'eastoutside';
            cbpos = cb.Position;
            cbpos(1) = 0.82;
            cbpos(2) = cbpos(2) + 0.2; %2.6*cbpos(2);
            cbpos(3) = 0.02;
            cbpos(4) = 0.4; %0.6*cbpos(4);
            cb.Position = cbpos;
            

            annotation('rectangle',[cb.Position(1) cb.Position(2)-0.06 cb.Position(3) 0.02],'FaceColor',color_not_sig0);
            annotation('textbox',[cb.Position(1)+0.02 cb.Position(2)-0.06 0.2 0.02],'String','Not Significant','FitBoxToText','on','EdgeColor','none','VerticalAlignment','middle');
            annotation('rectangle',[cb.Position(1) cb.Position(2)-0.09 cb.Position(3) 0.02],'FaceColor',color_nocov);
            annotation('textbox',[cb.Position(1)+0.02 cb.Position(2)-0.09 0.2 0.02],'String','No Coverage','FitBoxToText','on','EdgeColor','none','VerticalAlignment','middle');
            %-------------------------------------------------


            %export_fig(sprintf('figures/T14/Adj_%s_%s_2',metric,ss),'-eps');
            if (iii2 == 2)
                ss = [ss,'_var'];
            end
            print(h,sprintf('figures/T14_150/Adj_%s_%s',metric,ss),fig_fmt);
            if (trig_eps)
                print(h,sprintf('figures/T14_150/Adj_%s_%s',metric,ss),'-depsc');
            end
            close(h);
            
            
            if (iii2 == 1)
                
                % Save for json export
                S = struct();
                S.Img = Im;
                S.A = Adj_plt2;
                S.Astd = sqrt(Adj_plt2_var);
                S.A_npairs = Adj_plt2_npairs;
                S.A_npairs_sig = Adj_plt2_npairs_sig;
                S.A_nusubs = Adj_plt2_nusubs;
                S.A_nusubs_sig = Adj_plt2_nusubs_sig;
                S.A_min = nanmin(Adj_plt2(:));
                S.A_max = nanmax(Adj_plt2(:));
                S.dendro_Z = Z;
                S.dendro_reorder = fliplr(cluster_i);
                S.labels = rois_plt;
                S.atl_name = 'Custom Parcellation';
                S.atl_chans = 1:length(Adj_plt2);
                save(sprintf('./brainexport/controls_data_sub-0_freq-%i_atl-0',iM),'S');
                
                % Adj_plt2 is the master matrix
                save(sprintf('./cache/figure_t14_%i_150',iM),'Adj_plt','Adj_plt2','Adj_plt2_cl','cluster_i','rois_plt','rois_plt_all_parcellation','rois_plt_all_parcellation_cluster_i');
            end
            
            
            
            
            
            % Dendrogram
            h = figure('visible','off');
            set(h,'PaperUnits','Inches');
            set(h,'PaperPosition',[0 0 1.5 5]);
    %                 cluster_i2 = 1:length(cluster_i);
    %                 cluster_i2 = cluster_i2(cluster_i);
    %                 cluster_i2 = cluster_i2(cluster_i_manual);
            hD = dendrogram(Z,0,'Reorder',fliplr(cluster_i),'Orientation','right','ColorThreshold',0); %[a,b,c]
            for ihD = 1:length(hD)
                hD(ihD).Color = 0*[1 1 1];
            end
            axis off;
            print(h,sprintf('figures/T14_150/Adj_%s_%s_dendro',metric,ss),fig_fmt);
            if (trig_eps)
                print(h,sprintf('figures/T14_150/Adj_%s_%s_dendro',metric,ss),'-depsc');
            end
            close(h);
        
        
        
            % Get info
            if ((iii2 == 1) && ((iM == 1) || (iM == 5)))
                rois_plt2 = rois_plt;
                
                % calculate fraction sig
                Adj_plt = Adj_plt2;
                n_sig2 = 0;
                n_all2 = 0;
                for isig1 = 1:(length(Adj_plt)-1)
                    for isig2 = (isig1+1):length(Adj_plt)
                        cond_count = ~isnan(Adj_plt(isig1,isig2));
                        if (cond_count)
                            if (Adj_plt(isig1,isig2) ~= 0)
                                n_sig2 = n_sig2 + 1;
                            end
                            n_all2 = n_all2 + 1;
                        end
                    end
                end
                n_nosig2 = n_all2 - n_sig2;
                fprintf('[!] Significant in : %i of %i (%.2f%%), total possible: %i\n',n_sig2,(n_sig2 + n_nosig2),100*n_sig2/(n_sig2 + n_nosig2),nchoosek(length(Adj_plt),2));
            
                
                Adiag = diag(Adj_plt2);
                fprintf('[!!] %i of %i in diagonal are nans.\n',sum(isnan(Adiag)),length(Adiag) );
                % 6 of 31 in diagonal are nan, 25 are covered
                AdiagC = Adiag(~isnan(Adiag));
                fprintf('[!!] %i of %i covered are significant\n',sum(AdiagC~=0),length(AdiagC~=0));
                % 17 of 25 are significant (68%)

                
                
                AM = Adj_plt2; %(cluster_i,cluster_i);
                AMvar = Adj_plt2_var; %(cluster_i,cluster_i);
                AM(AM==0) = NaN;
                [min_vals,min_idxs] = min(AM);
                [minv,mini] = min(min_vals);
                minv_std = sqrt(AMvar(mini,min_idxs(mini)));
                fprintf('[!!] minimum at: %i, %i, coherence = %.2f +- %.2d\n',mini,min_idxs(mini),minv,minv_std);
                fprintf('\t%s, %s\n',rois_plt2{mini},rois_plt2{min_idxs(mini)})

                AM = Adj_plt2; %(cluster_i,cluster_i);
                AM(AM==0) = NaN;
                [max_vals,max_idxs] = max(AM);
                [maxv,maxi] = max(max_vals);
                maxv_std = sqrt(AMvar(maxi,max_idxs(maxi)));
                fprintf('[!!] maximum at: %i, %i, coherence = %.2f +- %.4f\n',maxi,max_idxs(maxi),maxv,maxv_std);
                fprintf('\t%s, %s\n',rois_plt2{maxi},rois_plt2{max_idxs(maxi)})

                % min and max regions
                AM2 = Adj_plt2; %(cluster_i,cluster_i);
                [a,b] = min(sum(~isnan(AM)));
                fprintf('[!!] Min region: %s, total pairs sig: %.2f\n',rois_plt2{b},a);
                am = AM2(b,:);
                am = am(~isnan(am));
                fprintf('\t%i of %i sig\n',sum(am~=0),length(am));

                AM2 = Adj_plt2; %(cluster_i,cluster_i);
                [a,b] = max(sum(~isnan(AM)));
                fprintf('[!!] Max region: %s, total pairs sig: %.2f\n',rois_plt2{b},a);
                am = AM2(b,:);
                am = am(~isnan(am));
                fprintf('\t%i of %i sig\n',sum(am~=0),length(am));

                
%                 
%                 AM = Adj_plt2; %Adj_plt2(cluster_i,cluster_i);
%                 AM(AM==0) = NaN;
%                 [min_vals,min_idxs] = min(AM);
%                 [minv,mini] = min(min_vals);
%                 fprintf('[!!] minimum at: %i, %i, coherence = %.2f\n',mini,min_idxs(mini),minv);
%                 fprintf('\t%s, %s\n',rois_plt2{mini},rois_plt2{min_idxs(mini)})
% 
%                 AM = Adj_plt2; %Adj_plt2(cluster_i,cluster_i);
%                 AM(AM==0) = NaN;
%                 [max_vals,max_idxs] = max(AM);
%                 [maxv,maxi] = max(max_vals);
%                 fprintf('[!!] maximum at: %i, %i, coherence = %.2f\n',maxi,max_idxs(maxi),maxv);
%                 fprintf('\t%s, %s\n',rois_plt2{maxi},rois_plt2{max_idxs(maxi)})
% 
%                 AM2 = Adj_plt2; % Adj_plt2(cluster_i,cluster_i);
%                 [a,b] = max(sum(~isnan(AM)));
%                 fprintf('[!!] Max region: %s, total pairs sig: %.2f\n',rois_plt2{b},a);
%                 am = AM2(b,:);
%                 am = am(~isnan(am));
%                 fprintf('\t%i of %i sig\n',sum(am~=0),length(am));
                %return
            end
            

            % print info
            frac_nocov = sum(sum(isnan(Adj_plt2)))/numel(Adj_plt2);
            frac_nosig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 == 0)))/numel(Adj_plt2);
            frac_sig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 ~= 0)))/numel(Adj_plt2);
            fprintf('[*] saved file: %s\n',sprintf('figures/T14_150/Adj_%s_%s',metric,ss));
            fprintf('[*] no coverage: %.6f (%i of %i)\n',frac_nocov,sum(sum(isnan(Adj_plt2))),numel(Adj_plt2));
            n_nosig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 == 0)));
            fprintf('[*] not significant: %.6f (%i of %i)\n',frac_nosig,n_nosig,numel(Adj_plt2));
            n_sig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 ~= 0)));
            fprintf('[*] significant: %.6f (%i of %i)\n',frac_sig,n_sig,numel(Adj_plt2));
            fprintf('[*] significance fraction: %.6f (%i of %i)\n',100*frac_sig/(frac_sig + frac_nosig),n_sig,(n_sig + n_nosig));
            
            
            if (iii2 == 1)
                
                % Calculate number of subjects per 150 area for web interface
                % Load previous 150 cache (used by web interface)
                CaRloc = load(sprintf('%s/fig_cluster2_reduce_%i_new.mat','cache',iM));
                Es2 = cell(1,n_rois_cl);
                Es2_new = cell(1,n_rois_cl);
                for iA = 1:n_rois_cl
                    % For each region
                    atmp = AdjMagEs2(iA,:);
                    
                    % Subjects list
                    es2 = [];
                    for jA = 1:length(atmp)
                        es2 = [es2, atmp{jA}];
                    end
                    es2_sid_int = unique(es2);
                    for kA = 1:length(es2)
                        sstr = Ca_hum.Subjects{es2(kA)};
                        es2(kA) = str2num(sstr(2:end));
                    end
                    
                    % Bchan list
                    btmp = AdjMagEs2B(iA,:);
                    es2b = [];
                    for jA = 1:length(btmp)
                        bchans = btmp{jA}; %double(reshape(btmp{jA},2,[])'+1);
                        es2b = [es2b; bchans];
                    end
                    
                    % Get bchans in each subject
%                     roi_150 = rois_plt_all_parcellation{iA};
%                     for jA = 1:length(es2_sid_int)
%                         sstr = Ca_hum.Subjects{es2_sid_int(jA)};
%                         CaAP = load(sprintf('%s/%s/label/all_parcellation_150.mat',dir_corL,sstr));
%                         CaS = load(sprintf('./cache/xsub_out_%s_1.mat',sstr));
%                         sint = str2num(sstr(2:end));
%                         
% %                         % get bips
% %                         ubchans = unique(es2b(es2 == sint,:));
% %                         
% %                         % Convert raw ROI labels to bip labels
% %                         bchan_atllabels = cell(CaS.ecog.n_bchan,1);
% %                         for kA = 1:CaS.ecog.n_bchan
% %                             b1c1 = CaS.ecog.bip(kA,1);
% %                             bchan_atllabels{kA} = CaAP.AtlLabels{1}{b1c1};
% %                         end
%                         
%                         return
%                         %CaAP;
%                     end
                    
                    orig_Es2 = CaRloc.Es2{iA};
                    [ne2,~] = size(orig_Es2);
                    Es2_end = zeros(ne2,1);
                    if (~isempty(es2b))
                        for iB = 1:ne2
                            sidint = orig_Es2(iB,1);
                            bchan = orig_Es2(iB,2);

                            % Is this sid-bchan combo significant?
%                             for jB = 1:length(es2)
%                                 sidint_sig = es2(jB);
%                             end
                            cond_combo = (es2 == sidint) & ((es2b(:,1)' == bchan) | (es2b(:,2)' == bchan));
                            if any(cond_combo)
                                Es2_end(iB) = sum(cond_combo);
                            end
                        end
                    end
                    new_Es2 = orig_Es2;
                    new_Es2(:,end) = Es2_end;
                    Es2_new{iA} = new_Es2;
                    %return

                    Es2{iA} = [es2']; % es2b
                end
                
                Es2lens = nan(size(Es2));
                for i = 1:length(Es2)
                    Es2lens(i) = length(Es2{i});
                end

                fprintf('[!] Final check, all significant interactions were significant in at least:\n')
                fprintf('\t%i subjects.\n',min(Es2lens(Es2lens~=0)));
                
                
                % --- Final stats --------------------------------------
                % calculate fraction sig
                Adj_plt = Adj_plt2;
                n_sig2 = 0;
                n_all2 = 0;

                % Vectorized Adj_plt
                Adj_plt_vec = nan(nchoosek(length(Adj_plt),2),1);
                c_tmp = 1;

                for isig1 = 1:(length(Adj_plt)-1)
                    for isig2 = (isig1+1):length(Adj_plt)
                        cond_count = ~isnan(Adj_plt(isig1,isig2));
                        if (cond_count)
                            if (Adj_plt(isig1,isig2) ~= 0)
                                n_sig2 = n_sig2 + 1;
                            end
                            n_all2 = n_all2 + 1;
                        end

                        % Vectorized Adj_plt
                        Adj_plt_vec(c_tmp) = Adj_plt(isig1,isig2);
                        c_tmp = c_tmp + 1;
                    end
                end
                n_nosig2 = n_all2 - n_sig2;
                fprintf('[!] Significant in : %i of %i (%.2f%%), total possible (%i choose 2): %i\n',n_sig2,(n_sig2 + n_nosig2),100*n_sig2/(n_sig2 + n_nosig2),length(Adj_plt),nchoosek(length(Adj_plt),2));

                
                
                %return
                %save('./cache/figure_t14_150-Es2','Es2_new');
                % Convert to pair vector
%                 Es2_new = cell(1,nchoosek(n_rois_cl,2));
%                 c3 = 1;
%                 for c1 = 1:(n_rois_cl-1)
%                     for c2 = (c1+1):n_rois_cl
%                         Es2_new{c3} = AdjMagEs2{c1,c2};
%                         c3 = c3 + 1;
%                     end
%                 end
                %return
            end
        end
        %##################################################################
        % End mean/variance loop
        %return
                    
                    
                    
        
        % all unique subjects for region 1: unique(CaRloc.Es2{1}(:,end))
        %return
%         CaRloc.Es2{1}
% 
%         ans =
% 
%             1.0000   43.0000   19.6195  -21.3235  -25.3364   21.0000
%             1.0000   44.0000   19.8943  -17.1689  -28.9102   34.0000
%            19.0000   15.0000   18.3579  -15.3166  -26.0573    3.0000
%            19.0000   23.0000   21.4229  -22.6023  -27.3582    7.0000

        % Regions mapping onto original bip electrodes
        Subjectsp = {'m00005'};
        roi_1 = 'superiortemporal';
        roi_2 = 'parsopercularis';
        bchan1_const = 35;
        bchan2_const = 84;
        r_samp_const = 72056321;
        Cm5 = load(sprintf('%s/%s/label/all_parcellation_150.mat',dir_corL,Subjectsp{1}));
        Ca5 = load(sprintf('./cache/xsub_out_%s_%i.mat',Subjectsp{1},iM));
        
        % Load 150 area labels from coreg directory
        roi150s = Cm5.AtlLabels{1};
        
        % Get 150 areas of the two bips
        roi_1_150 = roi150s{Ca5.ecog.bip(bchan1_const,1)};
        roi_2_150 = roi150s{Ca5.ecog.bip(bchan2_const,1)};
        
        CaC3 = load('cache/fig_cluster3_cluster_i.mat');
        idx1 = find(strcmp(rois_plt_all_parcellation, roi_1_150));
        idx2 = find(strcmp(rois_plt_all_parcellation, roi_2_150));
        fprintf('Original %s:\n',Subjectsp{1});
        fprintf('[*] Bipolar electrode 1:\n\tDK: %s\t\n\t150: %i\n',roi_1,idx1);
        fprintf('[*] Bipolar electrode 2:\n\tDK: %s\t\n\t150: %i\n',roi_2,idx2);
        
        
        if ((iii == 3) && trig_plot_mag_cov)
            
            n_comb_atl = Ca_hum.n_comb_atl;
            Mag = nan(n_comb_atl,1);
            MagVar = nan(n_comb_atl,1);
            CP = nan(n_comb_atl,1);
            Npairs = nan(n_comb_atl,1);
            Nsubs = nan(n_comb_atl,1);
            %DistsAtl = [];%nan(n_comb_atl,1);
            
            
            count = 1;
            for i1 = 1:(n_rois-1)
                for i2 = (i1+1):n_rois
                    
                    AA = Ca_hum.AdjAtl_sid{i1,i2};
                    Mag(count) = AdjMag(i1,i2);
                    MagVar(count) = AdjMagVar(i1,i2);
                    CP(count) = AdjCP(i1,i2);
                    Npairs(count) = length(AA);
                    Nsubs(count) = length(unique(AA));
                    %DistsAtl(count) = mean(adjct_dist{i1,i2});
                    %AA_dist = adjct_dist{i1,i2};
                    
                    
                    count = count + 1;
                end
            end
            
            msize = 3;
            
            % Number of pairs vs magnitude
            h = figure('visible','off');
            set(h,'Position',round(fig_size_scale*[0 0 1080 1080]))
            set(h,'PaperUnits','inches')
            set(h,'PaperPosition',[0 0 4 4])
            is_sig_mag = ((Mag ~= 0)&(~isnan(Mag)));
            
            
            
            % Display coherence info
            fprintf('[*] Coherence mean: %.2f\n',nanmean(Mag(is_sig_mag)));
            fprintf('[*] Coherence std: %.2f\n',nanstd(Mag(is_sig_mag)));
            fprintf('[*] Coherence n: %i\n',numel((Mag(is_sig_mag))));

            %MagVar
            MagVar = sqrt(MagVar);
            fprintf('[*] Coherence Std mean: %.4f\n',nanmean(MagVar(is_sig_mag)));
            fprintf('[*] Coherence Std std: %.4f\n',nanstd(MagVar(is_sig_mag)));
            fprintf('[*] Coherence Std min: %.4f\n',nanmin(MagVar(is_sig_mag)));
            fprintf('[*] Coherence Std max: %.4f\n',nanmax(MagVar(is_sig_mag)));
            fprintf('[*] Coherence Std n: %i\n',numel((MagVar(is_sig_mag))));

            % Number pairs
            fprintf('[*] N pairs mean: %.2f\n',nanmean(Npairs(is_sig_mag)));
            fprintf('[*] N pairs median: %.2f\n',nanmedian(Npairs(is_sig_mag)));
            fprintf('[*] N pairs std: %.2f\n',nanstd(Npairs(is_sig_mag)));
            fprintf('[*] N pairs min: %.2f\n',nanmin(Npairs(is_sig_mag)));
            fprintf('[*] N pairs max: %.2f\n',nanmax(Npairs(is_sig_mag)));
            fprintf('[*] N pairs n: %i\n',numel((Npairs(is_sig_mag))));

            % Number subs
            fprintf('[*] N subs mean: %.2f\n',nanmean(Nsubs(is_sig_mag)));
            fprintf('[*] N subs median: %.2f\n',nanmedian(Nsubs(is_sig_mag)));
            fprintf('[*] N subs std: %.2f\n',nanstd(Nsubs(is_sig_mag)));
            fprintf('[*] N subs min: %.2f\n',nanmin(Nsubs(is_sig_mag)));
            fprintf('[*] N subs max: %.2f\n',nanmax(Nsubs(is_sig_mag)));
            fprintf('[*] N subs n: %i\n',numel((Nsubs(is_sig_mag))));
            
            
            
            subplot(2,1,1);
            plot(Npairs(is_sig_mag),Mag(is_sig_mag),'black.','MarkerSize',msize);
            xlabel('Number of bipolar pairs');
            ylabel(sprintf('%s Coherence',metric(3:end)));
            axis tight;
            box off;
            set(gca,'TickDir','out');
            set(gca,'FontSize',10);
            hold all;
            [corr_val_mag, corr_pval_mag] = corr(Npairs(is_sig_mag),Mag(is_sig_mag));
            fprintf('[%s] corr number of pairs, mag: %.5f (p = %.5f), n=%i\n',metric(3:end),corr_val_mag,corr_pval_mag,sum(is_sig_mag));
            
%             subplot(3,1,2)
%             [f,x] = hist(Npairs(is_sig_mag),50); plot(x,f,'color',[1 1 1]*0.5);
%             xlabel('Number of bipolar pairs');
%             ylabel('Number of occurences');
%             box off;
%             set(gca,'TickDir','out');
%             axis tight;
            
            subplot(2,1,2);
            plot(Npairs(is_sig_mag),CP(is_sig_mag),'black.','MarkerSize',msize);
            xlabel('Number of bipolar pairs');
            ylabel(sprintf('%s Consistency',metric(3:end)));
            axis tight;
            box off;
            set(gca,'TickDir','out');
            set(gca,'FontSize',10);
            [corr_val_cp, corr_pval_cp] = corr(Npairs(is_sig_mag),CP(is_sig_mag));
            fprintf('[%s] corr number of pairs, cp: %.5f (p = %.5f)\n',metric(3:end),corr_val_cp,corr_pval_cp);
            
            ofname = sprintf('figures/T14_150/Npairs_%s_mag_%i_p_%d_cp_%i_p_%d',metric,...
                round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp);
            ofname = replace(ofname,'.','p');
            print(h,ofname,'-depsc');
            
            close(h);
            
            % Number of unique patients vs magnitude
            h = figure('visible','off');
            set(h,'Position',round(fig_size_scale*[0 0 1080 1080]))
            set(h,'PaperUnits','inches')
            set(h,'PaperPosition',[0 0 4 4])
            is_sig_mag = ((Mag ~= 0)&(~isnan(Mag)));
            subplot(2,1,1);
            plot(Nsubs(is_sig_mag),Mag(is_sig_mag),'black.','MarkerSize',msize);
            xlabel('Number of Subjects');
            ylabel(sprintf('%s Coherence',metric(3:end)));
            axis tight;
            box off;
            set(gca,'TickDir','out');
            set(gca,'FontSize',10);
            hold all;
            [corr_val_mag, corr_pval_mag] = corr(Nsubs(is_sig_mag),Mag(is_sig_mag));
            fprintf('[%s] corr number of patients, mag: %.5f (p = %.5f), n=%i\n',metric(3:end),corr_val_mag,corr_pval_mag,sum(is_sig_mag));
            
%             subplot(3,1,2)
%             [f,x] = hist(Nsubs(is_sig_mag),1:48); plot(x,f,'color',[1 1 1]*0.5);
%             xlabel('Number of patients');
%             ylabel('Number of occurences');
%             box off;
%             set(gca,'TickDir','out');
%             axis tight;
            
            
            subplot(2,1,2);
            plot(Nsubs(is_sig_mag),CP(is_sig_mag),'black.','MarkerSize',msize);
            xlabel('Number of patients');
            ylabel(sprintf('%s Consistency',metric(3:end)));
            axis tight;
            box off;
            set(gca,'TickDir','out');
            set(gca,'FontSize',10);
            [corr_val_cp, corr_pval_cp] = corr(Nsubs(is_sig_mag),CP(is_sig_mag));
            fprintf('[%s] corr number of patients, cp: %.5f (p = %.5f)\n',metric(3:end),corr_val_cp,corr_pval_cp);
            ofname = sprintf('figures/T14_150/Nsubs_%s_mag_%i_p_%d_cp_%i_p_%d',metric,...
                round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp);
            
            ofname = replace(ofname,'.','p');
            print(h,ofname,'-depsc');
            
            close(h);
            
            
            % Print diagnostics
            f_roi_cov = sum((Npairs > 0) & (Nsubs > 0)) / length(Nsubs);
            f_roi_cov_t = sum((Npairs >= n_pairs_thresh) & (Nsubs >= n_subs_thresh)) / length(Nsubs);
            fprintf('[*] Original ROI coverage: %.4f %%\n',100*f_roi_cov);
            fprintf('[*] Thresholded ROI coverage: %.4f %%\n',100*f_roi_cov_t);
            %return
            
            
            % --- REPEAT for variance ---------------------------------
            
            % Number of pairs vs magnitude
            h = figure('visible','off');
            set(h,'Position',round(fig_size_scale*[0 0 0.95*1080 0.8*1080]))
            is_sig_mag = ((MagVar ~= 0)&(~isnan(MagVar)));
            subplot(3,1,1);
            plot(Npairs(is_sig_mag),MagVar(is_sig_mag),'black.');
            xlabel('Number of bipolar pairs');
            ylabel(sprintf('%s Coherence Variance',metric(3:end)));
            axis tight;
            box off;
            set(gca,'TickDir','out');
            hold all;
            [corr_val_mag, corr_pval_mag] = corr(Npairs(is_sig_mag),MagVar(is_sig_mag));
            fprintf('[%s] corr number of pairs, magVar: %.5f (p = %.5f)\n',metric(3:end),corr_val_mag,corr_pval_mag);
            
            subplot(3,1,2)
            [f,x] = hist(Npairs(is_sig_mag),50); plot(x,f,'color',[1 1 1]*0.5);
            xlabel('Number of bipolar pairs');
            ylabel('Number of occurences');
            box off;
            set(gca,'TickDir','out');
            axis tight;
            
            subplot(3,1,3);
            plot(Npairs(is_sig_mag),CP(is_sig_mag),'black.');
            xlabel('Number of bipolar pairs');
            ylabel(sprintf('%s Consistency',metric(3:end)));
            axis tight;
            box off;
            set(gca,'TickDir','out');
            [corr_val_cp, corr_pval_cp] = corr(Npairs(is_sig_mag),CP(is_sig_mag));
            fprintf('[%s] corr number of pairs, cp: %.5f (p = %.5f)\n',metric(3:end),corr_val_cp,corr_pval_cp);
            ofname = sprintf('figures/T14_150/Npairs_%s_magVar_%i_p_%d_cp_%i_p_%d',...
                metric,round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp);
            ofname = replace(ofname,'.','p');
            %return
            print(h,ofname,'-depsc');
            
            close(h);
            
            % Number of unique patients vs magnitude
            h = figure('visible','off');
            set(h,'Position',round(fig_size_scale*[0 0 0.95*1080 0.8*1080]))
            is_sig_mag = ((MagVar ~= 0)&(~isnan(MagVar)));
            subplot(3,1,1);
            plot(Nsubs(is_sig_mag),MagVar(is_sig_mag),'black.');
            xlabel('Number of patients');
            ylabel(sprintf('%s Coherence Variance',metric(3:end)));
            axis tight;
            box off;
            set(gca,'TickDir','out');
            hold all;
            [corr_val_mag, corr_pval_mag] = corr(Nsubs(is_sig_mag),MagVar(is_sig_mag));
            fprintf('[%s] corr number of patients, magVar: %.5f (p = %.5f)\n',metric(3:end),corr_val_mag,corr_pval_mag);
            
            subplot(3,1,2)
            [f,x] = hist(Nsubs(is_sig_mag),1:48); plot(x,f,'color',[1 1 1]*0.5);
            xlabel('Number of patients');
            ylabel('Number of occurences');
            box off;
            set(gca,'TickDir','out');
            axis tight;
            
            
            subplot(3,1,3);
            plot(Nsubs(is_sig_mag),CP(is_sig_mag),'black.');
            xlabel('Number of patients');
            ylabel(sprintf('%s Consistency',metric(3:end)));
            axis tight;
            box off;
            set(gca,'TickDir','out');
            [corr_val_cp, corr_pval_cp] = corr(Nsubs(is_sig_mag),CP(is_sig_mag));
            fprintf('[%s] corr number of patients, cp: %.5f (p = %.5f)\n',metric(3:end),corr_val_cp,corr_pval_cp);
            ofname = sprintf('figures/T14_150/Nsubs_%s_magVar_%i_p_%d_cp_%i_p_%d',metric,...
                round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp);
            ofname = replace(ofname,'.','p');
            print(h,ofname,'-depsc');
            
            close(h);
            
            % --- END REPEAT ------------------------------------------
            
            
            % signal dependent noise
            h = figure('visible','off');
            plot(Mag(is_sig_mag),MagVar(is_sig_mag),'black.')
            [corr_val_cp, corr_pval_cp] = corr(Mag(is_sig_mag),MagVar(is_sig_mag));
            xlabel(sprintf('%s Coherence Mean',metric(3:end)))
            ylabel(sprintf('%s Coherence Variance',metric(3:end)))
            box off;
            set(gca,'TickDir','out');
            fprintf('[%s] corr coherence mean, variance: %.5f (p = %.5f)\n',metric(3:end),corr_val_cp,corr_pval_cp);
            ofname = sprintf('figures/T14_150/mag_%s_magVar_%i_p_%d',metric,...
                round(1000*corr_val_cp),corr_pval_cp);
            ofname = replace(ofname,'.','p');
            print(h,ofname,'-depsc');
            
            close(h);
            % plot(Npairs(is_sig_mag),MagVar(is_sig_mag),'black.')
            % plot(Nsubs(is_sig_mag),MagVar(is_sig_mag),'black.')
            
            
            % distance and connection strength
            h = figure('visible','off');
            plot(DistsAtl(:,1),DistsAtl(:,2),'black.')
            [corr_val_cp, corr_pval_cp] = corr(DistsAtl(:,1),DistsAtl(:,2));
            xlabel('Distance (mm)')
            ylabel(sprintf('%s Coherence',metric(3:end)))
            box off;
            set(gca,'TickDir','out');
            fprintf('[%s] corr dist, mag: %.5f (p = %.5f)\n',metric(3:end),corr_val_cp,corr_pval_cp);
            axis tight;
            %return
            ofname = sprintf('figures/T14_150/dist_%s_mag_%i_p_%d',metric,...
                round(1000*corr_val_cp),corr_pval_cp);
            ofname = replace(ofname,'.','p');
            print(h,ofname,'-depsc');
            close(h);
            
            
            % Diagonal vs k-th diagonal
            DiagMag = diag(AdjMag,0);
            TriuMag = zeros(nchoosek(n_rois,2),1);
            c = 1;
            for ii = 1:(n_rois - 1)
                for jj = (ii + 1):n_rois
                    TriuMag(c) = AdjMag(ii,jj);
                    c = c + 1;
                end
            end
            % Remove no-coverage
            DiagMag = DiagMag(~isnan(DiagMag));
            TriuMag = TriuMag(~isnan(TriuMag));
            % Remove not-connected
            DiagMag = DiagMag(DiagMag ~= 0);
            TriuMag = TriuMag(TriuMag ~= 0);
            n_bins = 50;
            x_bins = linspace(0,1,n_bins+1);
            [n_diag,~] = hist(DiagMag,x_bins);
            [n_triu,~] = hist(TriuMag,x_bins);
            h = figure('visible','off');
            plot(x_bins,n_diag/trapz(x_bins,n_diag),'black-');
            hold all;
            plot(x_bins,n_triu/trapz(x_bins,n_triu),'black--');
            xlabel(sprintf('%s Coherence',metric(3:end)));
            ylabel('Probability Density');
            legend({'Intra-region','Inter-region'})
            try
                [p_rs,~,stats] = ranksum(DiagMag,TriuMag);
            catch
                p_rs = NaN;
                stats = NaN;
            end
            fprintf('\tmean non-diag - mean diag: %.6f\n',mean(TriuMag) - mean(DiagMag))
            fprintf('\tmedian non-diag - median diag: %.6f\n',median(TriuMag) - median(DiagMag))
            fprintf('[*] ranksum test:\n')
            fprintf('\tp-val: %.6d\n',p_rs)
            try
                [~,p_tt2] = ttest2(DiagMag,TriuMag);
            catch
                p_tt2 = NaN;
            end
            fprintf('[*] unpaired t test:\n')
            fprintf('\tp-val: %.6d\n',p_tt2)
            % p  = 0.0033, z-val = -2.9403, ranksum = 2574
            box off;
            set(gca,'TickDir','out');
            print(h,sprintf('figures/T14_150/intra_inter_%s_ranksum-%i_ttest-%i',metric,...
                round(1000*p_rs),round(1000*p_tt2)),'-depsc');
            close(h);
            
            %return
        end
    end
end

% 29 janvier, 2020

% mkdir: impossible de créer le répertoire «figures»: Le fichier existe
% mkdir: impossible de créer le répertoire «figures/T14_150»: Le fichier existe
% pcBroadband - human fraction of ROIs significant: 0.5140
% Warning: ward's linkage specified with non-Euclidean distance metric. 
% > In linkage (line 204)
%   In figure_t14_150 (line 385) 
% cluster clash: 6.476126172322 mm
% [!] Significant in : 2387 of 4644 (51.40%), total possible: 11175
% [!!] 149 of 150 in diagonal are nans.
% [!!] 1 of 1 covered are significant
% [!!] minimum at: 121, 143, coherence = 0.13
% 	121, 143
% [!!] maximum at: 137, 140, coherence = 0.79
% 	137, 140
% [!!] Max region: 19, total pairs sig: 82.00
% 	82 of 111 sig
% [*] saved file: figures/T14_150/Adj_pcBroadband_Mag-bin
% [*] no coverage: 0.587156 (13211 of 22500)
% [*] not significant: 0.200622 (4514 of 22500)
% [*] significant: 0.212222 (4775 of 22500)
% [*] significance fraction: 51.404888 (4775 of 9289)
% [!] Final check, all significant interactions were significant in at least:
% 	1 subjects.
% [!] Significant in : 2387 of 4644 (51.40%), total possible (150 choose 2): 11175
% [*] saved file: figures/T14_150/Adj_pcBroadband_Mag-bin_var
% [*] no coverage: 0.587156 (13211 of 22500)
% [*] not significant: 0.242711 (5461 of 22500)
% [*] significant: 0.170133 (3828 of 22500)
% [*] significance fraction: 41.210033 (3828 of 9289)
% Original m00005:
% [*] Bipolar electrode 1:
% 	DK: superiortemporal	
% 	150: 111
% [*] Bipolar electrode 2:
% 	DK: parsopercularis	
% 	150: 148
% [*] Coherence mean: 0.26
% [*] Coherence std: 0.12
% [*] Coherence n: 2387
% [*] Coherence Std mean: 0.0547
% [*] Coherence Std std: 0.0662
% [*] Coherence Std min: 0.0000
% [*] Coherence Std max: 0.3880
% [*] Coherence Std n: 2387
% [*] N pairs mean: 23.10
% [*] N pairs median: 18.00
% [*] N pairs std: 13.43
% [*] N pairs min: 10.00
% [*] N pairs max: 93.00
% [*] N pairs n: 2387
% [*] N subs mean: 7.36
% [*] N subs median: 6.00
% [*] N subs std: 3.87
% [*] N subs min: 2.00
% [*] N subs max: 24.00
% [*] N subs n: 2387
% [Broadband] corr number of pairs, mag: -0.09116 (p = 0.00001), n=2387
% [Broadband] corr number of pairs, cp: -0.01847 (p = 0.36711)
% [Broadband] corr number of patients, mag: 0.00513 (p = 0.80217), n=2387
% [Broadband] corr number of patients, cp: -0.12851 (p = 0.00000)
% [*] Original ROI coverage: 82.9978 %
% [*] Thresholded ROI coverage: 41.5570 %
% [Broadband] corr number of pairs, magVar: 0.10136 (p = 0.00001)
% [Broadband] corr number of pairs, cp: -0.17812 (p = 0.00000)
% [Broadband] corr number of patients, magVar: 0.24018 (p = 0.00000)
% [Broadband] corr number of patients, cp: -0.27647 (p = 0.00000)
% [Broadband] corr coherence mean, variance: 0.53516 (p = 0.00000)
% [Broadband] corr dist, mag: 0.14725 (p = 0.00000)
% 	mean non-diag - mean diag: 0.047244
% 	median non-diag - median diag: 0.005105
% [*] ranksum test:
% 	p-val: 9.410246e-01
% [*] unpaired t test:
% 	p-val: 6.857151e-01
% pcGamma - human fraction of ROIs significant: 0.4544
% [!] Significant in : 2110 of 4644 (45.43%), total possible: 11175
% [!!] 149 of 150 in diagonal are nans.
% [!!] 1 of 1 covered are significant
% [!!] minimum at: 10, 51, coherence = 0.11
% 	10, 51
% [!!] maximum at: 4, 149, coherence = 0.78
% 	4, 149
% [!!] Max region: 19, total pairs sig: 74.00
% 	74 of 111 sig
% [*] saved file: figures/T14_150/Adj_pcGamma_Mag-bin
% [*] no coverage: 0.587156 (13211 of 22500)
% [*] not significant: 0.225244 (5068 of 22500)
% [*] significant: 0.187600 (4221 of 22500)
% [*] significance fraction: 45.440844 (4221 of 9289)
% [!] Final check, all significant interactions were significant in at least:
% 	2 subjects.
% [!] Significant in : 2110 of 4644 (45.43%), total possible (150 choose 2): 11175
% [*] saved file: figures/T14_150/Adj_pcGamma_Mag-bin_var
% [*] no coverage: 0.587156 (13211 of 22500)
% [*] not significant: 0.266800 (6003 of 22500)
% [*] significant: 0.146044 (3286 of 22500)
% [*] significance fraction: 35.375175 (3286 of 9289)
% Original m00005:
% [*] Bipolar electrode 1:
% 	DK: superiortemporal	
% 	150: 111
% [*] Bipolar electrode 2:
% 	DK: parsopercularis	
% 	150: 148
% [*] Coherence mean: 0.27
% [*] Coherence std: 0.12
% [*] Coherence n: 2110
% [*] Coherence Std mean: 0.0554
% [*] Coherence Std std: 0.0691
% [*] Coherence Std min: 0.0000
% [*] Coherence Std max: 0.4724
% [*] Coherence Std n: 2110
% [*] N pairs mean: 23.16
% [*] N pairs median: 18.00
% [*] N pairs std: 13.69
% [*] N pairs min: 10.00
% [*] N pairs max: 93.00
% [*] N pairs n: 2110
% [*] N subs mean: 7.43
% [*] N subs median: 6.00
% [*] N subs std: 3.98
% [*] N subs min: 2.00
% [*] N subs max: 24.00
% [*] N subs n: 2110
% [Gamma] corr number of pairs, mag: -0.08859 (p = 0.00005), n=2110
% [Gamma] corr number of pairs, cp: -0.01849 (p = 0.39591)
% [Gamma] corr number of patients, mag: -0.00739 (p = 0.73435), n=2110
% [Gamma] corr number of patients, cp: -0.14109 (p = 0.00000)
% [*] Original ROI coverage: 82.9978 %
% [*] Thresholded ROI coverage: 41.5570 %
% [Gamma] corr number of pairs, magVar: 0.08793 (p = 0.00036)
% [Gamma] corr number of pairs, cp: -0.18286 (p = 0.00000)
% [Gamma] corr number of patients, magVar: 0.20778 (p = 0.00000)
% [Gamma] corr number of patients, cp: -0.29834 (p = 0.00000)
% [Gamma] corr coherence mean, variance: 0.55439 (p = 0.00000)
% [Gamma] corr dist, mag: 0.16483 (p = 0.00000)
% 	mean non-diag - mean diag: 0.106526
% 	median non-diag - median diag: 0.061703
% [*] ranksum test:
% 	p-val: 1.386115e-01
% [*] unpaired t test:
% 	p-val: 3.630462e-01

