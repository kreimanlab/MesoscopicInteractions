close all;
clear;
clc;
rng shuffle;

n_perm = 10000;
perm_alpha_sec = 20;
cp_thresh_override = 0.05;
n_pairs_thresh = 10; % at least this many electrode pairs to be considered
n_subs_thresh = 2; % at least this many subjects to be considered
n_subs_ct_thresh = 2; % significant CTs in region pair must be from at least this many subjects
fig_size_scale = 0.75;
%n_resample = 20; % number of times to resample with n_pairs_thresh pairs and n_subs_thresh subjects

% n_pairs_thresh = 10 * 4; % at least this many electrode pairs to be considered
% n_subs_thresh = 2; % at least this many subjects to be considered
% [*] Original ROI coverage: 78.0952 %
% [*] Thresholded ROI coverage: 61.4286 %

fig_fmt = '-depsc';
trig_eps = false;
trig_mag_no_cp = true;
trig_plot_mag_cov = true;
system('mkdir figures');
system('mkdir figures/T14');


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
%metrics = {'pcBroadband'};

for iM = [1:5] %1:length(metrics) % [1 5] %
    metric = metrics{iM};
    
    % Load human cache
    Ca_hum = load(sprintf('%s/xsub_out_all_%i.mat',dir_cacheL,iM));
    n_rois = Ca_hum.n_rois;
    
    % Calculate final functional interaction matrix
    Adj = nan(n_rois,n_rois);
    AdjMag = nan(n_rois,n_rois);
    AdjCP = nan(n_rois,n_rois);
    AdjMagVar = nan(n_rois,n_rois);
    AdjMagReS = nan(n_rois,n_rois);
    AdjMagL = cell(n_rois,n_rois);
    %AdjCPVar = nan(n_rois,n_rois);
    Adj_dist = inf(n_rois,n_rois);
    N_bchan = nan(n_rois,n_rois);
    cp_thresh = cp_thresh_override;
    Dmat = Inf(n_rois,n_rois); % set to inf to avoid removing from every instance
    DistsAtl = [];
    
    for i1 = 1:n_rois
        for i2 = 1:n_rois
            AA = Ca_hum.AdjAtl{i1,i2};
            AA_dist = Ca_hum.adjct_dist{i1,i2};
            AA_sub = Ca_hum.AdjAtl_sid{i1,i2};
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
                    %Adj_dist(i1,i2) = mean(AA_dist(AA ~= 0));
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
                    AdjMagL{i1,i2} = [];
                   % DistsAtl = [DistsAtl; mean(AA_dist)];
                end
                %DistsAtl = [DistsAtl; [mean(AA_dist(AA ~= 0)), mean(AA(AA ~= 0))]];
            else
                %DistsAtl = [DistsAtl; Inf];
            end
            if (~isempty(AA))
                Adj_dist(i1,i2) = mean(AA_dist);
            end
        end
    end

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

        % Make figure
        h = figure;
        set(h,'Position',round(fig_size_scale*[0 0 0.95*1080 0.8*1080]))

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
                Adj_plt = AdjMag;
            else
                Adj_plt = AdjCP;
            end
            %Adj_plt(Adj == 0) = 0;
            color_not_sig = color_not_sig0;
        end
        
        rois = Ca_hum.rois;
        Dmat_plt = Dmat;
        rois_plt = rois;
        
        % checkpoint:
        % Adjacency matrix Adj_plt: 36 by 36
        % Adjacency rois_plt: 36 by 1
        % Distance matrix Dmat_plt: 
        
        % Filter out unknown from adjacency to plot
        known_idx = (~ strcmp(rois,'unknown'));
        ind_isnan_master(~known_idx) = true;
        Adj_plt = Adj_plt(known_idx,known_idx);
        rois_plt = rois_plt(known_idx);
        Dmat_plt = Dmat_plt(known_idx,known_idx);
        Adj_dist = Adj_dist(known_idx,known_idx);
        %DistsAtl = DistsAtl(known_idx,known_idx);
        
        for i = 1:length(rois_plt)
            rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '_','/');
            rois_plt{i} = convertRoiDK(rois_plt{i});
        end
        
        
        
        % Distance threshold: add nans to Adj_plt
        %Adj_plt(Dmat_plt <= dist_thresh) = nanmean(Adj_plt(:));
        dist_thresh = Ca_hum.dist_thresh;
        Adj_plt2 = Adj_plt;
        Adj_plt(Dmat_plt <= dist_thresh) = nan;
        
        %return
        % Filter out nans nodes
        cov_idx = ~ all(isnan(Adj_plt));
        Adj_plt = Adj_plt(cov_idx,cov_idx);
        Adj_plt2 = Adj_plt2(cov_idx,cov_idx);
        rois_plt = rois_plt(cov_idx);
        Dmat_plt = Dmat_plt(cov_idx,cov_idx);
        Adj_dist = Adj_dist(cov_idx,cov_idx);
        %DistsAtl = DistsAtl(cov_idx,cov_idx);
        ind_isnan_master(~cov_idx) = true;
        
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
        
        %rois_plt2 = rois_plt;
        
        
        
        
        % --- Clustering (output: cluster_i) --------------------------------------
        
        if (iM == 1)
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

            Z = linkage(Y,'average'); %,'centroid'
            cluster_i = optimalleaforder(Z,Y); % ,'transformation','inverse'
            roi_dist = roi_dist(cluster_i,cluster_i);
            %roi_dist(isinf(roi_dist)) = ;
            clash = nansum(nansum(triu(roi_dist,1) - triu(roi_dist,2)));
            fprintf('cluster clash: %.12f mm\n',clash)
            
        %cluster_i = optimalleaforder(Z,Y);
        end

        % -------------------------------------------------------------------------

        save(sprintf('./cache/figure_t14_%i',iM),'Adj_plt','Adj_plt2','cluster_i','rois_plt','Adj_dist');
        
        % Apply clustering
        Im = Im(cluster_i,cluster_i,:);
        rois_plt = rois_plt(cluster_i);
        
        % Manual adjustment
        n_wrap = 0; % 9 - Translate the whole clustering and wrap around
        cluster_i_manual = [(n_rois_cl - (n_wrap-1)):n_rois_cl, 1:(n_rois_cl - (n_wrap))];
        Im = Im(cluster_i_manual,cluster_i_manual,:);
        rois_plt = rois_plt(cluster_i_manual);
        
        
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
        yticks(1:length(rois_plt));
        yticklabels(rois_plt);
        xticks(1:length(rois_plt));
        xticklabels(rois_plt);
        xtickangle(90);
        set(gca,'tickdir','out');
        set(gca,'fontsize',fontsz);
        set(gca,'TickLength',[0.001, 0.001])
        if (strcmp(metric,'pcBroadband'))
            title('Functional Interactions')
        else
            title(sprintf('Functional Interactions - %s',metric(3:end)))
        end
        
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
        
        colormap(map);
        if (((~isnan(minv)) && (~isnan(maxv))) && (minv ~= maxv))
            ctick = linspace(minv,maxv,5);
            ctickT = cell(size(ctick));
            for ict = 1:length(ctickT)
                ctickT{ict} = sprintf('%.2f',ctick(ict));
            end
            cb = colorbar('Ytick',ctick,'TickLabels',ctickT);
            set(cb,'TickLength',0);
            caxis([minv maxv]);
        else
            cb = colorbar;
            set(cb,'TickLength',0);
        end

        %export_fig(sprintf('figures/T14/Adj_%s_%s_2',metric,ss),'-eps');
        print(h,sprintf('figures/T14/Adj_%s_%s',metric,ss),fig_fmt);
        if (trig_eps)
            print(h,sprintf('figures/T14/Adj_%s_%s',metric,ss),'-depsc');
        end
        close(h);
        
        % print info
        %Adj_plt = Adj_plt2;
        frac_nocov = sum(sum(isnan(Adj_plt2)))/numel(Adj_plt2);
        frac_nosig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 == 0)))/numel(Adj_plt2);
        frac_sig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 ~= 0)))/numel(Adj_plt2);
        fprintf('[*] saved file: %s\n',sprintf('figures/T14/Adj_%s_%s',metric,ss));
        fprintf('[*] no coverage: %.6f (%i of %i)\n',frac_nocov,sum(sum(isnan(Adj_plt2))),numel(Adj_plt2));
        n_nosig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 == 0)));
        fprintf('[*] not significant: %.6f (%i of %i)\n',frac_nosig,n_nosig,numel(Adj_plt2));
        n_sig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 ~= 0)));
        fprintf('[*] significant: %.6f (%i of %i)\n',frac_sig,n_sig,numel(Adj_plt2));
        fprintf('[*] significance fraction: %.6f (%i of %i)\n',100*frac_sig/(frac_sig + frac_nosig),n_sig,(n_sig + n_nosig));

        
        
        
        % Dendrogram
        h = figure;
        set(h,'Position',round(fig_size_scale*[0 0 0.95*1080 0.8*1080]));
        cluster_i2 = 1:length(cluster_i);
        cluster_i2 = cluster_i2(cluster_i);
        cluster_i2 = cluster_i2(cluster_i_manual);
        hD = dendrogram(Z,0,'Reorder',cluster_i2,'Orientation','right','ColorThreshold',0); %[a,b,c]
        for ihD = 1:length(hD)
            hD(ihD).Color = 0*[1 1 1];
        end
        axis off;
        print(h,sprintf('figures/T14/Adj_%s_%s_dendro',metric,ss),fig_fmt);
        if (trig_eps)
            print(h,sprintf('figures/T14/Adj_%s_%s_dendro',metric,ss),'-depsc');
        end
        %return
        close(h);
        
        %return
        
        %##################################################################
        % End mean/variance loop
        
        
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
            h = figure;
            set(h,'Position',round(fig_size_scale*[0 0 1080 1080]))
            set(h,'PaperUnits','inches')
            set(h,'PaperPosition',[0 0 4 4])
            is_sig_mag = ((Mag ~= 0)&(~isnan(Mag)));
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
            
            ofname = sprintf('figures/T14/Npairs_%s_mag_%i_p_%d_cp_%i_p_%d',metric,...
                round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp);
            ofname = replace(ofname,'.','p');
            print(h,ofname,'-depsc');
            
            close(h);
            
            % Number of unique patients vs magnitude
            h = figure;
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
            ofname = sprintf('figures/T14/Nsubs_%s_mag_%i_p_%d_cp_%i_p_%d',metric,...
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
            h = figure;
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
            ofname = sprintf('figures/T14/Npairs_%s_magVar_%i_p_%d_cp_%i_p_%d',...
                metric,round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp);
            ofname = replace(ofname,'.','p');
            %return
            print(h,ofname,'-depsc');
            
            close(h);
            
            % Number of unique patients vs magnitude
            h = figure;
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
            ofname = sprintf('figures/T14/Nsubs_%s_magVar_%i_p_%d_cp_%i_p_%d',metric,...
                round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp);
            ofname = replace(ofname,'.','p');
            print(h,ofname,'-depsc');
            
            close(h);
            
            % --- END REPEAT ------------------------------------------
            
            
            % signal dependent noise
            h = figure;
            plot(Mag(is_sig_mag),MagVar(is_sig_mag),'black.')
            [corr_val_cp, corr_pval_cp] = corr(Mag(is_sig_mag),MagVar(is_sig_mag));
            xlabel(sprintf('%s Coherence Mean',metric(3:end)))
            ylabel(sprintf('%s Coherence Variance',metric(3:end)))
            box off;
            set(gca,'TickDir','out');
            fprintf('[%s] corr coherence mean, variance: %.5f (p = %.5f)\n',metric(3:end),corr_val_cp,corr_pval_cp);
            ofname = sprintf('figures/T14/mag_%s_magVar_%i_p_%d',metric,...
                round(1000*corr_val_cp),corr_pval_cp);
            ofname = replace(ofname,'.','p');
            print(h,ofname,'-depsc');
            
            close(h);
            % plot(Npairs(is_sig_mag),MagVar(is_sig_mag),'black.')
            % plot(Nsubs(is_sig_mag),MagVar(is_sig_mag),'black.')
            
            
            % distance and connection strength
            h = figure;
            plot(DistsAtl(:,1),DistsAtl(:,2),'black.')
            [corr_val_cp, corr_pval_cp] = corr(DistsAtl(:,1),DistsAtl(:,2));
            xlabel('Distance (mm)')
            ylabel(sprintf('%s Coherence',metric(3:end)))
            box off;
            set(gca,'TickDir','out');
            fprintf('[%s] corr dist, mag: %.5f (p = %.5f)\n',metric(3:end),corr_val_cp,corr_pval_cp);
            axis tight;
            %return
            ofname = sprintf('figures/T14/dist_%s_mag_%i_p_%d',metric,...
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
            h = figure;
            plot(x_bins,n_diag/trapz(x_bins,n_diag),'black-');
            hold all;
            plot(x_bins,n_triu/trapz(x_bins,n_triu),'black--');
            xlabel(sprintf('%s Coherence',metric(3:end)));
            ylabel('Probability Density');
            legend({'Intra-region','Inter-region'})
            [p_rs,~,stats] = ranksum(DiagMag,TriuMag);
            fprintf('\tmean non-diag - mean diag: %.6f\n',mean(TriuMag) - mean(DiagMag))
            fprintf('\tmedian non-diag - median diag: %.6f\n',median(TriuMag) - median(DiagMag))
            fprintf('[*] ranksum test:\n')
            fprintf('\tp-val: %.6d\n',p_rs)
            [~,p_tt2] = ttest2(DiagMag,TriuMag);
            fprintf('[*] unpaired t test:\n')
            fprintf('\tp-val: %.6d\n',p_tt2)
            % p  = 0.0033, z-val = -2.9403, ranksum = 2574
            box off;
            set(gca,'TickDir','out');
            print(h,sprintf('figures/T14/intra_inter_%s_ranksum-%i_ttest-%i',metric,...
                round(1000*p_rs),round(1000*p_tt2)),'-depsc');
            close(h);
            
            %return
        end
        
        
        
        
        
    end

    
end


% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T14’: File exists
% pcBroadband - human fraction of ROIs significant: 0.4356
% cluster clash: 1.242633921848 mm
% [*] saved file: figures/T14/Adj_pcBroadband_Mag-bin
% [*] no coverage: 0.187305 (180 of 961)
% [*] not significant: 0.453694 (436 of 961)
% [*] significant: 0.359001 (345 of 961)
% [*] significance fraction: 44.174136 (345 of 781)
% [Broadband] corr number of pairs, mag: -0.17659 (p = 0.02163), n=169
% [Broadband] corr number of pairs, cp: -0.12822 (p = 0.09664)
% [Broadband] corr number of patients, mag: 0.12866 (p = 0.09550), n=169
% [Broadband] corr number of patients, cp: -0.21987 (p = 0.00407)
% [*] Original ROI coverage: 77.7778 %
% [*] Thresholded ROI coverage: 62.6984 %
% [Broadband] corr number of pairs, magVar: 0.11514 (p = 0.13605)
% [Broadband] corr number of pairs, cp: -0.12822 (p = 0.09664)
% [Broadband] corr number of patients, magVar: 0.40485 (p = 0.00000)
% [Broadband] corr number of patients, cp: -0.21987 (p = 0.00407)
% [Broadband] corr coherence mean, variance: 0.63511 (p = 0.00000)
% [Broadband] corr dist, mag: 0.29957 (p = 0.00000)
% 	mean non-diag - mean diag: 0.053433
% 	median non-diag - median diag: 0.064875
% [*] ranksum test:
% 	p-val: 7.062741e-03
% [*] unpaired t test:
% 	p-val: 1.101700e-02
% pcGamma - human fraction of ROIs significant: 0.4012
% [*] saved file: figures/T14/Adj_pcGamma_Mag-bin
% [*] no coverage: 0.187305 (180 of 961)
% [*] not significant: 0.478668 (460 of 961)
% [*] significant: 0.334027 (321 of 961)
% [*] significance fraction: 41.101152 (321 of 781)
% [Gamma] corr number of pairs, mag: -0.06865 (p = 0.39758), n=154
% [Gamma] corr number of pairs, cp: -0.15657 (p = 0.05248)
% [Gamma] corr number of patients, mag: 0.25520 (p = 0.00140), n=154
% [Gamma] corr number of patients, cp: -0.23034 (p = 0.00405)
% [*] Original ROI coverage: 77.7778 %
% [*] Thresholded ROI coverage: 62.6984 %
% [Gamma] corr number of pairs, magVar: 0.14729 (p = 0.06833)
% [Gamma] corr number of pairs, cp: -0.15657 (p = 0.05248)
% [Gamma] corr number of patients, magVar: 0.44706 (p = 0.00000)
% [Gamma] corr number of patients, cp: -0.23034 (p = 0.00405)
% [Gamma] corr coherence mean, variance: 0.71115 (p = 0.00000)
% [Gamma] corr dist, mag: 0.45336 (p = 0.00000)
% 	mean non-diag - mean diag: 0.060755
% 	median non-diag - median diag: 0.073388
% [*] ranksum test:
% 	p-val: 6.036360e-04
% [*] unpaired t test:
% 	p-val: 1.195240e-03




% ========== Run 2019 06 04 ===============================================
% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T14’: File exists
% pcBroadband - human fraction of ROIs significant: 0.4356
% cluster clash: 1.242633921848 mm
% [*] saved file: figures/T14/Adj_pcBroadband_Mag-bin
% [*] no coverage: 0.187305 (180 of 961)
% [*] not significant: 0.453694 (436 of 961)
% [*] significant: 0.359001 (345 of 961)
% [*] significance fraction: 44.174136 (345 of 781)
% [Broadband] corr number of pairs, mag: -0.17659 (p = 0.02163)
% [Broadband] corr number of pairs, cp: -0.12822 (p = 0.09664)
% [Broadband] corr number of patients, mag: 0.12866 (p = 0.09550)
% [Broadband] corr number of patients, cp: -0.21987 (p = 0.00407)
% [*] Original ROI coverage: 77.7778 %
% [*] Thresholded ROI coverage: 62.6984 %
% [Broadband] corr number of pairs, magVar: 0.11514 (p = 0.13605)
% [Broadband] corr number of pairs, cp: -0.12822 (p = 0.09664)
% [Broadband] corr number of patients, magVar: 0.40485 (p = 0.00000)
% [Broadband] corr number of patients, cp: -0.21987 (p = 0.00407)
% [Broadband] corr coherence mean, variance: 0.63511 (p = 0.00000)
% [Broadband] corr dist, mag: 0.29957 (p = 0.00000)
% 	mean non-diag - mean diag: 0.053433
% 	median non-diag - median diag: 0.064875
% [*] ranksum test:
% 	p-val: 7.062741e-03
% [*] unpaired t test:
% 	p-val: 1.101700e-02
% pcTheta - human fraction of ROIs significant: 0.5804
% [*] saved file: figures/T14/Adj_pcTheta_Mag-bin
% [*] no coverage: 0.187305 (180 of 961)
% [*] not significant: 0.324662 (312 of 961)
% [*] significant: 0.488033 (469 of 961)
% [*] significance fraction: 60.051216 (469 of 781)
% [Theta] corr number of pairs, mag: 0.24424 (p = 0.00022)
% [Theta] corr number of pairs, cp: -0.06629 (p = 0.32333)
% [Theta] corr number of patients, mag: 0.35372 (p = 0.00000)
% [Theta] corr number of patients, cp: -0.08357 (p = 0.21278)
% [*] Original ROI coverage: 77.7778 %
% [*] Thresholded ROI coverage: 62.6984 %
% [Theta] corr number of pairs, magVar: 0.18390 (p = 0.00577)
% [Theta] corr number of pairs, cp: -0.06629 (p = 0.32333)
% [Theta] corr number of patients, magVar: 0.33150 (p = 0.00000)
% [Theta] corr number of patients, cp: -0.08357 (p = 0.21278)
% [Theta] corr coherence mean, variance: 0.86074 (p = 0.00000)
% [Theta] corr dist, mag: 0.07430 (p = 0.10657)
% 	mean non-diag - mean diag: 0.000056
% 	median non-diag - median diag: -0.001991
% [*] ranksum test:
% 	p-val: 5.261806e-01
% [*] unpaired t test:
% 	p-val: 9.912088e-01
% pcAlpha - human fraction of ROIs significant: 0.3337
% [*] saved file: figures/T14/Adj_pcAlpha_Mag-bin
% [*] no coverage: 0.187305 (180 of 961)
% [*] not significant: 0.529657 (509 of 961)
% [*] significant: 0.283039 (272 of 961)
% [*] significance fraction: 34.827145 (272 of 781)
% [Alpha] corr number of pairs, mag: -0.00060 (p = 0.99470)
% [Alpha] corr number of pairs, cp: -0.16987 (p = 0.05623)
% [Alpha] corr number of patients, mag: 0.11274 (p = 0.20694)
% [Alpha] corr number of patients, cp: -0.19631 (p = 0.02697)
% [*] Original ROI coverage: 77.7778 %
% [*] Thresholded ROI coverage: 62.6984 %
% [Alpha] corr number of pairs, magVar: -0.02119 (p = 0.81304)
% [Alpha] corr number of pairs, cp: -0.16987 (p = 0.05623)
% [Alpha] corr number of patients, magVar: 0.13646 (p = 0.12608)
% [Alpha] corr number of patients, cp: -0.19631 (p = 0.02697)
% [Alpha] corr coherence mean, variance: 0.84817 (p = 0.00000)
% [Alpha] corr dist, mag: 0.07268 (p = 0.23219)
% 	mean non-diag - mean diag: 0.002467
% 	median non-diag - median diag: -0.002454
% [*] ranksum test:
% 	p-val: 4.885851e-01
% [*] unpaired t test:
% 	p-val: 7.120552e-01
% pcBeta - human fraction of ROIs significant: 0.3558
% [*] saved file: figures/T14/Adj_pcBeta_Mag-bin
% [*] no coverage: 0.187305 (180 of 961)
% [*] not significant: 0.510926 (491 of 961)
% [*] significant: 0.301769 (290 of 961)
% [*] significance fraction: 37.131882 (290 of 781)
% [Beta] corr number of pairs, mag: 0.04438 (p = 0.60659)
% [Beta] corr number of pairs, cp: -0.34025 (p = 0.00005)
% [Beta] corr number of patients, mag: 0.20307 (p = 0.01731)
% [Beta] corr number of patients, cp: -0.35493 (p = 0.00002)
% [*] Original ROI coverage: 77.7778 %
% [*] Thresholded ROI coverage: 62.6984 %
% [Beta] corr number of pairs, magVar: 0.00774 (p = 0.92850)
% [Beta] corr number of pairs, cp: -0.34025 (p = 0.00005)
% [Beta] corr number of patients, magVar: 0.12334 (p = 0.15102)
% [Beta] corr number of patients, cp: -0.35493 (p = 0.00002)
% [Beta] corr coherence mean, variance: 0.71462 (p = 0.00000)
% [Beta] corr dist, mag: 0.18952 (p = 0.00118)
% 	mean non-diag - mean diag: 0.012071
% 	median non-diag - median diag: 0.002717
% [*] ranksum test:
% 	p-val: 1.550274e-01
% [*] unpaired t test:
% 	p-val: 1.220723e-01
% pcGamma - human fraction of ROIs significant: 0.4012
% [*] saved file: figures/T14/Adj_pcGamma_Mag-bin
% [*] no coverage: 0.187305 (180 of 961)
% [*] not significant: 0.478668 (460 of 961)
% [*] significant: 0.334027 (321 of 961)
% [*] significance fraction: 41.101152 (321 of 781)
% [Gamma] corr number of pairs, mag: -0.06865 (p = 0.39758)
% [Gamma] corr number of pairs, cp: -0.15657 (p = 0.05248)
% [Gamma] corr number of patients, mag: 0.25520 (p = 0.00140)
% [Gamma] corr number of patients, cp: -0.23034 (p = 0.00405)
% [*] Original ROI coverage: 77.7778 %
% [*] Thresholded ROI coverage: 62.6984 %
% [Gamma] corr number of pairs, magVar: 0.14729 (p = 0.06833)
% [Gamma] corr number of pairs, cp: -0.15657 (p = 0.05248)
% [Gamma] corr number of patients, magVar: 0.44706 (p = 0.00000)
% [Gamma] corr number of patients, cp: -0.23034 (p = 0.00405)
% [Gamma] corr coherence mean, variance: 0.71115 (p = 0.00000)
% [Gamma] corr dist, mag: 0.45336 (p = 0.00000)
% 	mean non-diag - mean diag: 0.060755
% 	median non-diag - median diag: 0.073388
% [*] ranksum test:
% 	p-val: 6.036360e-04
% [*] unpaired t test:
% 	p-val: 1.195240e-03

% ========== Run 2019 05 23 ===============================================
% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T14’: File exists
% pcBroadband - human fraction of ROIs significant: 0.2739
% cluster clash: 1.242633921848 mm
% [*] saved file: figures/T14/Adj_pcBroadband_Mag-bin
% [*] no coverage: 0.187305 (180 of 961)
% [*] not significant: 0.453694 (436 of 961)
% [*] significant: 0.359001 (345 of 961)
% [*] significance fraction: 44.174136 (345 of 781)
% [Broadband] corr number of pairs, mag: -0.17659 (p = 0.02163)
% [Broadband] corr number of pairs, cp: -0.12822 (p = 0.09664)
% [Broadband] corr number of patients, mag: 0.12866 (p = 0.09550)
% [Broadband] corr number of patients, cp: -0.21987 (p = 0.00407)
% [*] Original ROI coverage: 77.7778 %
% [*] Thresholded ROI coverage: 62.6984 %
% [Broadband] corr number of pairs, magVar: 0.11514 (p = 0.13605)
% [Broadband] corr number of pairs, cp: -0.12822 (p = 0.09664)
% [Broadband] corr number of patients, magVar: 0.40485 (p = 0.00000)
% [Broadband] corr number of patients, cp: -0.21987 (p = 0.00407)
% [Broadband] corr coherence mean, variance: 0.63511 (p = 0.00000)
% [Broadband] corr dist, mag: 0.29957 (p = 0.00000)
% 	mean non-diag - mean diag: 0.053433
% 	median non-diag - median diag: 0.064875
% [*] ranksum test:
% 	p-val: 7.062741e-03
% [*] unpaired t test:
% 	p-val: 1.101700e-02
% pcTheta - human fraction of ROIs significant: 0.3650
% [*] saved file: figures/T14/Adj_pcTheta_Mag-bin
% [*] no coverage: 0.187305 (180 of 961)
% [*] not significant: 0.324662 (312 of 961)
% [*] significant: 0.488033 (469 of 961)
% [*] significance fraction: 60.051216 (469 of 781)
% [Theta] corr number of pairs, mag: 0.24424 (p = 0.00022)
% [Theta] corr number of pairs, cp: -0.06629 (p = 0.32333)
% [Theta] corr number of patients, mag: 0.35372 (p = 0.00000)
% [Theta] corr number of patients, cp: -0.08357 (p = 0.21278)
% [*] Original ROI coverage: 77.7778 %
% [*] Thresholded ROI coverage: 62.6984 %
% [Theta] corr number of pairs, magVar: 0.18390 (p = 0.00577)
% [Theta] corr number of pairs, cp: -0.06629 (p = 0.32333)
% [Theta] corr number of patients, magVar: 0.33150 (p = 0.00000)
% [Theta] corr number of patients, cp: -0.08357 (p = 0.21278)
% [Theta] corr coherence mean, variance: 0.86074 (p = 0.00000)
% [Theta] corr dist, mag: 0.07430 (p = 0.10657)
% 	mean non-diag - mean diag: 0.000056
% 	median non-diag - median diag: -0.001991
% [*] ranksum test:
% 	p-val: 5.261806e-01
% [*] unpaired t test:
% 	p-val: 9.912088e-01
% pcAlpha - human fraction of ROIs significant: 0.2099
% [*] saved file: figures/T14/Adj_pcAlpha_Mag-bin
% [*] no coverage: 0.187305 (180 of 961)
% [*] not significant: 0.529657 (509 of 961)
% [*] significant: 0.283039 (272 of 961)
% [*] significance fraction: 34.827145 (272 of 781)
% [Alpha] corr number of pairs, mag: -0.00060 (p = 0.99470)
% [Alpha] corr number of pairs, cp: -0.16987 (p = 0.05623)
% [Alpha] corr number of patients, mag: 0.11274 (p = 0.20694)
% [Alpha] corr number of patients, cp: -0.19631 (p = 0.02697)
% [*] Original ROI coverage: 77.7778 %
% [*] Thresholded ROI coverage: 62.6984 %
% [Alpha] corr number of pairs, magVar: -0.02119 (p = 0.81304)
% [Alpha] corr number of pairs, cp: -0.16987 (p = 0.05623)
% [Alpha] corr number of patients, magVar: 0.13646 (p = 0.12608)
% [Alpha] corr number of patients, cp: -0.19631 (p = 0.02697)
% [Alpha] corr coherence mean, variance: 0.84817 (p = 0.00000)
% [Alpha] corr dist, mag: 0.07268 (p = 0.23219)
% 	mean non-diag - mean diag: 0.002467
% 	median non-diag - median diag: -0.002454
% [*] ranksum test:
% 	p-val: 4.885851e-01
% [*] unpaired t test:
% 	p-val: 7.120552e-01
% pcBeta - human fraction of ROIs significant: 0.2238
% [*] saved file: figures/T14/Adj_pcBeta_Mag-bin
% [*] no coverage: 0.187305 (180 of 961)
% [*] not significant: 0.510926 (491 of 961)
% [*] significant: 0.301769 (290 of 961)
% [*] significance fraction: 37.131882 (290 of 781)
% [Beta] corr number of pairs, mag: 0.04438 (p = 0.60659)
% [Beta] corr number of pairs, cp: -0.34025 (p = 0.00005)
% [Beta] corr number of patients, mag: 0.20307 (p = 0.01731)
% [Beta] corr number of patients, cp: -0.35493 (p = 0.00002)
% [*] Original ROI coverage: 77.7778 %
% [*] Thresholded ROI coverage: 62.6984 %
% [Beta] corr number of pairs, magVar: 0.00774 (p = 0.92850)
% [Beta] corr number of pairs, cp: -0.34025 (p = 0.00005)
% [Beta] corr number of patients, magVar: 0.12334 (p = 0.15102)
% [Beta] corr number of patients, cp: -0.35493 (p = 0.00002)
% [Beta] corr coherence mean, variance: 0.71462 (p = 0.00000)
% [Beta] corr dist, mag: 0.18952 (p = 0.00118)
% 	mean non-diag - mean diag: 0.012071
% 	median non-diag - median diag: 0.002717
% [*] ranksum test:
% 	p-val: 1.550274e-01
% [*] unpaired t test:
% 	p-val: 1.220723e-01
% pcGamma - human fraction of ROIs significant: 0.2523
% [*] saved file: figures/T14/Adj_pcGamma_Mag-bin
% [*] no coverage: 0.187305 (180 of 961)
% [*] not significant: 0.478668 (460 of 961)
% [*] significant: 0.334027 (321 of 961)
% [*] significance fraction: 41.101152 (321 of 781)
% [Gamma] corr number of pairs, mag: -0.06865 (p = 0.39758)
% [Gamma] corr number of pairs, cp: -0.15657 (p = 0.05248)
% [Gamma] corr number of patients, mag: 0.25520 (p = 0.00140)
% [Gamma] corr number of patients, cp: -0.23034 (p = 0.00405)
% [*] Original ROI coverage: 77.7778 %
% [*] Thresholded ROI coverage: 62.6984 %
% [Gamma] corr number of pairs, magVar: 0.14729 (p = 0.06833)
% [Gamma] corr number of pairs, cp: -0.15657 (p = 0.05248)
% [Gamma] corr number of patients, magVar: 0.44706 (p = 0.00000)
% [Gamma] corr number of patients, cp: -0.23034 (p = 0.00405)
% [Gamma] corr coherence mean, variance: 0.71115 (p = 0.00000)
% [Gamma] corr dist, mag: 0.45336 (p = 0.00000)
% 	mean non-diag - mean diag: 0.060755
% 	median non-diag - median diag: 0.073388
% [*] ranksum test:
% 	p-val: 6.036360e-04
% [*] unpaired t test:
% 	p-val: 1.195240e-03
