close all;
clc;
rng shuffle;
clear;

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
%metrics = {'pcBroadband'};

makeMat = false;
for imm = 1:length(metrics)
    metric = metrics{imm};
    f_L = sprintf('figures/T14L/%s_Mag-bin_Adj.mat',metrics{imm});
    f_R = sprintf('figures/T14R/%s_Mag-bin_Adj.mat',metrics{imm});
    if (exist(f_L,'file') && exist(f_R,'file'))
        h = figure;
        fig_size_scale = 0.75;
        color_distt = 1*[1 1 1]; %0.6*[1 1 1];
        color_isnan = 1*[1 1 1]; %0.6*[1 1 1];
        %color_not_sig = 0.999*[1 1 1];
        color_not_siglr = 0.75*[1 1 1];
        color_not_covlr = [1 1 1]; %0.6*[1 1 1];
        fontsz = 10;
        fsize = fontsz;
        trig_mag_no_cp = true;
        fig_fmt = '-depsc';
        trig_eps = false;
        trig_mag_no_cp = true;
        trig_plot_mag_cov = true;
        
        set(h,'Position',round(fig_size_scale*[0 0 0.95*1080 0.8*1080]))
            
        L = load(f_L);
        R = load(f_R);
        
        % Statistics left vs. right
        p_val = 0.01;
        issig = zeros(size(L.AdjMagL));
        for iiL = 1:(length(issig)-1)
            for jjL = (iiL+1):length(issig)
                x = L.AdjMagL{iiL,jjL};
                y = R.AdjMagL{iiL,jjL};
                if (~isempty(x)) && (~isempty(y))
                    p = ranksum(x,y);
                    if (p < p_val)
                        is = 1;
                    else
                        is = 0;
                    end
                else
                    is = NaN;
                end
                issig(iiL,jjL) = is;
                issig(jjL,iiL) = is;
            end
        end
        
        Adj_plt = L.AdjMag - R.AdjMag;
        Adj = L.Adj - R.Adj;
        
         % Load human cache
        dir_cacheL = './cache';
        load(sprintf('%s/xsub_out_all_%i.mat',dir_cacheL,imm));
        iii = 3;
        
        % Calculate statistics
        
        
        % Calculate rows of nans
        Adj2 = Adj;
        ind_isnan = false(length(Adj2),1);
        for i = 1:length(ind_isnan)
            ind_isnan(i) = (sum(isnan(Adj2(i,:))) == length(ind_isnan));
        end
    
        Dmat = Inf(n_rois,n_rois);
        Dmat_plt = Dmat;
        rois_plt = rois;
        for i = 1:length(rois_plt)
            rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '_','/');
        end
        [n_roi_p,~] = size(Adj_plt);
        % generate image from matrix
        Adj_plt(Dmat_plt <= dist_thresh) = nanmean(Adj_plt(:));
        v = Adj_plt;
        [v_n, v_m] = size(v);
        
        % ------------------------ colormap
        %map = corrcmap(100);
        %map = viridis(100);
        map = redblue(100);
        %map = [viridis(50);flipud(inferno(50))];
        %map = [([0 1 0] + [1 0 0] .* mean(bone(50),2)); flipud([1 0 0] + [0 1 0] .* mean(bone(50),2))];
        
        if (iii == 3)
            %minv = min(v(v~=0));
            %maxv = max(v(v~=0));
            minv = min(v(~isnan(issig)));
            maxv = max(v(~isnan(issig)));
            maxabs = max(abs(minv),abs(maxv));
            minv = (-1)*maxabs;
            maxv = maxabs;
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
%                 if ((iii == 3) && (Adj(i,j) == 0))
%                     %Im(i,j,:) = color_not_sig;
%                 end
                if (issig(i,j) == 0)
                    Im(i,j,:) = color_not_siglr;
                end
                if (isnan(issig(i,j)))
                    Im(i,j,:) = color_not_covlr;
                end
            end
        end
        Im = Im(~ind_isnan,~ind_isnan,:); % remove rows of nans
        rois_plt = rois_plt(~ind_isnan);

        % --- Clustering (output: cluster_i) --------------------------------------
        %adjct_dist_cl = adjct_dist(ind_hum2mac,ind_hum2mac);

        % TODO: RERUN xsub_out_stats to enable adjct_dist clustering
        adjct_dist_cl = adjct_dist;
        adjct_dist_cl = adjct_dist_cl(~ind_isnan,~ind_isnan);
        %adjct_dist = adjct_dist(~ind_isnan)
        n_rois_cl = length(rois_plt);
        Y = zeros(1,nchoosek(n_rois_cl,2));
        yc = 1;
        for m1 = 1:(n_rois_cl-1)
            for m2 = (m1+1):n_rois_cl
                ad = adjct_dist_cl{m2,m1};
                if (isempty(ad))
                    Y(yc) = Inf;%Inf; 600;
                else
                    Y(yc) = mean(ad);
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

        Z = linkage(Y); %,'centroid'
        cluster_i = optimalleaforder(Z,Y,'transformation','inverse'); % ,'transformation','inverse'
        roi_dist = roi_dist(cluster_i,cluster_i);
        %roi_dist(isinf(roi_dist)) = ;
        clash = nansum(nansum(triu(roi_dist,1) - triu(roi_dist,2)));
        fprintf('cluster clash: %.12f mm\n',clash)
        %cluster_i = optimalleaforder(Z,Y);

        % -------------------------------------------------------------------------

        % Apply clustering
        Im = Im(cluster_i,cluster_i,:);
        rois_plt = rois_plt(cluster_i);
        
        % Remove unknown
        idx_remove = strcmp(rois_plt,'unknown');
        
        % Apply ROI filter
        Im = Im(~idx_remove,~idx_remove,:);
        rois_plt = rois_plt(~idx_remove);
        roi_dist = roi_dist(~idx_remove,~idx_remove);
        adjct_dist_cl = adjct_dist_cl(~idx_remove);

        % Convert ROI names
        for i = 1:length(rois_plt)
            rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '_','/');
            rois_plt{i} = convertRoiDK(rois_plt{i});
        end
        
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
        %return
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
            cb = colorbar('Ytick',linspace(minv,maxv,5));
            caxis([minv maxv]);
            set(cb,'TickLength',0);
        else
            cb = colorbar;
            set(cb,'TickLength',0);
        end

        %export_fig(sprintf('figures/T14/Adj_%s_%s_2',metric,ss),'-eps');
        print(h,sprintf('figures/T14%s/Adj_%s_%s','LR',metric,ss),fig_fmt);
        if (trig_eps)
            print(h,sprintf('figures/T14%s/Adj_%s_%s','LR',metric,ss),'-depsc');
        end
        close(h);
        
    end
end

if (makeMat)

Hem = {'R','L'};

for hL = 1:length(Hem)
    iii = 3;
    hemiL = Hem{hL};

    n_perm = 10000;
    perm_alpha_sec = 12;
    cp_thresh_override = 0.05; 
    n_pairs_thresh = 20; % at least this many electrode pairs to be considered
    n_subs_thresh = 4; % at least this many subjects to be considered
    n_subs_ct_thresh = 4; % significant CTs in region pair must be from at least this many subjects
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
    system(['mkdir figures/T14',hemiL]);
    system(['mkdir figures/T14',hemiL]);


    [~,host] = system('hostname');

    if contains(host,'ubuntu_1604')
        dir_artL = '/nas_share/RawData/data/h5_notch20/art';
        dir_resL = '/nas_share/RawData/data/results/coh_w10';
        dir_corL = '/nas_share/RawData/data/coreg';
        dir_cacheL = ['./cache',hemiL];
        subjects_dirL = '/mnt/cuenap_ssd/coregistration';
        dir_h5L = '/nas_share/RawData/data/h5_notch20';
    else
        dir_artL = '/media/klab/internal/data/h5_notch20/art';
        dir_resL = '/media/klab/internal/data/results/coh_w10';
        dir_corL = '/media/klab/internal/data/coreg';
        dir_cacheL = ['./cache_',hemiL];
        subjects_dirL = '/mnt/cuenap_ssd/coregistration';
        dir_h5L = '/media/klab/KLAB101/h5_notch20';
    end

    %metrics = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};
    %metrics = {'pcBroadband'};

    for iM = 1:length(metrics)
        metric = metrics{iM};

        % Load human cache
        load(sprintf('%s/xsub_out_all_%i.mat',dir_cacheL,iM));


        % Calculate final functional interaction matrix
        Adj = nan(n_rois,n_rois);
        AdjMag = nan(n_rois,n_rois);
        AdjCP = nan(n_rois,n_rois);
        AdjMagVar = nan(n_rois,n_rois);
        AdjMagReS = nan(n_rois,n_rois);
        AdjMagL = cell(n_rois,n_rois);
        %AdjCPVar = nan(n_rois,n_rois);
        N_bchan = nan(n_rois,n_rois);
        cp_thresh = cp_thresh_override;
        Dmat = Inf(n_rois,n_rois); % set to inf to avoid removing from every instance
        DistsAtl = [];

        for i1 = 1:n_rois
            for i2 = 1:n_rois
                AA = AdjAtl{i1,i2};
                AA_dist = adjct_dist{i1,i2};
                AA_sub = AdjAtl_sid{i1,i2};
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
                    end
                end
            end
        end

        frac_xsub_msu = sum(Adj(~isnan(Adj))) /numel((~isnan(Adj)));
        fprintf('%s - human fraction of ROIs significant: %.4f\n',metric,frac_xsub_msu)

        % --- PLOT  --------------------------------------------------
        color_distt = 0.6*[1 1 1];
        color_isnan = 0.6*[1 1 1];
        fontsz = 10;
        fsize = fontsz;

        % Calculate rows of nans
        Adj2 = Adj;
        ind_isnan = false(length(Adj2),1);
        for i = 1:length(ind_isnan)
            ind_isnan(i) = (sum(isnan(Adj2(i,:))) == length(ind_isnan));
        end

        % Save Adj
        ss = 'Mag-bin';
        save(sprintf('figures/T14%s/%s_%s_Adj',hemiL,metric,ss),'Adj','AdjMag','AdjCP','AdjMagVar','AdjMagL');
        
        for iii = 3:3

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
                color_not_sig = 0.999*[1 1 1];
            end
            Dmat_plt = Dmat;
            rois_plt = rois;
            for i = 1:length(rois_plt)
                rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '_','/');
            end
            [n_roi_p,~] = size(Adj_plt);
            % generate image from matrix
            Adj_plt(Dmat_plt <= dist_thresh) = nanmean(Adj_plt(:));
            v = Adj_plt;
            [v_n, v_m] = size(v);
            map = corrcmap(100);
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
                    if ((iii == 3) && (Adj(i,j) == 0))
                        Im(i,j,:) = color_not_sig;
                    end
                end
            end
            Im = Im(~ind_isnan,~ind_isnan,:); % remove rows of nans
            rois_plt = rois_plt(~ind_isnan);

            % --- Clustering (output: cluster_i) --------------------------------------
            %adjct_dist_cl = adjct_dist(ind_hum2mac,ind_hum2mac);

            % TODO: RERUN xsub_out_stats to enable adjct_dist clustering
            adjct_dist_cl = adjct_dist;
            adjct_dist_cl = adjct_dist_cl(~ind_isnan,~ind_isnan);
            %adjct_dist = adjct_dist(~ind_isnan)
            n_rois_cl = length(rois_plt);
            Y = zeros(1,nchoosek(n_rois_cl,2));
            yc = 1;
            for m1 = 1:(n_rois_cl-1)
                for m2 = (m1+1):n_rois_cl
                    ad = adjct_dist_cl{m2,m1};
                    if (isempty(ad))
                        Y(yc) = Inf;%Inf; 600;
                    else
                        Y(yc) = mean(ad);
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

            Z = linkage(Y); %,'centroid'
            cluster_i = optimalleaforder(Z,Y,'transformation','inverse'); % ,'transformation','inverse'
            roi_dist = roi_dist(cluster_i,cluster_i);
            %roi_dist(isinf(roi_dist)) = ;
            clash = nansum(nansum(triu(roi_dist,1) - triu(roi_dist,2)));
            fprintf('cluster clash: %.12f mm\n',clash)
            %cluster_i = optimalleaforder(Z,Y);

            % -------------------------------------------------------------------------

            % Apply clustering
            Im = Im(cluster_i,cluster_i,:);
            rois_plt = rois_plt(cluster_i);

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
            %return
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
                cb = colorbar('Ytick',linspace(minv,maxv,5));
                caxis([minv maxv]);
                set(cb,'TickLength',0);
            else
                cb = colorbar;
                set(cb,'TickLength',0);
            end

            %export_fig(sprintf('figures/T14/Adj_%s_%s_2',metric,ss),'-eps');
            
            print(h,sprintf('figures/T14%s/Adj_%s_%s',hemiL,metric,ss),fig_fmt);
            if (trig_eps)
                print(h,sprintf('figures/T14%s/Adj_%s_%s',hemiL,metric,ss),'-depsc');
            end
            close(h);

            % print info
            frac_nocov = sum(sum(isnan(Adj_plt)))/numel(Adj_plt);
            frac_nosig = sum(sum(~isnan(Adj_plt) & (Adj_plt == 0)))/numel(Adj_plt);
            frac_sig = sum(sum(~isnan(Adj_plt) & (Adj_plt ~= 0)))/numel(Adj_plt);
            fprintf('[*] saved file: %s\n',sprintf('figures/T14%s/Adj_%s_%s',hemiL,metric,ss));
            fprintf('[*] no coverage: %.6f\n',frac_nocov);
            fprintf('[*] not significant: %.6f\n',frac_nosig);
            fprintf('[*] significant: %.6f\n',frac_sig);
            fprintf('[*] significance fraction: %.6f\n',100*frac_sig/(frac_sig + frac_nosig));


            if ((iii == 3) && trig_plot_mag_cov)

                Mag = nan(n_comb_atl,1);
                MagVar = nan(n_comb_atl,1);
                CP = nan(n_comb_atl,1);
                Npairs = nan(n_comb_atl,1);
                Nsubs = nan(n_comb_atl,1);
                %DistsAtl = [];%nan(n_comb_atl,1);
                count = 1;
                for i1 = 1:(n_rois-1)
                    for i2 = (i1+1):n_rois

                        AA = AdjAtl_sid{i1,i2};
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


                % Number of pairs vs magnitude
                h = figure;
                set(h,'Position',round(fig_size_scale*[0 0 0.95*1080 0.8*1080]))
                is_sig_mag = ((Mag ~= 0)&(~isnan(Mag)));
                subplot(3,1,1);
                plot(Npairs(is_sig_mag),Mag(is_sig_mag),'black.');
                xlabel('Number of bipolar pairs');
                ylabel(sprintf('%s Coherence',metric(3:end)));
                axis tight;
                box off;
                set(gca,'TickDir','out');
                hold all;
                [corr_val_mag, corr_pval_mag] = corr(Npairs(is_sig_mag),Mag(is_sig_mag));
                fprintf('[%s] corr number of pairs, mag: %.5f (p = %.5f)\n',metric(3:end),corr_val_mag,corr_pval_mag);

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
                print(h,sprintf('figures/T14%s/Npairs_%s_mag_%i_p_%d_cp_%i_p_%d',hemiL,metric,...
                    round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp),'-depsc');

                close(h);

                % Number of unique patients vs magnitude
                h = figure;
                set(h,'Position',round(fig_size_scale*[0 0 0.95*1080 0.8*1080]))
                is_sig_mag = ((Mag ~= 0)&(~isnan(Mag)));
                subplot(3,1,1);
                plot(Nsubs(is_sig_mag),Mag(is_sig_mag),'black.');
                xlabel('Number of patients');
                ylabel(sprintf('%s Coherence',metric(3:end)));
                axis tight;
                box off;
                set(gca,'TickDir','out');
                hold all;
                [corr_val_mag, corr_pval_mag] = corr(Nsubs(is_sig_mag),Mag(is_sig_mag));
                fprintf('[%s] corr number of patients, mag: %.5f (p = %.5f)\n',metric(3:end),corr_val_mag,corr_pval_mag);

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
                print(h,sprintf('figures/T14%s/Nsubs_%s_mag_%i_p_%d_cp_%i_p_%d',hemiL,metric,...
                    round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp),'-depsc');

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
                print(h,sprintf('figures/T14%s/Npairs_%s_magVar_%i_p_%d_cp_%i_p_%d',hemiL,metric,...
                    round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp),'-depsc');

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
                print(h,sprintf('figures/T14%s/Nsubs_%s_magVar_%i_p_%d_cp_%i_p_%d',hemiL,metric,...
                    round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp),'-depsc');

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

                print(h,sprintf('figures/T14%s/mag_%s_magVar_%i_p_%d',hemiL,metric,...
                    round(1000*corr_val_cp),corr_pval_cp),'-depsc');

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
                print(h,sprintf('figures/T14%s/dist_%s_mag_%i_p_%d',hemiL,metric,...
                    round(1000*corr_val_cp),corr_pval_cp),'-depsc');
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
                try
                    [p_rs,~,stats] = ranksum(DiagMag,TriuMag);
                catch
                    p_rs = NaN;
                end
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
                print(h,sprintf('figures/T14%s/intra_inter_%s_ranksum-%i_ttest-%i',hemiL,metric,...
                    round(1000*p_rs),round(1000*p_tt2)),'-depsc');
                close(h);

                %return
            end





        end


    end

end


end




% Run May 25, 2019
% cluster clash: 817.567870037500 mm
% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T14R’: File exists
% mkdir: cannot create directory ‘figures/T14R’: File exists
% pcBroadband - human fraction of ROIs significant: 0.0764
% cluster clash: 922.660333772597 mm
% [*] saved file: figures/T14R/Adj_pcBroadband_Mag-bin
% [*] no coverage: 0.648920
% [*] not significant: 0.274691
% [*] significant: 0.076389
% [*] significance fraction: 21.758242
% [Broadband] corr number of pairs, mag: -0.09955 (p = 0.52534)
% [Broadband] corr number of pairs, cp: -0.04018 (p = 0.79809)
% [Broadband] corr number of patients, mag: 0.08444 (p = 0.59035)
% [Broadband] corr number of patients, cp: -0.25133 (p = 0.10401)
% [*] Original ROI coverage: 72.0635 %
% [*] Thresholded ROI coverage: 34.7619 %
% [Broadband] corr number of pairs, magVar: 0.10240 (p = 0.51350)
% [Broadband] corr number of pairs, cp: -0.04018 (p = 0.79809)
% [Broadband] corr number of patients, magVar: 0.24040 (p = 0.12047)
% [Broadband] corr number of patients, cp: -0.25133 (p = 0.10401)
% [Broadband] corr coherence mean, variance: 0.89072 (p = 0.00000)
% [Broadband] corr dist, mag: 0.24644 (p = 0.01393)
% 	mean non-diag - mean diag: 0.013356
% 	median non-diag - median diag: 0.012230
% [*] ranksum test:
% 	p-val: 6.694208e-01
% [*] unpaired t test:
% 	p-val: 3.946955e-01
% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T14L’: File exists
% mkdir: cannot create directory ‘figures/T14L’: File exists
% pcBroadband - human fraction of ROIs significant: 0.1211
% cluster clash: 802.969655754041 mm
% [*] saved file: figures/T14L/Adj_pcBroadband_Mag-bin
% [*] no coverage: 0.685957
% [*] not significant: 0.192901
% [*] significant: 0.121142
% [*] significance fraction: 38.574939
% [Broadband] corr number of pairs, mag: -0.34516 (p = 0.00260)
% [Broadband] corr number of pairs, cp: -0.11993 (p = 0.30880)
% [Broadband] corr number of patients, mag: -0.22245 (p = 0.05679)
% [Broadband] corr number of patients, cp: -0.20204 (p = 0.08430)
% [*] Original ROI coverage: 60.1587 %
% [*] Thresholded ROI coverage: 31.1111 %
% [Broadband] corr number of pairs, magVar: -0.15496 (p = 0.18741)
% [Broadband] corr number of pairs, cp: -0.11993 (p = 0.30880)
% [Broadband] corr number of patients, magVar: -0.07900 (p = 0.50348)
% [Broadband] corr number of patients, cp: -0.20204 (p = 0.08430)
% [Broadband] corr coherence mean, variance: 0.68291 (p = 0.00000)
% [Broadband] corr dist, mag: 0.13050 (p = 0.10329)
% 	mean non-diag - mean diag: 0.044696
% 	median non-diag - median diag: 0.059757
% [*] ranksum test:
% 	p-val: 9.645239e-02
% [*] unpaired t test:
% 	p-val: 1.034861e-01
