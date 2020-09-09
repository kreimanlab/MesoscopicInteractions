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
system('mkdir figures/T14_allatl_resampTv2');



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
for atl = [2] %1:20

    system(sprintf('mkdir figures/T14_allatl_resampTv2/atl%i_%s',atl,AtlNames{atl}));
    
    for iM = [1] %1:length(metrics) % [1 5] %
        metric = metrics{iM};

        % Temporal resampling
        AdjMagRep = {};
        AdjMagVarRep = {};
        BpairsRep = {};
        for repL = [5 6]
            % Load human cache
            Ca_hum = load(sprintf('%s/xsub_out_all_%i_atl%i_resampTv2-%i-1.mat',dir_cacheL,iM,atl,repL));
            n_rois = Ca_hum.n_rois;

            % Calculate final functional interaction matrix
            Adj = nan(n_rois,n_rois);
            AdjMag = nan(n_rois,n_rois);
            AdjMag4cl = nan(n_rois,n_rois);
            AdjCP = nan(n_rois,n_rois);
            AdjMagVar = nan(n_rois,n_rois);
            AdjMagReS = nan(n_rois,n_rois);
            AdjMagL = cell(n_rois,n_rois);
            %AdjCPVar = nan(n_rois,n_rois);
            N_bchan = nan(n_rois,n_rois);
            cp_thresh = cp_thresh_override;
            Dmat = Inf(n_rois,n_rois); % set to inf to avoid removing from every instance
            DistsAtl = [];
            %Bpairs = [];

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
                            %Bpairs = [Bpairs; AA];
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

                    % AdjMag for clustering
                    AdjMag4cl(i1,i2) = mean(AA(AA~=0));

                end
            end
            
            % Save for temporal
            AdjMagRep{repL} = AdjMag;
            AdjMagVarRep{repL} = AdjMagVar;
            

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

                    %h = figure('visible','on');
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
                    end

                    Adj_plt2_cl = AdjMag4cl;

                    rois = Ca_hum.rois;
                    Dmat_plt = Dmat;
                    rois_plt = rois;

                    % checkpoint:
                    % Adjacency matrix Adj_plt: 36 by 36
                    % Adjacency rois_plt: 36 by 1
                    % Distance matrix Dmat_plt: 


                    % Filter out unknown from adjacency to plot



                    cond_empty = false(length(rois),1);
                    for ce = 1:length(cond_empty)
                        elem = rois{ce};
                        cond_empty(ce) = (all(isspace(elem)) | isempty(elem));
                    end

                    cond_contains = (contains(lower(rois),'unknown')) | (contains(lower(rois),'???')) | (contains(lower(rois),'wall')) | (cond_empty);
                    known_idx = (~ cond_contains);
                    ind_isnan_master(~known_idx) = true;
                    Adj_plt = Adj_plt(known_idx,known_idx);
                    Adj_plt2_cl = Adj_plt2_cl(known_idx,known_idx);
                    rois_plt = rois_plt(known_idx);
                    Dmat_plt = Dmat_plt(known_idx,known_idx);
                    
                    % Save bipolar pairs
                    AdjAA = Ca_hum.AdjAtl;
                    AdjAA = AdjAA(known_idx,known_idx);

                    % Clean ROI labels for human readability
                    for i = 1:length(rois_plt)
                        % Single atlas exceptions
                        if (startsWith(rois_plt{i},'L_'))
                            rpt = strsplit(rois_plt{i},'L_');
                            rpt = strsplit(rpt{end},'_ROI');
                            rois_plt{i} = rpt{1};
                        end
                        if (endsWith(rois_plt{i},'_M132'))
                            rpt = strsplit(rois_plt{i},'_M132');
                            rois_plt{i} = rpt{1};
                        end
                        if (startsWith(rois_plt{i},'FVE_all_'))
                            rpt = strsplit(rois_plt{i},'FVE_all_');
                            rois_plt{i} = rpt{end};
                        elseif (startsWith(rois_plt{i},'FVE_'))
                            rpt = strsplit(rois_plt{i},'FVE_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'LVE00_'))
                            rpt = strsplit(rois_plt{i},'LVE00_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'PHT00_'))
                            rpt = strsplit(rois_plt{i},'PHT00_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'BoninBailey_'))
                            rpt = strsplit(rois_plt{i},'BoninBailey_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'FerryEtAl_00'))
                            rpt = strsplit(rois_plt{i},'FerryEtAl_00');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'UD86_'))
                            rpt = strsplit(rois_plt{i},'UD86_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'PGR91_'))
                            rpt = strsplit(rois_plt{i},'PGR91_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'LyonKaas02_'))
                            rpt = strsplit(rois_plt{i},'LyonKaas02_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'BRL87_'))
                            rpt = strsplit(rois_plt{i},'BRL87_');
                            rois_plt{i} = rpt{end};
                        end

                        % Remove underscores
                        rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '_',' ');
                        rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '.',' ');

                        if (atl == 2)
                            rois_plt{i} = convertRoiDK(rois_plt{i});
                        end
                    end



                    % Distance threshold: add nans to Adj_plt
                    %Adj_plt(Dmat_plt <= dist_thresh) = nanmean(Adj_plt(:));
                    Adj_plt2 = Adj_plt;
                    dist_thresh = Ca_hum.dist_thresh;
                    Adj_plt(Dmat_plt <= dist_thresh) = nan;
                    Adj_plt2_cl(Dmat_plt <= dist_thresh) = nan;

                    %return
                    % Filter out nans nodes
                    cov_idx = ~ all(isnan(Adj_plt));
                    % Bipolar pairs
                    AdjAA = AdjAA(cov_idx,cov_idx);
                    if (iii2 == 1)
                        BpairsRep{repL} = AdjAA;
                    end
                    
                    % Turn off filter
                    cov_idx = true(1,length(Adj_plt));
                    
                    Adj_plt = Adj_plt(cov_idx,cov_idx);
                    Adj_plt2 = Adj_plt2(cov_idx,cov_idx);
                    Adj_plt2_cl = Adj_plt2_cl(cov_idx,cov_idx);
                    rois_plt = rois_plt(cov_idx);
                    Dmat_plt = Dmat_plt(cov_idx,cov_idx);
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
                    %Adj_plt2 = Adj_plt;
                    %rois_plt2 = rois_plt;




                    % --- Clustering (output: cluster_i) --------------------------------------
    %adjct_dist_cl = adjct_dist(ind_hum2mac,ind_hum2mac);



                    %if ((iM == 1) && (iii2 == 1))
                    if ((iii2 == 1))
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
                                %ad = Adj_plt(m2,m1);
                                ad = Adj_plt(m1,m2);
                                if (isnan(ad))
                                    Y(yc) = 0;
                                else
                                    Y(yc) = ad;
                                end
                                yc = yc + 1;


                            end
                        end
                        %X = Y;
                        %Y = pdist(X);
                        roi_dist = squareform(Y);
                        % average   2778.832643576550 mm
                        % centroid  2773.791496333409 mm
                        % complete  2808.894068525726 mm
                        % median    2767.246370358684 mm
                        % single    2836.305543514634 mm
                        % ward      2778.832643576550 mm
                        % weighted  2778.918494966365 mm

                        %

                        % Imputed
                        %Z = linkage(Y,'complete'); %,'centroid'

                        % Marginalized
                        %M = Adj_plt2_cl;
                        M = Adj_plt2;
                        M = (M + M.')/2;

                        %Z = linkage(M,'ward',m);
                        %return
                        
                        % Turn off clustering
                        M = eye(size(M));
                        
                        Z = linkage(M,'ward',@nanDist4clustering);

                        cluster_i = optimalleaforder(Z,Y); % ,'transformation','inverse'
                        
                        % Override clustering
                        cluster_i = 1:length(M);
                        roi_dist = roi_dist(cluster_i,cluster_i);
                        %roi_dist(isinf(roi_dist)) = ;
                        clash = nansum(nansum(triu(roi_dist,1) - triu(roi_dist,2)));
                        fprintf('cluster clash: %.12f mm\n',clash)
                        save(sprintf('./cache/figure_t14_%i_atl%i_%s_resampTv2_rep%i',iM,atl,AtlNames{atl},repL),'AdjMag','AdjMagVar','Adj_plt','Adj_plt2','Adj_plt2_cl','Y','Z','M','cluster_i','rois_plt','atl');
                    %cluster_i = optimalleaforder(Z,Y);
                    end

                    % -------------------------------------------------------------------------

                    % moved saving up
    %                 if (iii2 == 1)
    %                     save(sprintf('./cache/figure_t14_%i_atl%i_%s',iM,atl,AtlNames{atl}),'Adj_plt','Adj_plt2','Adj_plt2_cl','Y','Z','M','cluster_i','rois_plt','atl');
    %                 end

                    % bypass clustering
                    %cluster_i = 1:length(cluster_i);
                    %Cai = load('cache/fig_cluster3_cluster_i.mat');
                    %cluster_i = Cai.cluster_i;

                    % Apply clustering
                    Im = Im(cluster_i,cluster_i,:);
                    %rois_plt = rois_plt(1:length(cluster_i));
                    rois_plt = rois_plt(cluster_i);

                    % Print Figure 2A coherence
                    if (atl==2)
                        if (iii2 == 1)
                            coh_mtxt = 'mean';
                            Adj_tmp = Adj_plt2(cluster_i,cluster_i);
                            %close all;
                        elseif (iii2 == 2)
                            coh_mtxt = 'std';
                            Adj_tmp = sqrt(Adj_plt2(cluster_i,cluster_i));
                            %imagesc(Adj_tmp); colormap inferno; colorbar;
                            %return
                        end

                        idx_1 = strcmp(rois_plt,'Pars Opercularis');
                        idx_2 = strcmp(rois_plt,'Superior Temporal');
                        fprintf('[*] %s - %s, coh %s: %.2f\n',...
                            rois_plt{idx_1},rois_plt{idx_2},coh_mtxt,Adj_tmp(idx_1,idx_2));
                    end


                    % Manual adjustment
                    n_wrap = 0;
                    %n_wrap = round(0.2258 * n_rois_cl); % 9 - Translate the whole clustering and wrap around
                    cluster_i_manual = [(n_rois_cl - (n_wrap-1)):n_rois_cl, 1:(n_rois_cl - (n_wrap))];
                    Im = Im(cluster_i_manual,cluster_i_manual,:);
                    rois_plt = rois_plt(cluster_i_manual);

                    rois_plt_sav = rois_plt;


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

                    % numbering
                    stride_num = 1;
                    rois_plt2 = cell(size(rois_plt));
                    for inum = 1:length(rois_plt)
                        if (mod(inum-1,stride_num) == 0)
                            n_suf = sprintf('%i',inum);
                            while (length(n_suf) < 4)
                                n_suf = [' ',n_suf];
                            end
                        else
                            n_suf = '    ';
                        end
                        rois_plt{inum} = [rois_plt{inum},'\bf{',n_suf,'}'];
                        rois_plt2{inum} = rois_plt{inum}; %['\bf{',n_suf,'}'];
                        %rois_plt{inum} = ['\bf{',rois_plt{inum},'} \t'];
                    end

                    sub_idx = stride_axislabel:stride_axislabel:length(rois_plt);
                    rty = 1:length(rois_plt);
                    yticks(rty(sub_idx));
                    yticklabels(rois_plt(sub_idx));

                    % Horizontal axis
                    rtx = 1:length(rois_plt2);
                    xticks(rtx(sub_idx));
                    xticklabels(rois_plt2(sub_idx));
                    xtickangle(90);

        %             yticks(1:length(rois_plt));
        %             yticklabels(rois_plt);
        %             xticks(1:length(rois_plt));
        %             xticklabels(rois_plt);
        %             xtickangle(90);

                    set(gca,'tickdir','out');
                    set(gca,'fontsize',fontsz*(31/n_rois_cl)^(0.5));
                    set(gca,'TickLength',[0.001, 0.001])
                    daspect([1 1 1]);
                    set(gca,'FontName','Arial');

                    % Move axis if needed
                    ax = gca;
                    tp = ax.Position;
                    ax1_xoffset = 0.05; % 0.05
                    tp(1) = tp(1) - ax1_xoffset;
                    ax.Position = tp;
                    if(ax.FontSize > 10)
                        ax.FontSize = 10;
                    end
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
                            ss = ['Mag-bin_',AtlNames{atl}];
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
                    cbpos(1) = ax.Position(1)+ax.Position(3)+0.11; %0.86;
                    cbpos(2) = ax.Position(2)+ax.Position(4)/3; %cbpos(2) + 0.2; %2.6*cbpos(2);
                    cbpos(3) = 0.02;
                    cbpos(4) = ax.Position(4)*0.5; %0.4; %0.6*cbpos(4);
                    cb.Position = cbpos;


                    annotation('rectangle',[cb.Position(1) cb.Position(2)-0.06 cb.Position(3) 0.02],'FaceColor',color_not_sig0);
                    annotation('textbox',[cb.Position(1)+0.02 cb.Position(2)-0.06 0.2 0.02],...
                        'String','Not Significant','FitBoxToText','on','EdgeColor','none',...
                        'VerticalAlignment','middle','FontName','Arial','FontSize',9);
                    annotation('rectangle',[cb.Position(1) cb.Position(2)-0.09 cb.Position(3) 0.02],'FaceColor',color_nocov);
                    annotation('textbox',[cb.Position(1)+0.02 cb.Position(2)-0.09 0.2 0.02],...
                        'String','No Coverage','FitBoxToText','on','EdgeColor','none',...
                        'VerticalAlignment','middle','FontName','Arial','FontSize',9);
                    %-------------------------------------------------



                    % Dendrogram
                    if (iii2 == 1)

                        ax2_offset = 0; %(-1)*(ax1.Position(1)+ax1.Position(3))*0.004; %-0.016;
                        switch (atl)
                            case (1)
                                ax2_offset = -0.016;
                            case (2)
                                ax2_offset = -0.022;
                            case (3)
                                ax2_offset = 0;
                            case (4)
                                ax2_offset = -0.02;
                            case (5)
                                ax2_offset = -0.02;
                            case (6)
                                ax2_offset = -0.003;
                                ax.Position(1) = ax.Position(1)-0.01;
                            case (7)
                                ax2_offset = -0.025;
                                ax.Position(1) = ax.Position(1)-0.0001;
                            case (8)
                                ax2_offset = -0.033;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (9)
                                ax2_offset = -0.003;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (10)
                                ax2_offset = -0.019;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (11)
                                ax2_offset = -0.003;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (12)
                                ax2_offset = -0.010;
                            case (13)
                                ax2_offset = -0.009;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (15)
                                ax2_offset = -0.033;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (16)
                                ax2_offset = -0.015;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (17)
                                ax2_offset = -0.08;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (18)
                                ax2_offset = -0.035;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (19)
                                ax2_offset = -0.05;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (20)
                                ax2_offset = -0.075;
                                ax.Position(1) = ax.Position(1)-0.00001;
                        end

                        %ax2 = axes('Position',[ax1.Position(1)+ax1.Position(3)+ax2_offset,ax1.Position(2),0.1,ax1.Position(4)]);
                        ax2 = axes('Position',[ax.Position(1)+ax.Position(3)+ax2_offset,ax.Position(2),0.1,ax.Position(4)]);
                        hD = dendrogram(Z,0,'Reorder',fliplr(cluster_i),'Orientation','right','ColorThreshold',0); %[a,b,c]



                        for ihD = 1:length(hD)
                            hD(ihD).Color = 0*[1 1 1];
                        end
                        axis tight;
                        ax2.YLim = [0.5 ax2.YLim(2)+0.5];
                        axis off;
                    end


                    % Get info
                    if ((iii2 == 1) && ((iM == 1) ||(iM == 5)) && (atl == 2))
                        Adiag = diag(Adj_plt2);
                        fprintf('[!!] %i of %i in diagonal are nans.\n',sum(isnan(Adiag)),length(Adiag) );
                        % 6 of 31 in diagonal are nan, 25 are covered
                        AdiagC = Adiag(~isnan(Adiag));
                        fprintf('[!!] %i of %i covered are significant\n',sum(AdiagC~=0),length(AdiagC~=0));
                        % 17 of 25 are significant (68%)

                        AM = Adj_plt2(cluster_i,cluster_i);
                        AM(AM==0) = NaN;
                        [min_vals,min_idxs] = min(AM);
                        [minv,mini] = min(min_vals);
                        fprintf('[!!] minimum at: %i, %i, coherence = %.2f\n',mini,min_idxs(mini),minv);
                        fprintf('\t%s, %s\n',rois_plt2{mini},rois_plt2{min_idxs(mini)})

                        AM = Adj_plt2(cluster_i,cluster_i);
                        AM(AM==0) = NaN;
                        [max_vals,max_idxs] = max(AM);
                        [maxv,maxi] = max(max_vals);
                        fprintf('[!!] maximum at: %i, %i, coherence = %.2f\n',maxi,max_idxs(maxi),maxv);
                        fprintf('\t%s, %s\n',rois_plt2{maxi},rois_plt2{max_idxs(maxi)})

                        AM2 = Adj_plt2(cluster_i,cluster_i);
                        [a,b] = max(sum(~isnan(AM)));
                        fprintf('[!!] Max region: %s, total pairs sig: %.2f\n',rois_plt2{b},a);
                        am = AM2(b,:);
                        am = am(~isnan(am));
                        fprintf('\t%i of %i sig\n',sum(am~=0),length(am));
                        %return
                    end


                    %export_fig(sprintf('figures/T14/Adj_%s_%s_2',metric,ss),'-eps');
                    if (iii2 == 2)
                        ss = [ss,'_var'];
                    end
                    %
                    print(h,sprintf('figures/T14_allatl_resampTv2/atl%i_%s/Adj_%s_%s_rep%i',atl,AtlNames{atl},metric,ss,repL),fig_fmt);
                    if (trig_eps)
                        print(h,sprintf('figures/T14_allatl_resampTv2/atl%i_%s/Adj_%s_%s_rep%i',atl,AtlNames{atl},metric,ss,repL),'-depsc');
                    end

                    % save to web figure
%                     if (iii2 == 1)
%                         print(h,sprintf('/home/jerry/Nextcloud2/Wangetal_PutativeInteractome1_nc/_figures_web/Figure_W3/Figure_W3_atl%i_%s_Adj_%i_%s',atl,AtlNames{atl},iM,metric(3:end)),'-dpng','-r300');
%                         print(h,sprintf('/home/jerry/Nextcloud2/Wangetal_PutativeInteractome1_nc/_figures_web/Figure_W3/Figure_W3_atl%i_%s_Adj_%i_%s',atl,AtlNames{atl},iM,metric(3:end)),'-depsc','-r300');
% 
%                         % Load atlas maps
%                         AC = load('./cache/custom_parcellation_atlas_overlap.mat');
% 
%                         %return
%                         % Save for json export
%                         S = struct();
%                         S.Img = Im;
%                         S.A = Adj_plt2(cluster_i,cluster_i);
%                         S.dendro_Z = Z;
%                         S.dendro_reorder = fliplr(cluster_i);
%                         S.labels = rois_plt_sav;
%                         S.atl_name = AtlNames{atl};
%                         S.atl_chans = AC.AChans{atl}; %1:length(Adj_plt2);
%                         S.atl_chansR = AC.AChansR{atl};
%                         save(sprintf('./brainexport/controls_data_sub-0_freq-%i_atl-%i',iM,atl),'S');
% 
%                     end

                    %return
                    close(h);



                    % Dendrogram
    %                 if (iii2 == 1)
    %                     h = figure('visible','off');
    %                     set(h,'PaperUnits','Inches');
    %                     set(h,'PaperPosition',[0 0 1.5 5]);
    %     %                 cluster_i2 = 1:length(cluster_i);
    %     %                 cluster_i2 = cluster_i2(cluster_i);
    %     %                 cluster_i2 = cluster_i2(cluster_i_manual);
    %                     hD = dendrogram(Z,0,'Reorder',fliplr(cluster_i),'Orientation','right','ColorThreshold',0); %[a,b,c]
    %                     for ihD = 1:length(hD)
    %                         hD(ihD).Color = 0*[1 1 1];
    %                     end
    %                     axis off;
    %                     print(h,sprintf('figures/T14_allatl_resampTv2/atl%i_%s/Adj_%s_%s_dendro',atl,AtlNames{atl},metric,ss),fig_fmt);
    %                     if (trig_eps)
    %                         print(h,sprintf('figures/T14_allatl_resampTv2/atl%i_%s/Adj_%s_%s_dendro',atl,AtlNames{atl},metric,ss),'-depsc');
    %                     end
    %                     close(h);
    %                 end


                    % print info
                    if (iii2 == 1)
                        frac_nocov = sum(sum(isnan(Adj_plt2)))/numel(Adj_plt2);
                        frac_nosig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 == 0)))/numel(Adj_plt2);
                        frac_sig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 ~= 0)))/numel(Adj_plt2);
                        fprintf('[*] saved file: %s\n',sprintf('figures/T14_allatl_resampTv2/atl%i_%s/Adj_%s_%s',atl,AtlNames{atl},metric,ss));

                        if (false)
                            fprintf('[*] no coverage: %.6f (%i of %i)\n',frac_nocov,sum(sum(isnan(Adj_plt2))),numel(Adj_plt2));
                            n_nosig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 == 0)));
                            fprintf('[*] not significant: %.6f (%i of %i)\n',frac_nosig,n_nosig,numel(Adj_plt2));
                            n_sig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 ~= 0)));
                            fprintf('[*] significant: %.6f (%i of %i)\n',frac_sig,n_sig,numel(Adj_plt2));
                            fprintf('[*] significance fraction: %.6f (%i of %i)\n',100*frac_sig/(frac_sig + frac_nosig),n_sig,(n_sig + n_nosig));
                        end

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

                        % Save values from mean

                        Adj_plt_mu = Adj_plt;
                        Adj_plt_mu_vec = Adj_plt_vec;
                    end

                    % GAMMA:
                    % [!] Significant in : 150 of 378 (39.68%), total possible: 465
                end
                %##################################################################
                % End mean/variance loop
                %return

                if ((iii == 3) && trig_plot_mag_cov)
                    n_roisO = sum(cov_idx);
                    n_comb_atl = nchoosek(n_roisO,2);%Ca_hum.n_comb_atl;
                    Mag = nan(n_comb_atl,1);
                    MagVar = nan(n_comb_atl,1);
                    CP = nan(n_comb_atl,1);
                    Npairs = nan(n_comb_atl,1);
                    Nsubs = nan(n_comb_atl,1);
                    %DistsAtl = [];%nan(n_comb_atl,1);

                    % Filter for known, covered rois
                    AdjMagO = AdjMag(known_idx,known_idx);
                    AdjMagVarO = AdjMagVar(known_idx,known_idx);
                    AdjCPO = AdjCP(known_idx,known_idx);
                    AAO = Ca_hum.AdjAtl_sid(known_idx,known_idx);
                    AA2 = Ca_hum.AdjAtl(known_idx,known_idx);
                    AdjMagO = AdjMagO(cov_idx,cov_idx);
                    AdjMagVarO = AdjMagVarO(cov_idx,cov_idx);
                    AdjCPO = AdjCPO(cov_idx,cov_idx);
                    AAO = AAO(cov_idx,cov_idx);
                    AA2 = AA2(cov_idx,cov_idx);

                    count = 1;
                    for i1 = 1:(n_roisO-1)
                        for i2 = (i1+1):n_roisO


                            AA = AAO{i1,i2}; %Ca_hum.AdjAtl_sid{i1,i2};
                            Aa = AA2{i1,i2};
                            Mag(count) = AdjMagO(i1,i2);
                            %Mag(count) = Adj_plt(i1,i2);
                            MagVar(count) = AdjMagVarO(i1,i2);
                            CP(count) = AdjCPO(i1,i2);
                            Npairs(count) = length(AA);
                            Nsubs(count) = length(unique(AA));
                            %DistsAtl(count) = mean(adjct_dist{i1,i2});
                            %AA_dist = adjct_dist{i1,i2};

    %                         if (Npairs(count) == 10) && (Nsubs(count) == 4)
    %                             return
    %                         end

                            count = count + 1;
                        end
                    end

                    msize = 3;


                    % Number of pairs vs magnitude
                    h = figure('visible','off');
                    set(h,'Position',round(fig_size_scale*[0 0 1080 1080]))
                    set(h,'PaperUnits','inches')
                    set(h,'PaperPosition',[0 0 4 4])

                    % Correlation condition
                    %is_sig_mag = (~isnan(Mag));
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

                    %[Mag(is_sig_mag),MagVar(is_sig_mag),Npairs(is_sig_mag)/10,Nsubs(is_sig_mag)]

    %                 % Display variance
    %                 V_std = nan(nchoosek(length(AdjMagVar),2),1);
    %                 c = 1;
    %                 for i1 = 1:(length(AdjMagVar))
    %                     for i2 = (i1+1):length(AdjMagVar)
    %                         V_std(c) = sqrt(AdjMagVar(i1,i2));
    %                         c = c + 1;
    %                     end
    %                 end
    %                 pass_idx = (V_std ~= 0 ) & (~isnan(V_std));
    %                 V_std = V_std(pass_idx);
    %                 fprintf('[*] Coherence Stdev mean: %.2f\n',mean(V_std));
    %                 fprintf('[*] Coherence Stdev median: %.2f\n',median(V_std));
    %                 fprintf('[*] Coherence Stdev min: %.2f\n',min(V_std));
    %                 fprintf('[*] Coherence Stdev max: %.2f\n',max(V_std));
    %                 fprintf('[*] Coherence Stdev n: %.2f\n',numel(V_std));

                    subplot(2,1,1);
                    plot(Npairs(is_sig_mag),Mag(is_sig_mag),'black.','MarkerSize',msize);
                    xlabel('Number of bipolar pairs');
                    ylabel(sprintf('%s Coherence',metric(3:end)));
                    axis tight;
                    box off;
                    set(gca,'TickDir','out');
                    set(gca,'FontSize',10);
                    hold all;
                    try
                        [corr_val_mag, corr_pval_mag] = corr(Npairs(is_sig_mag),Mag(is_sig_mag));                    
                        fprintf('[%s] corr number of pairs, mag: %.5f (p = %.5f), n=%i\n',metric(3:end),corr_val_mag,corr_pval_mag,length(Mag(is_sig_mag)));
                    catch
                        corr_val_mag = NaN;
                        corr_pval_mag = NaN;
                    end
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
                    try
                        [corr_val_cp, corr_pval_cp] = corr(Npairs(is_sig_mag),CP(is_sig_mag));
                    catch
                        corr_val_cp = NaN;
                        corr_pval_cp = NaN;
                    end
                    fprintf('[%s] corr number of pairs, cp: %.5f (p = %.5f), n=%i\n',metric(3:end),corr_val_cp,corr_pval_cp,length(CP(is_sig_mag)));

                    ofname = sprintf('figures/T14_allatl_resampTv2/atl%i_%s/Npairs_%s_mag_%i_p_%d_cp_%i_p_%d_rep%i',atl,AtlNames{atl},metric,...
                        round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp,repL);
                    ofname = replace(ofname,'.','p');
                    print(h,ofname,'-depsc');

                    close(h);

                    % Number of unique patients vs magnitude
                    h = figure('visible','off');
                    set(h,'Position',round(fig_size_scale*[0 0 1080 1080]))
                    set(h,'PaperUnits','inches')
                    set(h,'PaperPosition',[0 0 4 4])



                    subplot(2,1,1);
                    plot(Nsubs(is_sig_mag),Mag(is_sig_mag),'black.','MarkerSize',msize);
                    xlabel('Number of Subjects');
                    ylabel(sprintf('%s Coherence',metric(3:end)));
                    axis tight;
                    box off;
                    set(gca,'TickDir','out');
                    set(gca,'FontSize',10);
                    hold all;
                    try
                        [corr_val_mag, corr_pval_mag] = corr(Nsubs(is_sig_mag),Mag(is_sig_mag));
                    catch
                        corr_val_mag = NaN;
                        corr_pval_mag = NaN;
                    end
                    fprintf('[%s] corr number of patients, mag: %.5f (p = %.5f), n=%i\n',metric(3:end),corr_val_mag,corr_pval_mag,length(Mag(is_sig_mag)));

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
                    try
                        [corr_val_cp, corr_pval_cp] = corr(Nsubs(is_sig_mag),CP(is_sig_mag));
                    catch
                        corr_val_cp = NaN;
                        corr_pval_cp = NaN;
                    end
                    fprintf('[%s] corr number of patients, cp: %.5f (p = %.5f), n=%i\n',metric(3:end),corr_val_cp,corr_pval_cp,length(CP(is_sig_mag)));
                    ofname = sprintf('figures/T14_allatl_resampTv2/atl%i_%s/Nsubs_%s_mag_%i_p_%d_cp_%i_p_%d_rep%i',atl,AtlNames{atl},metric,...
                        round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp,repL);

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

                    %is_sig_mag = ((Mag ~= 0)&(~isnan(Mag)));
                    %is_sig_mag = ((MagVar ~= 0)&(~isnan(MagVar)));

                    subplot(3,1,1);
                    plot(Npairs(is_sig_mag),MagVar(is_sig_mag),'black.');
                    xlabel('Number of bipolar pairs');
                    ylabel(sprintf('%s Coherence Variance',metric(3:end)));
                    axis tight;
                    box off;
                    set(gca,'TickDir','out');
                    hold all;
                    try
                        [corr_val_mag, corr_pval_mag] = corr(Npairs(is_sig_mag),MagVar(is_sig_mag));
                    catch
                        corr_val_mag = NaN;
                        corr_pval_mag = NaN;
                    end
                    fprintf('[%s] corr number of pairs, magVar: %.5f (p = %.5f), n=%i\n',metric(3:end),corr_val_mag,corr_pval_mag,length(MagVar(is_sig_mag)));

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
                    try
                        [corr_val_cp, corr_pval_cp] = corr(Npairs(is_sig_mag),CP(is_sig_mag));
                    catch
                        corr_val_mag = NaN;
                        corr_pval_mag = NaN;
                    end
                    fprintf('[%s] corr number of pairs, cp: %.5f (p = %.5f), n=%i\n',metric(3:end),corr_val_cp,corr_pval_cp,length(CP(is_sig_mag)));
                    ofname = sprintf('figures/T14_allatl_resampTv2/atl%i_%s/Npairs_%s_magVar_%i_p_%d_cp_%i_p_%d_rep%i',...
                        atl,AtlNames{atl},metric,round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp,repL);
                    ofname = replace(ofname,'.','p');
                    %return
                    print(h,ofname,'-depsc');

                    close(h);

                    % Number of unique patients vs magnitude
                    h = figure('visible','off');
                    set(h,'Position',round(fig_size_scale*[0 0 0.95*1080 0.8*1080]))
                    %is_sig_mag = ((MagVar ~= 0)&(~isnan(MagVar)));
                    subplot(3,1,1);
                    plot(Nsubs(is_sig_mag),MagVar(is_sig_mag),'black.');
                    xlabel('Number of patients');
                    ylabel(sprintf('%s Coherence Variance',metric(3:end)));
                    axis tight;
                    box off;
                    set(gca,'TickDir','out');
                    hold all;
                    try
                        [corr_val_mag, corr_pval_mag] = corr(Nsubs(is_sig_mag),MagVar(is_sig_mag));
                    catch
                        corr_val_mag = NaN;
                        corr_pval_mag = NaN;
                    end
                    fprintf('[%s] corr number of patients, magVar: %.5f (p = %.5f), n=%i\n',metric(3:end),corr_val_mag,corr_pval_mag,length(MagVar(is_sig_mag)));

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
                    try
                        [corr_val_cp, corr_pval_cp] = corr(Nsubs(is_sig_mag),CP(is_sig_mag));
                    catch
                        corr_val_cp = NaN;
                        corr_pval_cp = NaN;
                    end
                    fprintf('[%s] corr number of patients, cp: %.5f (p = %.5f), n=%i\n',metric(3:end),corr_val_cp,corr_pval_cp,length(CP(is_sig_mag)));
                    ofname = sprintf('figures/T14_allatl_resampTv2/atl%i_%s/Nsubs_%s_magVar_%i_p_%d_cp_%i_p_%d_rep%i',atl,AtlNames{atl},metric,...
                        round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp,repL);
                    ofname = replace(ofname,'.','p');
                    print(h,ofname,'-depsc');

                    close(h);

                    % --- END REPEAT ------------------------------------------


                    % signal dependent noise
                    h = figure('visible','off');
                    plot(Mag(is_sig_mag),MagVar(is_sig_mag),'black.')
                    try
                        [corr_val_cp, corr_pval_cp] = corr(Mag(is_sig_mag),MagVar(is_sig_mag));
                    catch
                        corr_val_cp = NaN;
                        corr_pval_cp = NaN;
                    end
                    xlabel(sprintf('%s Coherence Mean',metric(3:end)))
                    ylabel(sprintf('%s Coherence Variance',metric(3:end)))
                    box off;
                    set(gca,'TickDir','out');
                    fprintf('[%s] corr coherence mean, variance: %.5f (p = %.5f), n=%i\n',metric(3:end),corr_val_cp,corr_pval_cp,length(MagVar(is_sig_mag)));
                    ofname = sprintf('figures/T14_allatl_resampTv2/atl%i_%s/mag_%s_magVar_%i_p_%d_rep%i',atl,AtlNames{atl},metric,...
                        round(1000*corr_val_cp),corr_pval_cp,repL);
                    ofname = replace(ofname,'.','p');
                    print(h,ofname,'-depsc');

                    close(h);
                    % plot(Npairs(is_sig_mag),MagVar(is_sig_mag),'black.')
                    % plot(Nsubs(is_sig_mag),MagVar(is_sig_mag),'black.')


                    % distance and connection strength
                    h = figure('visible','off');
                    plot(DistsAtl(:,1),DistsAtl(:,2),'black.')
                    try
                        [corr_val_cp, corr_pval_cp] = corr(DistsAtl(:,1),DistsAtl(:,2));
                    catch
                        corr_val_cp = NaN;
                        corr_pval_cp = NaN;
                    end
                    xlabel('Distance (mm)')
                    ylabel(sprintf('%s Coherence',metric(3:end)))
                    box off;
                    set(gca,'TickDir','out');
                    fprintf('[%s] corr dist, mag: %.5f (p = %.5f), n=%i\n',metric(3:end),corr_val_cp,corr_pval_cp,length(DistsAtl(:,2)));
                    axis tight;
                    %return
                    ofname = sprintf('figures/T14_allatl_resampTv2/atl%i_%s/dist_%s_mag_%i_p_%d_rep%i',atl,AtlNames{atl},metric,...
                        round(1000*corr_val_cp),corr_pval_cp,repL);
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
                    [~,p_tt2] = ttest2(DiagMag,TriuMag);
                    fprintf('[*] unpaired t test:\n')
                    fprintf('\tp-val: %.6d\n',p_tt2)
                    % p  = 0.0033, z-val = -2.9403, ranksum = 2574
                    box off;
                    set(gca,'TickDir','out');
                    print(h,sprintf('figures/T14_allatl_resampTv2/atl%i_%s/intra_inter_%s_ranksum-%i_ttest-%i_rep%i',atl,AtlNames{atl},metric,...
                        round(1000*p_rs),round(1000*p_tt2),repL),'-depsc');
                    close(h);

                    %return
                end

            end
        end
        
        %%
        Ca14 = load(sprintf('./cache/figure_t14_%i_atl%i_Desikan-Killiany',iM,atl));
        AdjMagFull = Ca14.AdjMag(known_idx,known_idx);
        AdjMagVarFull = Ca14.AdjMagVar(known_idx,known_idx);
        
        
        % Robustness of the time-averaged coherence
        AdjAA_1 = BpairsRep{5};
        AdjAA_2 = BpairsRep{6};
        vAA_1 = [];
        vAA_2 = [];
        for ir1 = 1:(length(AdjAA_1))
            for ir2 = (ir1):length(AdjAA_1)
                vAA_1 = [vAA_1; AdjAA_1{ir1,ir2}];
                vAA_2 = [vAA_2; AdjAA_2{ir1,ir2}];
            end
        end
%         idx_both = (vAA_1 ~= 0) | (vAA_2 ~= 0);
%         vAA_1 = vAA_1(idx_both);
%         vAA_2 = vAA_2(idx_both);
        
        [p,~,~] = ranksum(vAA_1,vAA_2);
        [~,pt,~] = ttest(vAA_1,vAA_2);
        fprintf('[*] Robustness of time-averaged coherence\n\tRanskum p=%.6d, T p=%.6d, n=%i\n\tDifference mean=%.8f, std=%.8f\n',...
            p,pt,length(vAA_1),mean(vAA_1-vAA_2),std(vAA_1-vAA_2));
        
        % Rep1 : keep second half
        % Rep2 : keep first half
        AdjMagRep1 = AdjMagRep{5};
        AdjMagRep1 = AdjMagRep1(known_idx,known_idx);
        AdjMagRep2 = AdjMagRep{6};
        AdjMagRep2 = AdjMagRep2(known_idx,known_idx);
        
        diff1 = AdjMagFull - AdjMagRep1;
        diff2 = AdjMagFull - AdjMagRep2;
        diff3 = AdjMagRep1 - AdjMagRep2;
        v_d1 = zeros(nchoosek(length(diff1),2),1);
        v_d2 = zeros(nchoosek(length(diff1),2),1);
        v_d3 = zeros(nchoosek(length(diff1),2),1);
        v_A = zeros(nchoosek(length(diff1),2),1);
        v_A1 = zeros(nchoosek(length(diff1),2),1);
        v_A2 = zeros(nchoosek(length(diff1),2),1);
        c = 1;
        for i1 = 1:(length(diff1)-1)
            for i2 = (i1+1):length(diff1)
                v_d1(c) = diff1(i1,i2);
                v_d2(c) = diff2(i1,i2);
                v_d3(c) = diff3(i1,i2);
                v_A(c) = AdjMagFull(i1,i2);
                v_A1(c) = AdjMagRep1(i1,i2);
                v_A2(c) = AdjMagRep2(i1,i2);
                c = c + 1;
            end
        end
        is_sig = (~isnan(v_A)) & (v_A ~= 0);
        v_d1 = v_d1(is_sig);
        v_d2 = v_d2(is_sig);
        v_d3 = v_d3(is_sig);
        v_A = v_A(is_sig);
        v_A1 = v_A1(is_sig);
        v_A2 = v_A2(is_sig);
        
        fprintf('Remove first half of data.\n');
        fprintf('[*] mean diff: %0.6f\n',mean(v_d1));
        fprintf('[*] stdev diff: %0.6f\n',std(v_d1));
        f_d = 1 - (v_A1./v_A);
        fprintf('[*] mean diff: %0.2f%%\n',mean(100*f_d));
        fprintf('[*] stdev diff: %0.2f%%\n',std(100*f_d));
        
        
        ind_both = (v_A1 > 0) & (v_A2 > 0);
        fprintf('Compare first half with second.\n');
        fprintf('[*] mean diff: %0.6f\n',mean(v_d3(ind_both)));
        fprintf('[*] stdev diff: %0.6f\n',std(v_d3(ind_both)));
        f_d = 1 - (v_A1./v_A2);
        f_d(isinf(f_d)) = NaN;
        fprintf('[*] mean diff: %0.2f%%\n',nanmean(100*f_d(ind_both)));
        fprintf('[*] stdev diff: %0.2f%%\n',nanstd(100*f_d(ind_both)));
        fprintf('[*] n: %i\n',length(f_d(ind_both)))
        
        
        [p,~,~] = ranksum(v_A1(ind_both),v_A2(ind_both));
        fprintf('\tRanksum\tp=%.4d\n',p);
        [~,p,~] = ttest(v_A1(ind_both),v_A2(ind_both),'tail','both');
        fprintf('\tT-test\tp=%.4d\n',p);
        fprintf('\tn=%i\n',length(v_A1))
        n_both = sum(ind_both);
        fprintf('\tboth n=%i\n',n_both);
        
%         [p,h,stats] = ranksum(v_A,v_A1);
%         fprintf('p-val: %.1d\n',p);
%         fprintf('Remove second half of data.\n');
%         fprintf('[*] mean diff: %0.6f\n',mean(v_d2));
%         fprintf('[*] stdev diff: %0.6f\n',std(v_d2));
        
        %return

    end

end

% 19 fevr. 2020
% (Gamma)

% Compare first half with second.
% [*] mean diff: 0.002852
% [*] stdev diff: 0.014583
% [*] mean diff: -0.95%
% [*] stdev diff: 4.88%
% [*] n: 179
% 	Ranksum	p=7.8508e-01
% 	T-test	p=9.6347e-03
% 	n=183
% 	both n=179


% 18 fevr. 2020

% Compare first half with second.
% [*] mean diff: -0.000672
% [*] stdev diff: 0.025493
% [*] mean diff: -0.04%
% [*] stdev diff: 8.13%
% [*] n: 186
% 	Ranksum	p=8.8309e-01
% 	T-test	p=7.1961e-01
% 	n=193
% 	both n=186