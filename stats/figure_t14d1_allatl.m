close all;
clc;
rng shuffle;

n_perm = 10000;
perm_alpha_sec = 12;
cp_thresh_override = 0.05;
n_pairs_thresh = 0; % 20 at least this many electrode pairs to be considered
n_subs_thresh = 0; % 4 at least this many subjects to be considered
font_size = 10;

fig_fmt = '-depsc'; %'-depsc';
trig_write_figs = true; % whether to print figs or not
trig_eps = true;
trig_mag_no_cp = true;
system('mkdir figures');
system('mkdir figures/T14d1');

% Fast i/o definitions
%dir_artL = '/media/klab/internal/data/h5_notch20/art';
%dir_resL = '/media/klab/internal/data/results/coh_w10';
%dir_corL = '/media/klab/internal/data/coreg';

dir_artL = '/media/klab/internal/data/h5_notch20/art';
dir_resL = '/media/klab/internal/data/results/coh_w10';
dir_corL = '/media/klab/internal/data/coreg';
dir_cacheL = './cache';
    
%dir_cacheL = './cache';
subjects_dirL = '/mnt/cuenap_ssd/coregistration';

SubjectsL = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124',...
    'mSu'};



% SubjectsL = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
%     'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
%     'm00035','m00037','m00038','m00039','m00043','m00044'         ,'m00047',... % ,'m00045'
%     'm00048','m00049','m00052','m00053','m00055',         'm00058','m00059',... % ,'m00056'
%     'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
%     'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124',...
%     'mSu'};

% Exclude monkey
SubjectsL = SubjectsL(1:(end-1));

%SubjectsL = {'m00001'};

% Slow i/o definitions
%dir_h5L = '/media/klab/KLAB101/h5_notch20';

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'}; % 'pcDelta'
metrics_suffix = {'0.5-125 Hz','3-8 Hz','8-12 Hz','12-30 Hz','30-100 Hz'};

for iM = 1:length(metrics)
    for iSub = 1:length(SubjectsL)
        
        for atlM = 1:20 %1:20
        
            sid = SubjectsL{iSub};
            sid_int = iSub;

            metric = metrics{iM};

            % Load human cache
            sid_sav = sid;
            iSub_sav = iSub;
            %load(sprintf('%s/xsub_out_all_%i.mat',dir_cacheL,iM));
            load(sprintf('%s/xsub_out_all_%i_atl%i.mat',dir_cacheL,iM,atlM));
            sid = sid_sav;
            iSub = iSub_sav;

            % Calculate final functional interaction matrix
            Adj = nan(n_rois,n_rois);
            AdjMag = nan(n_rois,n_rois);
            AdjVar = nan(n_rois,n_rois);
            AdjNpairs = nan(n_rois,n_rois);
            AdjNpairs_sig = nan(n_rois,n_rois);
            AdjMag4cl = nan(n_rois,n_rois);
            AdjCP = nan(n_rois,n_rois);
            N_bchan = nan(n_rois,n_rois);
            cp_thresh = cp_thresh_override;
            Dmat = Inf(n_rois,n_rois); % set to inf to avoid removing from every instance

            for i1 = 1:n_rois
                for i2 = 1:n_rois
                    AA = AdjAtl{i1,i2};
                    AA_dist = adjct_dist{i1,i2};
                    AA_sub = AdjAtl_sid{i1,i2};
                    %n_pairs = length(AA_sub);
                    n_pairs = length(AA(AA_sub == iSub));
                    n_subs = length(unique(AA_sub));

                    % index by subject
                    AA = AA(AdjAtl_sid{i1,i2}==iSub);

                    % ROI pair coverage condition (grey)
                    if ( (~ isempty(AA))  && (n_pairs >= n_pairs_thresh) && (n_subs >= n_subs_thresh))
                    %if (~ isempty(AA))
                        N_bchan(i1,i2) = length(AA);
                        frac_cp = sum(AA ~= 0)/length(AA);
                        AdjCP(i1,i2) = frac_cp;
                        % if ( (frac_cp > cp_thresh) && (n_subs_ct >= n_subs_ct_thresh) )
                        if (frac_cp > cp_thresh)
                            % ROI pair significance condition (white)
                            Adj(i1,i2) = 1;
                            AdjMag(i1,i2) = mean(AA(AA ~= 0));
                            AdjVar(i1,i2) = var(AA(AA ~= 0));
                            AdjNpairs(i1,i2) = n_pairs;
                            AdjNpairs_sig(i1,i2) = length(AA(AA ~= 0));
                        else
                            Adj(i1,i2) = 0;
                            AdjMag(i1,i2) = 0;
                            AdjVar(i1,i2) = 0;
                        end
                    end

                    AdjMag4cl(i1,i2) = mean(AA(AA~=0));
                end
            end

            frac_xsub_msu = sum(Adj(~isnan(Adj))) /numel(Adj(~isnan(Adj)));
            fprintf('%s - human fraction of ROIs significant: %.4f\n',metric,frac_xsub_msu)

            % --- PLOT  --------------------------------------------------
            colordef_adj;
            color_nocov = COLOR_ADJ_NOCOV;
            color_distt = color_nocov;
            color_isnan = color_nocov;
            color_not_sig0 = COLOR_ADJ_NOSIG;

            %color_distt = 0.6*[1 1 1];
            %color_isnan = 0.6*[1 1 1];
            fontsz = font_size; %10;
            %fsize = fontsz;

            % Calculate rows of nans
            Adj2 = Adj;
            ind_isnan = false(length(Adj2),1);
            for i = 1:length(ind_isnan)
                ind_isnan(i) = (sum(isnan(Adj2(i,:))) == length(ind_isnan));
            end

            for iii = 3 %1:3

                % Make figure
                h = figure('visible','off');
                %set(h,'Position',round(1*[0 0 0.95*1080 0.8*1080]));
                set(h,'PaperUnits','Inches');
                set(h,'PaperPosition',[0 0 8.5 6.8]);
                %subplot(2,1,1);

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
                    %Adj_plt = AdjCP;
                    %Adj_plt(Adj == 0) = 0;
                    color_not_sig = color_not_sig0; %0.999*[1 1 1];
                end
                Dmat_plt = Dmat;
                rois_plt = rois;

                % for hierarchichal clustering
                Adj_plt2_cl = AdjMag4cl;
                Adj_plt_var = AdjVar;
                Adj_plt_npairs = AdjNpairs;
                Adj_plt_npairs_sig = AdjNpairs_sig;
                
                % remove unknown areas
                cond_empty = false(length(rois),1);
                for ce = 1:length(cond_empty)
                    elem = rois{ce};
                    cond_empty(ce) = (all(isspace(elem)) | isempty(elem));
                end
                cond_contains = (contains(lower(rois),'unknown')) | (contains(lower(rois),'???')) | (contains(lower(rois),'wall')) | (cond_empty);
                known_idx = (~ cond_contains);
                Adj_plt = Adj_plt(known_idx,known_idx);
                Adj_plt2_cl = Adj_plt2_cl(known_idx,known_idx);
                Adj_plt_var = Adj_plt_var(known_idx,known_idx);
                Adj_plt_npairs = Adj_plt_npairs(known_idx,known_idx);
                Adj_plt_npairs_sig = Adj_plt_npairs_sig(known_idx,known_idx);
                rois_plt = rois_plt(known_idx);
                Dmat_plt = Dmat_plt(known_idx,known_idx);


                % Filter out nans nodes
                cov_idx = ~ all(isnan(Adj_plt));
                Adj_plt = Adj_plt(cov_idx,cov_idx);
                Adj_plt2_cl = Adj_plt2_cl(cov_idx,cov_idx);
                Adj_plt_var = Adj_plt_var(cov_idx,cov_idx);
                Adj_plt_npairs = Adj_plt_npairs(cov_idx,cov_idx);
                Adj_plt_npairs_sig = Adj_plt_npairs_sig(cov_idx,cov_idx);
                rois_plt = rois_plt(cov_idx);
                Dmat_plt = Dmat_plt(cov_idx,cov_idx);

                if ((length(rois_plt) <= 1) || (nanmin(Adj_plt(:)) == nanmax(Adj_plt(:))))
                    fprintf(2,'[!] No matrix to plot, skipped:\n');
                    fprintf(2,'\t%s',sprintf('./brainexport/controls_data_sub-%i_freq-%i_atl-%i\n',iSub,iM,atl))
                else

%                     for i = 1:length(rois_plt)
%                         rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '_','/');
%                     end
                    [n_roi_p,~] = size(Adj_plt);
                    Adj_plt_sav = Adj_plt;
                    % generate image from matrix
                    Adj_plt(Dmat_plt <= dist_thresh) = nan; %nanmean(Adj_plt(:));
                    Adj_plt2_cl(Dmat_plt <= dist_thresh) = nan;
                    Adj_plt_var(Dmat_plt <= dist_thresh) = nan;
                    Adj_plt_npairs(Dmat_plt <= dist_thresh) = nan;
                    Adj_plt_npairs_sig(Dmat_plt <= dist_thresh) = nan;
                    v = Adj_plt;
                    [v_n, v_m] = size(v);
                    %map = corrcmap(100);
                    map = inferno(100);
                    if (iii == 3)
                        minv = min(v(v~=0));
                        maxv = max(v(v~=0));
                    else
                        minv = min(v(:));
                        maxv = max(v(:));
                    end
                    %minv = min(v(:));
                    %maxv = max(v(:));
                    ncol = size(map,1);
                    % breaks for iSub=18, iM=1, atl=4
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
            %         Im = Im(~ind_isnan,~ind_isnan,:); % remove rows of nans
            %         rois_plt = rois_plt(~ind_isnan);
            %         Adj_plt2_cl = Adj_plt2_cl(~ind_isnan,~ind_isnan);

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
                            %ad = adjct_dist_cl{m2,m1};
                            ad = Adj_plt(m1,m2);
                            if (isnan(ad))
                                Y(yc) = 0; %Inf;%Inf; 600;
                            else
                                Y(yc) = ad; %mean(ad);
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

                    %Z = linkage(Y); %,'centroid'
                    M = Adj_plt2_cl;
                    M = (M + M.')/2;
                    Z = linkage(M,'ward',@nanDist4clustering);

                    try
                        cluster_i = optimalleaforder(Z,Y); %,'transformation','inverse'); % ,'transformation','inverse'
                    catch
                        % Fails for iSub=6, atl=7, iM=3
                        cluster_i = 1:length(Adj_plt);
                    end
                    roi_dist = roi_dist(cluster_i,cluster_i);
                    %roi_dist(isinf(roi_dist)) = ;
                    clash = nansum(nansum(triu(roi_dist,1) - triu(roi_dist,2)));
                    fprintf('cluster clash: %.12f mm\n',clash)
                    %cluster_i = optimalleaforder(Z,Y);

                    % -------------------------------------------------------------------------

                    % Apply clustering
                    Im = Im(cluster_i,cluster_i,:);
                    rois_plt = rois_plt(cluster_i);
                    rois_plt_sav = rois_plt;

                    imagesc(Im);
                    daspect([1 1 1]);
                    set(gca,'FontSize', font_size);
                    set(gca,'FontName','Arial');
                    %imagesc(Im, 'Parent', ax);
            %         colormap(map);
            %         colorbar;
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

                    % Pad ROI names
                    for i = 1:length(rois_plt)
                        %rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '_','/');
                        rois_plt{i} = convertRoiDK(rois_plt{i});
                    end

            %         maxlen = 0;
            %         for jj = 1:length(rois_plt)
            %             l = length(rois_plt{jj});
            %             if (l > maxlen)
            %                 maxlen = l;
            %             end
            %         end
            %         for jj = 1:length(rois_plt)
            %             l = length(rois_plt{jj});
            %             sat = (l < maxlen);
            %             while (sat)
            %                 rois_plt{jj} = [' ',rois_plt{jj}];
            %                 l = length(rois_plt{jj});
            %                 sat = (l < maxlen);
            %             end
            %         end


                    % numbering
                    ax = gca;
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

                    if (atlM ~= 2)
                       rois_plt = rois_plt_sav; 
                    end
                    
                    yticks(1:length(rois_plt));
                    yticklabels(rois_plt);
                    xticks(1:length(rois_plt));
                    xticklabels(rois_plt);
                    xtickangle(90);
                    set(gca,'tickdir','out');
                    set(gca,'fontsize',fontsz);
                    set(gca,'TickLength',[0.001, 0.001])
                    daspect([1 1 1]);
                    set(gca,'FontName','Arial');

                    % Move axis if needed
                    ax = gca;
                    tp = ax.Position;
                    ax1_xoffset = 0.05;
                    tp(1) = tp(1) - ax1_xoffset;
                    ax.Position = tp;
                    if(ax.FontSize > 10)
                        ax.FontSize = 10;
                    end

            %         if (strcmp(metric,'pcBroadband'))
            %             title(sprintf('Functional Interactions - %s',sid))
            %         else
            %             title(sprintf('Functsional Interactions - %s - %s',sid,metric(3:end)))
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
                            ss = 'Mag';
                        else
                            ss = 'CP';
                        end
                    end


                    % New colormap handling
                    %-----------------------------------------
                    colormap(map);
                    if (((~isnan(minv)) && (~isnan(maxv))) && (minv ~= maxv))
                        cytick = linspace(minv,maxv,3);
                        cytick_str = cell(1,length(cytick));
                        for kk = 1:length(cytick)
                            cytick_str{kk} = sprintf('%.2f',cytick(kk));
                        end
                        cb = colorbar('ytick',cytick,'yticklabel',cytick_str,'FontSize',font_size);
                        caxis([minv maxv]);
                    else
                        minv = minv - 0.0001;
                        maxv = maxv + 0.0001;
                        cytick = linspace(minv,maxv,2);
                        cytick_str = cell(1,length(cytick));
                        for kk = 1:length(cytick)
                            cytick_str{kk} = sprintf('%.2f',cytick(kk));
                        end
                        if (any(isnan(cytick)))
                            cb = colorbar('ytick',[0 1],'yticklabel',{'NaN','NaN'},'FontSize',font_size);
                            caxis([0 1]);
                        else
                            cb = colorbar('ytick',cytick,'yticklabel',cytick_str,'FontSize',font_size);
                            caxis([minv maxv]);
                        end
                    end
                    set(cb,'TickDir','out');
                    set(cb,'TickLength',0);

                    %--- display frequency band range ----
                    metricTxt = metric(3:end);
                    iii2 = 1;
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

                    %ylabel(cb,sprintf('%s Coherence',metric(3:end)));

                    colormap(map);
                    cbpos = cb.Position;
            %         if (strcmp(sid,'m00124') || strcmp(sid,'m00024') || strcmp(sid,'m00038') || strcmp(sid,'m00048') || strcmp(sid,'m00061')|| strcmp(sid,'m00095'))
            %             cbpos(2) = 2*cbpos(2);
            %             cbpos(4) = 0.5*cbpos(4);
            %         else
            %             cbpos(2) = 2.6*cbpos(2);
            %             cbpos(4) = 0.6*cbpos(4);
            %         end
                    cbpos(1) = ax.Position(1)+ax.Position(3)+0.11; %0.86;
                    cbpos(2) = ax.Position(2)+ax.Position(4)/3; %cbpos(2) + 0.2; %2.6*cbpos(2);
                    cbpos(3) = 0.02;
                    cbpos(4) = ax.Position(4)*0.5; %0.4; %0.6*cbpos(4);

                    cb.Position = cbpos;

                    annotation('rectangle',[cb.Position(1) cb.Position(2)-0.06 cb.Position(3) 0.02],'FaceColor',color_not_sig0);
                    annotation('textbox',[cb.Position(1)+0.02 cb.Position(2)-0.06 0.2 0.02],'String','Not Significant','FitBoxToText','on','EdgeColor','none','VerticalAlignment','middle','FontName','Arial','FontSize',9);
                    annotation('rectangle',[cb.Position(1) cb.Position(2)-0.09 cb.Position(3) 0.02],'FaceColor',color_nocov);
                    annotation('textbox',[cb.Position(1)+0.02 cb.Position(2)-0.09 0.2 0.02],'String','No Coverage','FitBoxToText','on','EdgeColor','none','VerticalAlignment','middle','FontName','Arial','FontSize',9);
                    %-------------------------------------------------


                    % Dendrogram
                    if (iii2 == 1)

                        ax2_offset = -0.02; %(-1)*(ax1.Position(1)+ax1.Position(3))*0.004; %-0.016;
                        %ax.Position(1) = ax.Position(1)-0.01;
            %             switch (atl)
            %                 case (1)
            %                     ax2_offset = -0.016;
            %                 case (2)
            %                     ax2_offset = -0.022;
            %                 case (3)
            %                     ax2_offset = 0;
            %                 case (4)
            %                     ax2_offset = -0.02;
            %                 case (5)
            %                     ax2_offset = -0.02;
            %                 case (6)
            %                     ax2_offset = -0.003;
            %                     ax.Position(1) = ax.Position(1)-0.01;
            %                 case (7)
            %                     ax2_offset = -0.025;
            %                     ax.Position(1) = ax.Position(1)-0.0001;
            %                 case (8)
            %                     ax2_offset = -0.033;
            %                     ax.Position(1) = ax.Position(1)-0.00001;
            %                 case (9)
            %                     ax2_offset = -0.003;
            %                     ax.Position(1) = ax.Position(1)-0.00001;
            %                 case (10)
            %                     ax2_offset = -0.019;
            %                     ax.Position(1) = ax.Position(1)-0.00001;
            %                 case (11)
            %                     ax2_offset = -0.003;
            %                     ax.Position(1) = ax.Position(1)-0.00001;
            %                 case (12)
            %                     ax2_offset = -0.010;
            %                 case (13)
            %                     ax2_offset = -0.009;
            %                     ax.Position(1) = ax.Position(1)-0.00001;
            %                 case (15)
            %                     ax2_offset = -0.033;
            %                     ax.Position(1) = ax.Position(1)-0.00001;
            %                 case (16)
            %                     ax2_offset = -0.015;
            %                     ax.Position(1) = ax.Position(1)-0.00001;
            %                 case (17)
            %                     ax2_offset = -0.08;
            %                     ax.Position(1) = ax.Position(1)-0.00001;
            %                 case (18)
            %                     ax2_offset = -0.035;
            %                     ax.Position(1) = ax.Position(1)-0.00001;
            %                 case (19)
            %                     ax2_offset = -0.05;
            %                     ax.Position(1) = ax.Position(1)-0.00001;
            %                 case (20)
            %                     ax2_offset = -0.075;
            %                     ax.Position(1) = ax.Position(1)-0.00001;
            %             end

                        %ax2 = axes('Position',[ax1.Position(1)+ax1.Position(3)+ax2_offset,ax1.Position(2),0.1,ax1.Position(4)]);
                        ax2 = axes('Position',[ax.Position(1)+ax.Position(3)+ax2_offset,ax.Position(2),0.1,ax.Position(4)]);

                        try
                            hD = dendrogram(Z,0,'Reorder',fliplr(cluster_i),'Orientation','right','ColorThreshold',0); %[a,b,c]
                            for ihD = 1:length(hD)
                                hD(ihD).Color = 0*[1 1 1];
                            end
                            cond_dendro = true;
                        catch
                            cond_dendro = false;
                        end
                        axis tight;
                        ax2.YLim = [0.5 ax2.YLim(2)+0.5];
                        axis off;
                    end

                    % calculate fraction sig
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

                    % print info
                    Adj_plt = Adj_plt_sav;
                    frac_nocov = sum(sum(isnan(Adj_plt)))/numel(Adj_plt);
                    frac_nosig = sum(sum(~isnan(Adj_plt) & (Adj_plt == 0)))/numel(Adj_plt);
                    frac_sig = sum(sum(~isnan(Adj_plt) & (Adj_plt ~= 0)))/numel(Adj_plt);
                    fprintf('[*] saved file: %s\n',sprintf('figures/T14/Adj_%s_%s',metric,ss));
                    fprintf('[*] no coverage: %.6f (%i of %i)\n',frac_nocov,sum(sum(isnan(Adj_plt))),numel(Adj_plt));
                    n_nosig = sum(sum(~isnan(Adj_plt) & (Adj_plt == 0)));
                    fprintf('[*] not significant: %.6f (%i of %i)\n',frac_nosig,n_nosig,numel(Adj_plt));
                    n_sig = sum(sum(~isnan(Adj_plt) & (Adj_plt ~= 0)));
                    fprintf('[*] significant: %.6f (%i of %i)\n',frac_sig,n_sig,numel(Adj_plt));
                    fprintf('[*] significance fraction: %.6f (%i of %i)\n',100*frac_sig/(frac_sig + frac_nosig),n_sig,(n_sig + n_nosig));


                    % --- Output legend ---
                    Ao = nan(1,nchoosek(length(Adj_plt),2));
                    ca = 1;
                    for iao = 1:(length(Adj_plt)-1)
                        for jao = (iao+1):length(Adj_plt)
                            Ao(ca) = Adj_plt(iao,jao);
                            ca = ca + 1;
                        end
                    end

                    if (strcmp(metric,'pcBroadband'))
                        mtxt = '';
                        mtxt2 = '';
                    else
                        mtxt = metric(3:end);
                        mtxt2 = [mtxt,' '];
                    end
                    legend_txt = sprintf('Figure W2B%s Significant time-averaged interactions in one subject according to Desikan areas (Desikan et al., 2006). Pairwise time-average %scoherence measurements between cortical areas in subject %i over %.1f days. There were %i bipolar electrodes with coverage across %i region pairs, out of %i total. Of the covered pairs of regions, %i were found to be significant (%.0f%%).',...
                        mtxt,mtxt2,sid_int,ecog.n_samples/(ecog.fs*3600*24),ecog.n_bchan,sum(~isnan(Ao)),nchoosek(ecog.n_bchan,2),sum(Ao(~isnan(Ao)) ~= 0),100*(sum(Ao(~isnan(Ao)) ~= 0))/(sum(~isnan(Ao))));

                    sid = 'mX';
                    ofi = fopen(sprintf('figures/T14d1/Figure_W2-%iB_Adj_%s_%s_%s_atl-%i_sig-%i-of-%i.txt',...
                        sid_int,metric,sid,ss,atlM,n_sig2,(n_sig2 + n_nosig2)),'w');
                    fprintf(ofi,'%s',legend_txt);
                    fclose(ofi);

                    % ---------------------

                    if (trig_write_figs)
                        print(h,sprintf('figures/T14d1/Figure_W2-%iB_Adj_%s_%s_%s_atl-%i_sig-%i-of-%i',...
                            sid_int,metric,sid,ss,atlM,n_sig2,(n_sig2 + n_nosig2)),fig_fmt,'-r100');
                        if (trig_eps)
                            %print(h,sprintf('figures/T14d1/Figure_W2-%iB_Adj_%s_%s_%s',sid_int,metric,sid,ss),'-depsc');
                            print(h,sprintf('figures/T14d1/Figure_W2-%iB_Adj_%s_%s_%s_atl-%i_sig-%i-of-%i',...
                                sid_int,metric,sid,ss,atlM,n_sig2,(n_sig2 + n_nosig2)),'-depsc','-r500');
                        end
                    end

                   % return
                    close(h);



                    if (atl ~= atlM)
                        fprintf(2,'[!] Error, loaded atlas index does not match.\n');
                    end


                    % channel list
                    Cta = load(sprintf('./cache/xsub_out_%s_%i_atl%i',SubjectsL{iSub},iM,atl));
                    atl_chans = nan(1,Cta.ecog.n_bchan);


                    Ct = Cta.C;
                    if (strcmp('R',Ct.EleHemi{1}))
                        reg_dict = Ct.AtlROIs{atl}.RH.struct_names;
                    else
                        reg_dict = Ct.AtlROIs{atl}.LH.struct_names;
                    end

                    bip_label = Ct.AtlLabels{atl}(Cta.ecog.bip(:,1));
                    for iac = 1:length(atl_chans)
                        reg_str = Ct.AtlLabels{atl}{iac};
                        reg_idx = find(strcmp(reg_dict,reg_str));
                        if (~isempty(reg_idx))
                            atl_chans(iac) = reg_idx;
                        end
                    end

                    %return
                    % Bipolar electrode matrix, atl=0
                    %atl = 0;
                    S = struct();
                    S.Img = Im;
                    S.A = Adj_plt(cluster_i,cluster_i);
                    S.Astd = sqrt(Adj_plt_var(cluster_i,cluster_i));
                    S.A_npairs = Adj_plt_npairs(cluster_i,cluster_i);
                    S.A_npairs_sig = Adj_plt_npairs_sig(cluster_i,cluster_i);
                    S.A_min = nanmin(S.A(S.A ~= 0));
                    S.A_max = nanmax(S.A);
                    S.dendro_Z = Z;
                    S.dendro_reorder = fliplr(cluster_i);
                    S.labels = rois_plt_sav;
                    S.atl_name = C.AtlNames{atl};
                    S.atl_chans = atl_chans; %1:length(rois_plt); %1:length(Adj_plt2);
                    %S.atl_chansR = AC.AChansR{atl};
                    save(sprintf('./brainexport/controls_data_sub-%i_freq-%i_atl-%i',iSub,iM,atl),'S');

                end
            end
        end

    end

end
