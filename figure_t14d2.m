close all;
clear;
clc;
rng shuffle;

n_perm = 10000;
perm_alpha_sec = 12;
cp_thresh_override = 0.05;
n_pairs_thresh = 40; %exactly this number is resampled  %5*12 * 4; % at least this many electrode pairs to be considered
n_subs_thresh = 1; % at least this many subjects to be considered
n_subs_ct_thresh = 1; % significant CTs in region pair must be from at least this many subjects
n_resample = 10*400; % number of times to resample with n_pairs_thresh pairs and n_subs_thresh subjects
trig_mag0_is_nan = false; % only include 0 in average if all resampled is zero

% n_pairs_thresh = 10 * 4; % at least this many electrode pairs to be considered
% n_subs_thresh = 2; % at least this many subjects to be considered
% [*] Original ROI coverage: 78.0952 %
% [*] Thresholded ROI coverage: 61.4286 %

fig_fmt = '-depsc';
trig_eps = false;
trig_mag_no_cp = true;
trig_plot_mag_cov = true;
trig_sym = true;
system('mkdir figures');
system('mkdir figures/T14d2');

% Fast i/o definitions
dir_artL = '/media/klab/internal/data/h5_notch20/art';
dir_resL = '/media/klab/internal/data/results/coh_w10';
dir_corL = '/media/klab/internal/data/coreg';
dir_cacheL = './cache'; %'./cache/old_coh-stats_ct-100_cp-100';
subjects_dirL = '/mnt/cuenap_ssd/coregistration';

% Slow i/o definitions
dir_h5L = '/media/klab/KLAB101/h5_notch20';

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

for iM = 1:length(metrics)
    metric = metrics{iM};
    
    % Load human cache
    load(sprintf('%s/xsub_out_all_%i.mat',dir_cacheL,iM));
    
    
    % Calculate final functional interaction matrix
    Adj = nan(n_rois,n_rois);
    AdjMag = nan(n_rois,n_rois);
    AdjCP = nan(n_rois,n_rois);
    AdjMagVar = nan(n_rois,n_rois);
    AdjMagReS = nan(n_rois,n_rois,n_resample);
    %AdjCPVar = nan(n_rois,n_rois);
    N_bchan = nan(n_rois,n_rois);
    cp_thresh = cp_thresh_override;
    Dmat = Inf(n_rois,n_rois); % set to inf to avoid removing from every instance
    DistsAtl = [];
    
    for i3 = 1:n_resample
        parfor i1 = 1:n_rois
            for i2 = 1:n_rois
                AA = AdjAtl{i1,i2};
                AA_dist = adjct_dist{i1,i2};
                AA_sub = AdjAtl_sid{i1,i2};

                % resample
                rs_idx = 1:length(AA);
                rs_idx = rs_idx(randperm(length(rs_idx)));
                if (length(rs_idx) >= n_pairs_thresh)
                    rs_idx = rs_idx(1:n_pairs_thresh);
                end
                AA = AA(rs_idx);
                AA_dist = AA_dist(rs_idx);
                AA_sub = AA_sub(rs_idx);

                n_pairs = length(AA_sub);
                n_subs = length(unique(AA_sub));
                n_subs_ct = length(unique(AA_sub(AA ~= 0)));
                % ROI pair coverage condition (grey)
                if ( (~ isempty(AA))  && (n_pairs >= n_pairs_thresh) && (n_subs >= n_subs_thresh))
                    N_bchan(i1,i2) = length(AA);
                    frac_cp = sum(AA ~= 0)/length(AA);
                    AdjCP(i1,i2) = frac_cp;
                    % ROI pair significance condition (white)
                    %if ( (frac_cp > cp_thresh) )
                    if ( (frac_cp > cp_thresh) && (n_subs_ct >= n_subs_ct_thresh) )
                        Adj(i1,i2) = 1;
                        AdjMag(i1,i2) = mean(AA(AA ~= 0));
                        AdjMagVar(i1,i2) = var(AA(AA ~= 0));
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
                    end
                end
            end
        end
        
        AdjMagReS(:,:,i3) = AdjMag;
    end
    %return
    AdjMag1 = AdjMag;
    
    if (trig_mag0_is_nan)
        %AdjMagReS_sav = AdjMagReS;
        AdjMagReS_avg = nan(n_rois,n_rois);
        for i2 = 1:n_rois
            for j2 = 1:n_rois
                AdjMagReS_t = squeeze(AdjMagReS(i2,j2,:));
                if (all( squeeze(AdjMagReS(12,13,:)) == 0 ))
                    AdjMagReS_avg(i2,j2) = 0;
                else
                    AdjMagReS_avg(i2,j2) = nanmean(AdjMagReS_t(AdjMagReS_t~=0));
                end
            end
        end
        AdjMag = AdjMagReS_avg;
        %AdjMagReS(AdjMagReS==0) = NaN;
    else
        AdjMag = nanmean(AdjMagReS,3);
    end
    
    
    
    
    if (trig_sym)
        AdjMagSym = nan(size(AdjMag));
        AdjSym = nan(size(Adj));
        for ii = 1:length(AdjMag)
            for jj = (ii):length(AdjMag)
                AdjMagSym(ii,jj) = AdjMag(ii,jj);
                AdjMagSym(jj,ii) = AdjMag(ii,jj);
                AdjSym(ii,jj) = Adj(ii,jj);
                AdjSym(jj,ii) = Adj(ii,jj);
            end
        end
        AdjMag = AdjMagSym;
        Adj = AdjSym;
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
    
    for iii = 1:3

        % Make figure
        h = figure;
        set(h,'Position',round(1*[0 0 0.95*1080 0.8*1080]))

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
%                 mag_min = min(AdjMag(AdjMag ~= 0 ));
%                 mag_max = max(AdjMag(AdjMag ~= 0 ));
                
            else
                ss = 'CP-bin';
            end
        end
        
        colormap(map);
        if (((~isnan(minv)) && (~isnan(maxv))) && (minv ~= maxv))
            cb = colorbar('Ytick',linspace(minv,maxv,5));
            caxis([minv maxv]);
            set(cb,'TickLength',0);
        else
            cb = colorbar;
            set(cb,'TickLength',0);
        end

        print(h,sprintf('figures/T14d2/Adj_%s_%s',metric,ss),fig_fmt);
        if (trig_eps)
            print(h,sprintf('figures/T14d2/Adj_%s_%s',metric,ss),'-depsc');
        end
        close(h);

        
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
            set(h,'Position',round(1*[0 0 0.95*1080 0.8*1080]))
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
            print(h,sprintf('figures/T14d2/Npairs_%s_mag_%i_p_%d_cp_%i_p_%d',metric,...
                round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp),'-depsc');
            
            close(h);
            
            % Number of unique patients vs magnitude
            h = figure;
            set(h,'Position',round(1*[0 0 0.95*1080 0.8*1080]))
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
            print(h,sprintf('figures/T14d2/Nsubs_%s_mag_%i_p_%d_cp_%i_p_%d',metric,...
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
            set(h,'Position',round(1*[0 0 0.95*1080 0.8*1080]))
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
            print(h,sprintf('figures/T14d2/Npairs_%s_magVar_%i_p_%d_cp_%i_p_%d',metric,...
                round(1000*corr_val_mag),corr_pval_mag,round(1000*corr_val_cp),corr_pval_cp),'-depsc');
            
            close(h);
            
            % Number of unique patients vs magnitude
            h = figure;
            set(h,'Position',round(1*[0 0 0.95*1080 0.8*1080]))
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
            print(h,sprintf('figures/T14d2/Nsubs_%s_magVar_%i_p_%d_cp_%i_p_%d',metric,...
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
            
            print(h,sprintf('figures/T14d2/mag_%s_magVar_%i_p_%d',metric,...
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
            print(h,sprintf('figures/T14d2/dist_%s_mag_%i_p_%d',metric,...
                round(1000*corr_val_cp),corr_pval_cp),'-depsc');
            
            close(h);
            %return
        end
        
        
        
        
        
    end

    
end