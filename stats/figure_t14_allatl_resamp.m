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
system('mkdir figures/T14_allatl_resamp');



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

    system(sprintf('mkdir figures/T14_allatl_resamp/atl%i_%s',atl,AtlNames{atl}));
    
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
        
        
        
        % Permutation parameters <-- important
        DropoutN = 0:round(length(Ca_hum.Subjects)/2); %[0 8]; %0:length(Ca_hum.Subjects);
        n_perm = 50000;
        
        
        
        % Drop out subjects
        n_sub = length(Ca_hum.Subjects);
        i_nSub = 1;
        AdjMag_Hi = nan(n_roi_pairs,length(DropoutN));
        AdjMag_Lo = nan(n_roi_pairs,length(DropoutN));
        AdjMag_Med = nan(n_roi_pairs,length(DropoutN));
        AdjMag_Avg = nan(n_roi_pairs,length(DropoutN));
        AdjMag_Std = nan(n_roi_pairs,length(DropoutN));
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
            tic;
            parfor ip = 1:n_perm_ac
                % Pick nSub subjects to drop out by randperm
                sIdx = 1:n_sub;
                sIdx = sIdx(randperm(n_sub));
                sub_ex = sIdx(1:nSub);
                
                % Rebuild AdjMag excluding these subjects
                Adj = nan(n_rois,n_rois);
                AdjMag = nan(n_rois,n_rois);
                AdjMag4cl = nan(n_rois,n_rois);
                AdjCP = nan(n_rois,n_rois);
                AdjMagVar = nan(n_rois,n_rois);
                AdjMagReS = nan(n_rois,n_rois);
                AdjMagL = cell(n_rois,n_rois);
                N_bchan = nan(n_rois,n_rois);
                cp_thresh = cp_thresh_override;
                Dmat = Inf(n_rois,n_rois); % set to inf to avoid removing from every instance
                DistsAtl = [];
                for i1 = 1:n_rois
                    for i2 = (i1):n_rois
                        AA_sub = Ca_hum.AdjAtl_sid{i1,i2};
                        AA = Ca_hum.AdjAtl{i1,i2};
                        AA_dist = Ca_hum.adjct_dist{i1,i2};
                        
                        % Filter sub
                        i_trim = false(size(AA_sub));
                        for i3 = 1:length(sub_ex)
                            i_trim = i_trim | (AA_sub == sub_ex(i3));
                        end
                        AA_sub = AA_sub(~i_trim);
                        AA = AA(~i_trim);
                        AA_dist = AA_dist(~i_trim);
                        
                        
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
                
                AdjMag3(:,:,ip) = AdjMag;
            end
            
            % 95%CI for AdjMag3
            CI_alpha = 0.05;
            AdjMag_hi = nan(n_rois,n_rois);
            AdjMag_lo = nan(n_rois,n_rois);
            AdjMag_med = nan(n_rois,n_rois);
            AdjMag_avg = nan(n_rois,n_rois);
            AdjMag_std = nan(n_rois,n_rois);
            for i1 = 1:n_rois
                for i2 = (i1):n_rois
                    mags = squeeze(AdjMag3(i1,i2,:));
                    s_mags = sort(mags,'descend');
                    sIdx = (~isnan(s_mags)) & (s_mags ~= 0);
                    s_mags = s_mags(sIdx);
                    
                    n_mags = length(s_mags);
                    idx_hi = round((CI_alpha/2) * n_mags);
                    idx_lo = round((1-(CI_alpha/2)) * n_mags);
                    idx_med = round((0.5) * n_mags);
                    
                    % Exception for no resampling (nSub = 0)
                    if (idx_hi == 0)
                        idx_hi = 1;
                    end
                    
                    if (isempty(s_mags))
                        hi = NaN;
                        lo = NaN;
                        med = NaN;
                        avg = NaN;
                        std = NaN;
                    else
                        hi = s_mags(idx_hi);
                        lo = s_mags(idx_lo);
                        med = s_mags(idx_med);
                        avg = nanmean(s_mags);
                        std = nanstd(s_mags);
                    end
                    
                    AdjMag_hi(i1,i2) = hi;
                    AdjMag_lo(i1,i2) = lo;
                    AdjMag_med(i1,i2) = med;
                    AdjMag_avg(i1,i2) = avg;
                    AdjMag_std(i1,i2) = std;
                    
%                     if (i1 == 4) && (i2 == 10)
%                         return
%                     end
                end
            end
            
            % Vectorize CI
            AdjMag_hi_v = nan(n_roi_pairs,1);
            AdjMag_lo_v = nan(n_roi_pairs,1);
            AdjMag_med_v = nan(n_roi_pairs,1);
            AdjMag_avg_v = nan(n_roi_pairs,1);
            AdjMag_std_v = nan(n_roi_pairs,1);
            for c = 1:n_roi_pairs
                i1 = Vec_ij(c,1);
                i2 = Vec_ij(c,2);
                AdjMag_hi_v(c) = AdjMag_hi(i1,i2);
                AdjMag_lo_v(c) = AdjMag_lo(i1,i2);
                AdjMag_med_v(c) = AdjMag_med(i1,i2);
                AdjMag_avg_v(c) = AdjMag_avg(i1,i2);
                AdjMag_std_v(c) = AdjMag_std(i1,i2);
            end
            
            % Store
            AdjMag_Hi(:,i_nSub) = AdjMag_hi_v;
            AdjMag_Lo(:,i_nSub) = AdjMag_lo_v;
            AdjMag_Med(:,i_nSub) = AdjMag_med_v;
            AdjMag_Avg(:,i_nSub) = AdjMag_avg_v;
            AdjMag_Std(:,i_nSub) = AdjMag_std_v;
            
            t_sing = toc;
            fprintf('per atl-metric ETA: %.1f mins\n',(t_sing * (length(DropoutN) - i_nSub))/(60))
            
%             if (i_nSub == 3)
%                 return
%             end
            
            i_nSub = i_nSub + 1;
        end
        
        save(sprintf('./cache/figure_t14_allatl_resamp_im-%i',iM));
        
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
        load(sprintf('./cache/figure_t14_allatl_resamp_im-%i',iM));
        
        % Unknown areas filter
        rois = Ca_hum.rois;
        cond_empty = false(length(rois),1);
        for ce = 1:length(cond_empty)
            elem = rois{ce};
            cond_empty(ce) = (all(isspace(elem)) | isempty(elem));
        end
        cond_contains = (contains(lower(rois),'unknown')) | (contains(lower(rois),'???')) | (contains(lower(rois),'wall')) | (cond_empty);
        known_idx = (~ cond_contains);
        for iz = 1:length(AdjMag_Hi)
            i1 = Vec_ij(iz,1);
            i2 = Vec_ij(iz,2);
            
            %if ( (~ known_idx(i1)) || (~ known_idx(i1)) )
            if ( (~ known_idx(i1)) || (~ known_idx(i2)) )
                AdjMag_Hi(iz) = NaN;
                AdjMag_Lo(iz) = NaN;
                AdjMag_Med(iz) = NaN;
                AdjMag_Avg(iz) = NaN;
                AdjMag_Std(iz) = NaN;
            end
        end
        
        
        
        pct_hi = nanmean(100 * ((AdjMag_Hi(:,2:end) ./ AdjMag_Hi(:,1)) - 1));
        pct_lo = nanmean(100 * (1 - (AdjMag_Lo(:,2:end) ./ AdjMag_Lo(:,1))));
        pct_hi_s = nanstd(100 * ((AdjMag_Hi(:,2:end) ./ AdjMag_Hi(:,1)) - 1));
        pct_lo_s = nanstd(100 * (1 - (AdjMag_Lo(:,2:end) ./ AdjMag_Lo(:,1))));
        diff_avg = ((AdjMag_Hi+AdjMag_Lo)/2);
        diff_avg = abs(diff_avg(:,1) - diff_avg(:,2:end));
        
        fprintf('[*] Number of permutations: %i\n',n_perm);
        for i = 1:length(pct_hi)
%             fprintf('# subject drop: %i, 95%% CI Hi: %.2f%% (+- %.2f%%), Lo:%.2f%% (+- %.2f%%)\n',...
%                 DropoutN(i+1),pct_hi(i),pct_hi_s(i),pct_lo(i),pct_lo_s(i));
            fprintf('# subject drop: %i, absolute diff mu: %.3f, std: %.3f, min: %.3f, max: %.3f\n',...
                DropoutN(i+1), nanmean(diff_avg(:,i)), nanstd(diff_avg(:,i)), nanmin(diff_avg(:,i)), nanmax(diff_avg(:,i)) );
            
            %[p,H,stats] = ranksum(v_A,v_A1);
        end
        
        
        % Remove zeros and nans
        key = AdjMag_Hi(:,1);
        idx_sig = (~isnan(key)) & ((key) ~= 0);
        AdjMag_Hi = AdjMag_Hi(idx_sig,:);
        AdjMag_Lo = AdjMag_Lo(idx_sig,:);
        AdjMag_Med = AdjMag_Med(idx_sig,:);
        AdjMag_Avg = AdjMag_Avg(idx_sig,:);
        AdjMag_Std = AdjMag_Std(idx_sig,:);
        Vec_ij = Vec_ij(idx_sig,:);
        
        % Sort
        key = AdjMag_Hi(:,1);
        [~,idx_sort] = sort(key,'descend');
        AdjMag_Hi = AdjMag_Hi(idx_sort,:);
        AdjMag_Lo = AdjMag_Lo(idx_sort,:);
        AdjMag_Med = AdjMag_Med(idx_sort,:);
        AdjMag_Avg = AdjMag_Avg(idx_sort,:);
        AdjMag_Std = AdjMag_Std(idx_sort,:);
        Vec_ij = Vec_ij(idx_sort,:);
        
        fprintf('[*] Number of edges compared: %i\n',length(AdjMag_Hi(:,1)));
        
        
        
        % trim number of dropout subjects
        %DropoutN = DropoutN(1:(1+6));
        
        
        
        % X-axis
        X = 1:length(AdjMag_Hi);
        
        close all;
        h = figure('Position',[0 0 800 200],'visible','off');
        hold all;
        
        % Colormap
        %cmap = inferno(length(DropoutN) - 1);
        cmap = gray(length(DropoutN) - 1 + 1);
        cmap = cmap(2:end,:);
        cmap = [(0.5*[0 0 1]); cmap];
        
        % Plot
        DropoutNp = DropoutN(1:end);
        for iR = 1:length(DropoutNp)
            i = length(DropoutNp) - iR + 1;
            if (i == 1)
                mk_size = 4;
            else
                mk_size = 2;
            end
%             Y = AdjMag_Hi(:,i);
%             plot(X(Y~=0),Y(Y~=0),'.','Color',cmap(i,:),'MarkerSize',mk_size);
%             Y = AdjMag_Lo(:,i);
%             plot(X(Y~=0),Y(Y~=0),'.','Color',cmap(i,:),'MarkerSize',mk_size);
            %Y = AdjMag_Avg(:,i);
            %plot(X(Y~=0),Y(Y~=0),'.','Color',cmap(i,:),'MarkerSize',mk_size);
            Y = AdjMag_Avg(:,i) + AdjMag_Std(:,i);
            plot(X(Y~=0),Y(Y~=0),'.','Color',cmap(i,:),'MarkerSize',mk_size);
            Y = AdjMag_Avg(:,i) - AdjMag_Std(:,i);
            plot(X(Y~=0),Y(Y~=0),'.','Color',cmap(i,:),'MarkerSize',mk_size);
        end
        
        % X labels
        xlab_txt = cell(length(X),1);
        roiname = Ca_hum.C.AtlROIs{2}.RH.struct_names;
        %return
%         def_rois_short;
%         roiname = rois_short;
        for iL = 1:length(X)
            %xlab_txt{iL} = sprintf('%2i,%2i',Vec_ij(iL,1),Vec_ij(iL,2));
            xlab_txt{iL} = sprintf('%s,%s',...
                convertRoiDK(roiname{Vec_ij(iL,1)}),...
                convertRoiDK(roiname{Vec_ij(iL,2)}));
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
        
        print(h,sprintf('figures/figure_t14_allatl_resamp_metric%i_atl%i',iM,atl),'-depsc');
        close(h);
        
        % correlation
        corr_avg_r = nan(1,nSub);
        corr_avg_p = nan(1,nSub);
        corr_avg_n = length(AdjMag_Avg(:,1));
        for i = 1:nSub
            [r,p] = corr(AdjMag_Avg(:,1),AdjMag_Avg(:,1+i));
            corr_avg_r(i) = r;
            corr_avg_p(i) = p;
        end
        
        % Show difference in means
        diff_avg = (AdjMag_Avg(:,1) - AdjMag_Avg(:,2:end));
        h = figure('visible','off','Position',[0 0 320 320]);
        plot(0:(length(AdjMag_Avg(1,:))-1),[0 ,max(abs(diff_avg))],'black.');
        xlabel('Number of Dropout Subjects')
        ylabel(sprintf('Absolute Difference in %s Coherence',metric(3:end)));
        %axis tight;
        set(gca,'TickDir','Out');
        box off;
        print(h,sprintf('figures/figure_t14_allatl_resamp_metric%i_atl%i_diff_avg',iM,atl),'-depsc');
        close(h);
        
        % Show difference in std
        diff_std = (AdjMag_Std(:,1) - AdjMag_Std(:,2:end));
        h = figure('visible','off','Position',[0 0 320 320]);
        plot(0:(length(AdjMag_Std(1,:))-1),[0 ,max(abs(diff_std))],'black.');
        xlabel('Number of Dropout Subjects')
        ylabel(sprintf('Absolute Difference in %s Coherence SD',metric(3:end)));
        %axis tight;
        set(gca,'TickDir','Out');
        box off;
        print(h,sprintf('figures/figure_t14_allatl_resamp_metric%i_atl%i_diff_std',iM,atl),'-depsc');
        close(h);
        
        [maxdiffavg,maxdiffavg_mIdx] = max(abs(diff_avg));
        [maxdiffstd,maxdiffstd_mIdx] = max(abs(diff_std));
        fprintf('[*] metric: %s\n',metric(3:end));
        fprintf('[*] Mean\n');
        for i = 1:length(maxdiffavg)
            fprintf('sub drop %i: max diff avg = %.6f (%.6f, %.2f%%), r=%.3d, p=%.3d, n=%i, %s, max diff std = %.6f (%.6f, %.2f%%), %s\n',...
                i,maxdiffavg(i),abs(AdjMag_Avg(maxdiffavg_mIdx(i),1)),100*maxdiffavg(i)/(abs(AdjMag_Avg(maxdiffavg_mIdx(i),1))),...
                corr_avg_r(i),corr_avg_p(i),corr_avg_n,xlab_txt{maxdiffavg_mIdx(i)},...
                maxdiffstd(i),...
                abs(AdjMag_Std(maxdiffstd_mIdx(i),1)),100*maxdiffstd(i)/(abs(AdjMag_Std(maxdiffstd_mIdx(i),1))),...
                xlab_txt{maxdiffstd_mIdx(i)});
        end
        
        diff_med = (AdjMag_Med(:,1) - AdjMag_Med(:,2:end));
        [maxdiffmed,maxdiffmed_mIdx] = max(abs(diff_med));
        fprintf('[*] Median\n');
        for i = 1:length(maxdiffmed)
            fprintf('sub drop %i: max diff med = %.6f, %s\n',...
                i,maxdiffmed(i),xlab_txt{maxdiffmed_mIdx(i)});
        end
        
        % log-normal test
        clear std;
        for i = 1:length(maxdiffmed)
            x = AdjMag_Avg(:,i+1);
            [h_lin,pval_lin,~] = kstest((x-mean(x))/std(x));
            x = log10(AdjMag_Avg(:,i+1));
            [h_log,pval_log,~] = kstest((x-mean(x))/std(x));
            fprintf('sub drop %i, KS normal p=%.6d, KS lognormal p=%.6d\n',...
                i,pval_lin,pval_log);
        end
    end

end




% 20 oct. 2020
% [*] metric: Broadband
% sub drop 1: max diff avg = 0.003356, IPA,POR, max diff std = 0.031170, CUN,PHC
% sub drop 2: max diff avg = 0.006915, IPA,POR, max diff std = 0.042752, CUN,PHC
% sub drop 3: max diff avg = 0.009913, IPA,POR, max diff std = 0.051655, CUN,PHC
% sub drop 4: max diff avg = 0.012209, IPA,POR, max diff std = 0.058656, CUN,PHC
% sub drop 5: max diff avg = 0.014904, IPA,POR, max diff std = 0.064652, CUN,PHC
% sub drop 6: max diff avg = 0.016997, IPA,POR, max diff std = 0.069949, CUN,PHC
% sub drop 7: max diff avg = 0.019201, IPA,POR, max diff std = 0.074339, CUN,PHC
% sub drop 8: max diff avg = 0.021115, IPA,POR, max diff std = 0.078033, CUN,PHC
% sub drop 9: max diff avg = 0.023748, IPA,POR, max diff std = 0.081650, CUN,PHC
% sub drop 10: max diff avg = 0.025291, CUN,LIN, max diff std = 0.085375, CUN,PHC
% sub drop 11: max diff avg = 0.027502, CUN,LIN, max diff std = 0.088479, CUN,PHC
% sub drop 12: max diff avg = 0.029597, CUN,LIN, max diff std = 0.090091, CUN,PHC
% sub drop 13: max diff avg = 0.031170, CUN,PHC, max diff std = 0.092596, POR,FRP
% sub drop 14: max diff avg = 0.033547, CUN,PHC, max diff std = 0.094574, CUN,PHC
% sub drop 15: max diff avg = 0.036147, CUN,PHC, max diff std = 0.097049, POR,FRP
% sub drop 16: max diff avg = 0.039936, CUN,PHC, max diff std = 0.102227, PHC,POP
% sub drop 17: max diff avg = 0.042716, CUN,PHC, max diff std = 0.107529, PHC,POP
% sub drop 18: max diff avg = 0.043654, CUN,PHC, max diff std = 0.113541, PHC,POP
% sub drop 19: max diff avg = 0.047074, CUN,PHC, max diff std = 0.120822, PHC,POP
% sub drop 20: max diff avg = 0.050303, CUN,PHC, max diff std = 0.126321, PHC,POP
% sub drop 21: max diff avg = 0.053644, CUN,PHC, max diff std = 0.131736, PHC,POP
% sub drop 22: max diff avg = 0.056039, CUN,PHC, max diff std = 0.138033, PHC,POP
% sub drop 23: max diff avg = 0.059405, CUN,PHC, max diff std = 0.142310, PHC,POP
% sub drop 24: max diff avg = 0.061773, CUN,PHC, max diff std = 0.148020, PHC,POP
% [*] Median
% sub drop 1: max diff med = 0.000000, CMF,PHC
% sub drop 2: max diff med = 0.000003, ITP,STP
% sub drop 3: max diff med = 0.000286, ITP,STP
% sub drop 4: max diff med = 0.000635, ITP,STP
% sub drop 5: max diff med = 0.001775, STP,TPP
% sub drop 6: max diff med = 0.003856, PRC,RMF
% sub drop 7: max diff med = 0.008713, MTP,TPP
% sub drop 8: max diff med = 0.008713, MTP,TPP
% sub drop 9: max diff med = 0.008871, MTP,TPP
% sub drop 10: max diff med = 0.011107, STP,TPP
% sub drop 11: max diff med = 0.011107, STP,TPP
% sub drop 12: max diff med = 0.011253, STP,TPP
% sub drop 13: max diff med = 0.012215, STP,TPP
% sub drop 14: max diff med = 0.014991, POR,TPP
% sub drop 15: max diff med = 0.030133, CMF,FUS
% sub drop 16: max diff med = 0.030293, PHC,POP
% sub drop 17: max diff med = 0.030293, PHC,POP
% sub drop 18: max diff med = 0.030293, PHC,POP
% sub drop 19: max diff med = 0.030293, PHC,POP
% sub drop 20: max diff med = 0.030293, PHC,POP
% sub drop 21: max diff med = 0.030293, PHC,POP
% sub drop 22: max diff med = 0.030293, PHC,POP
% sub drop 23: max diff med = 0.030293, PHC,POP
% sub drop 24: max diff med = 0.061305, FUS,SPA
% sub drop 1, KS normal p=3.310498e-02, KS lognormal p=5.887445e-01
% sub drop 2, KS normal p=3.093798e-02, KS lognormal p=5.527295e-01
% sub drop 3, KS normal p=3.015473e-02, KS lognormal p=4.867535e-01
% sub drop 4, KS normal p=2.605114e-02, KS lognormal p=4.650428e-01
% sub drop 5, KS normal p=2.322974e-02, KS lognormal p=4.478843e-01
% sub drop 6, KS normal p=2.867303e-02, KS lognormal p=4.114829e-01
% sub drop 7, KS normal p=2.764575e-02, KS lognormal p=3.884929e-01
% sub drop 8, KS normal p=2.828790e-02, KS lognormal p=4.803641e-01
% sub drop 9, KS normal p=2.820058e-02, KS lognormal p=4.856742e-01
% sub drop 10, KS normal p=3.449964e-02, KS lognormal p=4.786667e-01
% sub drop 11, KS normal p=3.835140e-02, KS lognormal p=4.275613e-01
% sub drop 12, KS normal p=4.080417e-02, KS lognormal p=3.720453e-01
% sub drop 13, KS normal p=4.380047e-02, KS lognormal p=3.467069e-01
% sub drop 14, KS normal p=4.456362e-02, KS lognormal p=3.845371e-01
% sub drop 15, KS normal p=5.004605e-02, KS lognormal p=5.028521e-01
% sub drop 16, KS normal p=5.054548e-02, KS lognormal p=5.106763e-01
% sub drop 17, KS normal p=5.277681e-02, KS lognormal p=4.941947e-01
% sub drop 18, KS normal p=5.086797e-02, KS lognormal p=4.882962e-01
% sub drop 19, KS normal p=5.399647e-02, KS lognormal p=4.996563e-01
% sub drop 20, KS normal p=5.151908e-02, KS lognormal p=5.032079e-01
% sub drop 21, KS normal p=5.576742e-02, KS lognormal p=5.060208e-01
% sub drop 22, KS normal p=6.555997e-02, KS lognormal p=5.318636e-01
% sub drop 23, KS normal p=5.981421e-02, KS lognormal p=5.583690e-01
% sub drop 24, KS normal p=6.397453e-02, KS lognormal p=5.555515e-01

% Gamma

% [*] Number of permutations: 50000
% # subject drop: 1, absolute diff mu: 0.002, std: 0.004, min: 0.000, max: 0.031
% # subject drop: 2, absolute diff mu: 0.009, std: 0.011, min: 0.000, max: 0.065
% # subject drop: 3, absolute diff mu: 0.009, std: 0.011, min: 0.000, max: 0.065
% # subject drop: 4, absolute diff mu: 0.009, std: 0.011, min: 0.000, max: 0.065
% # subject drop: 5, absolute diff mu: 0.009, std: 0.011, min: 0.000, max: 0.065
% # subject drop: 6, absolute diff mu: 0.009, std: 0.011, min: 0.000, max: 0.065
% # subject drop: 7, absolute diff mu: 0.009, std: 0.011, min: 0.000, max: 0.065
% # subject drop: 8, absolute diff mu: 0.009, std: 0.012, min: 0.000, max: 0.065
% # subject drop: 9, absolute diff mu: 0.010, std: 0.012, min: 0.000, max: 0.063
% # subject drop: 10, absolute diff mu: 0.011, std: 0.012, min: 0.000, max: 0.063
% # subject drop: 11, absolute diff mu: 0.011, std: 0.012, min: 0.000, max: 0.056
% # subject drop: 12, absolute diff mu: 0.010, std: 0.012, min: 0.000, max: 0.065
% # subject drop: 13, absolute diff mu: 0.011, std: 0.012, min: 0.000, max: 0.065
% # subject drop: 14, absolute diff mu: 0.011, std: 0.013, min: 0.000, max: 0.065
% # subject drop: 15, absolute diff mu: 0.011, std: 0.013, min: 0.000, max: 0.065
% # subject drop: 16, absolute diff mu: 0.011, std: 0.012, min: 0.000, max: 0.065
% # subject drop: 17, absolute diff mu: 0.012, std: 0.013, min: 0.000, max: 0.065
% # subject drop: 18, absolute diff mu: 0.013, std: 0.014, min: 0.000, max: 0.065
% # subject drop: 19, absolute diff mu: 0.014, std: 0.018, min: 0.000, max: 0.177
% # subject drop: 20, absolute diff mu: 0.015, std: 0.020, min: 0.000, max: 0.177
% # subject drop: 21, absolute diff mu: 0.015, std: 0.020, min: 0.000, max: 0.177
% # subject drop: 22, absolute diff mu: 0.016, std: 0.021, min: 0.000, max: 0.177
% # subject drop: 23, absolute diff mu: 0.016, std: 0.021, min: 0.000, max: 0.177
% # subject drop: 24, absolute diff mu: 0.017, std: 0.022, min: 0.000, max: 0.177
% [*] Number of edges compared: 184
% [*] metric: Gamma
% sub drop 1: max diff avg = 0.002021, POR,FRP, max diff std = 0.032496, PHC,POP
% sub drop 2: max diff avg = 0.003677, POR,SMA, max diff std = 0.046279, PHC,POP
% sub drop 3: max diff avg = 0.005974, POR,FRP, max diff std = 0.055268, PHC,POP
% sub drop 4: max diff avg = 0.007684, POR,FRP, max diff std = 0.063395, PHC,POP
% sub drop 5: max diff avg = 0.008983, POR,FRP, max diff std = 0.070283, PHC,POP
% sub drop 6: max diff avg = 0.010404, POR,FRP, max diff std = 0.076279, PHC,POP
% sub drop 7: max diff avg = 0.012482, POR,FRP, max diff std = 0.082280, PHC,POP
% sub drop 8: max diff avg = 0.013613, POR,FRP, max diff std = 0.087477, PHC,POP
% sub drop 9: max diff avg = 0.015261, POR,FRP, max diff std = 0.092077, PHC,POP
% sub drop 10: max diff avg = 0.016750, POR,FRP, max diff std = 0.096986, PHC,POP
% sub drop 11: max diff avg = 0.017929, POR,SMA, max diff std = 0.100775, PHC,POP
% sub drop 12: max diff avg = 0.019222, POR,FRP, max diff std = 0.104629, PHC,POP
% sub drop 13: max diff avg = 0.020880, POR,SMA, max diff std = 0.108342, PHC,POP
% sub drop 14: max diff avg = 0.022072, POR,SMA, max diff std = 0.112271, PHC,POP
% sub drop 15: max diff avg = 0.022945, POR,SMA, max diff std = 0.115492, PHC,POP
% sub drop 16: max diff avg = 0.024630, POR,SMA, max diff std = 0.119421, PHC,POP
% sub drop 17: max diff avg = 0.026156, POR,SMA, max diff std = 0.122948, PHC,POP
% sub drop 18: max diff avg = 0.027421, POR,SMA, max diff std = 0.126011, PHC,POP
% sub drop 19: max diff avg = 0.029397, POR,SMA, max diff std = 0.130054, PHC,POP
% sub drop 20: max diff avg = 0.030685, POR,SMA, max diff std = 0.134086, PHC,POP
% sub drop 21: max diff avg = 0.031970, POR,SMA, max diff std = 0.137509, PHC,POP
% sub drop 22: max diff avg = 0.033569, POR,SMA, max diff std = 0.141865, PHC,POP
% sub drop 23: max diff avg = 0.034257, POR,SMA, max diff std = 0.144823, PHC,POP
% sub drop 24: max diff avg = 0.036002, POR,SMA, max diff std = 0.149323, PHC,POP
% [*] Median
% sub drop 1: max diff med = 0.000000, CMF,PHC
% sub drop 2: max diff med = 0.000251, MTP,STP
% sub drop 3: max diff med = 0.000410, MTP,STP
% sub drop 4: max diff med = 0.002274, PSC,SMA
% sub drop 5: max diff med = 0.002535, PSC,SMA
% sub drop 6: max diff med = 0.003001, POP,PRC
% sub drop 7: max diff med = 0.005793, FUS,PSC
% sub drop 8: max diff med = 0.006446, FUS,PSC
% sub drop 9: max diff med = 0.006446, FUS,PSC
% sub drop 10: max diff med = 0.008316, MTP,POP
% sub drop 11: max diff med = 0.009403, PSC,PRC
% sub drop 12: max diff med = 0.009769, FUS,PSC
% sub drop 13: max diff med = 0.013450, FUS,PSC
% sub drop 14: max diff med = 0.014270, FUS,PSC
% sub drop 15: max diff med = 0.014545, FUS,PSC
% sub drop 16: max diff med = 0.015274, FUS,PSC
% sub drop 17: max diff med = 0.015274, FUS,PSC
% sub drop 18: max diff med = 0.016759, FUS,PSC
% sub drop 19: max diff med = 0.017903, FUS,PSC
% sub drop 20: max diff med = 0.021903, FUS,PSC
% sub drop 21: max diff med = 0.024170, ENT,IPA
% sub drop 22: max diff med = 0.024170, ENT,IPA
% sub drop 23: max diff med = 0.062361, PHC,POP
% sub drop 24: max diff med = 0.062361, PHC,POP
% sub drop 1, KS normal p=5.225995e-02, KS lognormal p=1.195340e-01
% sub drop 2, KS normal p=5.264930e-02, KS lognormal p=1.197721e-01
% sub drop 3, KS normal p=5.182415e-02, KS lognormal p=1.162815e-01
% sub drop 4, KS normal p=5.159939e-02, KS lognormal p=1.352306e-01
% sub drop 5, KS normal p=4.969681e-02, KS lognormal p=1.126918e-01
% sub drop 6, KS normal p=4.833563e-02, KS lognormal p=1.285611e-01
% sub drop 7, KS normal p=4.629876e-02, KS lognormal p=1.135094e-01
% sub drop 8, KS normal p=4.375997e-02, KS lognormal p=1.118237e-01
% sub drop 9, KS normal p=4.054326e-02, KS lognormal p=1.092150e-01
% sub drop 10, KS normal p=4.034073e-02, KS lognormal p=1.107737e-01
% sub drop 11, KS normal p=4.880195e-02, KS lognormal p=1.387363e-01
% sub drop 12, KS normal p=5.612308e-02, KS lognormal p=1.587270e-01
% sub drop 13, KS normal p=6.455514e-02, KS lognormal p=1.929610e-01
% sub drop 14, KS normal p=6.867791e-02, KS lognormal p=2.245001e-01
% sub drop 15, KS normal p=7.225465e-02, KS lognormal p=2.125440e-01
% sub drop 16, KS normal p=6.880120e-02, KS lognormal p=2.294363e-01
% sub drop 17, KS normal p=7.174513e-02, KS lognormal p=2.310078e-01
% sub drop 18, KS normal p=7.609327e-02, KS lognormal p=2.288450e-01
% sub drop 19, KS normal p=7.641032e-02, KS lognormal p=2.329055e-01
% sub drop 20, KS normal p=7.711432e-02, KS lognormal p=2.366074e-01
% sub drop 21, KS normal p=8.224344e-02, KS lognormal p=2.358662e-01
% sub drop 22, KS normal p=8.694217e-02, KS lognormal p=2.315267e-01
% sub drop 23, KS normal p=9.190033e-02, KS lognormal p=2.282869e-01
% sub drop 24, KS normal p=9.092420e-02, KS lognormal p=2.269723e-01







% 29 janv. 2020

% Gamma n=1
% [n_perm: 1] atl: 2	metric: 5	#drop: 0	of 1
% per atl-metric ETA: 0.0 mins
% [n_perm: 100000] atl: 2	metric: 5	#drop: 1	of 1
% per atl-metric ETA: 0.0 mins
% [*] Number of permutations: 100000
% # subject drop: 1, 95% CI Hi: 1.74% (+- 2.38%), Lo:1.37% (+- 2.49%)
% [*] Number of edges compared: 183

% Gamma n=8
% [n_perm: 1] atl: 2	metric: 5	#drop: 0	of 8
% Starting parallel pool (parpool) using the 'local' profile ...
% Connected to the parallel pool (number of workers: 6).
% per atl-metric ETA: 0.5 mins
% [n_perm: 100000] atl: 2	metric: 5	#drop: 8	of 8
% per atl-metric ETA: 0.0 mins
% [*] Number of permutations: 100000
% # subject drop: 8, 95% CI Hi: 12.00% (+- 12.25%), Lo:10.58% (+- 10.46%)
% [*] Number of edges compared: 183


% Broadband all
% [*] Number of permutations: 100000
% # subject drop: 1, 95% CI Hi: 2.18% (+- 3.27%), Lo:1.60% (+- 2.60%)
% # subject drop: 2, 95% CI Hi: 11.05% (+- 13.65%), Lo:10.18% (+- 10.92%)
% # subject drop: 3, 95% CI Hi: 11.10% (+- 13.64%), Lo:10.24% (+- 10.89%)
% # subject drop: 4, 95% CI Hi: 11.36% (+- 13.68%), Lo:10.37% (+- 10.90%)
% # subject drop: 5, 95% CI Hi: 11.80% (+- 13.71%), Lo:10.79% (+- 10.96%)
% # subject drop: 6, 95% CI Hi: 12.40% (+- 13.87%), Lo:11.30% (+- 11.11%)
% # subject drop: 7, 95% CI Hi: 13.32% (+- 14.31%), Lo:12.04% (+- 11.40%)
% # subject drop: 8, 95% CI Hi: 14.28% (+- 15.21%), Lo:12.37% (+- 11.47%)
% # subject drop: 9, 95% CI Hi: 16.58% (+- 17.38%), Lo:14.77% (+- 13.49%)
% # subject drop: 10, 95% CI Hi: 17.48% (+- 17.97%), Lo:15.38% (+- 13.68%)
% # subject drop: 11, 95% CI Hi: 17.90% (+- 18.24%), Lo:15.70% (+- 13.76%)
% # subject drop: 12, 95% CI Hi: 18.40% (+- 18.63%), Lo:15.91% (+- 13.84%)
% # subject drop: 13, 95% CI Hi: 19.05% (+- 18.90%), Lo:16.41% (+- 14.11%)
% # subject drop: 14, 95% CI Hi: 19.78% (+- 19.42%), Lo:16.63% (+- 14.23%)
% # subject drop: 15, 95% CI Hi: 20.82% (+- 21.36%), Lo:16.98% (+- 14.40%)
% # subject drop: 16, 95% CI Hi: 21.95% (+- 21.93%), Lo:17.50% (+- 14.57%)
% [*] Number of edges compared: 193



