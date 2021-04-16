close all;
clear;

addpath(genpath('SWP'));


[~,host] = system('hostname');
if contains(host,'hopperu')
    % Path definitions
    dir_art = '/mnt/cuenap/data/scripts/v2/art_szs';%'/media/klab/KLAB101/h5_notch20/art_nosz';%'/mnt/cuenap/data/h5_notch20/art_nosz';
    dir_res = '/mnt/cuenap/data/results/coh_w20';
    dir_cor = '/mnt/cuenap_ssd/coregistration';
    dir_h5 = '/mnt/cuenap/data/h5_notch20';
    %dir_out = '/';
elseif contains(host,'o2.rc.hms.harvard.edu')
    % Path definitions
    dir_art = '/n/groups/kreiman/jerry/data/h5/art';
    dir_res = '/n/scratch3/users/j/jw324/opencl2/coh_w10';
    dir_cor = '/n/scratch2/jw324/data/coreg';
    dir_h5 = '/n/groups/kreiman/jerry/data/h5';
elseif contains(host,'kraken')
    dir_art = '/media/klab/untitled/h5_notch20/art_nosz'; %'/media/klab/internal/data/h5_notch20/art';
    dir_res = '/media/klab/untitled/results/coh_w10';
    dir_cor = '/media/klab/internal/data/coreg'; %'/mnt/cuenap_ssd/coregistration';
    dir_h5 = '/media/klab/untitled/h5_notch20';
elseif contains(host,'Leibniz')
    dir_art = '/Volumes/KLAB101/h5_notch20/art_nosz'; %'/media/klab/internal/data/h5_notch20/art';
    dir_res = '/Volumes/KLAB101/results/coh_w10';
    dir_cor = '/Volumes/KLAB101/coreg'; %'/mnt/cuenap_ssd/coregistration';
    dir_h5 = '/Volumes/KLAB101/h5_notch20';
elseif contains(host,'ubuntu_1604')
    dir_art = '/nas_share/RawData/data/h5_notch20/art_nosz';
    dir_res = '/nas_share/RawData/scripts/opencl/results_res5hz';
    dir_cor = '/mnt/cuenap_ssd/coregistration';
    %dir_h5 = '/nas_share/RawData/data/h5_notch20';
    dir_h5 = '/mnt/cuenap/data/h5_notch20';
end

SubjectsSleep = {'m00019','m00023','m00024','m00026','m00030','m00032','m00035',...
    'm00039','m00043','m00044','m00049','m00079','m00083','m00084'};

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

% classification params
trig_svm = false;
trig_loCoh = false;
%n_resample = 10;
n_resample = 1;

dir_out_fig = './figures/sleep_sw';
mkdir(dir_out_fig);
n_perm = 10000;
atl = 2;
trig_plot_filter = false;
for metrici = [1]
    metric = metrics{metrici};
    %features = {};
    
    Acc_mean = [];
    Acc_null = [];
    Acc_std = [];
    Acc_n = [];
    Acc_nwake = [];
    Acc_nsleep = [];
    Sigma_diff_sleep = [];
    Sigma_ranksum_p = [];
    Omega_diff_sleep = [];
    Omega_ranksum_p = [];
    Coh_diff_sleep = [];
    Coh_ranksum_sig = [];
    Sub_i = [];
    i_feat = 1;
    %for i = [17 20 22 39] %1:length(Subjects) % 15
    for i = 1:length(Subjects) % 15
        isleep = strcmp(Subjects{i},SubjectsSleep);
        sub_has_sleep = ~isempty(find(isleep));
        
%         % manually exclude subjects with no usable annotations
%         sub_has_sleep = (sub_has_sleep && (i ~= 17)); 
%         sub_has_sleep = (sub_has_sleep && (i ~= 20));
%         sub_has_sleep = (sub_has_sleep && (i ~= 22));
%         sub_has_sleep = (sub_has_sleep && (i ~= 39));
        if (sub_has_sleep)
            sid = Subjects{i};
            Ca = load(sprintf('./cache/%s_sleepsig_metric-%i_atl-%i.mat',sid,metrici,atl));
            
            pvals = Ca.pvals;
            pvals2 = (-1)*log10(pvals);
            iidx = ~isinf(pvals2);
            pvals2 = pvals2(iidx);
            pvals = pvals(iidx);
            [pvals2s,sIdx] = sort(pvals2);
            
            % Load
            fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
            fn_dist = sprintf('%s/%s_dists-%s-%i.mat',dir_res,sid,metric,n_perm);
            fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
            art_idx = h5read(fn_art,'/art_idx');
            art_idxB = (art_idx == 1);
            R = h5read(fn_graph,'/R',[1 1],size(art_idx));
            
            % plot sig coherence
%             for j = 13
%                 R1 = R(j,:);
%                 R1(art_idxB(j,:)) = NaN; % remove artifacts
%                 thresh = Ca.Ca.coh_thresh(j);
%                 t = linspace(0,(length(R1)*Ca.w)/(3600),length(R1));
%                 h = figure('Position',[0 0 1200 1200],'visible','on'); hold all;
%                 plot(t,R1,'black.'); hold on;
%                 plot([t(1) t(end)],thresh*[1 1],'--','color',[0 0 1]*0.8);
%                 return;
%             end

            % plot swi
            
            
            % compute small world
            % cache_plot_sleep_sw_%s - only consider significant
            % cache_plot_sleep_sw_all_%s - consider all
            fn_cache_sw = sprintf('./cache/plot_sleep_sw_all_%s_metric-%i',sid,metrici);
            if (exist([fn_cache_sw,'.mat'],'file'))
                fprintf('[*] found cache, loading: %s\n',fn_cache_sw);
                Ca2 = load(fn_cache_sw);
                Sigma = Ca2.Sigma;
                Omega = Ca2.Omega;
            else
                Sigma = nan(1,Ca.n_graph);
                Omega = nan(1,Ca.n_graph);
                parfor j = 1:Ca.n_graph
                    % count
                    if (mod(j-1,round(Ca.n_graph/100)) == 0)
                        fprintf('[%s] window %i of %i (%.1f %%)\n',sid,j,Ca.n_graph,100*(j/Ca.n_graph));
                    end

                    Ahs = nan(Ca.ecog.n_bchan,Ca.ecog.n_bchan);
                    for c = 1:Ca.n_pairs
                        b1 = Ca.chan1(c) + 1;
                        b2 = Ca.chan2(c) + 1;
                        cond_dist = (Ca.Ca.Dmats(c) > Ca.Ca.dist_thresh);
                        cond_sig = (R(c,j) > Ca.Ca.coh_thresh(c));
                        if (cond_dist)
                            if (true) %(cond_sig)
                                Ahs(b1,b2) = R(c,j);
                                Ahs(b2,b1) = R(c,j);
                            else
                                Ahs(b1,b2) = 0;
                                Ahs(b2,b1) = 0;
                            end
                        else
                            Ahs(b1,b2) = 0;
                            Ahs(b2,b1) = 0;
                        end
                    end

                    % ----------------------------------------------------------------
                    % Remove NaNs
                    %return
                    Ahs_nodiag = Ahs;
                    for inod = 1:length(Ahs)
                        Ahs_nodiag(inod,inod) = 0;
                    end
                    A = Ahs_nodiag;
                    %A = Ahs;
                    cond_pass = true;
                    comp_idx = 1:length(A);
                    while (cond_pass)
                        [~,A_b] = max(sum(isnan(A)));
                        A(A_b,:) = [];
                        A(:,A_b) = [];
                        comp_idx(A_b) = [];
                        cond_pass = (sum(isnan(A(:))) ~= 0);
                    end
                    Ahs = A;
                    %rois2 = rois2(comp_idx);
                    %Ahs_dist = Ahs_dist(comp_idx,comp_idx);
                    % ----------------------------------------------------------------

                    if (~isempty(Ahs))
                        % Compute small world metric
                        %swp = nan;

                        if (~all(Ahs(:)==0))
                            %fprintf('[!]\n')

                            %try
                                [swp,dc,dl,net_clus,rand_clus,net_path,rand_path,reg_clus,reg_path,gamma,lambda,sigma] = small_world_propensity_raw(Ahs);
                                Sigma(j) = sigma;
                                Omega(j) = (rand_path/net_path) - (net_clus/reg_clus);
                            %catch
                                % cannot compute for whatever reason
                            %end 
                        end
                    end
                end
                save(fn_cache_sw,'Sigma','Omega');
            end
            
            %return
            
            % plot coherence
            n_plot = 2;
            h = figure('Position',[0 0 1600 800],'visible','off'); hold all;
            pIdx = round(linspace(1,length(sIdx),n_plot));
            for j = 1:n_plot
                subplot(n_plot,1,j);
                rIdx = sIdx(pIdx(j));
                
                % variable to plot
                %R1 = R(rIdx,:);
                if (j == 1)
                    R1 = Sigma;
                    yltxt = 'Sigma';
                else
                    R1 = Omega;
                    yltxt = 'Omega';
                end
                
                R1(art_idxB(rIdx,:)) = NaN;
                % remove ambiguous label
                R1(Ca.label_sleep==2) = NaN;
                L1 = Ca.label_sleep;
                % remove ambiguous label
                L1(L1==2) = NaN;
                
                % harmonize length
                if (length(R1) > length(L1))
                    R1 = R1(1:length(L1));
                elseif (length(R1) < length(L1))
                    L1 = L1(1:length(R1));
                end
                
                % Cut ambiguous label and shrink
                cond_same_length = ((length(Ca.label_sleep) == length(L1)) && (length(R1) == length(L1)));
                if (~cond_same_length)
                    fprintf(2,'[!] Warn: variables R1, L1, Ca.label_sleep not the same length.\n');
                end
                state_begin = 1;
                state_prev = Ca.label_sleep(1);
                amb_start = [];
                amb_stop = [];
                % check for starting labels
                if (state_prev == 2)
                    amb_start = [amb_start; 1];
                end
                for l = 1:length(Ca.label_sleep)
                    state_change = (Ca.label_sleep(l) ~= state_prev);
                    
                    % start of ambiguous
                    if (state_change && (Ca.label_sleep(l) == 2))
                        amb_start = [amb_start; l];
                    % end of ambiguous 
                    elseif (state_change && (state_prev == 2))
                        amb_stop = [amb_stop; l];
                    end
                    
                    state_prev = Ca.label_sleep(l);
                end
                % check for and close any hanging ambiguous labels
                if (state_prev == 2)
                    amb_stop = [amb_stop; l];
                end
                
%                 % Pop ambiguous
%                 if (~isempty(amb_start))
%                     popped = nan(size(amb_start));
%                     popped_idx = false(size(amb_start));
%                     offset = 0;
%                     for p = 1:length(amb_start)
%                         astart = amb_start(p);
%                         astop = amb_stop(p);
%                         popped_idx(astart:astop) = true;
%                         popped(p) = amb_start(p) - offset;
%                         offset = offset + (astop - astart + 1);
%                     end
%                     R1 = R1(~popped_idx);
%                     L1 = L1(~popped_idx);
%                     % last cut compensation
%                     if (popped(end) > length(R1))
%                         popped(end) = length(R1);
%                     end
%                 end
                if (~isempty(amb_start))
                    popped = nan(size(amb_start));
                    popped_idx = false(size(L1));
                    offset = 0;
                    for p = 1:length(amb_start)
                        astart = amb_start(p);
                        astop = amb_stop(p);
                        popped_idx(astart:astop) = true;
                        popped(p) = amb_start(p) - offset;
                        offset = offset + (astop - astart + 1);
                    end
                    R1 = R1(~popped_idx);
                    L1 = L1(~popped_idx);
                    % last cut compensation
                    if (popped(end) > length(R1))
                        popped(end) = length(R1);
                    end
                end
                
                
                % show day/night
                t = linspace(0,(length(R1)*Ca.w)/(3600),length(L1));
                
                Rplt = R1(~isinf(R1));
                state_prev = L1(1);
                state_begin = t(1);
                col_day = [253 223 161]/255;
                col_night = [161 207 253]/255;
                for l = 1:length(L1)
                    if (~isnan(L1(l)))
                        state_change = (L1(l) ~= state_prev);

                        % shade rectangle on state change
                        if (state_change)
                            % --- plot rectangle ---
                            x = state_begin;
                            y = min(Rplt);
                            w = t(l) - x;
                            he = max(Rplt) - min(Rplt);
                            if (L1(state_begin==t) == 1)
                                col_shade = col_night;
                            else
                                col_shade = col_day;
                            end
                            rectangle('Position',[x y w he],'FaceColor',col_shade,'LineStyle','none'); hold on;
                            state_begin = t(l);
                            % --- stop plot rectangle ---
                        end

                        % state update
                        state_prev = L1(l);
                    end
                end
                % plot last rectangle
                % --- plot rectangle ---
                x = state_begin;
                y = min(Rplt);
                w = t(l) - x;
                he = max(Rplt) - min(Rplt);
                if (L1(state_begin==t) == 1)
                    col_shade = col_night;
                else
                    col_shade = col_day;
                end
                rectangle('Position',[x y w he],'FaceColor',col_shade,'LineStyle','none'); hold on;
                state_begin = t(l);
                % --- stop plot rectangle ---
                            
                
                % Plot small world metric
                t = linspace(0,(length(R1)*Ca.w)/(3600),length(R1));
                plot(t,R1,'black.','MarkerSize',1); hold on;
                
                % Plot filtered
                
                
                % interp nans
                R2 = R1;
                nanr = isnan(R1);
                xtmp = 1:numel(R1);
                R2(nanr) = interp1(xtmp(~nanr), R1(~nanr), xtmp(nanr));
                
                % imperfect final interp
                if (isnan(R2(end)))
                    R2(isnan(R2)) = nanmean(R2);
                end
                if (isnan(R2(1)))
                    R2(isnan(R2)) = nanmean(R2);
                end
                
                %[F1,~] = envelope(R2);
                per_min_center = 90; % minute
                pmw = 30; % minute bandwidth
                F1 = bandpass(R2,[(1/((per_min_center+(pmw/2))*60)) (1/((per_min_center-(pmw/2))*60))],(1/w));
                F1 = F1 - (nanmean(F1(1:100) - nanmean(R1(1:100))));
                %F1 = F1 - mean(F1) + mean(R1);
                if (trig_plot_filter)
                    plot(t,F1,'-','Color',[0 1 0]*0.7); hold on;
                end
                %plot(t,F1,'-','Color',[0 1 0]*0.7); hold on;
                
                % Plot discontinuities due to ambiguous
                if (~isempty(amb_start))
                    for p = 1:length(popped)
                        plot(t(popped(p))*[1 1],[min(Rplt) max(Rplt)],'red-'); hold on;
                    end
                end
                
                box off;
                ax = gca;
                set(ax,'TickDir','Out');
                xlabel('Time (hours)');
                xticks(t(1):24:t(end));
                ylabel(yltxt);
                %ylabel(sprintf('%s Coherence',metric(3:end)));
                axis tight;
                set(gca, 'Layer', 'top')
                %title(sprintf('-log_1_0(p)=%.1d',pvals2(rIdx)));
                
%                 if (~isempty(find(Ca.label_sleep == 2))) %(j == n_plot)
%                     return
%                 end

                if (j == 1)
                    X = R2;
                    Xlabel = yltxt;
                elseif (j > 1)
                    X = [X; R2];
                    Xlabel = [Xlabel; {yltxt}];
                end
            end
            set(gcf,'renderer','Painters');
            print(h,sprintf('%s/coh_sub-%i_metric-%i',dir_out_fig,i,metrici),'-depsc');
            print(h,sprintf('%s/coh_sub-%i_metric-%i',dir_out_fig,i,metrici),'-dpng');
            %return
            close(h)
            
            h = figure('Position',[0 0 1600 400],'visible','off');
            cwt(R2,minutes(1/60));
            print(h,sprintf('%s/coh_sub-%i_metric-%i_cwt',dir_out_fig,i,metrici),'-dpng');
            close(h)
            
            
            
            % manhattan plot
            h = figure('Position',[0,0,600,300],'visible','off'); hold all;
            plot(1:length(pvals2s),pvals2s,'.black','MarkerSize',1);
            p_thresh = (-1)*log10(0.05/(Ca.n_pairs*length(SubjectsSleep)));
            plot([1 length(pvals2s)],p_thresh*[1 1],'red--');
            ax = gca;
            set(ax,'TickDir','out');
            xlabel('Coherence Rank');
            ylabel('Significance, -log_1_0(p-value)')
            axis tight;
            title(sprintf('Subject %i',i))
            legend({'Significance','Threshold'},'Location','NorthWest');
            set(gcf,'renderer','Painters');
            print(h,sprintf('%s/manhattan_sub-%i_metric-%i',dir_out_fig,i,metrici),'-depsc');
            print(h,sprintf('%s/manhattan_sub-%i_metric-%i',dir_out_fig,i,metrici),'-dpng');
            close(h);
            %return
            
            
            % Save variables for classification
%             F = struct();
%             F.X = X;
%             F.Y = L1;
%             F.R = R;
%             features{i_feat} = F;
%             i_feat = i_feat + 1;

            
            % Train nearest neighbors
            n_obs = length(L1);
            
            % Define K-fold crossval K range
            %K = 2:4:40; %2:11; % full run
            K = 2; % testing only
            
            if (trig_svm)
                fprintf('[*] training SVM classifier..\n');
            else
                fprintf('[*] training KNN classifier..\n');
            end
            
            Loss = nan(1,length(K)*n_resample);
            cc = 1;
            for r = 1:n_resample
                [L_ri] = resample_label(L1);
                if (trig_svm)
                    mdl = fitcsvm(X(:,L_ri)',L1(L_ri)');
                    %mdl = fitcsvm(X(:,L_ri)',L1(L_ri)','OptimizeHyperparameters','all');
                else
                    %mdl = fitcknn(X(:,L_ri)',L1(L_ri)');
                    %mdl = fitcknn(X(:,L_ri)',L1(L_ri)','OptimizeHyperparameters','all');
                    mdl = fitcknn(X(:,L_ri)',L1(L_ri)','NumNeighbors',length(L1(L_ri))-1,'Distance','chebychev','DistanceWeight','squaredinverse');
                end
                for k = K
                    cvmdl = crossval(mdl,'KFold',k);
                    cvmdlloss = kfoldLoss(cvmdl);
                    Loss(cc) = cvmdlloss;
                    cc = cc + 1;
                end
            end
            fprintf('\tSmall-world k-fold accuracy: %.5f +- %.5f, n=%i\n',nanmean(1-Loss),nanstd(1-Loss),length(Loss))
            acc_sw_mean = nanmean(1-Loss);
            acc_sw_std = nanstd(1-Loss);
            n_k_sw = length(K);
            
            % coherence feature selection
            n_feat = 10;
            if (trig_loCoh)
                fIdx = sIdx(1:n_feat);
            else
                fIdx = sIdx((end-n_feat+1):end);
            end
            coh_p = pvals2(fIdx); % pvals2 - magnitude, pvals - raw pvalue
            % classify
            Loss = nan(1,length(K)*n_resample);
            cc = 1;
            for r = 1:n_resample
                [L_ri] = resample_label(L1);
                if (trig_svm)
                    mdl = fitcsvm(R(fIdx,L_ri)',L1(L_ri)');
                else
                    mdl = fitcknn(R(fIdx,L_ri)',L1(L_ri)');
                end
                for k = K
                    cvmdl = crossval(mdl,'KFold',k);
                    cvmdlloss = kfoldLoss(cvmdl);
                    Loss(cc) = cvmdlloss;
                    cc = cc + 1;
                end
            end
            
            %%
            % ------------------
            % null
            %
            Kn = 2;
            n_perm = 400;
            fprintf('[*] Permutation test started n=%i\n',n_perm);
            %Loss_null = nan(1,n_perm*length(Kn)*n_resample);
            Loss_null = nan(1,n_perm);
            %cc = 1;
            for ip = 1:n_perm
                %for r = 1:n_resample
                [L_ri] = resample_label(L1);
                X_n = R(fIdx,L_ri)';
                Y_n = L1(L_ri)';
                permIdx = randperm(length(Y_n));
                %X_n = X_n(permIdx,:);
                % shuffle labels only
                Y_n = Y_n(permIdx);
                if (trig_svm)
                    mdl_null = fitcsvm(X_n,Y_n);
                else
                    mdl_null = fitcknn(X_n,Y_n,'NumNeighbors',length(Y_n)-1,'Distance','chebychev','DistanceWeight','squaredinverse');
                end
                %for k = Kn
                cvmdl_null = crossval(mdl_null,'KFold',Kn);
                cvmdlloss_n = kfoldLoss(cvmdl_null);
                %Loss_null(cc) = cvmdlloss_n;
                %cc = cc + 1;
                %cc = (ip-1)*n_resample + r;
                Loss_null(ip) = cvmdlloss_n;
                %end
                %end
            end
            % ------------------
            
            fprintf('\tCoherence k-fold accuracy: %.5f +- %.5f, n=%i\n',nanmean(1-Loss),nanstd(1-Loss),length(Loss))
            acc_coh_mean = nanmean(1-Loss);
            acc_coh_std = nanstd(1-Loss);
            n_k_coh = length(K);
            
            
            
            %both
            Loss = nan(1,length(K)*n_resample);
            cc = 1;
            for r = 1:n_resample
                [L_ri] = resample_label(L1);
                if (trig_svm)
                    mdl = fitcsvm([X(:,L_ri); R(fIdx,L_ri)]',L1(L_ri)');
                else
                    mdl = fitcknn([X(:,L_ri); R(fIdx,L_ri)]',L1(L_ri)');
                end
                for k = K
                    cvmdl = crossval(mdl,'KFold',k);
                    cvmdlloss = kfoldLoss(cvmdl);
                    Loss(cc) = cvmdlloss;
                    cc = cc + 1;
                end
            end
            fprintf('\tBoth k-fold accuracy: %.5f +- %.5f, n=%i\n',nanmean(1-Loss),nanstd(1-Loss),length(Loss))
            acc_both_mean = nanmean(1-Loss);
            acc_both_std = nanstd(1-Loss);
            n_k_both = length(K);
            
            % show dataset numbers
            fprintf('\t#wake: %i, #sleep: %i\n',nansum(L1==0),nansum(L1==1));
            
            
            % Save vars
            Sub_i = [Sub_i, i];
            Acc_mean = [Acc_mean,[acc_sw_mean;acc_coh_mean;acc_both_mean]];
            Acc_null = [Acc_null; (1-Loss_null)];
            Acc_std = [Acc_std,[acc_sw_std;acc_coh_std;acc_both_std]];
            Acc_n = [Acc_n,[n_k_sw;n_k_coh;n_k_both]];
            Acc_nwake = [Acc_nwake,nansum(L1==0)];
            Acc_nsleep = [Acc_nsleep,nansum(L1==1)];
            Sigma_diff_sleep = [Sigma_diff_sleep,(mean(X(1,L1==1),2)-mean(X(1,L1==0)))/mean(X(1,L1==0),2)];
            [pSigma,~] = ranksum(X(1,L1==1),X(1,L1==0));
            Sigma_ranksum_p = [Sigma_ranksum_p,pSigma];
            Omega_diff_sleep = [Omega_diff_sleep,(mean(abs(X(2,L1==1)),2)-mean(abs(X(2,L1==0)),2))/mean(abs(X(2,L1==0)),2)];
            [pOmega,~] = ranksum(X(2,L1==1),X(2,L1==0));
            Omega_ranksum_p = [Omega_ranksum_p,pOmega];
            Rd = R(fIdx,1:length(L1));
            Coh_diff_sleep = [Coh_diff_sleep,(mean(Rd(:,L1==1),2)-mean(Rd(:,L1==0),2))./mean(Rd(:,L1==0),2)];
            Coh_ranksum_sig = [Coh_ranksum_sig,coh_p];
            i_feat = i_feat + 1;
%             K = 2:5;
%             for k = K
%                 %h = figure;
%                 % k-fold cross validation
%                 for k2 = 1:k
%                     % build test and train indices
%                     idx_test = false(1,n_obs);
%                     width = round((1/k)*n_obs);
%                     start = (k2-1)*width + 1;
%                     stop = start + width - 1;
%                     idx_test(start:stop) = true;
%                     X_train_coh = 1;
%                     idx_train = ~idx_test;
%                     
%                     n_roc = 10;
%                     NNeigh = linspace(1,sum(idx_train),n_roc);
%                     for n = NNeigh
%                         % train coherence
%                         Mdl_knn_coh = fitcknn(R(:,idx_train)',L1(idx_train)','NumNeighbors',n,'Standardize',1);
%                         % train small world
%                         Mdl_knn_sw = fitcknn(X(:,idx_train)',L1(idx_train)','NumNeighbors',n,'Standardize',1);
%                         % train both
%                         Mdl_knn_both = fitcknn([X(:,idx_train);R(:,idx_train)]',L1(idx_train)','NumNeighbors',n,'Standardize',1);
%                     
%                         % test coherence
%                         [label,score,cost] = predict(Mdl_knn_coh,R(:,idx_test)');
%                         return
%                     end
%                     
%                     %subplot(k,1,k2)
%                     %imagesc(idx_train)
%                 end
%                 
%                 return
%             end
            
        end
    end
    fprintf('[*] %i-subject Small-world accuracy: %.5f +- %.5f\n',i_feat-1,nanmean(Acc_mean(1,:)),nanstd(Acc_mean(1,:)));
    fprintf('[*] %i-subject Coherence accuracy: %.5f +- %.5f\n',i_feat-1,nanmean(Acc_mean(2,:)),nanstd(Acc_mean(2,:)));
    fprintf('[*] %i-subject Both accuracy: %.5f +- %.5f\n',i_feat-1,nanmean(Acc_mean(3,:)),nanstd(Acc_mean(3,:)));

    [p12,~] = ranksum(Acc_mean(1,:),Acc_mean(2,:));
    [p13,~] = ranksum(Acc_mean(1,:),Acc_mean(3,:));
    [p23,~] = ranksum(Acc_mean(2,:),Acc_mean(3,:));
    fprintf('[*] p-value Small-world vs Coherence means: %.4d\n',p12);
    fprintf('[*] p-value Small-world vs Both means: %.4d\n',p13);
    fprintf('[*] p-value Coherence vs Both means: %.4d\n',p23);
    
    if (trig_svm)
        ofname = sprintf('./cache/plot_sleep_sw_trainsvm_metric-%i_perm',metrici);
    else
        ofname = sprintf('./cache/plot_sleep_sw_train_metric-%i_perm',metrici);
    end
    
    if (trig_loCoh)
        ofname = sprintf('%s_lowCoh_perm',ofname);
    end
    
    save(ofname,'Acc_null','Acc_mean','Acc_std','Acc_n','Acc_nwake','Acc_nsleep','Sub_i','K','Sigma_diff_sleep','Sigma_ranksum_p','Omega_diff_sleep','Omega_ranksum_p','Coh_diff_sleep','Coh_ranksum_sig');
    
    if (~(contains(host,'Leibniz')))
        writematrix([Sub_i;Acc_mean;Acc_std;Sigma_diff_sleep;Sigma_ranksum_p;Omega_diff_sleep;Omega_ranksum_p;...
            Coh_diff_sleep;Coh_ranksum_sig;Acc_nsleep;Acc_nwake]',sprintf('./plot_sleep_sw_train_metric-%i.xlsx',metrici),'Sheet',1);
    end
    
end

function [LtIdx] = resample_label(L1)
    Lt = L1;
    LtIdx = 1:length(Lt);
    n1 = sum(Lt==1);
    n0 = sum(Lt==0);
    %cond_pass = (n1 ~= n0);
    if (n1 > n0)
        rlabel = 1;
    elseif (n1 < n0)
        rlabel = 0;
    end
    n_gap = abs(n1 - n0);
    rp = find(Lt==rlabel);
    rp_seed = randperm(length(rp)-n_gap+1);
    LtIdx(rp(rp_seed:(rp_seed+n_gap-1)))=[];
    %Lt(rp(rp_seed:(rp_seed+n_gap-1)))=[];
end


% janv. 30 2021 optimized hyperparams

%averages
% BoxConstraint: 12.4224
% KernelScale: 0.0104
% KernelFunction: gaussian
% PolynomialOrder: NaN
% Standardize: false

% Best estimated feasible point (according to models):
%     BoxConstraint    KernelScale    KernelFunction    PolynomialOrder    Standardize
%     _____________    ___________    ______________    _______________    ___________
% 
%        4.9569         0.0054154        gaussian             NaN             false 

% Best estimated feasible point (according to models):
%     BoxConstraint    KernelScale    KernelFunction    PolynomialOrder    Standardize
%     _____________    ___________    ______________    _______________    ___________
% 
%        41.919         0.016224         gaussian             NaN             false   

% Best estimated feasible point (according to models):
%     BoxConstraint    KernelScale    KernelFunction    PolynomialOrder    Standardize
%     _____________    ___________    ______________    _______________    ___________
% 
%        2.3843         0.0044148        gaussian             NaN             false  

% Best estimated feasible point (according to models):
%     BoxConstraint    KernelScale    KernelFunction    PolynomialOrder    Standardize
%     _____________    ___________    ______________    _______________    ___________
% 
%        0.42957        0.015655         gaussian             NaN             false   



% janvier 26 2021

% >> plot_sleep_sw_train
% Warning: Directory already exists. 
% > In plot_sleep_sw_train (line 48) 
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00019_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.57996 +- 0.00913, n=400
% 	Coherence k-fold accuracy: 0.63966 +- 0.01214, n=400
% 	Both k-fold accuracy: 0.67379 +- 0.00648, n=400
% 	#wake: 17102, #sleep: 13129
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00023_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.64921 +- 0.02396, n=400
% 	Coherence k-fold accuracy: 0.63766 +- 0.01092, n=400
% 	Both k-fold accuracy: 0.73143 +- 0.01238, n=400
% 	#wake: 51763, #sleep: 8752
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00024_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.64998 +- 0.04396, n=400
% 	Coherence k-fold accuracy: 0.72099 +- 0.06732, n=400
% 	Both k-fold accuracy: 0.78566 +- 0.05853, n=400
% 	#wake: 34210, #sleep: 13847
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00026_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.71767 +- 0.11392, n=400
% 	Coherence k-fold accuracy: 0.84925 +- 0.03391, n=400
% 	Both k-fold accuracy: 0.85222 +- 0.03008, n=400
% 	#wake: 13528, #sleep: 4047
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00030_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50676 +- 0.00319, n=400
% 	Coherence k-fold accuracy: 0.54658 +- 0.00698, n=400
% 	Both k-fold accuracy: 0.54892 +- 0.00839, n=400
% 	#wake: 18544, #sleep: 13970
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00032_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.56135 +- 0.02532, n=400
% 	Coherence k-fold accuracy: 0.55927 +- 0.02629, n=400
% 	Both k-fold accuracy: 0.57286 +- 0.03696, n=400
% 	#wake: 6052, #sleep: 3416
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00043_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.52493 +- 0.00444, n=400
% 	Coherence k-fold accuracy: 0.54620 +- 0.01372, n=400
% 	Both k-fold accuracy: 0.56169 +- 0.01321, n=400
% 	#wake: 21722, #sleep: 17903
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00049_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.58610 +- 0.05430, n=400
% 	Coherence k-fold accuracy: 0.57070 +- 0.02736, n=400
% 	Both k-fold accuracy: 0.62124 +- 0.05498, n=400
% 	#wake: 29832, #sleep: 11752
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00083_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.54909 +- 0.02479, n=400
% 	Coherence k-fold accuracy: 0.59466 +- 0.02234, n=400
% 	Both k-fold accuracy: 0.60956 +- 0.02577, n=400
% 	#wake: 31113, #sleep: 13808
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00084_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.63664 +- 0.03439, n=400
% 	Coherence k-fold accuracy: 0.83673 +- 0.02916, n=400
% 	Both k-fold accuracy: 0.84984 +- 0.02528, n=400
% 	#wake: 27474, #sleep: 12265
% [*] 10-subject Small-world accuracy: 0.59617 +- 0.06580
% [*] 10-subject Coherence accuracy: 0.65017 +- 0.11495
% [*] 10-subject Both accuracy: 0.68072 +- 0.11710
% [*] p-value Small-world vs Coherence means: 4.2736e-01
% [*] p-value Small-world vs Both means: 1.2122e-01
% [*] p-value Coherence vs Both means: 3.8467e-01

% Gamma

% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00019_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.57056 +- 0.00509, n=400
% 	Coherence k-fold accuracy: 0.60834 +- 0.00728, n=400
% 	Both k-fold accuracy: 0.62015 +- 0.00984, n=400
% 	#wake: 17102, #sleep: 13129
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00023_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.60433 +- 0.04688, n=400
% 	Coherence k-fold accuracy: 0.61042 +- 0.02133, n=400
% 	Both k-fold accuracy: 0.64698 +- 0.03935, n=400
% 	#wake: 51763, #sleep: 8752
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00024_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.67405 +- 0.04297, n=400
% 	Coherence k-fold accuracy: 0.77252 +- 0.05206, n=400
% 	Both k-fold accuracy: 0.79100 +- 0.05827, n=400
% 	#wake: 34210, #sleep: 13847
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00026_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.72490 +- 0.11147, n=400
% 	Coherence k-fold accuracy: 0.75679 +- 0.08508, n=400
% 	Both k-fold accuracy: 0.76031 +- 0.09087, n=400
% 	#wake: 13528, #sleep: 4047
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00030_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50622 +- 0.00252, n=400
% 	Coherence k-fold accuracy: 0.53030 +- 0.01126, n=400
% 	Both k-fold accuracy: 0.53297 +- 0.01027, n=400
% 	#wake: 18544, #sleep: 13970
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00032_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.55963 +- 0.01328, n=400
% 	Coherence k-fold accuracy: 0.53243 +- 0.01806, n=400
% 	Both k-fold accuracy: 0.54125 +- 0.01615, n=400
% 	#wake: 6052, #sleep: 3416
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00043_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.54382 +- 0.01087, n=400
% 	Coherence k-fold accuracy: 0.54209 +- 0.01104, n=400
% 	Both k-fold accuracy: 0.56049 +- 0.00744, n=400
% 	#wake: 21722, #sleep: 17903
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.60865 +- 0.04601, n=400
% 	Coherence k-fold accuracy: 0.53898 +- 0.02214, n=400
% 	Both k-fold accuracy: 0.63071 +- 0.05245, n=400
% 	#wake: 29832, #sleep: 11752
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.51963 +- 0.00290, n=400
% 	Coherence k-fold accuracy: 0.58767 +- 0.02323, n=400
% 	Both k-fold accuracy: 0.59170 +- 0.02312, n=400
% 	#wake: 31113, #sleep: 13808
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.82612 +- 0.06902, n=400
% 	Coherence k-fold accuracy: 0.83062 +- 0.02365, n=400
% 	Both k-fold accuracy: 0.87342 +- 0.02320, n=400
% 	#wake: 27474, #sleep: 12265
% [*] 10-subject Small-world accuracy: 0.61379 +- 0.10073
% [*] 10-subject Coherence accuracy: 0.63102 +- 0.11291
% [*] 10-subject Both accuracy: 0.65490 +- 0.11544
% [*] p-value Small-world vs Coherence means: 7.9134e-01
% [*] p-value Small-world vs Both means: 3.8467e-01
% [*] p-value Coherence vs Both means: 3.8467e-01



% jeudi janvier 21 2021

% Gamma

% >> plot_sleep_sw_train
% Warning: Directory already exists. 
% > In plot_sleep_sw_train (line 47) 
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00019_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.58392 +- 0.00924, n=100
% 	Coherence k-fold accuracy: 0.64172 +- 0.01345, n=100
% 	Both k-fold accuracy: 0.67164 +- 0.00643, n=100
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00023_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.64097 +- 0.02493, n=100
% 	Coherence k-fold accuracy: 0.64031 +- 0.00903, n=100
% 	Both k-fold accuracy: 0.72518 +- 0.00956, n=100
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00024_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.62935 +- 0.03845, n=100
% 	Coherence k-fold accuracy: 0.73451 +- 0.06696, n=100
% 	Both k-fold accuracy: 0.78606 +- 0.05548, n=100
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00026_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.71423 +- 0.14627, n=100
% 	Coherence k-fold accuracy: 0.85323 +- 0.04287, n=100
% 	Both k-fold accuracy: 0.85734 +- 0.04227, n=100
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00030_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50519 +- 0.00289, n=100
% 	Coherence k-fold accuracy: 0.54537 +- 0.00686, n=100
% 	Both k-fold accuracy: 0.54932 +- 0.00677, n=100
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00032_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.56211 +- 0.02623, n=100
% 	Coherence k-fold accuracy: 0.55466 +- 0.02596, n=100
% 	Both k-fold accuracy: 0.58565 +- 0.02262, n=100
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00043_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.52231 +- 0.00255, n=100
% 	Coherence k-fold accuracy: 0.54076 +- 0.00447, n=100
% 	Both k-fold accuracy: 0.56486 +- 0.01860, n=100
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00049_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.59405 +- 0.06148, n=100
% 	Coherence k-fold accuracy: 0.56682 +- 0.01999, n=100
% 	Both k-fold accuracy: 0.63757 +- 0.06221, n=100
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00083_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.53772 +- 0.02278, n=100
% 	Coherence k-fold accuracy: 0.58498 +- 0.01344, n=100
% 	Both k-fold accuracy: 0.61282 +- 0.02732, n=100
% [*] found cache, loading: ./cache/plot_sleep_sw_all_m00084_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.64770 +- 0.05124, n=100
% 	Coherence k-fold accuracy: 0.84616 +- 0.02829, n=100
% 	Both k-fold accuracy: 0.85051 +- 0.02741, n=100
% [*] 10-subject Small-world accuracy: 0.59375 +- 0.06507
% [*] 10-subject Coherence accuracy: 0.65085 +- 0.12027
% [*] 10-subject Both accuracy: 0.68410 +- 0.11512
% [*] p-value Small-world vs Coherence means: 3.4470e-01
% [*] p-value Small-world vs Both means: 7.5662e-02
% [*] p-value Coherence vs Both means: 4.2736e-01


% Broadband

% 
% > In plot_24h_sw_train (line 54) 
% [*] Behavior: ARM MOVEMENT
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00030
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.76072 +- 0.00122, n=10
% 	Coherence k-fold accuracy: 0.81899 +- 0.00146, n=10
% 	Both k-fold accuracy: 0.81954 +- 0.00064, n=10
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00043
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00083
% [*] 1-subject Small-world accuracy: 0.76072 +- 0.00000
% [*] 1-subject Coherence accuracy: 0.81899 +- 0.00000
% [*] 1-subject Both accuracy: 0.81954 +- 0.00000
% [*] p-value Small-world vs Coherence means: 0001
% [*] p-value Small-world vs Both means: 0001
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: CONTACT
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00030
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.76116 +- 0.00136, n=10
% 	Coherence k-fold accuracy: 0.81961 +- 0.00111, n=10
% 	Both k-fold accuracy: 0.81876 +- 0.00191, n=10
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00043
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00083
% [*] 1-subject Small-world accuracy: 0.76116 +- 0.00000
% [*] 1-subject Coherence accuracy: 0.81961 +- 0.00000
% [*] 1-subject Both accuracy: 0.81876 +- 0.00000
% [*] p-value Small-world vs Coherence means: 0001
% [*] p-value Small-world vs Both means: 0001
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: EATING
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00030
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.76135 +- 0.00106, n=10
% 	Coherence k-fold accuracy: 0.81923 +- 0.00084, n=10
% 	Both k-fold accuracy: 0.81896 +- 0.00126, n=10
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00043
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00083
% [*] 1-subject Small-world accuracy: 0.76135 +- 0.00000
% [*] 1-subject Coherence accuracy: 0.81923 +- 0.00000
% [*] 1-subject Both accuracy: 0.81896 +- 0.00000
% [*] p-value Small-world vs Coherence means: 0001
% [*] p-value Small-world vs Both means: 0001
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: HEAD MOVEMENT
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00030
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.76086 +- 0.00159, n=10
% 	Coherence k-fold accuracy: 0.81910 +- 0.00130, n=10
% 	Both k-fold accuracy: 0.81819 +- 0.00187, n=10
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00043
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00083
% [*] 1-subject Small-world accuracy: 0.76086 +- 0.00000
% [*] 1-subject Coherence accuracy: 0.81910 +- 0.00000
% [*] 1-subject Both accuracy: 0.81819 +- 0.00000
% [*] p-value Small-world vs Coherence means: 0001
% [*] p-value Small-world vs Both means: 0001
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: PATIENT IS TALKING
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00030
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.76082 +- 0.00161, n=10
% 	Coherence k-fold accuracy: 0.81851 +- 0.00127, n=10
% 	Both k-fold accuracy: 0.81956 +- 0.00085, n=10
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00043
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00083
% [*] 1-subject Small-world accuracy: 0.76082 +- 0.00000
% [*] 1-subject Coherence accuracy: 0.81851 +- 0.00000
% [*] 1-subject Both accuracy: 0.81956 +- 0.00000
% [*] p-value Small-world vs Coherence means: 0001
% [*] p-value Small-world vs Both means: 0001
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: QUIET
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00030
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.76124 +- 0.00079, n=10
% 	Coherence k-fold accuracy: 0.81905 +- 0.00142, n=10
% 	Both k-fold accuracy: 0.81960 +- 0.00089, n=10
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00043
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00083
% [*] 1-subject Small-world accuracy: 0.76124 +- 0.00000
% [*] 1-subject Coherence accuracy: 0.81905 +- 0.00000
% [*] 1-subject Both accuracy: 0.81960 +- 0.00000
% [*] p-value Small-world vs Coherence means: 0001
% [*] p-value Small-world vs Both means: 0001
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: SLEEP
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00030
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.76133 +- 0.00099, n=10
% 	Coherence k-fold accuracy: 0.81877 +- 0.00107, n=10
% 	Both k-fold accuracy: 0.81893 +- 0.00144, n=10
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00043
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00083
% [*] 1-subject Small-world accuracy: 0.76133 +- 0.00000
% [*] 1-subject Coherence accuracy: 0.81877 +- 0.00000
% [*] 1-subject Both accuracy: 0.81893 +- 0.00000
% [*] p-value Small-world vs Coherence means: 0001
% [*] p-value Small-world vs Both means: 0001
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: SOMEONE IS TALKING
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00030
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.76072 +- 0.00079, n=10
% 	Coherence k-fold accuracy: 0.81934 +- 0.00102, n=10
% 	Both k-fold accuracy: 0.81859 +- 0.00197, n=10
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00043
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00083
% [*] 1-subject Small-world accuracy: 0.76072 +- 0.00000
% [*] 1-subject Coherence accuracy: 0.81934 +- 0.00000
% [*] 1-subject Both accuracy: 0.81859 +- 0.00000
% [*] p-value Small-world vs Coherence means: 0001
% [*] p-value Small-world vs Both means: 0001
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: VIDEO GAMES
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00030
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.76139 +- 0.00101, n=10
% 	Coherence k-fold accuracy: 0.81895 +- 0.00178, n=10
% 	Both k-fold accuracy: 0.81924 +- 0.00094, n=10
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00043
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00083
% [*] 1-subject Small-world accuracy: 0.76139 +- 0.00000
% [*] 1-subject Coherence accuracy: 0.81895 +- 0.00000
% [*] 1-subject Both accuracy: 0.81924 +- 0.00000
% [*] p-value Small-world vs Coherence means: 0001
% [*] p-value Small-world vs Both means: 0001
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: WATCH TV
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00030
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.76088 +- 0.00135, n=10
% 	Coherence k-fold accuracy: 0.81881 +- 0.00230, n=10
% 	Both k-fold accuracy: 0.81902 +- 0.00149, n=10
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00043
% [*] found cache, loading: ./cache_plot_sleep_sw_all_m00083
% [*] 1-subject Small-world accuracy: 0.76088 +- 0.00000
% [*] 1-subject Coherence accuracy: 0.81881 +- 0.00000
% [*] 1-subject Both accuracy: 0.81902 +- 0.00000
% [*] p-value Small-world vs Coherence means: 0001
% [*] p-value Small-world vs Both means: 0001
% [*] p-value Coherence vs Both means: 0001
