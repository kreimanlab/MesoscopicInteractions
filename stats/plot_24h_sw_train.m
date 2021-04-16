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
    dir_art = '/n/scratch3/users/j/jw324/opencl2/h5_notch20/art';
    dir_res = '/n/groups/kreiman/jerry/opencl/results_w2';
    dir_cor = '/n/scratch2/jw324/data/coreg';
    dir_h5 = '/n/scratch3/users/j/jw324/opencl2/h5_notch20';
elseif contains(host,'kraken')
    dir_art = '/media/jerry/internal/data/h5_notch20/art';
    dir_res = '/media/jerry/internal/data/MesoscopicInteractions-master/opencl/results/w2';
    dir_cor = '/media/klab/internal/data/coreg'; %'/mnt/cuenap_ssd/coregistration';
    dir_h5 = '/media/jerry/internal/data/h5_notch20';
elseif contains(host,'ubuntu_1604')
    dir_art = '/nas_share/RawData/data/h5_notch20/art_nosz';
    dir_res = '/nas_share/RawData/scripts/opencl/results_res5hz';
    dir_cor = '/mnt/cuenap_ssd/coregistration';
    %dir_h5 = '/nas_share/RawData/data/h5_notch20';
    dir_h5 = '/mnt/cuenap/data/h5_notch20';
end

SubjectsSleep = {'m00030','m00043','m00083'};

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

% classification params
trig_svm = false;
trig_loCoh = true;
%n_resample = 40;
n_resample = 40;

dir_out_fig = './figures/24h_sw';
mkdir(dir_out_fig);
n_perm = 10000;
atl = 2;
for metrici = 1 %[1 5]
    metric = metrics{metrici};
    %features = {};
    
    Ca = load(sprintf('./cache/%s_24h_metric-%i_atl-%i.mat',SubjectsSleep{1},metrici,atl));
    for j2 = [1 2 4 10]%1:Ca.n_beh
        state_skip = false;
        fprintf('[*] Behavior: %s\n',Ca.T.AnnotLabels{j2});
        
        Acc_mean = [];
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
                %Ca = load(sprintf('./cache/%s_sleepsig_metric-%i_atl-%i.mat',sid,metrici,atl));

                Ca = load(sprintf('./cache/%s_24h_metric-%i_atl-%i.mat',sid,metrici,atl));

                % TEMPORARY
                %Ca.w = 2;

                % Loop through behaviors
                %for j2 = 1:Ca.n_beh
                Ca_pvals = Ca.Pvals{j2};
                Ca_label_sleep = Ca.Label_sleep{j2};

                if (~isempty(Ca_label_sleep))

                    pvals = Ca_pvals;
                    pvals2 = (-1)*log10(pvals);
                    iidx = ~isinf(pvals2);
                    pvals2 = pvals2(iidx);
                    pvals = pvals(iidx);
                    [pvals2s,sIdx] = sort(pvals2);
                    
                    % bip pair labels
                    labelBP = cell(1,nchoosek(Ca.ecog.n_bchan,2));
                    ibp = 1;
                    for ibp1 = 1:(Ca.ecog.n_bchan-1)
                        for ibp2 = (ibp1 + 1):Ca.ecog.n_bchan
                            labelBP{ibp} = sprintf('%s-%s',Ca.C.AtlLabels{atl}{Ca.ecog.bip(ibp1,1)},Ca.C.AtlLabels{atl}{Ca.ecog.bip(ibp2,1)});
                            ibp = ibp + 1;
                        end
                    end

                    % Load
                    fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
                    fn_dist = sprintf('%s/%s_dists-%s-%i.mat',dir_res,sid,metric,n_perm);
                    fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
%                     art_idx = h5read(fn_art,'/art_idx');
%                     art_idxB = (art_idx == 1);
%                     R = h5read(fn_graph,'/R',[1 1],size(art_idx));
                    fn_cache = sprintf('%s/%s-art_idxB-w%i',dir_art,sid,Ca.w);
                    if (exist([fn_cache,'.mat'],'file'))
                        fprintf('[*] Reading art_idxB from: %s\n',fn_cache);
                        load(fn_cache);
                    else
                        fprintf(2,'[*] Missing file: %s.mat\n',fn_cache);
                    end
                    R = h5read(fn_graph,'/R',[1 1],size(art_idxB));

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
                    fn_cache_sw = sprintf('./cache/plot_24h_sw_all_%s_metric-%i',sid,metrici);
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
                        R1(Ca_label_sleep==2) = NaN;
                        L1 = Ca_label_sleep;
                        % remove ambiguous label
                        L1(L1==2) = NaN;

                        % harmonize length
                        if (length(R1) > length(L1))
                            R1 = R1(1:length(L1));
                        elseif (length(R1) < length(L1))
                            L1 = L1(1:length(R1));
                        end

                        % Cut ambiguous label and shrink
                        cond_same_length = ((length(Ca_label_sleep) == length(L1)) && (length(R1) == length(L1)));
                        if (~cond_same_length)
                            fprintf(2,'[!] Warn: variables R1, L1, Ca_label_sleep not the same length.\n');
                        end
                        state_begin = 1;
                        state_prev = Ca_label_sleep(1);
                        amb_start = [];
                        amb_stop = [];
                        % check for starting labels
                        if (state_prev == 2)
                            amb_start = [amb_start; 1];
                        end
                        for l = 1:length(Ca_label_sleep)
                            state_change = (Ca_label_sleep(l) ~= state_prev);

                            % start of ambiguous
                            if (state_change && (Ca_label_sleep(l) == 2))
                                amb_start = [amb_start; l];
                            % end of ambiguous 
                            elseif (state_change && (state_prev == 2))
                                amb_stop = [amb_stop; l];
                            end

                            state_prev = Ca_label_sleep(l);
                        end
                        % check for and close any hanging ambiguous labels
                        if (state_prev == 2)
                            amb_stop = [amb_stop; l];
                        end

                        % Pop ambiguous
                        if (~isempty(amb_start))
                            popped = nan(size(amb_start));
                            popped_idx = false(size(amb_start));
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


%                         % show day/night
%                         t = linspace(0,(length(R1)*Ca.w)/(3600),length(L1));
% 
%                         Rplt = R1(~isinf(R1));
%                         state_prev = L1(1);
%                         state_begin = t(1);
%                         col_day = [253 223 161]/255;
%                         col_night = [161 207 253]/255;
%                         for l = 1:length(L1)
%                             if (~isnan(L1(l)))
%                                 state_change = (L1(l) ~= state_prev);
% 
%                                 % shade rectangle on state change
%                                 if (state_change)
%                                     % --- plot rectangle ---
%                                     x = state_begin;
%                                     y = min(Rplt);
%                                     w = t(l) - x;
%                                     he = max(Rplt) - min(Rplt);
%                                     if (L1(state_begin==t) == 1)
%                                         col_shade = col_night;
%                                     else
%                                         col_shade = col_day;
%                                     end
%                                     rectangle('Position',[x y w he],'FaceColor',col_shade,'LineStyle','none'); hold on;
%                                     state_begin = t(l);
%                                     % --- stop plot rectangle ---
%                                 end
% 
%                                 % state update
%                                 state_prev = L1(l);
%                             end
%                         end
%                         % plot last rectangle
%                         % --- plot rectangle ---
%                         x = state_begin;
%                         y = min(Rplt);
%                         w = t(l) - x;
%                         he = max(Rplt) - min(Rplt);
%                         if (L1(state_begin==t) == 1)
%                             col_shade = col_night;
%                         else
%                             col_shade = col_day;
%                         end
%                         rectangle('Position',[x y w he],'FaceColor',col_shade,'LineStyle','none'); hold on;
%                         state_begin = t(l);
%                         % --- stop plot rectangle ---
% 
% 
%                         % Plot small world metric
%                         t = linspace(0,(length(R1)*Ca.w)/(3600),length(R1));
%                         plot(t,R1,'black.','MarkerSize',1); hold on;
% 
%                         % Plot filtered
% 
% 

                        if (isempty(R1))
                            state_skip = true;
                        else
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
    % 
    %                         %[F1,~] = envelope(R2);
    %                         per_min_center = 90; % minute
    %                         pmw = 30; % minute bandwidth
    %                         F1 = bandpass(R2,[(1/((per_min_center+(pmw/2))*60)) (1/((per_min_center-(pmw/2))*60))],(1/w));
    %                         F1 = F1 - (nanmean(F1(1:100) - nanmean(R1(1:100))));
    %                         %F1 = F1 - mean(F1) + mean(R1);
    %                         plot(t,F1,'-','Color',[0 1 0]*0.7); hold on;
    % 
    %                         % Plot discontinuities due to ambiguous
    %                         if (~isempty(amb_start))
    %                             for p = 1:length(popped)
    %                                 plot(t(popped(p))*[1 1],[min(Rplt) max(Rplt)],'red-'); hold on;
    %                             end
    %                         end
    % 
    %                         box off;
    %                         ax = gca;
    %                         set(ax,'TickDir','Out');
    %                         xlabel('Time (hours)');
    %                         xticks(t(1):24:t(end));
    %                         ylabel(yltxt);
    %                         %ylabel(sprintf('%s Coherence',metric(3:end)));
    %                         axis tight;
    %                         set(gca, 'Layer', 'top')
    %                         %title(sprintf('-log_1_0(p)=%.1d',pvals2(rIdx)));
    % 
    %         %                 if (~isempty(find(Ca_label_sleep == 2))) %(j == n_plot)
    %         %                     return
    %         %                 end
    % 
                            if (j == 1)
                                X = R2;
                                Xlabel = yltxt;
                            elseif (j > 1)
                                X = [X; R2];
                                Xlabel = [Xlabel; {yltxt}];
                            end
                        end
                    end

                    if (~state_skip)
                        bstr = replace(Ca.T.AnnotLabels{j2},' ','_');
    %                     print(h,sprintf('%s/coh_sub-%i_metric-%i_%s',dir_out_fig,i,metrici,bstr),'-depsc');
    %                     print(h,sprintf('%s/coh_sub-%i_metric-%i_%s',dir_out_fig,i,metrici,bstr),'-dpng');
    %                     %return
    %                     close(h)
    % 
    %                     h = figure('Position',[0 0 1600 400],'visible','off');
    %                     cwt(R2,minutes(1/60));
    %                     print(h,sprintf('%s/coh_sub-%i_metric-%i_cwt',dir_out_fig,i,metrici),'-dpng');
    %                     close(h)
    % 
    % 
    % 
    %                     % manhattan plot
    %                     h = figure('Position',[0,0,600,300],'visible','off'); hold all;
    %                     plot(1:length(pvals2s),pvals2s,'.black','MarkerSize',1);
    %                     p_thresh = (-1)*log10(0.05/(Ca.n_pairs*length(SubjectsSleep)));
    %                     plot([1 length(pvals2s)],p_thresh*[1 1],'red--');
    %                     ax = gca;
    %                     set(ax,'TickDir','out');
    %                     xlabel('Coherence Rank');
    %                     ylabel('Significance, -log_1_0(p-value)')
    %                     axis tight;
    %                     title(sprintf('Subject %i',i))
    %                     legend({'Significance','Threshold'},'Location','NorthWest')
    %                     print(h,sprintf('%s/manhattan_sub-%i_metric-%i_%s',dir_out_fig,i,metrici,bstr),'-depsc');
    %                     print(h,sprintf('%s/manhattan_sub-%i_metric-%i_%s',dir_out_fig,i,metrici,bstr),'-dpng');
    %                     close(h);
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
                        %K = 2:4:40;
                        K = 2:2:20;
                        %K = 2:6:20;
                        %K = [2,12];
                        %K = [2];
                        
                        % small-world
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
                                
                                mdl = fitcknn(X(:,L_ri)',L1(L_ri)','NumNeighbors',round(length(L1(L_ri)) * 0.01),'Distance','cityblock');
                                
                                %mdl = fitcknn(X(:,L_ri)',L1(L_ri)','NumNeighbors',length(L1(L_ri))-1,'Distance','chebychev','DistanceWeight','inverse');
%                                 mdl = fitcknn(X(:,L_ri)',L1(L_ri)','OptimizeHyperparameters','auto',...
%                                     'HyperparameterOptimizationOptions',...
%                                     struct('AcquisitionFunctionName','expected-improvement-plus'));
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
                        
                        % print area pairs
                        for ibp = 1:n_feat
                            fprintf('\t%s\n',labelBP{fIdx(ibp)});
                        end
                        
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
                        fprintf('\tCoherence k-fold accuracy: %.5f +- %.5f, n=%i\n',nanmean(1-Loss),nanstd(1-Loss),length(Loss))
                        acc_coh_mean = nanmean(1-Loss);
                        acc_coh_std = nanstd(1-Loss);
                        n_k_coh = length(K);

                        % both
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
            %end

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
            ofname = sprintf('./cache/plot_24h_sw_trainsvm_metric-%i_%s',metrici,bstr);
        else
            ofname = sprintf('./cache/plot_24h_sw_train_metric-%i_%s',metrici,bstr);
        end
        
        if (trig_loCoh)
            ofname = sprintf('%s_lowCoh',ofname);
        end
        save(ofname,'Acc_mean','Acc_std','Acc_n','Acc_nwake','Acc_nsleep','Sub_i','K','Sigma_diff_sleep','Sigma_ranksum_p','Omega_diff_sleep','Omega_ranksum_p','Coh_diff_sleep','Coh_ranksum_sig');
        writematrix([Sub_i;Acc_mean;Acc_std;Sigma_diff_sleep;Sigma_ranksum_p;Omega_diff_sleep;Omega_ranksum_p;...
            Coh_diff_sleep;Coh_ranksum_sig;Acc_nsleep;Acc_nwake]',sprintf('./plot_24h_sw_train_metric-%i_%s.xlsx',metrici,bstr),'Sheet',1);

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

% Best estimated feasible point (according to models):
%     BoxConstraint    KernelScale    KernelFunction    PolynomialOrder    Standardize
%     _____________    ___________    ______________    _______________    ___________
% 
%        187.95            NaN            linear              NaN             true 

% Best estimated feasible point (according to models):
%     BoxConstraint    KernelScale    KernelFunction    PolynomialOrder    Standardize
%     _____________    ___________    ______________    _______________    ___________
% 
%        6.2052          0.46669         gaussian             NaN             true  

% Best estimated feasible point (according to models):
%     BoxConstraint    KernelScale    KernelFunction    PolynomialOrder    Standardize
%     _____________    ___________    ______________    _______________    ___________
% 
%        949.44         0.055834         gaussian             NaN             false   

% Best estimated feasible point (according to models):
%     BoxConstraint    KernelScale    KernelFunction    PolynomialOrder    Standardize
%     _____________    ___________    ______________    _______________    ___________
% 
%        998.31            NaN          polynomial             4              false  

% Best estimated feasible point (according to models):
%     BoxConstraint    KernelScale    KernelFunction    PolynomialOrder    Standardize
%     _____________    ___________    ______________    _______________    ___________
% 
%        1.1706         0.0044479        gaussian             NaN             false 

% > In plot_24h_sw_train (line 47) 
% [*] Behavior: ARM MOVEMENT
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50469 +- 0.00472, n=400
% 	Coherence k-fold accuracy: 0.66374 +- 0.02378, n=400
% 	Both k-fold accuracy: 0.66189 +- 0.02419, n=400
% 	#wake: 38045, #sleep: 4848
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.53771 +- 0.02324, n=400
% 	Coherence k-fold accuracy: 0.71801 +- 0.03125, n=400
% 	Both k-fold accuracy: 0.71693 +- 0.03476, n=400
% 	#wake: 36505, #sleep: 6650
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.52335 +- 0.00790, n=400
% 	Coherence k-fold accuracy: 0.55562 +- 0.02862, n=400
% 	Both k-fold accuracy: 0.55501 +- 0.03132, n=400
% 	#wake: 18103, #sleep: 5851
% [*] 3-subject Small-world accuracy: 0.52191 +- 0.01655
% [*] 3-subject Coherence accuracy: 0.64579 +- 0.08267
% [*] 3-subject Both accuracy: 0.64461 +- 0.08234
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 7.0000e-01
% [*] Behavior: CONTACT
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.51190 +- 0.01483, n=400
% 	Coherence k-fold accuracy: 0.56388 +- 0.01522, n=400
% 	Both k-fold accuracy: 0.56406 +- 0.01595, n=400
% 	#wake: 42467, #sleep: 426
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.55382 +- 0.02415, n=400
% 	Coherence k-fold accuracy: 0.85272 +- 0.01059, n=400
% 	Both k-fold accuracy: 0.85433 +- 0.01132, n=400
% 	#wake: 40192, #sleep: 2963
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.51078 +- 0.00796, n=400
% 	Coherence k-fold accuracy: 0.56569 +- 0.04811, n=400
% 	Both k-fold accuracy: 0.57557 +- 0.05692, n=400
% 	#wake: 19345, #sleep: 4609
% [*] 3-subject Small-world accuracy: 0.52550 +- 0.02453
% [*] 3-subject Coherence accuracy: 0.66076 +- 0.16624
% [*] 3-subject Both accuracy: 0.66465 +- 0.16437
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 7.0000e-01
% [*] Behavior: EATING
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.49519 +- 0.00987, n=400
% 	Coherence k-fold accuracy: 0.67767 +- 0.03262, n=400
% 	Both k-fold accuracy: 0.68719 +- 0.03278, n=400
% 	#wake: 41855, #sleep: 1038
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.52982 +- 0.01741, n=400
% 	Coherence k-fold accuracy: 0.97535 +- 0.00872, n=400
% 	Both k-fold accuracy: 0.97664 +- 0.00936, n=400
% 	#wake: 42654, #sleep: 501
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.57468 +- 0.00756, n=400
% 	Coherence k-fold accuracy: 0.57876 +- 0.00720, n=400
% 	Both k-fold accuracy: 0.58016 +- 0.00819, n=400
% 	#wake: 22486, #sleep: 1468
% [*] 3-subject Small-world accuracy: 0.53323 +- 0.03985
% [*] 3-subject Coherence accuracy: 0.74393 +- 0.20643
% [*] 3-subject Both accuracy: 0.74800 +- 0.20512
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 7.0000e-01
% [*] Behavior: HEAD MOVEMENT
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50739 +- 0.00473, n=400
% 	Coherence k-fold accuracy: 0.64199 +- 0.01808, n=400
% 	Both k-fold accuracy: 0.63136 +- 0.01959, n=400
% 	#wake: 38538, #sleep: 4355
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.53111 +- 0.01826, n=400
% 	Coherence k-fold accuracy: 0.68672 +- 0.03564, n=400
% 	Both k-fold accuracy: 0.69395 +- 0.03631, n=400
% 	#wake: 39117, #sleep: 4038
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.55575 +- 0.02021, n=400
% 	Coherence k-fold accuracy: 0.63686 +- 0.02299, n=400
% 	Both k-fold accuracy: 0.63695 +- 0.02187, n=400
% 	#wake: 23779, #sleep: 175
% [*] 3-subject Small-world accuracy: 0.53142 +- 0.02418
% [*] 3-subject Coherence accuracy: 0.65519 +- 0.02743
% [*] 3-subject Both accuracy: 0.65409 +- 0.03463
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: PATIENT IS TALKING
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.49966 +- 0.01831, n=400
% 	Coherence k-fold accuracy: 0.68196 +- 0.07383, n=400
% 	Both k-fold accuracy: 0.70033 +- 0.07369, n=400
% 	#wake: 42287, #sleep: 606
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.53306 +- 0.01240, n=400
% 	Coherence k-fold accuracy: 0.87928 +- 0.06242, n=400
% 	Both k-fold accuracy: 0.87583 +- 0.06213, n=400
% 	#wake: 41583, #sleep: 1572
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.54591 +- 0.01167, n=400
% 	Coherence k-fold accuracy: 0.56671 +- 0.02818, n=400
% 	Both k-fold accuracy: 0.57025 +- 0.03070, n=400
% 	#wake: 23000, #sleep: 954
% [*] 3-subject Small-world accuracy: 0.52621 +- 0.02387
% [*] 3-subject Coherence accuracy: 0.70932 +- 0.15807
% [*] 3-subject Both accuracy: 0.71547 +- 0.15335
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: QUIET
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50390 +- 0.00328, n=400
% 	Coherence k-fold accuracy: 0.60554 +- 0.01577, n=400
% 	Both k-fold accuracy: 0.60201 +- 0.01686, n=400
% 	#wake: 34170, #sleep: 8723
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.55303 +- 0.03383, n=400
% 	Coherence k-fold accuracy: 0.68984 +- 0.10623, n=400
% 	Both k-fold accuracy: 0.70262 +- 0.10935, n=400
% 	#wake: 35426, #sleep: 7729
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.61057 +- 0.02220, n=400
% 	Coherence k-fold accuracy: 0.57422 +- 0.03004, n=400
% 	Both k-fold accuracy: 0.57353 +- 0.02518, n=400
% 	#wake: 20236, #sleep: 3718
% [*] 3-subject Small-world accuracy: 0.55583 +- 0.05339
% [*] 3-subject Coherence accuracy: 0.62320 +- 0.05980
% [*] 3-subject Both accuracy: 0.62606 +- 0.06782
% [*] p-value Small-world vs Coherence means: 4.0000e-01
% [*] p-value Small-world vs Both means: 4.0000e-01
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: SLEEP
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50469 +- 0.00304, n=400
% 	Coherence k-fold accuracy: 0.52713 +- 0.00942, n=400
% 	Both k-fold accuracy: 0.53130 +- 0.01249, n=400
% 	#wake: 31712, #sleep: 11181
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.55604 +- 0.02283, n=400
% 	Coherence k-fold accuracy: 0.57710 +- 0.04009, n=400
% 	Both k-fold accuracy: 0.58867 +- 0.04167, n=400
% 	#wake: 31602, #sleep: 11553
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.67138 +- 0.02315, n=400
% 	Coherence k-fold accuracy: 0.70669 +- 0.04723, n=400
% 	Both k-fold accuracy: 0.71532 +- 0.05356, n=400
% 	#wake: 22461, #sleep: 1493
% [*] 3-subject Small-world accuracy: 0.57737 +- 0.08537
% [*] 3-subject Coherence accuracy: 0.60364 +- 0.09268
% [*] 3-subject Both accuracy: 0.61176 +- 0.09416
% [*] p-value Small-world vs Coherence means: 7.0000e-01
% [*] p-value Small-world vs Both means: 7.0000e-01
% [*] p-value Coherence vs Both means: 7.0000e-01
% [*] Behavior: SOMEONE IS TALKING
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50452 +- 0.00322, n=400
% 	Coherence k-fold accuracy: 0.56286 +- 0.01253, n=400
% 	Both k-fold accuracy: 0.56830 +- 0.01383, n=400
% 	#wake: 32302, #sleep: 10591
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.54322 +- 0.01752, n=400
% 	Coherence k-fold accuracy: 0.73818 +- 0.04633, n=400
% 	Both k-fold accuracy: 0.74115 +- 0.04862, n=400
% 	#wake: 31925, #sleep: 11230
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.53454 +- 0.00991, n=400
% 	Coherence k-fold accuracy: 0.56967 +- 0.05817, n=400
% 	Both k-fold accuracy: 0.57294 +- 0.05103, n=400
% 	#wake: 7256, #sleep: 16698
% [*] 3-subject Small-world accuracy: 0.52743 +- 0.02031
% [*] 3-subject Coherence accuracy: 0.62357 +- 0.09931
% [*] 3-subject Both accuracy: 0.62746 +- 0.09848
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 7.0000e-01
% [*] Behavior: VIDEO GAMES
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50258 +- 0.01685, n=400
% 	Coherence k-fold accuracy: 0.78463 +- 0.08840, n=400
% 	Both k-fold accuracy: 0.80540 +- 0.08784, n=400
% 	#wake: 42349, #sleep: 544
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.53804 +- 0.02609, n=400
% 	Coherence k-fold accuracy: 0.93894 +- 0.00962, n=400
% 	Both k-fold accuracy: 0.93899 +- 0.00951, n=400
% 	#wake: 40650, #sleep: 2505
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.55098 +- 0.01886, n=400
% 	Coherence k-fold accuracy: 0.66203 +- 0.02910, n=400
% 	Both k-fold accuracy: 0.66632 +- 0.03326, n=400
% 	#wake: 23785, #sleep: 169
% [*] 3-subject Small-world accuracy: 0.53053 +- 0.02506
% [*] 3-subject Coherence accuracy: 0.79520 +- 0.13876
% [*] 3-subject Both accuracy: 0.80357 +- 0.13635
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 7.0000e-01
% [*] Behavior: WATCH TV
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50580 +- 0.00326, n=400
% 	Coherence k-fold accuracy: 0.57334 +- 0.01114, n=400
% 	Both k-fold accuracy: 0.57359 +- 0.01159, n=400
% 	#wake: 36982, #sleep: 5911
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.54778 +- 0.01982, n=400
% 	Coherence k-fold accuracy: 0.74678 +- 0.03371, n=400
% 	Both k-fold accuracy: 0.74790 +- 0.03310, n=400
% 	#wake: 33858, #sleep: 9297
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-1
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50197 +- 0.01707, n=400
% 	Coherence k-fold accuracy: 0.74110 +- 0.07885, n=400
% 	Both k-fold accuracy: 0.73730 +- 0.07824, n=400
% 	#wake: 23025, #sleep: 929
% [*] 3-subject Small-world accuracy: 0.51852 +- 0.02542
% [*] 3-subject Coherence accuracy: 0.68707 +- 0.09853
% [*] 3-subject Both accuracy: 0.68627 +- 0.09772
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: ARM MOVEMENT
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50619 +- 0.00381, n=400
% 	Coherence k-fold accuracy: 0.64013 +- 0.03200, n=400
% 	Both k-fold accuracy: 0.64145 +- 0.02705, n=400
% 	#wake: 38045, #sleep: 4848
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.53161 +- 0.02271, n=400
% 	Coherence k-fold accuracy: 0.71595 +- 0.04800, n=400
% 	Both k-fold accuracy: 0.74712 +- 0.03888, n=400
% 	#wake: 36505, #sleep: 6650
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.52785 +- 0.00685, n=400
% 	Coherence k-fold accuracy: 0.56048 +- 0.02701, n=400
% 	Both k-fold accuracy: 0.56451 +- 0.02920, n=400
% 	#wake: 18103, #sleep: 5851
% [*] 3-subject Small-world accuracy: 0.52188 +- 0.01372
% [*] 3-subject Coherence accuracy: 0.63885 +- 0.07775
% [*] 3-subject Both accuracy: 0.65103 +- 0.09168
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 7.0000e-01
% [*] Behavior: CONTACT
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.51348 +- 0.01270, n=400
% 	Coherence k-fold accuracy: 0.66214 +- 0.02177, n=400
% 	Both k-fold accuracy: 0.66266 +- 0.03192, n=400
% 	#wake: 42467, #sleep: 426
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.55896 +- 0.03173, n=400
% 	Coherence k-fold accuracy: 0.76989 +- 0.04997, n=400
% 	Both k-fold accuracy: 0.77732 +- 0.05237, n=400
% 	#wake: 40192, #sleep: 2963
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50654 +- 0.00653, n=400
% 	Coherence k-fold accuracy: 0.57080 +- 0.05039, n=400
% 	Both k-fold accuracy: 0.58469 +- 0.06438, n=400
% 	#wake: 19345, #sleep: 4609
% [*] 3-subject Small-world accuracy: 0.52632 +- 0.02847
% [*] 3-subject Coherence accuracy: 0.66761 +- 0.09966
% [*] 3-subject Both accuracy: 0.67489 +- 0.09690
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 7.0000e-01
% [*] Behavior: EATING
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.51561 +- 0.00985, n=400
% 	Coherence k-fold accuracy: 0.66889 +- 0.03877, n=400
% 	Both k-fold accuracy: 0.67125 +- 0.04675, n=400
% 	#wake: 41855, #sleep: 1038
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.53320 +- 0.03306, n=400
% 	Coherence k-fold accuracy: 0.97784 +- 0.00975, n=400
% 	Both k-fold accuracy: 0.97832 +- 0.00817, n=400
% 	#wake: 42654, #sleep: 501
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.55861 +- 0.00609, n=400
% 	Coherence k-fold accuracy: 0.57842 +- 0.00816, n=400
% 	Both k-fold accuracy: 0.58123 +- 0.00755, n=400
% 	#wake: 22486, #sleep: 1468
% [*] 3-subject Small-world accuracy: 0.53580 +- 0.02162
% [*] 3-subject Coherence accuracy: 0.74172 +- 0.20943
% [*] 3-subject Both accuracy: 0.74360 +- 0.20820
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 7.0000e-01
% [*] Behavior: HEAD MOVEMENT
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50690 +- 0.00391, n=400
% 	Coherence k-fold accuracy: 0.63998 +- 0.02962, n=400
% 	Both k-fold accuracy: 0.62893 +- 0.02314, n=400
% 	#wake: 38538, #sleep: 4355
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.53461 +- 0.01469, n=400
% 	Coherence k-fold accuracy: 0.69021 +- 0.04814, n=400
% 	Both k-fold accuracy: 0.69345 +- 0.04677, n=400
% 	#wake: 39117, #sleep: 4038
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.57924 +- 0.03179, n=400
% 	Coherence k-fold accuracy: 0.61051 +- 0.02478, n=400
% 	Both k-fold accuracy: 0.60509 +- 0.02114, n=400
% 	#wake: 23779, #sleep: 175
% [*] 3-subject Small-world accuracy: 0.54025 +- 0.03650
% [*] 3-subject Coherence accuracy: 0.64690 +- 0.04029
% [*] 3-subject Both accuracy: 0.64249 +- 0.04572
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: PATIENT IS TALKING
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.48277 +- 0.01146, n=400
% 	Coherence k-fold accuracy: 0.67398 +- 0.07243, n=400
% 	Both k-fold accuracy: 0.66441 +- 0.08370, n=400
% 	#wake: 42287, #sleep: 606
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.54100 +- 0.01665, n=400
% 	Coherence k-fold accuracy: 0.85388 +- 0.07492, n=400
% 	Both k-fold accuracy: 0.85823 +- 0.05703, n=400
% 	#wake: 41583, #sleep: 1572
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50618 +- 0.00857, n=400
% 	Coherence k-fold accuracy: 0.56081 +- 0.02619, n=400
% 	Both k-fold accuracy: 0.55901 +- 0.02579, n=400
% 	#wake: 23000, #sleep: 954
% [*] 3-subject Small-world accuracy: 0.50998 +- 0.02930
% [*] 3-subject Coherence accuracy: 0.69622 +- 0.14779
% [*] 3-subject Both accuracy: 0.69388 +- 0.15178
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: QUIET
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50547 +- 0.00324, n=400
% 	Coherence k-fold accuracy: 0.60020 +- 0.01997, n=400
% 	Both k-fold accuracy: 0.59514 +- 0.01913, n=400
% 	#wake: 34170, #sleep: 8723
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.55896 +- 0.03987, n=400
% 	Coherence k-fold accuracy: 0.69569 +- 0.08273, n=400
% 	Both k-fold accuracy: 0.69991 +- 0.09056, n=400
% 	#wake: 35426, #sleep: 7729
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.54874 +- 0.01445, n=400
% 	Coherence k-fold accuracy: 0.53418 +- 0.03386, n=400
% 	Both k-fold accuracy: 0.53614 +- 0.03281, n=400
% 	#wake: 20236, #sleep: 3718
% [*] 3-subject Small-world accuracy: 0.53772 +- 0.02840
% [*] 3-subject Coherence accuracy: 0.61003 +- 0.08120
% [*] 3-subject Both accuracy: 0.61040 +- 0.08294
% [*] p-value Small-world vs Coherence means: 4.0000e-01
% [*] p-value Small-world vs Both means: 4.0000e-01
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: SLEEP
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50039 +- 0.00335, n=400
% 	Coherence k-fold accuracy: 0.53853 +- 0.00827, n=400
% 	Both k-fold accuracy: 0.54051 +- 0.00760, n=400
% 	#wake: 31712, #sleep: 11181
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.55061 +- 0.02767, n=400
% 	Coherence k-fold accuracy: 0.57971 +- 0.04587, n=400
% 	Both k-fold accuracy: 0.59518 +- 0.04828, n=400
% 	#wake: 31602, #sleep: 11553
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.56352 +- 0.00680, n=400
% 	Coherence k-fold accuracy: 0.69610 +- 0.03499, n=400
% 	Both k-fold accuracy: 0.69188 +- 0.03757, n=400
% 	#wake: 22461, #sleep: 1493
% [*] 3-subject Small-world accuracy: 0.53817 +- 0.03335
% [*] 3-subject Coherence accuracy: 0.60478 +- 0.08172
% [*] 3-subject Both accuracy: 0.60919 +- 0.07665
% [*] p-value Small-world vs Coherence means: 4.0000e-01
% [*] p-value Small-world vs Both means: 4.0000e-01
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: SOMEONE IS TALKING
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50369 +- 0.00302, n=400
% 	Coherence k-fold accuracy: 0.62816 +- 0.02089, n=400
% 	Both k-fold accuracy: 0.62596 +- 0.01769, n=400
% 	#wake: 32302, #sleep: 10591
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.54589 +- 0.01594, n=400
% 	Coherence k-fold accuracy: 0.71811 +- 0.02035, n=400
% 	Both k-fold accuracy: 0.71437 +- 0.03149, n=400
% 	#wake: 31925, #sleep: 11230
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.52135 +- 0.00944, n=400
% 	Coherence k-fold accuracy: 0.53159 +- 0.03291, n=400
% 	Both k-fold accuracy: 0.52098 +- 0.02768, n=400
% 	#wake: 7256, #sleep: 16698
% [*] 3-subject Small-world accuracy: 0.52364 +- 0.02119
% [*] 3-subject Coherence accuracy: 0.62595 +- 0.09328
% [*] 3-subject Both accuracy: 0.62044 +- 0.09681
% [*] p-value Small-world vs Coherence means: 2.0000e-01
% [*] p-value Small-world vs Both means: 4.0000e-01
% [*] p-value Coherence vs Both means: 7.0000e-01
% [*] Behavior: VIDEO GAMES
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.48847 +- 0.01374, n=400
% 	Coherence k-fold accuracy: 0.79376 +- 0.08367, n=400
% 	Both k-fold accuracy: 0.78527 +- 0.08389, n=400
% 	#wake: 42349, #sleep: 544
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.53778 +- 0.02540, n=400
% 	Coherence k-fold accuracy: 0.95163 +- 0.01597, n=400
% 	Both k-fold accuracy: 0.94973 +- 0.01475, n=400
% 	#wake: 40650, #sleep: 2505
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.54332 +- 0.02171, n=400
% 	Coherence k-fold accuracy: 0.65101 +- 0.02813, n=400
% 	Both k-fold accuracy: 0.65208 +- 0.02744, n=400
% 	#wake: 23785, #sleep: 169
% [*] 3-subject Small-world accuracy: 0.52319 +- 0.03020
% [*] 3-subject Coherence accuracy: 0.79880 +- 0.15037
% [*] 3-subject Both accuracy: 0.79569 +- 0.14910
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 0001
% [*] Behavior: WATCH TV
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00030-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00030_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50397 +- 0.00390, n=400
% 	Coherence k-fold accuracy: 0.57335 +- 0.00780, n=400
% 	Both k-fold accuracy: 0.57121 +- 0.00789, n=400
% 	#wake: 36982, #sleep: 5911
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00043-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00043_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.53124 +- 0.02225, n=400
% 	Coherence k-fold accuracy: 0.68483 +- 0.03941, n=400
% 	Both k-fold accuracy: 0.68074 +- 0.04015, n=400
% 	#wake: 33858, #sleep: 9297
% [*] Reading art_idxB from: /media/jerry/internal/data/h5_notch20/art/m00083-art_idxB-w2
% [*] found cache, loading: ./cache/plot_24h_sw_all_m00083_metric-5
% [*] training KNN classifier..
% 	Small-world k-fold accuracy: 0.50903 +- 0.00769, n=400
% 	Coherence k-fold accuracy: 0.79470 +- 0.07307, n=400
% 	Both k-fold accuracy: 0.80645 +- 0.07251, n=400
% 	#wake: 23025, #sleep: 929
% [*] 3-subject Small-world accuracy: 0.51475 +- 0.01450
% [*] 3-subject Coherence accuracy: 0.68429 +- 0.11068
% [*] 3-subject Both accuracy: 0.68613 +- 0.11771
% [*] p-value Small-world vs Coherence means: 1.0000e-01
% [*] p-value Small-world vs Both means: 1.0000e-01
% [*] p-value Coherence vs Both means: 0001
