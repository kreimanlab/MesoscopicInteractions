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

%     dir_art = '/media/klab/untitled/h5_notch20/art_nosz'; %'/media/klab/internal/data/h5_notch20/art';
%     dir_res = '/media/klab/untitled/results/coh_w10';
%     dir_cor = '/media/klab/internal/data/coreg'; %'/mnt/cuenap_ssd/coregistration';
%     dir_h5 = '/media/klab/untitled/h5_notch20';
elseif contains(host,'ubuntu_1604')
    dir_art = '/nas_share/RawData/data/h5_notch20/art_nosz';
    dir_res = '/nas_share/RawData/scripts/opencl/results_res5hz';
    dir_cor = '/mnt/cuenap_ssd/coregistration';
    %dir_h5 = '/nas_share/RawData/data/h5_notch20';
    dir_h5 = '/mnt/cuenap/data/h5_notch20';
end



SubjectsSleep = {'m00030','m00043','m00083'};
%SubjectsSleep = {'m00043','m00083'};
%SubjectsSleep = {'m00030'};
%SubjectsSleep = {'m00083'};

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

trig_plot_filter = false;
dir_out_fig = './figures/24h_sw';
mkdir(dir_out_fig);
n_perm = 10000;
atl = 2;
for metrici = [1 5]
    metric = metrics{metrici};
    features = {};
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
            Ca2 = load(sprintf('./cache/%s_sleepsig_metric-%i_atl-%i.mat',sid,metrici,atl));
            
            Ca = load(sprintf('./cache/%s_24h_metric-%i_atl-%i.mat',sid,metrici,atl));
            Ca.Ca = Ca2.Ca;
            % TEMPORARY
            %Ca.w = 2;
            
            %return
            % Loop through behaviors
            for j2 = 1:Ca.n_beh
                state_skip = false; %used to skip unusable behavior labels
                
                Ca_pvals = Ca.Pvals{j2};
                Ca_label_sleep = Ca.Label_sleep{j2};
                
                if (~isempty(Ca_label_sleep))
                    fprintf('[*] plotting behavior: %s\n',Ca.T.AnnotLabels{j2});

                    pvals = Ca_pvals;
                    pvals2 = (-1)*log10(pvals);
                    pvals2 = pvals2(~isinf(pvals2));
                    [pvals2s,sIdx] = sort(pvals2);

                    % Load
                    fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
                    fn_dist = sprintf('%s/%s_dists-%s-%i.mat',dir_res,sid,metric,n_perm);
                    fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
                    %art_idx = h5read(fn_art,'/art_idx');
                    %art_idxB = (art_idx == 1);
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
                        n_graph = Ca.n_graph;
                        n_bchan = Ca.ecog.n_bchan;
                        Dmats = Ca.Ca.Dmats;
                        dist_thresh = Ca.Ca.dist_thresh;
                        coh_thresh = Ca.Ca.coh_thresh;
                        chan1 = Ca.chan1;
                        chan2 = Ca.chan2;
                        n_pairs = Ca.n_pairs;
                        for j = 1:n_graph
                            tic;
                            % count

                            Ahs = nan(n_bchan,n_bchan);
                            for c = 1:n_pairs
                                b1 = chan1(c) + 1;
                                b2 = chan2(c) + 1;
                                cond_dist = (Dmats(c) > dist_thresh);
                                cond_sig = (R(c,j) > coh_thresh(c));
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
                            
                            t_single = toc;
                            if (mod(j-1,round(n_graph/100)) == 0)
                                fprintf('[%s] window %i of %i (%.1f %%) ETA: %.2f hrs\n',...
                                    sid,j,n_graph,100*(j/n_graph),(n_graph-j)*t_single/3600);
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
%                         if (~isempty(amb_start))
%                             popped = nan(size(amb_start));
%                             popped_idx = false(size(amb_start));
%                             offset = 0;
%                             for p = 1:length(amb_start)
%                                 astart = amb_start(p);
%                                 astop = amb_stop(p);
%                                 popped_idx(astart:astop) = true;
%                                 popped(p) = amb_start(p) - offset;
%                                 offset = offset + (astop - astart + 1);
%                             end
%                             R1 = R1(~popped_idx);
%                             L1 = L1(~popped_idx);
%                             % last cut compensation
%                             if (popped(end) > length(R1))
%                                 popped(end) = length(R1);
%                             end
%                         end

                        % Pop ambiguous
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

                        if (isempty(R1))
                            state_skip = true;
                        else
                            
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
                            xticks(t(1):6:t(end));
                            ylabel(yltxt);
                            %ylabel(sprintf('%s Coherence',metric(3:end)));
                            axis tight;
                            set(gca, 'Layer', 'top')
                            %title(sprintf('-log_1_0(p)=%.1d',pvals2(rIdx)));

            %                 if (~isempty(find(Ca_label_sleep == 2))) %(j == n_plot)
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
                    end

                    if (~state_skip)
                        
                        bstr = replace(Ca.T.AnnotLabels{j2},' ','_');
                        set(gcf,'renderer','Painters');
                        print(h,sprintf('%s/coh_sub-%i_metric-%i_%s',dir_out_fig,i,metrici,bstr),'-depsc');
                        print(h,sprintf('%s/coh_sub-%i_metric-%i_%s',dir_out_fig,i,metrici,bstr),'-dpng');
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
                        legend({'Significance','Threshold'},'Location','NorthWest')
                        set(gcf,'renderer','Painters');
                        print(h,sprintf('%s/manhattan_sub-%i_metric-%i_%s',dir_out_fig,i,metrici,bstr),'-depsc');
                        print(h,sprintf('%s/manhattan_sub-%i_metric-%i_%s',dir_out_fig,i,metrici,bstr),'-dpng');
                        close(h);
                    end

                    %return
                end
            end
            
            % Save variables for classification
%             F = struct();
%             F.X = X;
%             F.Y = L1;
%             F.R = R;
%             features{i_feat} = F;
            i_feat = i_feat + 1;
        end
    end
end
