close all;
clear;

def_rois_short;

system('mkdir figures');
system('mkdir figures/T20_lognormal');

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
metrics_suffix = {'0.5-125 Hz','3-8 Hz','8-12 Hz','12-30 Hz','30-100 Hz'};
dir_cacheL = './cache';

perm_alpha_sec = 20;
cp_thresh_override = 0.05;
n_pairs_thresh = 10;% 10; % at least this many electrode pairs to be considered
n_subs_thresh = 2;% 2; % at least this many subjects to be considered
n_subs_ct_thresh = 0;% 2; % significant CTs in region pair must be from at least this many subjects

% Get atlas names
CaAtl = load(sprintf('./cache/xsub_out_all_%i',1));
AtlNames = CaAtl.C.AtlNames;
for atl = [2] %1:20

    %system(sprintf('mkdir figures/T14_allatl/atl%i_%s',atl,AtlNames{atl}));
    
    for iM = 1 %1:5 %1:length(metrics) % [1 5] %
        metric = metrics{iM};

        % Load human cache
        Ca_hum = load(sprintf('%s/xsub_out_all_%i_atl%i_sCorr.mat',dir_cacheL,iM,atl));
        n_rois = Ca_hum.n_rois;
        rois = Ca_hum.rois;

        % Calculate final functional interaction matrix
        Adj = nan(n_rois,n_rois);
        AdjMag = nan(n_rois,n_rois);
        AdjCT = nan(n_rois,n_rois);
        AdjMag4cl = nan(n_rois,n_rois);
        AdjCP = nan(n_rois,n_rois);
        AdjMagVar = nan(n_rois,n_rois);
        AdjCTVar = nan(n_rois,n_rois);
        AdjMagReS = nan(n_rois,n_rois);
        AdjMagL = cell(n_rois,n_rois);
        AdjMagL_sub = cell(n_rois,n_rois);
        %AdjCPVar = nan(n_rois,n_rois);
        N_bchan = nan(n_rois,n_rois);
        cp_thresh = cp_thresh_override;
        Dmat = Inf(n_rois,n_rois); % set to inf to avoid removing from every instance
        DistsAtl = [];

        for i1 = 1:n_rois
            for i2 = 1:n_rois
                AA = Ca_hum.AdjAtl{i1,i2};
                AAct = Ca_hum.AdjAtlCT{i1,i2};
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
                    %if ( ~isempty(AA(AA ~= 0)) )
                        %return;
                        Adj(i1,i2) = 1;
                        AdjMag(i1,i2) = nanmean(AA(AA ~= 0));
                        AdjMagVar(i1,i2) = nanvar(AA(AA ~= 0));
                        AdjCT(i1,i2) = nanmean(AAct(AA ~= 0));
                        AdjCTVar(i1,i2) = nanvar(AAct(AA ~= 0));
                        AdjMagL{i1,i2} = AA(AA~=0);
                        AdjMagL_sub{i1,i2} = AA_sub(AA~=0);
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
                        AdjCT(i1,i2) = 0;
                        AdjCTVar(i1,i2) = 0;
                        AdjMagL{i1,i2} = [];
                    end
                end
                
                % AdjMag for clustering
                AdjMag4cl(i1,i2) = mean(AA(AA~=0));
                
            end
        end
        
%         % Test Adj normalize
%         AdjMagNorm = nan(size(AdjMag));
%         for iN = 1:length(AdjMag)
%             
%         end
%         return
        
        % Histogram of log10(Coherence)
        n_rois = length(AdjMagL);
        AdjV = [];
        AdjVct = [];
        AdjVv = [];
        AdjVctv = [];
        for i1 = 1:(n_rois-1)
            for i2 = (i1+1):n_rois
                if ((i1 ~= 1) && (i2 ~= 1))
                    L = AdjMagL{i1,i2};
                    Ls = AdjMagL_sub{i1,i2}';
                    %L2 = L(L>0);
                    %L2 = L((L>0)&(Ls==21));
                    L2 = AdjMag(i1,i2);
                    L3 = AdjCT(i1,i2);
                    L4 = AdjMagVar(i1,i2);
                    L5 = AdjCTVar(i1,i2);
                    AdjV = [AdjV; L2];
                    AdjVct = [AdjVct; L3];
                    AdjVv = [AdjVv; L4];
                    AdjVctv = [AdjVctv; L5];
                end
            end
        end
        %histogram(log10(AdjV),50,'Normalization','pdf')
        AdjVct = AdjVct(AdjV>0);
        AdjVctv = AdjVctv(AdjV>0);
        AdjVv = AdjVv(AdjV>0);
        AdjV = AdjV(AdjV>0);
        %return
        
        % Show relationship between avg coherence and CT
        h = figure('visible','off');
        set(h,'Position',[0 0 200 200]);
        [r,pval] = corr(AdjV,AdjVct,'Type','Pearson');
        plot(AdjV,AdjVct,'black.');
        xlabel('Coherence');
        ylabel('Time Consistency');
        box off;
        set(gca,'TickDir','out');
        print(h,sprintf('figures/T20_lognormal/sCorr_metric%i_atl%i_mag_vs_ct_r-%i_p-%i',iM,atl,round(r*1000),round(pval*1000)),'-depsc');
        print(h,sprintf('figures/T20_lognormal/sCorr_metric%i_atl%i_mag_vs_ct_r-%i_p-%i',iM,atl,round(r*1000),round(pval*1000)),'-dpng');
        close(h);
        
        
        %return
        hh = figure('visible','off');
        set(hh,'Position',[0 0 200 200]);
        %h = histogram(log10(AdjV),'Normalization','pdf','DisplayStyle','bar','FaceColor',0.5*[1 1 1]);
        h = histogram(AdjV,'Normalization','pdf','DisplayStyle','bar','FaceColor',0.5*[1 1 1]);
        x_cen = linspace(h.BinLimits(1),h.BinLimits(2),100); %h.BinEdges + h.BinWidth/2;
        %y_norm = normpdf(x_cen,(nanmean(log10(AdjV))),(nanstd(log10(AdjV))));
        y_norm = lognpdf(x_cen,(nanmean(log(AdjV))),(nanstd(log(AdjV))));
        xticks([x_cen(1), 0.5*(x_cen(1)+x_cen(end)), x_cen(end)]);
        hold on;
        plot(x_cen,y_norm,'black-');
        
        % overlay linear fit
%         y_norm_linear = normpdf(x_cen,(nanmean((AdjV))),(nanstd((AdjV))));
%         hold on;
%         plot(x_cen,y_norm_linear,'blue--');

        
        hold on;
        y_norm2 = normpdf(x_cen,(nanmean((AdjV))),(nanstd((AdjV))));
        cc = plasma(2);
        plot(x_cen,y_norm2,':','color',cc(1,:));
        %return
        
        xlabel('Coherence');
        ylabel('pdf')
        set(gca,'TickDir','out');
        box off;
        
        % --- KS test ---
        Xt = (log10(AdjV));
        [h,p,ksstat,cv] = kstest((Xt-mean(Xt))/std(Xt));
        fprintf('[*] Log(C) kstest p=%.4d, n=%i\n',p,length(Xt));
        fprintf('\tmin:%.8f, max:%.8f\n',min(Xt),max(Xt));
        % --- Linear test ---
        Xt = (AdjV);
        [h,p,ksstat,cv] = kstest((Xt-mean(Xt))/std(Xt));
        fprintf('[*] Linear kstest p=%.4d, n=%i\n',p,length(Xt));
        fprintf('\tmin:%.8f, max:%.8f\n',min(Xt),max(Xt));
        
        
        % Broadband
        % [*] Log(C) kstest p=5.9133e-01, n=193
        % Gamma
        % [*] Log(C) kstest p=1.4780e-01, n=183
        %return
        print(hh,sprintf('figures/T20_lognormal/sCorr_metric%i_atl%i_hist',iM,atl),'-dsvg');
        %saveas(gcf,sprintf('figures/T20_lognormal/metric%i_atl%i_hist',iM,atl),'epsc');
        close(hh);
        
        %return
        
        % Cycle through ROIs
        for ir = 1:length(rois)
            %rois1 = rois;
            rois1 = rois_short;
            % area name
            if (false) %(atl == 2)
                roi_label = convertRoiDK(rois{ir});
                for irr = 1:length(rois1)
                    rois1{irr} = convertRoiDK(rois1{irr});
                end
            else
                roi_label = rois{ir};
            end
            
            % Connectivity profile of area
            Adj1 = AdjMagL(ir,:);
            
            % Get mean, SD
            Adj1_mu = zeros(size(Adj1));
            Adj1_sigma = zeros(size(Adj1));
            for ir2 = 1:length(Adj1_mu)
                Adj2 = Adj1{ir2};
                if (length(Adj2) < n_pairs_thresh)
                    Adj1_mu(ir2) = NaN;
                else
                    Adj1_mu(ir2) = nanmean(Adj2);
                end
                Adj1_sigma(ir2) = nanstd(Adj2);
            end

            % Sort
            [~,sIdx] = sort(Adj1_mu,'descend');
            Adj1 = Adj1(sIdx);
            Adj1_mu = Adj1_mu(sIdx);
            Adj1_sigma = Adj1_sigma(sIdx);
            rois1 = rois1(sIdx);

            % remove nans
            nIdx = isnan(Adj1_mu);
            Adj1 = Adj1(~nIdx);
            Adj1_mu = Adj1_mu(~nIdx);
            Adj1_sigma = Adj1_sigma(~nIdx);
            rois1 = rois1(~nIdx);

            % Plot
            if (length(Adj1) >= 10)
                close all;
                h = figure('visible','off');
                set(h,'Position',[0 0 400 200]);
%                 set(h,'PaperUnits','Inches');
%                 set(h,'PaperPosition',[0 0 6 4.5]);
                col_pt = 0.3*[1 1 1];
                col_ebar = 0.3*[1 1 1];
                hold all;
                Adj1_x = 1:length(Adj1);
                n_bip_pairs = 0;
                for ir2 = 1:length(Adj1_mu)
                    Adj2 = Adj1{ir2};
%                     for ir3 = 1:length(Adj2)
%                         plot([Adj1_x(ir2)],log10(Adj2(ir3)),'.','color',col_pt,'MarkerSize',1);
%                     end
                    %plot(Adj1_x(ir2),log10(nanmean(Adj2)),'o','color',col_pt);
                    %plot(Adj1_x(ir2)*[1 1],[log10(nanmean(Adj2)-nanstd(Adj2)) log10(nanmean(Adj2)+nanstd(Adj2))],'-','color',col_ebar);
                    plot(Adj1_x(ir2),(nanmean(Adj2)),'o','color',col_pt);
                    plot(Adj1_x(ir2)*[1 1],[(nanmean(Adj2)-nanstd(Adj2)) (nanmean(Adj2)+nanstd(Adj2))],'-','color',col_ebar);
                    n_bip_pairs = n_bip_pairs + length(Adj2);
                end

                % Expected normal profile
                n_rand = 100000;
                %AdjN = mean(sort(normrnd(mean(log10(Adj1_mu)),std(log10(Adj1_mu)),n_rand,length(Adj1_mu)),2,'descend'));
                AdjN = mean(sort(lognrnd(mean(log(Adj1_mu)),std(log(Adj1_mu)),n_rand,length(Adj1_mu)),2,'descend'));
                plot(Adj1_x,(AdjN),'black-','LineWidth',1.5);

                set(gca,'TickDir','out');
                axis tight;
                xlim([min(Adj1_x)-0.5, max(Adj1_x)+0.5]);
                ax = gca;
                yrange = ax.YLim(2) - ax.YLim(1);
                ymargin = 0.05*yrange;
                ylim([ax.YLim(1)-ymargin,ax.YLim(2)+ymargin]);
                xticks(Adj1_x);
                xticklabels(rois1);
                xtickangle(90);
                ylabel('Coherence');
                title(sprintf('Area %s',roi_label));

                print(h,sprintf('figures/T20_lognormal/sCorr_metric%i_atl%i_%s',iM,atl,rois{ir}),'-depsc');
                print(h,sprintf('figures/T20_lognormal/sCorr_metric%i_atl%i_%s',iM,atl,rois{ir}),'-dpng');
                close(h);
                
                fprintf('[%s] n_bip_pairs: %i\n',rois{ir},n_bip_pairs);
                if (strcmp(rois{ir},'superiortemporal'))
                    %return
                end
            end
            
            %return
        end
        
        

    end
end

%%
% --- KS test for 150 parcellation ---
atl = 0;
for iM = [1]
    Ca = load(sprintf('./cache/figure_t14_%i_150',iM));
    Adj = Ca.Adj_plt2;
    
    %Vectorize adjacency
    n_rois = length(Adj);
    AdjV = zeros(nchoosek(n_rois,2),1);
    count = 1;
    for i1 = 1:(n_rois-1)
        for i2 = (i1+1):n_rois
            AdjV(count) = Adj(i1,i2);
            count = count + 1;
        end
    end
    AdjV = AdjV(AdjV>0);
    
    hh = figure('visible','off');
    set(hh,'Position',[0 0 200 200]);
    %h = histogram(log10(AdjV),'Normalization','pdf','DisplayStyle','bar','FaceColor',0.5*[1 1 1]);
    h = histogram(AdjV,'Normalization','pdf','DisplayStyle','bar','FaceColor',0.5*[1 1 1]);
    x_cen = linspace(h.BinLimits(1),h.BinLimits(2),100); %h.BinEdges + h.BinWidth/2;
    
    y_norm = lognpdf(x_cen,(nanmean(log(AdjV))),(nanstd(log(AdjV))));
    xticks([x_cen(1), 0.5*(x_cen(1)+x_cen(end)), x_cen(end)]);
    hold on;
    plot(x_cen,y_norm,'black-');
    
    % overlay linear fit
%     y_norm_linear = normpdf(x_cen,(nanmean(log10(AdjV))),(nanstd(log10(AdjV))));
%     hold on;
%     plot(x_cen,y_norm_linear,'blue--');
    
    
    hold on;
    y_norm2 = normpdf(x_cen,(nanmean((AdjV))),(nanstd((AdjV))));
    cc = plasma(2);
    plot(x_cen,y_norm2,':','color',cc(1,:));
    
    xlabel('Coherence');
    ylabel('pdf')
    set(gca,'TickDir','out');
    box off;
    
    

    % --- KS test ---
    Xt = (log10(AdjV));
    [h,p,ksstat,cv] = kstest((Xt-mean(Xt))/std(Xt));
    fprintf('[*] Log(C) kstest p=%.4d, n=%i\n',p,length(Xt));
    fprintf('\tmin:%.8f, max:%.8f\n',min(Xt),max(Xt));
    
    
    % --- Linear KS test ---
    Xt = ((AdjV));
    [h,p,ksstat,cv] = kstest((Xt-mean(Xt))/std(Xt));
    fprintf('[*] Linear kstest p=%.4d, n=%i\n',p,length(Xt));
    fprintf('\tmin:%.8f, max:%.8f\n',min(Xt),max(Xt));
    % Broadband
    % [*] Log(C) kstest p=5.9133e-01, n=193
    % Gamma
    % [*] Log(C) kstest p=1.4780e-01, n=183
    %return
    print(hh,sprintf('figures/T20_lognormal/sCorr_150_metric%i_atl%i_hist',iM,atl),'-dsvg');
    %print(hh,sprintf('figures/T20_lognormal/sCorr_150_metric%i_atl%i_hist',iM,atl),'-depsc');
    %saveas(gcf,sprintf('figures/T20_lognormal/metric%i_atl%i_hist',iM,atl),'epsc');
    close(hh);
end



% 23 Oct 2020

% figure_t20_lognormal_sCorr
% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T20_lognormal’: File exists
% [*] Log(C) kstest p=7.4609e-02, n=294
% 	min:-0.48186548, max:-0.09345223
% [*] Linear kstest p=6.5008e-03, n=294
% 	min:0.32971182, max:0.80639489
% [fusiform] n_bip_pairs: 1553
% [inferiorparietal] n_bip_pairs: 702
% [inferiortemporal] n_bip_pairs: 2555
% [lateraloccipital] n_bip_pairs: 872
% [lateralorbitofrontal] n_bip_pairs: 537
% [lingual] n_bip_pairs: 726
% [middletemporal] n_bip_pairs: 2519
% [parsopercularis] n_bip_pairs: 309
% [parsorbitalis] n_bip_pairs: 285
% [parstriangularis] n_bip_pairs: 248
% [postcentral] n_bip_pairs: 540
% [precentral] n_bip_pairs: 433
% [rostralmiddlefrontal] n_bip_pairs: 976
% [superiorfrontal] n_bip_pairs: 671
% [superiorparietal] n_bip_pairs: 246
% [superiortemporal] n_bip_pairs: 1620
% [supramarginal] n_bip_pairs: 957
% [*] Log(C) kstest p=8.4181e-23, n=2387
% 	min:-0.87155125, max:-0.10260642
% [*] Linear kstest p=9.6007e-47, n=2387
% 	min:0.13441531, max:0.78957534