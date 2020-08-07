close all;
clear;

% threshold by number of unique subjects
thresh_nsub = 1; %2;
thresh_npair = 5; %10;

dir_cacheLp = './cache';
cmap = jet(5);
msize = 3;
title_txt = {};
legend_txt = {};
metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
for iM_local = [1 5]
    metric = metricsp{iM_local};
    fn_cache = [dir_cacheLp,'/xsub_out_all_',num2str(iM_local)];
    Ca = load(fn_cache);
    fn_ca = sprintf('%s/fig_cluster2_A_%i.mat',dir_cacheLp,iM_local);
    Ca1 = load(fn_ca);
    fn_ca2 = sprintf('%s/fig_cluster2_reduce_%i_new.mat',dir_cacheLp,iM_local);
    Ca2 = load(fn_ca2);

    [n_A,~] = size(Ca2.A);
    count = 1;
    n_Acomb = nchoosek(n_A,2);
    Mag = zeros(n_Acomb,1);
    Ct = zeros(n_Acomb,1);
    Dist = zeros(n_Acomb,1);
    Nusub = zeros(n_Acomb,1);
    Nelec = zeros(n_Acomb,1);
    Npair = zeros(n_Acomb,1);
    Median = zeros(n_Acomb,1);
    for i = 1:(n_A-1)
        for j = (i+1):n_A
            es1 = Ca2.Es{i};
            es2 = Ca2.Es{j};
            es = [es1; es2];
            n_usub = length(unique(es(:,1)));
            %n_usub = min([length(unique(es1(:,1))), length(unique(es2(:,1)))]);
            n_elec = length(es(:,1));
            %n_elec = min([length(es1(:,1)), length(es2(:,1))]);
            Mag(count) = Ca2.A(i,j);
            Ct(count) = Ca2.Act(i,j);
            Median(count) = Ca2.A_median(i,j);
            Nusub(count) = Ca2.A_nusubs(i,j); %n_usub;
            Nelec(count) = n_elec;
            Npair(count) = Ca2.A_npairs(i,j);
            Dist(count) = Ca2.Ad(i,j);
            count = count + 1;
        end
    end
    Median(Median == 0) = NaN;

%     
%     h = figure;
%     fig_size_scale = 0.5;
%     set(h,'Position',round(fig_size_scale*[0 0 1080 1080]))
%     set(h,'PaperUnits','inches','PaperPosition',[0 0 4 4]);
%     subplot(2,1,1);
%     plot(Nusub,Median,'black.','MarkerSize',msize);
%     [r,p] = corr(Nusub(~isnan(Median)),Median(~isnan(Median)));
%     %title( sprintf( 'Pearson r=%.3f, p=%.3d',r,p ) );
%     box off;
%     xlabel('Number of subjects');
%     ylabel('Median Coherence')
%     set(gca,'TickDir','out');
%     set(gca,'FontSize',10)
%     axis tight;
%     
%     subplot(2,1,2);
%     plot(Npair,Median,'black.','MarkerSize',msize);
%     [r2,p2] = corr(Npair(~isnan(Median)),Median(~isnan(Median)));
%     %title( sprintf( 'Pearson r = %.3f, p=%.3d',r,p ) );
%     box off;
%     xlabel('Number of pairs');
%     ylabel('Coherence')
%     set(gca,'TickDir','out');
%     set(gca,'FontSize',10)
%     axis tight;
    
    %return
    
    
    h = figure;
    fig_size_scale = 0.5;
    set(h,'Position',round(fig_size_scale*[0 0 1080 1080]))
    set(h,'PaperUnits','inches','PaperPosition',[0 0 4 4]);
    subplot(2,1,1);
    plot(Nusub(Nusub>=thresh_nsub),Mag(Nusub>=thresh_nsub),'black.','MarkerSize',msize);
    sIdx = (Nusub>=thresh_nsub) & (~isnan(Mag));
    [r,p] = corr(Nusub(sIdx),Mag(sIdx));
    %title( sprintf( 'Pearson r=%.3f, p=%.3d',r,p ) );
    box off;
    xlabel('Number of subjects');
    ylabel('Coherence')
    set(gca,'TickDir','out');
    set(gca,'FontSize',10)
    axis tight;
    
    %h = figure;
    subplot(2,1,2);
    plot(Npair(Npair>=thresh_npair),Mag(Npair>=thresh_npair),'black.','MarkerSize',msize);
    sIdx = (Npair>=thresh_npair) & (~isnan(Mag));
    [r2,p2] = corr(Npair(sIdx),Mag(sIdx));
    %title( sprintf( 'Pearson r = %.3f, p=%.3d',r,p ) );
    box off;
    xlabel('Number of pairs');
    ylabel('Coherence')
    set(gca,'TickDir','out');
    set(gca,'FontSize',10)
    axis tight;
    
    % save
    savet = sprintf('usub_r-%.3d_p-%.3d_npair_r-%.3d_p-%.3d',r,p,r2,p2);
    savet = replace(savet,'.','p');
    print(h,['figures/cluster3/mag_',savet,'_',num2str(iM_local)],'-dpng','-r400');
    print(h,['figures/cluster3/mag_',savet,'_',num2str(iM_local)],'-depsc','-r400');
    
    
%     h = figure;
%     plot(Dist,Mag,'black.');
%     [r,p] = corr(Dist(~isnan(Mag)),Mag(~isnan(Mag)));
%     title( sprintf( 'Pearson r = %.3f, p=%.3d',r,p ) );
%     box off;
%     xlabel('Distance');
%     ylabel('Coherence')

    x = Ct;
    y = Mag;
    h = figure;
    plot(x,y,'.','Color',[0 0 0]); % ,cmap(iM_local,:)
    hold on;
    [r,p] = corr(x(~isnan(y)),y(~isnan(y)));
    %title( sprintf( '%s Pearson r = %.3f, p=%.3d',metric(3:end),r,p ) );
    box off;
    xlabel('X');
    ylabel('Y')

    % 
    % h = figure;
%     plot(Ct,Mag,'.','Color',cmap(iM_local,:)); hold on;
%     [r,p] = corr(Ct(~isnan(Mag)),Mag(~isnan(Mag)),'type','Pearson');
%     title_txt = [title_txt; { sprintf( '%s r = %.3f, p=%.3d',metric(3:end),r,p )}];
%     legend_txt = [legend_txt, {metric(3:end)}];
end
% 
title(title_txt);
%legend(legend_txt);
box off;
xlabel('Consistency across time');
ylabel('Coherence')