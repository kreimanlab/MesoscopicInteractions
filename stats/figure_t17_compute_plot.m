close all;
clear;
rng('shuffle');

%/media/klab/internal/data/coreg/fsaverage_sym/label/all_surf_ielvis_m00037.label
dir_artLp = '/media/klab/internal/data/h5_notch20/art';
dir_corLp = '/media/klab/internal/data/coreg';
%dir_resLp = '/home/jerry/data/results/coh_w10';
setenv('SUBJECTS_DIR',dir_corLp);
dir_cacheLp = './cache';
dir_h5Lp = '/media/klab/internal/data/h5_notch20';
mkdir('figures/T17');

metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

load('cache/figure_t17_1_kmedoids_compute-all.mat');
% n_compute = 32;
% % 'K','Ssumd','Res','A2km';
% compute_K = cell(n_compute,1);
% compute_Ssumd = cell(n_compute,1);
% compute_Res = cell(n_compute,1);
% compute_A2km = cell(n_compute,1);

n_roi = length(cluster_i);
K = compute_K{1};
maj_mean = zeros(length(K),1);
maj_std = zeros(length(K),1);
maj_min = zeros(length(K),1);
maj_max = zeros(length(K),1);
for k = 1:length(K)
    cluster_mem = zeros(length(A2km),n_compute);
    cluster_sumd = zeros(1,n_compute);
    for compute = 1:n_compute
        %fprintf('[*] Compute stage %i of %i\n',compute,n_compute)
        cluster_mem(:,compute) = compute_Res{compute}{k,1};
        cluster_sumd(compute) = mean(compute_Res{compute}{k,3}); % mm
    end
    
%     majority = zeros(n_roi,1);
%     for j = 1:n_roi
%         majority(j) = sum(cluster_mem(j,:)==mode(cluster_mem(j,:)))/n_compute;
%     end
    majority = cluster_sumd;
    
    
    maj_mean(k) = mean(majority);
    maj_std(k) = std(majority);
    maj_min(k) = min(majority);
    maj_max(k) = max(majority);
    
    
end

h = figure('visible','on','PaperUnits','inches','PaperPosition',[0 0 6 6]);
set(gca,'TickDir','Out');
box off;
axis tight;
hold all;
trig_derivative = false;
if (trig_derivative)
    plot(K(1:(end-1)),diff(maj_mean),'black-');
    plot(K(1:(end-1)),diff(maj_mean + maj_std),'black--');
    plot(K(1:(end-1)),diff(maj_mean - maj_std),'black--');
else
    plot(K,maj_mean,'black-');
    plot(K,maj_mean + maj_std,'black--');
    plot(K,maj_mean - maj_std,'black--');
end

xlabel('K');
ylabel(sprintf('Mean within-cluster distance (unitless)'));
ax = gca;
XLim2 = ax.XLim;
YLim2 = ax.YLim;
XLim2 = sign(XLim2).*(ceil(abs(XLim2))); %round(ax.XLim);
%YLim2 = sign(YLim2).*(ceil(abs(YLim2))); %round(ax.YLim);
yvals = linspace(YLim2(1),YLim2(2),5);
ystrs = cell(size(yvals));
for i = 1:length(yvals)
    ystrs{i} = sprintf('%.0f',yvals(i));
end
try
    set(gca,'XTick', XLim2(1):2:XLim2(2) );
    set(gca,'YTick', yvals);
    set(gca,'YTickLabels', ystrs);
    set(gca,'XLim',XLim2);
    set(gca,'YLim',YLim2);
catch
end

set(gca,'FontName','Helvetica');
set(gca,'FontSize',10);
