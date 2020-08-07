% Helper function called by plot_xsub_coh

close all;
%clear;

% Set paths
if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        resultsDir = '/media/klab/44/data/results';
        h5Dir = '/media/klab/44/h5';
    elseif strcmp(strip(hname),'ubuntu_1604')
        resultsDir = '/nas_share/RawData/data/results';
        h5Dir = '/nas_share/RawData/scripts/synth/out';
    else
        resultsDir = '/mnt/cuenap2/data/results';
        h5Dir = '/mnt/cuenap2/scripts/synth/out';
    end
end

% Read xsub file
%load('xsub-s_Destrieux-2009.mat');
%load('xsub-sP_Destrieux-2009.mat');
%load('xsub-scBroadband_Destrieux-2009.mat');
%load('xsub-s_Desikan-Killiany.mat');
%load('xsub-sP_Desikan-Killiany.mat');
%load('xsub-scBroadband_Desikan-Killiany.mat');
%load('xsub-sP_HCP-MMP1.mat');

%load('xsub/xsub-pcBroadband_Destrieux-2009.mat');
%load('xsub/xsub-pcBroadband_Desikan-Killiany.mat');

%load('xsub/xsub2-pcAlpha_atl-2.mat')
roi_excl = [1,5];
%m_threshold = 0.2;

% Load
clear A;
if (~exist('A','var'))
    A = Analysis(resultsDir,h5Dir);
end

% alias
atl_i = ATL_I;
atl_ind = atl_i;

at_name = AT.P.AtlNames{atl_ind};


% Trim off extra rois
e_i = true(1,length(adjct_dist));
e_i(roi_excl) = false;
adjct_dist = adjct_dist(e_i,e_i);
AdjCS = AdjCS(e_i,e_i);
AdjCS_mag = AdjCS_mag(e_i,e_i);
adjct_bchan = adjct_bchan(e_i,e_i);
adjct_sub = adjct_sub(e_i,e_i);
adjct_xsub = adjct_xsub(e_i,e_i);
adjct_xsubM = adjct_xsubM(e_i,e_i);
AT.P.AtlROIs{atl_ind}.RH.struct_names = AT.P.AtlROIs{atl_ind}.RH.struct_names(e_i);
AT.P.AtlROIs{atl_ind}.LH.struct_names = AT.P.AtlROIs{atl_ind}.LH.struct_names(e_i);
n_rois = n_rois - length(roi_excl);

% --- Clustering (output: cluster_i) --------------------------------------
Y = zeros(1,nchoosek(n_rois,2));
yc = 1;
for m1 = 1:(n_rois-1)
    for m2 = (m1+1):n_rois
        ad = adjct_dist{m2,m1};
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

Z = linkage(Y,'median');
cluster_i = optimalleaforder(Z,Y,'transformation','inverse'); % ,'transformation','inverse'
roi_dist = roi_dist(cluster_i,cluster_i);
%roi_dist(isinf(roi_dist)) = ;
clash = nansum(nansum(triu(roi_dist,1) - triu(roi_dist,2)));
fprintf('cluster clash: %.12f mm\n',clash)
%cluster_i = optimalleaforder(Z,Y);
%return
% -------------------------------------------------------------------------

% Adjust nans
AdjCS_mag((AdjCS == -1)) = NaN;
AdjCS((AdjCS == -1)) = NaN;
for i = 1:length(AdjCS)
    AdjCS_mag(i,i) = NaN;
end

AdjM = AdjCS_mag;

load('xsub_coh_isI');
%AdjM_isI = false(size(AdjM));
% Threshold metric
for i = 1:length(AdjM)
    for j = 1:length(AdjM)
        %if ((~isnan(AdjM(i,j))) && (AdjM(i,j) < m_threshold))
        %    AdjCS(i,j) = 0;
        %    AdjM(i,j) = m_threshold;
            %AdjCS(j,i) = 0;
%        elseif ((~isnan(AdjM(i,j))) && (AdjM(i,j) >= m_threshold))
%            AdjM_isI(i,j) = true;
        if ((~isnan(AdjM(i,j))) && (~AdjM_isI(i,j)))
            %AdjM(i,j) = m_threshold;
            AdjCS(i,j) = 0;
        end
    end
end


AT_c = AT;

% % Apply clustering
AdjM = AdjM(cluster_i,cluster_i);
AT_c.P.AtlROIs{atl_ind}.RH.struct_names = AT_c.P.AtlROIs{atl_ind}.RH.struct_names(cluster_i);
AT_c.P.AtlROIs{atl_ind}.LH.struct_names = AT_c.P.AtlROIs{atl_ind}.LH.struct_names(cluster_i);
AdjCS = AdjCS(cluster_i,cluster_i);
adjct_bchan = adjct_bchan(cluster_i,cluster_i);
adjct_sub = adjct_sub(cluster_i,cluster_i);
adjct_xsub = adjct_xsub(cluster_i,cluster_i);
adjct_xsubM = adjct_xsubM(cluster_i,cluster_i);
adjct_dist = adjct_dist(cluster_i,cluster_i);
AdjM_isI = AdjM_isI(cluster_i,cluster_i);

% compact
com_i = true(1,length(AdjM));
for i2 = 1:length(AdjM)
    if (sum(isnan(AdjM(i2,:))) >= (length(AdjM)-1))
        com_i(i2) = false;
    end
end
AdjM = AdjM(com_i,com_i);
AT_c.P.AtlROIs{atl_ind}.RH.struct_names = AT_c.P.AtlROIs{atl_ind}.RH.struct_names(com_i);
AT_c.P.AtlROIs{atl_ind}.LH.struct_names = AT_c.P.AtlROIs{atl_ind}.LH.struct_names(com_i);
AdjCS = AdjCS(com_i,com_i);
roi_dist = roi_dist(com_i,com_i);
adjct_bchan = adjct_bchan(com_i,com_i);
adjct_sub = adjct_sub(com_i,com_i);
adjct_xsub = adjct_xsub(com_i,com_i);
adjct_xsubM = adjct_xsubM(com_i,com_i);
adjct_dist = adjct_dist(com_i,com_i);
AdjM_isI = AdjM_isI(com_i,com_i);

% --- Plot ---
h = figure;
fig_w = 10;
fig_h = 10;
set(h,'Position',[0 0 fig_w*100 fig_h*100])
set(h, 'PaperUnits', 'Inches')
set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
ax1 = axes('Parent',h,'Units','Normalized','Position',[0.21 0.28 0.5 0.5]);

ax2 = A.plotAdjElectrodeDistThresh(ax1, 1, AdjM, ...
    sprintf('Functional Interactivity, %s, %s',...
    metric,at_name), inf(size(AdjM)), AdjCS, 0.5, AT_c, atl_ind); %A.const.P_VALUE_CT

% colorbar
c = colorbar(ax1,'Location','manual');
n_cbarticks = 5;
cmarks = linspace((c.Limits(1)),(c.Limits(2)),n_cbarticks);
cmarklabs = cell(1,n_cbarticks);
for j = 1:n_cbarticks
    cmarklabs{j} = sprintf('%.2f',cmarks(j));
end

set(c,'Position',[0.9 0.28 0.02 0.5]);
set(c,'FontSize',8);
set(c,'TicksMode','manual');
set(c,'Ticks',cmarks);
set(c,'TickLabelsMode','manual');
set(c,'TickLabels',cmarklabs);
set(c,'TickLength',0);

% colorbar
c2 = colorbar(ax2,'Location','southoutside');
%cmarks = round(linspace(c2.Limits(1),c2.Limits(2),3));
set(c2,'Position',[0.21 0.08 0.5 0.014]);
set(c2,'FontSize',8);
set(c2,'TicksMode','manual');
set(c2,'Ticks',cmarks);
set(c2,'TickLabelsMode','manual');
set(c2,'TickLabels',cmarklabs);
set(c2,'TickLength',0);

A.saveFig(h,sprintf('./figures/xsub_coh/xsub_coh_%s_%s',...
    metric,at_name));
%print(h,sprintf('./figures/adjacency_matrix_atlas/xsub_%s_%s',metric,at_name),'-depsc');


%close(h);

h2 = figure;
Me = AdjM(AdjM > m_threshold);
Di = roi_dist(AdjM > m_threshold);

Me = Me(~isinf(Di));
Di = Di(~isinf(Di));
plot(Di,Me,'black.');
[r,p] = corr(Me,Di);
[r_s,p_s] = corr(Me,Di,'Type','Spearman');
title(sprintf('Pearson: %.2f (p = %.2d), Spearman: %.2f (p = %0.2d)',r,p,r_s,p_s),'FontSize',6);
set(gca,'TickDir','out');
axis([min(Di) max(Di) min(Me) max(Me)*1.05]);
xlabel('Distance (mm)');
ylabel(sprintf('Metric (%s)',metric));
A.saveFig(h2,sprintf('./figures/xsub_coh/xsub_coh_%s_%s_dist',...
  metric,at_name));
