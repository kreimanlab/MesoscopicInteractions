close all;

if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        resultsDir = '/media/klab/44/data/results';
        h5Dir = '/media/klab/44/h5';
    else
        resultsDir = '/mnt/cuenap2/data/results';
        h5Dir = '/mnt/cuenap2/scripts/synth/out';
    end
end

if (~exist('A','var'))
    A = Analysis(resultsDir,h5Dir);
end

% Read xsub file
%load('xsub-s_Destrieux-2009.mat');
%load('xsub-sP_Destrieux-2009.mat');
%load('xsub-scBroadband_Destrieux-2009.mat');
%load('xsub-s_Desikan-Killiany.mat');
%load('xsub-sP_Desikan-Killiany.mat');
%load('result_xsub-scBroadband_Desikan-Killiany.mat');
load('xsub-scGamma_M132.mat');

% alias
atl_ind = atl_i;

at_name = AT.P.AtlNames{atl_ind};
h = figure;
fig_w = 10.5;
fig_h = 10.5;
set(h,'Position',[0 0 fig_w*100 fig_h*100])
set(h, 'PaperUnits', 'Inches')
set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
ax1 = axes('Parent',h,'Units','Normalized','Position',[0.2 0.24 0.55 0.55]);
%ax2 = A.plotAdjElectrode(ax1, i, dmat, sprintf('Distance Matrix: %s',sid));
%AdjM = AdjMa{atl_ind};

% Adjust nans
AdjCS_mag((AdjCS == -1)) = NaN;
AdjCS((AdjCS == -1)) = NaN;
for i = 1:length(AdjCS)
    AdjCS_mag(i,i) = NaN;
end

% plot
AdjM = AdjCS_mag;

AT_c = AT;
ax2 = A.plotAdjElectrodeDistThresh(ax1, 1, AdjM, ...
    sprintf('Functional Interactivity %s, %s, %s',sid,...
    metric,at_name), inf(size(AdjM)), AdjCS, 0.5, AT_c, atl_ind); %A.const.P_VALUE_CT
% colorbar
c = colorbar(ax1,'Location','manual');
n_cbarticks = 5;
cmarks = linspace((c.Limits(1)),(c.Limits(2)),n_cbarticks);
cmarklabs = cell(1,n_cbarticks);
for j = 1:n_cbarticks
    cmarklabs{j} = sprintf('%.2f',cmarks(j));
end

set(c,'Position',[0.9 0.24 0.02 0.553]);
set(c,'FontSize',8);
set(c,'TicksMode','manual');
set(c,'Ticks',cmarks);
set(c,'TickLabelsMode','manual');
set(c,'TickLabels',cmarklabs);
set(c,'TickLength',0);

% colorbar
c2 = colorbar(ax2,'Location','southoutside');
%cmarks = round(linspace(c2.Limits(1),c2.Limits(2),3));
set(c2,'Position',[0.2 0.08 0.55 0.014]);
set(c2,'FontSize',8);
set(c2,'TicksMode','manual');
set(c2,'Ticks',cmarks);
set(c2,'TickLabelsMode','manual');
set(c2,'TickLabels',cmarklabs);
set(c2,'TickLength',0);
A.saveFig(h,sprintf('./figures/adjacency_matrix_atlas/xsub_%s_%s',...
  metric,at_name));
%close(h);

% compare to macaque
load('coreg/AdjMKFV')
AdjMKneurons = AdjMKneurons + AdjMKneurons';
AdjMKneurons(AdjMKneurons == 0) = NaN;
AdjMK = log10(AdjMKneurons);
AdjMK_bin = ~isnan(AdjMK);
AdjMK_bin = double(AdjMK_bin);

% Align macaque map to human map
MK_human_labels = AT.P.AtlROIs{atl_i}.LH.struct_names;
for i = 1:length(MK_human_labels)
    st = strsplit(MK_human_labels{i},'_M132');
    st = st{1};
    MK_human_labels{i} = st;
end
n_found = 0;
ind_hum2mac = zeros(1,length(AdjMK));
for i = 1:length(MK_labels)
    % Edits
    MK_labels{i} = replace(MK_labels{i},'/','_');
    if strcmp(MK_labels{i},'Core')
        MK_labels{i} = 'Aud-core';
    end
    if strcmp(MK_labels{i},'ENTO')
        MK_labels{i} = 'Entorhinal';
    end
    if strcmp(MK_labels{i},'INS')
        MK_labels{i} = 'Insula';
    end
    if strcmp(MK_labels{i},'PERI')
        MK_labels{i} = 'Perirhinal';
    end
    if strcmp(MK_labels{i},'PIR')
        MK_labels{i} = 'Piriform';
    end
    if strcmp(MK_labels{i},'POLE')
        MK_labels{i} = 'Temporal-pole';
    end
    if strcmp(MK_labels{i},'Pro.St.')
        MK_labels{i} = 'Prostriate';
    end
    if strcmp(MK_labels{i},'SII')
        MK_labels{i} = 'S2';
    end
    if strcmp(MK_labels{i},'SUB')
        MK_labels{i} = 'Subiculum';
    end
    if strcmp(MK_labels{i},'TEa_ma')
        MK_labels{i} = 'TEam-a';
    end
    if strcmp(MK_labels{i},'TEa_mp')
        MK_labels{i} = 'TEam-p';
    end
    
    la = find(strcmp(MK_human_labels,MK_labels{i}),1);
    if (~isempty(la))
        %AdjCS_mac(i,:) = AdjCS(la,:);
        ind_hum2mac(i) = la;
        n_found = n_found + 1;
    end
end

% test that indexing works
%strcmp(MK_human_labels(ind_hum2mac)' , MK_labels)
AdjCS_mac = AdjCS(ind_hum2mac,ind_hum2mac);
for i = 1:length(AdjCS_mac)
    %if (~isnan(AdjCS_mac(i,i)))
    AdjCS_mac(i,i) = NaN;
    %end
end
AdjMK_bin(isnan(AdjCS_mac)) = NaN;

% -------------- plot macaque ----------
h = figure;
fig_w = 10.5;
fig_h = 10.5;
set(h,'Position',[0 0 fig_w*100 fig_h*100])
set(h, 'PaperUnits', 'Inches')
set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
ax1 = axes('Parent',h,'Units','Normalized','Position',[0.2 0.24 0.55 0.55]);

AdjM = AdjMK_bin;

AT_c = AT;
AT_c.P.AtlROIs{atl_ind}.RH.struct_names = MK_labels;
ax2 = A.plotAdjElectrodeDistThresh(ax1, 1, AdjM, ...
    sprintf('Functional Interactivity %s, %s, %s',sid,...
    metric,at_name), inf(size(AdjM)), AdjM, 0.5, AT_c, atl_ind); %A.const.P_VALUE_CT
% colorbar
c = colorbar(ax1,'Location','manual');
n_cbarticks = 5;
cmarks = linspace((c.Limits(1)),(c.Limits(2)),n_cbarticks);
cmarklabs = cell(1,n_cbarticks);
for j = 1:n_cbarticks
    cmarklabs{j} = sprintf('%.2f',cmarks(j));
end

set(c,'Position',[0.9 0.24 0.02 0.553]);
set(c,'FontSize',8);
set(c,'TicksMode','manual');
set(c,'Ticks',cmarks);
set(c,'TickLabelsMode','manual');
set(c,'TickLabels',cmarklabs);
set(c,'TickLength',0);

% colorbar
c2 = colorbar(ax2,'Location','southoutside');
%cmarks = round(linspace(c2.Limits(1),c2.Limits(2),3));
set(c2,'Position',[0.2 0.08 0.55 0.014]);
set(c2,'FontSize',8);
set(c2,'TicksMode','manual');
set(c2,'Ticks',cmarks);
set(c2,'TickLabelsMode','manual');
set(c2,'TickLabels',cmarklabs);
set(c2,'TickLength',0);
A.saveFig(h,sprintf('./figures/adjacency_matrix_atlas/xsub_tract_%s_%s',...
  metric,at_name));
%close(h);

% -------------- plot human ----------
h = figure;
fig_w = 10.5;
fig_h = 10.5;
set(h,'Position',[0 0 fig_w*100 fig_h*100])
set(h, 'PaperUnits', 'Inches')
set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
ax1 = axes('Parent',h,'Units','Normalized','Position',[0.2 0.24 0.55 0.55]);

AdjM = AdjCS_mac;

AT_c = AT;
AT_c.P.AtlROIs{atl_ind}.RH.struct_names = MK_labels;
ax2 = A.plotAdjElectrodeDistThresh(ax1, 1, AdjM, ...
    sprintf('Functional Interactivity %s, %s, %s',sid,...
    metric,at_name), inf(size(AdjM)), AdjM, 0.5, AT_c, atl_ind); %A.const.P_VALUE_CT
% colorbar
c = colorbar(ax1,'Location','manual');
n_cbarticks = 5;
cmarks = linspace((c.Limits(1)),(c.Limits(2)),n_cbarticks);
cmarklabs = cell(1,n_cbarticks);
for j = 1:n_cbarticks
    cmarklabs{j} = sprintf('%.2f',cmarks(j));
end

set(c,'Position',[0.9 0.24 0.02 0.553]);
set(c,'FontSize',8);
set(c,'TicksMode','manual');
set(c,'Ticks',cmarks);
set(c,'TickLabelsMode','manual');
set(c,'TickLabels',cmarklabs);
set(c,'TickLength',0);

% colorbar
c2 = colorbar(ax2,'Location','southoutside');
%cmarks = round(linspace(c2.Limits(1),c2.Limits(2),3));
set(c2,'Position',[0.2 0.08 0.55 0.014]);
set(c2,'FontSize',8);
set(c2,'TicksMode','manual');
set(c2,'Ticks',cmarks);
set(c2,'TickLabelsMode','manual');
set(c2,'TickLabels',cmarklabs);
set(c2,'TickLength',0);
A.saveFig(h,sprintf('./figures/adjacency_matrix_atlas/xsub_hum_%s_%s',...
  metric,at_name));
%close(h);


% Correlation
%%
v_mac = AdjMK_bin(:);
v_mac = v_mac(~isnan(v_mac));
v_hum = AdjCS_mac(:);
v_hum = v_hum(~isnan(v_hum));

[r,pval] = corr(v_hum,v_mac);
fprintf('Metric: %s, M132 corr human to macaque: %.6f, p-val: %d\n',metric,r,pval)

% ROC
AdjCSM_mac = AdjCS_mag(ind_hum2mac,ind_hum2mac);
n_roc = 2^13;
Thetas = linspace(min(AdjCSM_mac(:)),max(AdjCSM_mac(:)),n_roc);
FP = zeros(1,n_roc);
TP = zeros(1,n_roc);
for i = 1:length(Thetas)
    theta = Thetas(i);
    
    % Thresholding
    AdjCS_bin = double(AdjCSM_mac > theta);
    AdjCS_bin(isnan(AdjMK_bin)) = NaN;
    h = AdjCS_bin(:);
    h = h(~isnan(h));
    m = AdjMK_bin(:);
    m = m(~isnan(m));
    
    tp = 0;
    cp = 0;
    fp = 0;
    cn = 0;
    for j = 1:length(h)
        if (m(j) == 1)
            cp = cp + 1;
        else
            cn = cn + 1;
        end
        if ((h(j) == 1) && (m(j) == 1))
            tp = tp + 1;
        elseif ((h(j) == 1) && (m(j) == 0))
            fp = fp + 1;
        end
    end
    FP(i) = fp/cn;
    TP(i) = tp/cp;
end
h = figure;
plot(FP,TP,'black-'); hold on;
plot([0 1],[0 1],'--','color',0.333*[1 1 1])
axis([0 1 0 1])
xlabel('False Positive Rate')
ylabel('True Positive Rate')
set(gca,'xtick',linspace(0,1,11),'Ticklength',[2e-3 1e-5],...
                'xticklabel',linspace(0,1,11),'TickDir','out');
set(gca,'ytick',linspace(0,1,11),'Ticklength',[2e-3 1e-5],...
                'yticklabel',linspace(0,1,11),'TickDir','out');
%title('ROC - Average metric threshold')
title(sprintf('ROC - Average metric threshold - AUC: %.2f',abs(trapz(FP,TP))))
fprintf('AUC: %.4f\n',abs(trapz(FP,TP)))
A.saveFig(h,sprintf('./figures/adjacency_matrix_atlas/xsub_ROC_%s_%s',...
  metric,at_name));


% --- Grab threshold at 50% "true" positive -------------------------------
[~,i50] = min(abs(TP - 0.5));
M_thresh_TP50 = Thetas(i50);


% -------------- plot human with new threshold ----------
h = figure;
fig_w = 10.5;
fig_h = 10.5;
set(h,'Position',[0 0 fig_w*100 fig_h*100])
set(h, 'PaperUnits', 'Inches')
set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
ax1 = axes('Parent',h,'Units','Normalized','Position',[0.2 0.24 0.55 0.55]);

AdjM = double(AdjCSM_mac > M_thresh_TP50);
AdjM(isnan(AdjMK_bin)) = NaN;

AT_c = AT;
AT_c.P.AtlROIs{atl_ind}.RH.struct_names = MK_labels;
ax2 = A.plotAdjElectrodeDistThresh(ax1, 1, AdjM, ...
    sprintf('%s, %s, TP50 thresh: %.6f',sid,...
    metric,M_thresh_TP50), inf(size(AdjM)), AdjM, 0.5, AT_c, atl_ind); %A.const.P_VALUE_CT
% colorbar
c = colorbar(ax1,'Location','manual');
n_cbarticks = 5;
cmarks = linspace((c.Limits(1)),(c.Limits(2)),n_cbarticks);
cmarklabs = cell(1,n_cbarticks);
for j = 1:n_cbarticks
    cmarklabs{j} = sprintf('%.2f',cmarks(j));
end

set(c,'Position',[0.9 0.24 0.02 0.553]);
set(c,'FontSize',8);
set(c,'TicksMode','manual');
set(c,'Ticks',cmarks);
set(c,'TickLabelsMode','manual');
set(c,'TickLabels',cmarklabs);
set(c,'TickLength',0);

% colorbar
c2 = colorbar(ax2,'Location','southoutside');
%cmarks = round(linspace(c2.Limits(1),c2.Limits(2),3));
set(c2,'Position',[0.2 0.08 0.55 0.014]);
set(c2,'FontSize',8);
set(c2,'TicksMode','manual');
set(c2,'Ticks',cmarks);
set(c2,'TickLabelsMode','manual');
set(c2,'TickLabels',cmarklabs);
set(c2,'TickLength',0);
A.saveFig(h,sprintf('./figures/adjacency_matrix_atlas/xsub_hum-TP50_%s_%s',...
  metric,at_name));
%close(h);


% Correlation after ROC thresholding
v_mac = AdjMK_bin(:);
v_mac = v_mac(~isnan(v_mac));
v_hum = AdjM(:);
v_hum = v_hum(~isnan(v_hum));

[r,pval] = corr(v_hum,v_mac);
fprintf('Metric: %s, M132 corr human TP50 to macaque: %.6f, p-val: %d\n',metric,r,pval)




% Threshold on fraction significant
AdjSig = zeros(size(adjct_xsub));
for m = 1:length(adjct_xsub)
    for n = 1:length(adjct_xsub)
        at = adjct_xsub{m,n};
        AdjSig(m,n) = sum(~isnan(at))/numel(at);
    end
end
AdjSig = AdjSig(ind_hum2mac,ind_hum2mac);
%AdjCSM_mac = AdjCS_mag(ind_hum2mac,ind_hum2mac);
n_roc = 2^13;
Thetas = linspace(min(AdjSig(:)),max(AdjSig(:)),n_roc);
FP = zeros(1,n_roc);
TP = zeros(1,n_roc);
for i = 1:length(Thetas)
    theta = Thetas(i);
    
    % Thresholding
    AdjCS_bin = double(AdjSig > theta);
    AdjCS_bin(isnan(AdjMK_bin)) = NaN;
    h = AdjCS_bin(:);
    h = h(~isnan(h));
    m = AdjMK_bin(:);
    m = m(~isnan(m));
    
    tp = 0;
    cp = 0;
    fp = 0;
    cn = 0;
    for j = 1:length(h)
        if (m(j) == 1)
            cp = cp + 1;
        else
            cn = cn + 1;
        end
        if ((h(j) == 1) && (m(j) == 1))
            tp = tp + 1;
        elseif ((h(j) == 1) && (m(j) == 0))
            fp = fp + 1;
        end
    end
    FP(i) = fp/cn;
    TP(i) = tp/cp;
end
FP = [1 FP];
TP = [1 TP];
h = figure;
plot(FP,TP,'black-'); hold on;
plot([0 1],[0 1],'--','color',0.333*[1 1 1])
axis([0 1 0 1])
xlabel('False Positive Rate')
ylabel('True Positive Rate')
set(gca,'xtick',linspace(0,1,11),'Ticklength',[2e-3 1e-5],...
                'xticklabel',linspace(0,1,11),'TickDir','out');
set(gca,'ytick',linspace(0,1,11),'Ticklength',[2e-3 1e-5],...
                'yticklabel',linspace(0,1,11),'TickDir','out');
title(sprintf('ROC - Fraction significant threshold - AUC: %.2f',abs(trapz(FP,TP))))
fprintf('AUC: %.4f\n',abs(trapz(FP,TP)))
A.saveFig(h,sprintf('./figures/adjacency_matrix_atlas/xsub_ROC-sig_%s_%s',...
  metric,at_name));

% --- mSu
% dmat = A.getDistanceMatrix(52);
% CT = A.getCT(52,metric,false);
%AT = A.getAT(52,CT,dmat);
