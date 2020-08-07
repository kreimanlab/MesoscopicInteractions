close all;
clear;

if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        %resultsDir = '/media/klab/44/data/results';
        resultsDir = '/media/klab/D0BA8713BA86F56E/data/results';
        h5Dir = '/media/klab/44/h5';
    else
        resultsDir = '/mnt/cuenap2/data/results';
        h5Dir = '/mnt/cuenap2/scripts/synth/out';
    end
end

if (~exist('A','var'))
    A = Analysis(resultsDir,h5Dir);
    % Correct for multiple comparisons - 60 minutes
    A.const.P_VALUE = A.const.P_VALUE/60;
end

%Metrics = {'p','s','pP','sP','sd','st','sa','sb','sg',...
%    'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma',...
%    'scBroadband','scDelta','scTheta','scAlpha','scBeta','scGamma'};
%Metrics = {'s'};
Metrics = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma','s','sP'};

NotEnvMetrics = {'p','s','pP','sP'};

for metric_i = 1:length(Metrics)
metric = Metrics{metric_i};


CT_all = cell(1,A.h5eeg.n_sub);
for i = 1:A.h5eeg.n_sub
    close all;
    
    try
        
        
    sid = A.h5eeg.subjects{i};
    dmat = A.getDistanceMatrix(i);
    %metric = Metrics{2};
    
    % Load CT
    CT = A.getCT(i,metric,false);
    CTnull = A.getCT(i,CT.metric,true);
    
    
    % Save for xsub
    CT_all{i} = CT;
    AdjCT_thresh = binoinv((1-(A.const.P_VALUE_CT/nchoosek(A.h5eeg.n_bchan{i},2))),CT.n_w,(A.const.P_VALUE))/CT.n_w;
    %AdjCT_thresh = binoinv((1-A.const.P_VALUE_CT),CT.n_w,A.const.P_VALUE)/CT.n_w;

    % --- Plot AdjCT ------------------------------------------------------
    h = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    set(h,'Position',[0 0 fig_w*100 fig_h*100])
    set(h, 'PaperUnits', 'Inches')
    set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    ax1 = axes('Parent',h,'Units','Normalized','Position',[0.18 0.35 0.55 0.4]);
    %ax2 = A.plotAdjElectrode(ax1, i, dmat, sprintf('Distance Matrix: %s',sid));
    ax2 = A.plotAdjElectrodeDistThresh(ax1, i, CT.AdjCT, ...
        sprintf('Consistency across time %s, %s',sid, metric), dmat, ...
        CT.AdjCT, AdjCT_thresh); %A.const.P_VALUE_CT
    
    % colorbar
    c = colorbar(ax1,'Location','manual');
    n_cbarticks = 5;
    cmarks = linspace((c.Limits(1)),(c.Limits(2)),n_cbarticks);
    cmarklabs = cell(1,n_cbarticks);
    for j = 1:n_cbarticks
        cmarklabs{j} = sprintf('%.2f',cmarks(j));
    end
    
    set(c,'Position',[0.86 0.35 0.02 0.4]);
    set(c,'FontSize',8);
    set(c,'TicksMode','manual');
    set(c,'Ticks',cmarks);
    set(c,'TickLabelsMode','manual');
    set(c,'TickLabels',cmarklabs);
    set(c,'TickLength',0);
    
    % colorbar
    c2 = colorbar(ax2,'Location','southoutside');
    %cmarks = round(linspace(c2.Limits(1),c2.Limits(2),3));
    set(c2,'Position',[0.18 0.24 0.55 0.014]);
    set(c2,'FontSize',8);
    set(c2,'TicksMode','manual');
    set(c2,'Ticks',cmarks);
    set(c2,'TickLabelsMode','manual');
    set(c2,'TickLabels',cmarklabs);
    set(c2,'TickLength',0);
    
    A.saveFig(h,sprintf('./figures/adjacency_matrix/%s_%s_AdjCT',...
        sid,metric));
    
    close(h);
    
    % --- Plot AdjCT Negative ---------------------------------------------
    h = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    set(h,'Position',[0 0 fig_w*100 fig_h*100])
    set(h, 'PaperUnits', 'Inches')
    set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    ax1 = axes('Parent',h,'Units','Normalized','Position',[0.18 0.35 0.55 0.4]);
    %ax2 = A.plotAdjElectrode(ax1, i, dmat, sprintf('Distance Matrix: %s',sid));
    ax2 = A.plotAdjElectrodeDistThresh(ax1, i, CTnull.AdjCT, ...
        sprintf('Consistency across time %s, %s',sid, metric), dmat, ...
        CTnull.AdjCT, AdjCT_thresh); %A.const.P_VALUE_CT
    
    % colorbar
    c = colorbar(ax1,'Location','manual');
    n_cbarticks = 5;
    cmarks = linspace((c.Limits(1)),(c.Limits(2)),n_cbarticks);
    cmarklabs = cell(1,n_cbarticks);
    for j = 1:n_cbarticks
        cmarklabs{j} = sprintf('%.2f',cmarks(j));
    end
    
    set(c,'Position',[0.86 0.35 0.02 0.4]);
    set(c,'FontSize',8);
    set(c,'TicksMode','manual');
    set(c,'Ticks',cmarks);
    set(c,'TickLabelsMode','manual');
    set(c,'TickLabels',cmarklabs);
    set(c,'TickLength',0);
    
    % colorbar
    c2 = colorbar(ax2,'Location','southoutside');
    %cmarks = round(linspace(c2.Limits(1),c2.Limits(2),3));
    set(c2,'Position',[0.18 0.24 0.55 0.014]);
    set(c2,'FontSize',8);
    set(c2,'TicksMode','manual');
    set(c2,'Ticks',cmarks);
    set(c2,'TickLabelsMode','manual');
    set(c2,'TickLabels',cmarklabs);
    set(c2,'TickLength',0);
    
    A.saveFig(h,sprintf('./figures/adjacency_matrix/%s_%s_AdjCTneg',...
        sid,metric));
    
    close(h);
    
    
    % --- Plot AdjMag -----------------------------------------------------
    h = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    set(h,'Position',[0 0 fig_w*100 fig_h*100])
    set(h, 'PaperUnits', 'Inches')
    set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    ax1 = axes('Parent',h,'Units','Normalized','Position',[0.18 0.35 0.55 0.4]);
    %ax2 = A.plotAdjElectrode(ax1, i, dmat, sprintf('Distance Matrix: %s',sid));
    ax2 = A.plotAdjElectrodeDistThresh(ax1, i, CT.AdjMag, ...
        sprintf('Average magnitude %s, %s',sid, metric), dmat, ...
        CT.AdjCT, AdjCT_thresh); %A.const.P_VALUE_CT
    
    % colorbar
    c = colorbar(ax1,'Location','manual');
    n_cbarticks = 5;
    cmarks = linspace((c.Limits(1)),(c.Limits(2)),n_cbarticks);
    cmarklabs = cell(1,n_cbarticks);
    for j = 1:n_cbarticks
        cmarklabs{j} = sprintf('%.2f',cmarks(j));
    end
    
    set(c,'Position',[0.86 0.35 0.02 0.4]);
    set(c,'FontSize',8);
    set(c,'TicksMode','manual');
    set(c,'Ticks',cmarks);
    set(c,'TickLabelsMode','manual');
    set(c,'TickLabels',cmarklabs);
    set(c,'TickLength',0);
    
    % colorbar
    c2 = colorbar(ax2,'Location','southoutside');
    %cmarks = round(linspace(c2.Limits(1),c2.Limits(2),3));
    set(c2,'Position',[0.18 0.24 0.55 0.014]);
    set(c2,'FontSize',8);
    set(c2,'TicksMode','manual');
    set(c2,'Ticks',cmarks);
    set(c2,'TickLabelsMode','manual');
    set(c2,'TickLabels',cmarklabs);
    set(c2,'TickLength',0);
    
    A.saveFig(h,sprintf('./figures/adjacency_matrix/%s_%s_AdjMag',...
        sid,metric));
    
    close(h);
    
    % --- Plot AdjMag Negative ---------------------------------------------
    h = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    set(h,'Position',[0 0 fig_w*100 fig_h*100])
    set(h, 'PaperUnits', 'Inches')
    set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    ax1 = axes('Parent',h,'Units','Normalized','Position',[0.18 0.35 0.55 0.4]);
    %ax2 = A.plotAdjElectrode(ax1, i, dmat, sprintf('Distance Matrix: %s',sid));
    ax2 = A.plotAdjElectrodeDistThresh(ax1, i, CTnull.AdjMag, ...
        sprintf('Consistency across time %s, %s',sid, metric), dmat, ...
        CTnull.AdjCT, AdjCT_thresh); %A.const.P_VALUE_CT
    
    % colorbar
    c = colorbar(ax1,'Location','manual');
    n_cbarticks = 5;
    cmarks = linspace((c.Limits(1)),(c.Limits(2)),n_cbarticks);
    cmarklabs = cell(1,n_cbarticks);
    for j = 1:n_cbarticks
        cmarklabs{j} = sprintf('%.2f',cmarks(j));
    end
    
    set(c,'Position',[0.86 0.35 0.02 0.4]);
    set(c,'FontSize',8);
    set(c,'TicksMode','manual');
    set(c,'Ticks',cmarks);
    set(c,'TickLabelsMode','manual');
    set(c,'TickLabels',cmarklabs);
    set(c,'TickLength',0);
    
    % colorbar
    c2 = colorbar(ax2,'Location','southoutside');
    %cmarks = round(linspace(c2.Limits(1),c2.Limits(2),3));
    set(c2,'Position',[0.18 0.24 0.55 0.014]);
    set(c2,'FontSize',8);
    set(c2,'TicksMode','manual');
    set(c2,'Ticks',cmarks);
    set(c2,'TickLabelsMode','manual');
    set(c2,'TickLabels',cmarklabs);
    set(c2,'TickLength',0);
    
    A.saveFig(h,sprintf('./figures/adjacency_matrix/%s_%s_AdjMagneg',...
        sid,metric));
    
    close(h);
    
    % --- Plot Significance across time -----------------------------------
    % Pick electrode pair to plot
    AdjCT = CT.AdjCT;
    AdjCT(AdjCT <= AdjCT_thresh) = NaN;
    AdjCT(dmat <= A.const.DIST_THRESH_MM) = NaN;
    AdjCT = abs(AdjCT);
    la = length(AdjCT);
    [~,mini] = min(AdjCT(:));
    min_i = mod(mini-1,la)+1;
    min_j = ceil(mini/la);
    [~,maxi] = max(AdjCT(:));
    max_i = mod(maxi-1,la)+1;
    max_j = ceil(maxi/la);
    med = nanmedian(AdjCT(:));
    medi = find(AdjCT(:)==med,1);
    if (isempty(medi))
        fprintf(2,'W: Did not find median index for plot significance across time.\n');
        medi = la*(round(la/2)) + 1;
    end
    med_i = mod(medi-1,la)+1;
    med_j = ceil(medi/la);
    
    %Pick non-significant electrode pair
    AdjCT = CT.AdjCT;
    AdjCT(AdjCT > AdjCT_thresh) = NaN;
    AdjCT = abs(AdjCT);
    [~,mini] = min(AdjCT(:));
    min_i_n = mod(mini-1,la)+1;
    min_j_n = ceil(mini/la);
    
    h2 = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    set(h2,'Position',[0 0 fig_w*100 fig_h*100])
    set(h2, 'PaperUnits', 'Inches')
    set(h2, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h2, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    
    bchan1 = max_j;%30;
    bchan2 = max_i;%32;
    ax_1 = axes('Parent',h2,'Units','Normalized','Position',[0.1 0.85 0.8 0.12]);
    A.plotSig(i, ax_1, CT, bchan1, bchan2, dmat);
    
    bchan1 = med_j;%26;
    bchan2 = med_i;%59;
    ax_2 = axes('Parent',h2,'Units','Normalized','Position',[0.1 0.65 0.8 0.12]);
    A.plotSig(i, ax_2, CT, bchan1, bchan2, dmat);
    
    bchan1 = min_j;%10;
    bchan2 = min_i;%61;
    ax_3 = axes('Parent',h2,'Units','Normalized','Position',[0.1 0.45 0.8 0.12]);
    A.plotSig(i, ax_3, CT, bchan1, bchan2, dmat);
    ax_3n = axes('Parent',h2,'Units','Normalized','Position',[0.1 0.25 0.8 0.12]);
    A.plotSig(i, ax_3n, CTnull, bchan1, bchan2, dmat);
    
    bchan1 = min_j_n;
    bchan2 = min_i_n;
    ax_3 = axes('Parent',h2,'Units','Normalized','Position',[0.1 0.05 0.8 0.12]);
    A.plotSig(i, ax_3, CT, bchan1, bchan2, dmat);
    
    A.saveFig(h2,sprintf('./figures/adjacency_matrix/%s_%s_Sig',sid,metric));
    close(h2);
    
    
    % --- Plot timeseries -------------------------------------------------
    
    %skip_env = sum(contains(NotEnvMetrics,metric)) > 0;
    h3 = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    set(h3,'Position',[0 0 fig_w*100 fig_h*100])
    set(h3, 'PaperUnits', 'Inches')
    set(h3, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h3, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    
    % --- maximum ---
    % Get index of the average segment
    chanBip = max_j;
    chan2Bip = max_i;
    [~,a2] = min(abs(abs(CT.Adj{chanBip,chan2Bip}(CT.SigIdx{chanBip,chan2Bip}))...
        - CT.AdjMag(chanBip,chan2Bip)));
    samp = (a2 - 1)*round(CT.w*A.h5eeg.fs{i}) + 1;
    if (~isempty(samp))
    atlas = '';
    bchan1 = A.h5eeg.bip{i}(chanBip,1);
    bchan2 = A.h5eeg.bip{i}(chanBip,2);
    vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan1 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
    vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan2 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
    v = vb1 - vb2;
    b2chan1 = A.h5eeg.bip{i}(chan2Bip,1);
    b2chan2 = A.h5eeg.bip{i}(chan2Bip,2);
    vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan1 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
    vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan2 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
    v2 = vb1 - vb2;
    
    ax1 = axes('Parent',h3,'Units','Normalized','Position',[0.1 0.85 0.65 0.09]);
    A.plotLine(ax1,v,i); 
    ax2 = axes('Parent',h3,'Units','Normalized','Position',[0.1 0.75 0.65 0.09]);
    A.plotLine(ax2,v2,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h3); close(hb);
    set(h3.Children(1),'Position',[0.8 0.85 0.12 0.09]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h3); close(hb);
    set(h3.Children(1),'Position',[0.8 0.75 0.12 0.09]);
    end
    
    % --- median ---
    % Get index of the average segment
    chanBip = med_j;
    chan2Bip = med_i;
    [~,a2] = min(abs(abs(CT.Adj{chanBip,chan2Bip}(CT.SigIdx{chanBip,chan2Bip}))...
        - CT.AdjMag(chanBip,chan2Bip)));
    samp = (a2 - 1)*round(CT.w*A.h5eeg.fs{i}) + 1;
    if (~isempty(samp))
        
    %atlas = '';
    bchan1 = A.h5eeg.bip{i}(chanBip,1);
    bchan2 = A.h5eeg.bip{i}(chanBip,2);
    vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan1 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
    vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan2 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
    v = vb1 - vb2;
    b2chan1 = A.h5eeg.bip{i}(chan2Bip,1);
    b2chan2 = A.h5eeg.bip{i}(chan2Bip,2);
    vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan1 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
    vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan2 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
    v2 = vb1 - vb2;
    
    ax1 = axes('Parent',h3,'Units','Normalized','Position',[0.1 0.65 0.65 0.09]);
    A.plotLine(ax1,v,i); 
    ax2 = axes('Parent',h3,'Units','Normalized','Position',[0.1 0.55 0.65 0.09]);
    A.plotLine(ax2,v2,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h3); close(hb);
    set(h3.Children(1),'Position',[0.8 0.65 0.12 0.09]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h3); close(hb);
    set(h3.Children(1),'Position',[0.8 0.55 0.12 0.09]);
    end
    
    % --- minimum ---
    % Get index of the average segment
    chanBip = min_j;
    chan2Bip = min_i;
    [~,a2] = min(abs(abs(CT.Adj{chanBip,chan2Bip}(CT.SigIdx{chanBip,chan2Bip}))...
        - CT.AdjMag(chanBip,chan2Bip)));
    samp = (a2 - 1)*round(CT.w*A.h5eeg.fs{i}) + 1;
    if (~isempty(samp))
    %atlas = '';
    bchan1 = A.h5eeg.bip{i}(chanBip,1);
    bchan2 = A.h5eeg.bip{i}(chanBip,2);
    vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan1 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
    vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan2 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
    v = vb1 - vb2;
    b2chan1 = A.h5eeg.bip{i}(chan2Bip,1);
    b2chan2 = A.h5eeg.bip{i}(chan2Bip,2);
    vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan1 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
    vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan2 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
    v2 = vb1 - vb2;
    
    ax1 = axes('Parent',h3,'Units','Normalized','Position',[0.1 0.45 0.65 0.09]);
    A.plotLine(ax1,v,i); 
    ax2 = axes('Parent',h3,'Units','Normalized','Position',[0.1 0.35 0.65 0.09]);
    A.plotLine(ax2,v2,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h3); close(hb);
    set(h3.Children(1),'Position',[0.8 0.45 0.12 0.09]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h3); close(hb);
    set(h3.Children(1),'Position',[0.8 0.35 0.12 0.09]);
    end
    
    % --- non significant minimum ---
    % Get index of the average segment
    chanBip = min_j_n;
    chan2Bip = min_i_n;
    [~,a2] = min(abs(abs(CT.Adj{chanBip,chan2Bip}(:))...
        - CT.AdjMag(chanBip,chan2Bip))); % CT.SigIdx{chanBip,chan2Bip}
    samp = (a2 - 1)*round(CT.w*A.h5eeg.fs{i}) + 1;
    if (~isempty(samp))
    
        %atlas = '';
        bchan1 = A.h5eeg.bip{i}(chanBip,1);
        bchan2 = A.h5eeg.bip{i}(chanBip,2);
        vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan1 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
        vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan2 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
        v = vb1 - vb2;
        b2chan1 = A.h5eeg.bip{i}(chan2Bip,1);
        b2chan2 = A.h5eeg.bip{i}(chan2Bip,2);
        vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan1 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
        vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan2 samp],[1 round(CT.w*A.h5eeg.fs{i})]);
        v2 = vb1 - vb2;

        ax1 = axes('Parent',h3,'Units','Normalized','Position',[0.1 0.15 0.65 0.09]);
        A.plotLine(ax1,v,i); 
        ax2 = axes('Parent',h3,'Units','Normalized','Position',[0.1 0.05 0.65 0.09]);
        A.plotLine(ax2,v2,i);
        hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h3); close(hb);
        set(h3.Children(1),'Position',[0.8 0.15 0.12 0.09]);
        hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h3); close(hb);
        set(h3.Children(1),'Position',[0.8 0.05 0.12 0.09]);
    end
    
    
    A.saveFig(h3,sprintf('./figures/adjacency_matrix/%s_%s_times',sid,metric));
    
    %return
    close(h3)
    
    % --- Print artifacts ---
    
    h4 = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    set(h4,'Position',[0 0 fig_w*100 fig_h*100])
    set(h4, 'PaperUnits', 'Inches')
    set(h4, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h4, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    
    
    %return
    close(h4)
    
    catch e
       fprintf('Skipped.\n')
       %rethrow(e)
    end
    
end
end