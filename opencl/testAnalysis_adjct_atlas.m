close all;
%clear;

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

Metrics = {'p','s','pP','sP','sd','st','sa','sb','sg',...
    'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma',...
    'scBroadband','scDelta','scTheta','scAlpha','scBeta','scGamma'};
%Metrics = {'s'};

NotEnvMetrics = {'p','s','pP','sP'};

for metric_i = 1:length(Metrics)
metric = Metrics{metric_i};


CT_all = cell(1,A.h5eeg.n_sub);
for i = 1:A.h5eeg.n_sub
    close all;
    
    %try
    
    sid = A.h5eeg.subjects{i};
    dmat = A.getDistanceMatrix(i);
    
    % Load CT
    CT = A.getCT(i,metric,false);
    %CTnull = A.getCT(i,CT.metric,true);
    
    % Save for xsub
    CT_all{i} = CT;
    
    % Map to atlas
    AT = A.getAT(i,CT,dmat);
    
    % CT statistical threshold
    AdjCT_thresh = binoinv((1-A.const.P_VALUE_CT),CT.n_w,A.const.P_VALUE)/CT.n_w;
    
    % average AdjCT then threshold
    AdjMa = cell(1,AT.n_atlas);
    for a = 1:AT.n_atlas
        AdjM = zeros(AT.n_rois(a),AT.n_rois(a));
        for l = 1:(AT.n_rois(a)-1)
            for k = (l+1):AT.n_rois(a)
                AdjM(l,k) = mean(AT.AdjCT{a}{l,k});
                AdjM(k,l) = mean(AT.AdjCT{a}{k,l});
            end
        end
        AdjMa{a} = AdjM;
    end

    %%
    % --- Plot AdjCT ------------------------------------------------------
    for atl_ind = 1:AT.n_atlas
    %atl_ind = 2;
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
    AdjM = AdjMa{atl_ind};
    ax2 = A.plotAdjElectrodeDistThresh(ax1, i, AdjM, ...
        sprintf('Consistency across time (avg then threshold) %s, %s, %s',sid,...
        metric,at_name), inf(size(AdjM)), AdjM, AdjCT_thresh, AT, atl_ind); %A.const.P_VALUE_CT
    
    % colorbar
    c = colorbar(ax1,'Location','manual');
    n_cbarticks = 5;
    cmarks = linspace((c.Limits(1)),(c.Limits(2)),n_cbarticks);
    cmarklabs = cell(1,n_cbarticks);
    for j = 1:n_cbarticks
        cmarklabs{j} = sprintf('%.2f',cmarks(j));
    end
    
    set(c,'Position',[0.93 0.24 0.02 0.553]);
    set(c,'FontSize',8);
    set(c,'TicksMode','manual');
    set(c,'Ticks',cmarks);
    set(c,'TickLabelsMode','manual');
    set(c,'TickLabels',cmarklabs);
    set(c,'TickLength',0);
    
    % colorbar
    c2 = colorbar(ax2,'Location','southoutside');
    %cmarks = round(linspace(c2.Limits(1),c2.Limits(2),3));
    set(c2,'Position',[0.2 0.04 0.55 0.014]);
    set(c2,'FontSize',8);
    set(c2,'TicksMode','manual');
    set(c2,'Ticks',cmarks);
    set(c2,'TickLabelsMode','manual');
    set(c2,'TickLabels',cmarklabs);
    set(c2,'TickLength',0);
    
    A.saveFig(h,sprintf('./figures/adjacency_matrix_atlas/%s_%s_%s_AdjCT',...
       sid,metric,at_name));
    %return
    
    close(h);
    end

    
%     catch
%        fprintf('Skipped.\n')
%     end
    
end
end