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

for i = 1:A.h5eeg.n_sub
    sid = A.h5eeg.subjects{i};
    dmat = A.getDistanceMatrix(i);

    h = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    set(h,'Position',[0 0 fig_w*100 fig_h*100])
    set(h, 'PaperUnits', 'Inches')
    set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    ax1 = axes('Parent',h,'Units','Normalized','Position',[0.18 0.35 0.55 0.4]);
    %ax2 = A.plotAdjElectrode(ax1, i, dmat, sprintf('Distance Matrix: %s',sid));
    ax2 = A.plotAdjElectrodeDistThresh(ax1, i, dmat, sprintf('Distance Matrix: %s',sid), dmat);
    
    % colorbar
    c = colorbar(ax1,'Location','manual');
    %cmarks = round(linspace(ceil(c.Limits(1)),floor(c.Limits(2)),2));
    %cmarks = [ceil(c.Limits(1)),15,floor(c.Limits(2))];
    n_cbarticks = 5;
    cmarks = linspace(c.Limits(1),c.Limits(2),n_cbarticks);
    cmarklabs = cell(1,n_cbarticks);
    for j = 1:n_cbarticks
        cmarklabs{j} = sprintf('%.0f',cmarks(j));
    end
    %cmarks = c.Limits(1):50:(c.Limits(2));
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
    
    A.saveFig(h,sprintf('./figures/distance_matrix_dthresh/%s_nchan-%i',...
        sid,A.h5eeg.n_chan{i}));
    
    %return
    close(h);
    %break
end