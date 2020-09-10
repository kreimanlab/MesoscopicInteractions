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
    %h = figure;
    sid = A.h5eeg.subjects{i};
    dmat = A.getDistanceMatrix(i);
    
    w = 5;
    
    % Channel selection
    chanBip = randi([1 A.h5eeg.n_bchan{i}]);
    chan2Bip = randi([1 A.h5eeg.n_bchan{i}]);

    % Sample selection
    samp = randi([1 (A.h5eeg.n_samples{i}-round(w*A.h5eeg.fs{i}))]);
    
    % Read
    bchan1 = A.h5eeg.bip{i}(chanBip,1);
    bchan2 = A.h5eeg.bip{i}(chanBip,2);
    vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan1 samp],[1 round(w*A.h5eeg.fs{i})]);
    vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan2 samp],[1 round(w*A.h5eeg.fs{i})]);
    v = vb1 - vb2;

    b2chan1 = A.h5eeg.bip{i}(chan2Bip,1);
    b2chan2 = A.h5eeg.bip{i}(chan2Bip,2);
    vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan1 samp],[1 round(w*A.h5eeg.fs{i})]);
    vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan2 samp],[1 round(w*A.h5eeg.fs{i})]);
    v2 = vb1 - vb2;

    fprintf('\n%s\n',A.h5eeg.subjects{i})
    fprintf('bchan1: %i\n',bchan1)
    fprintf('bchan2: %i\n',bchan2)
    fprintf('dist bchan1,bchan2: %.4f\n',A.h5eeg.bip{i}(chanBip,3))
    fprintf('b2chan1: %i\n',b2chan1)
    fprintf('b2chan2: %i\n',b2chan2)
    fprintf('dist b2chan1,b2chan2: %.4f\n',A.h5eeg.bip{i}(chan2Bip,3))
    
    % Plot Delta
    h = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    fig_n_w0 = 0.05; % width starting position
    fig_n_w = 0.4; % width
    fig_n_ws = 0.001; % space between brain
    fig_n_wbrain = 0.2; % width brain
    fig_n_h0 = 0.048; % height per axis
    atlas = 'd';
    set(h,'Position',[0 0 fig_w*100 fig_h*100])
    set(h, 'PaperUnits', 'Inches')
    set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    % --- raw pair ---
    ax1 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.90 fig_n_w fig_n_h0]);
    A.plotLine(ax1,v,i);    
    ax1b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.85 fig_n_w fig_n_h0]);
    A.plotLine(ax1b,v2,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.90 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.85 fig_n_wbrain fig_n_h0]);
    % --- filt pair ---
    vf = A.del(v,i);
    v2f = A.del(v2,i);
    ax2 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.80 fig_n_w fig_n_h0]);
    A.plotLine(ax2,vf,i);
    ax2b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.75 fig_n_w fig_n_h0]);
    A.plotLine(ax2b,v2f,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.80 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.75 fig_n_wbrain fig_n_h0]);
    % --- env pair ---
    vfe = A.env(vf);
    v2fe = A.env(v2f);
    ax3 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.70 fig_n_w fig_n_h0]);
    A.plotLine(ax3,vfe,i);
    ax3b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.65 fig_n_w fig_n_h0]);
    A.plotLine(ax3b,v2fe,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.70 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.65 fig_n_wbrain fig_n_h0]);
    r_s = corr(v', v2', 'Type','Spearman');
    r_senv = corr(vfe', v2fe', 'Type','Spearman');
    A.saveFig(h,sprintf('./figures/voltages/filtering_%s_chans%i-%i_samp-%i_w-%i_rs_%i_rsenv_%i_del',...
        sid,bchan1,b2chan1,samp,w,round(1000*r_s),round(1000*r_senv)))
    close(h);

    % Plot Theta
    h = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    fig_n_w0 = 0.05; % width starting position
    fig_n_w = 0.4; % width
    fig_n_ws = 0.001; % space between brain
    fig_n_wbrain = 0.2; % width brain
    fig_n_h0 = 0.048; % height per axis
    %atlas = 'd';
    set(h,'Position',[0 0 fig_w*100 fig_h*100])
    set(h, 'PaperUnits', 'Inches')
    set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    % --- raw pair ---
    ax1 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.90 fig_n_w fig_n_h0]);
    A.plotLine(ax1,v,i);    
    ax1b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.85 fig_n_w fig_n_h0]);
    A.plotLine(ax1b,v2,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.90 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.85 fig_n_wbrain fig_n_h0]);
    % --- filt pair ---
    vf = A.the(v,i);
    v2f = A.the(v2,i);
    ax2 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.80 fig_n_w fig_n_h0]);
    A.plotLine(ax2,vf,i);
    ax2b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.75 fig_n_w fig_n_h0]);
    A.plotLine(ax2b,v2f,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.80 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.75 fig_n_wbrain fig_n_h0]);
    % --- env pair ---
    vfe = A.env(vf);
    v2fe = A.env(v2f);
    ax3 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.70 fig_n_w fig_n_h0]);
    A.plotLine(ax3,vfe,i);
    ax3b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.65 fig_n_w fig_n_h0]);
    A.plotLine(ax3b,v2fe,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.70 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.65 fig_n_wbrain fig_n_h0]);
    r_s = corr(v', v2', 'Type','Spearman');
    r_senv = corr(vfe', v2fe', 'Type','Spearman');
    A.saveFig(h,sprintf('./figures/voltages/filtering_%s_chans%i-%i_samp-%i_w-%i_rs_%i_rsenv_%i_the',...
        sid,bchan1,b2chan1,samp,w,round(1000*r_s),round(1000*r_senv)))
    close(h);
    
    % Plot Alpha
    h = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    fig_n_w0 = 0.05; % width starting position
    fig_n_w = 0.4; % width
    fig_n_ws = 0.001; % space between brain
    fig_n_wbrain = 0.2; % width brain
    fig_n_h0 = 0.048; % height per axis
    %atlas = 'd';
    set(h,'Position',[0 0 fig_w*100 fig_h*100])
    set(h, 'PaperUnits', 'Inches')
    set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    % --- raw pair ---
    ax1 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.90 fig_n_w fig_n_h0]);
    A.plotLine(ax1,v,i);    
    ax1b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.85 fig_n_w fig_n_h0]);
    A.plotLine(ax1b,v2,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.90 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.85 fig_n_wbrain fig_n_h0]);
    % --- filt pair ---
    vf = A.alp(v,i);
    v2f = A.alp(v2,i);
    ax2 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.80 fig_n_w fig_n_h0]);
    A.plotLine(ax2,vf,i);
    ax2b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.75 fig_n_w fig_n_h0]);
    A.plotLine(ax2b,v2f,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.80 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.75 fig_n_wbrain fig_n_h0]);
    % --- env pair ---
    vfe = A.env(vf);
    v2fe = A.env(v2f);
    ax3 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.70 fig_n_w fig_n_h0]);
    A.plotLine(ax3,vfe,i);
    ax3b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.65 fig_n_w fig_n_h0]);
    A.plotLine(ax3b,v2fe,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.70 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.65 fig_n_wbrain fig_n_h0]);
    r_s = corr(v', v2', 'Type','Spearman');
    r_senv = corr(vfe', v2fe', 'Type','Spearman');
    A.saveFig(h,sprintf('./figures/voltages/filtering_%s_chans%i-%i_samp-%i_w-%i_rs_%i_rsenv_%i_alp',...
        sid,bchan1,b2chan1,samp,w,round(1000*r_s),round(1000*r_senv)))
    close(h);
    
    % Plot Beta
    h = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    fig_n_w0 = 0.05; % width starting position
    fig_n_w = 0.4; % width
    fig_n_ws = 0.001; % space between brain
    fig_n_wbrain = 0.2; % width brain
    fig_n_h0 = 0.048; % height per axis
    %atlas = 'd';
    set(h,'Position',[0 0 fig_w*100 fig_h*100])
    set(h, 'PaperUnits', 'Inches')
    set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    % --- raw pair ---
    ax1 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.90 fig_n_w fig_n_h0]);
    A.plotLine(ax1,v,i);    
    ax1b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.85 fig_n_w fig_n_h0]);
    A.plotLine(ax1b,v2,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.90 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.85 fig_n_wbrain fig_n_h0]);
    % --- filt pair ---
    vf = A.bet(v,i);
    v2f = A.bet(v2,i);
    ax2 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.80 fig_n_w fig_n_h0]);
    A.plotLine(ax2,vf,i);
    ax2b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.75 fig_n_w fig_n_h0]);
    A.plotLine(ax2b,v2f,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.80 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.75 fig_n_wbrain fig_n_h0]);
    % --- env pair ---
    vfe = A.env(vf);
    v2fe = A.env(v2f);
    ax3 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.70 fig_n_w fig_n_h0]);
    A.plotLine(ax3,vfe,i);
    ax3b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.65 fig_n_w fig_n_h0]);
    A.plotLine(ax3b,v2fe,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.70 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.65 fig_n_wbrain fig_n_h0]);
    r_s = corr(v', v2', 'Type','Spearman');
    r_senv = corr(vfe', v2fe', 'Type','Spearman');
    A.saveFig(h,sprintf('./figures/voltages/filtering_%s_chans%i-%i_samp-%i_w-%i_rs_%i_rsenv_%i_bet',...
        sid,bchan1,b2chan1,samp,w,round(1000*r_s),round(1000*r_senv)))
    close(h);
    
    % Plot Gamma
    h = figure;
    fig_w = 8.5;
    fig_h = 11.0;
    fig_n_w0 = 0.05; % width starting position
    fig_n_w = 0.4; % width
    fig_n_ws = 0.001; % space between brain
    fig_n_wbrain = 0.2; % width brain
    fig_n_h0 = 0.048; % height per axis
    %atlas = 'd';
    set(h,'Position',[0 0 fig_w*100 fig_h*100])
    set(h, 'PaperUnits', 'Inches')
    set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    % --- raw pair ---
    ax1 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.90 fig_n_w fig_n_h0]);
    A.plotLine(ax1,v,i);    
    ax1b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.85 fig_n_w fig_n_h0]);
    A.plotLine(ax1b,v2,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.90 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.85 fig_n_wbrain fig_n_h0]);
    % --- filt pair ---
    vf = A.gam(v,i);
    v2f = A.gam(v2,i);
    ax2 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.80 fig_n_w fig_n_h0]);
    A.plotLine(ax2,vf,i);
    ax2b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.75 fig_n_w fig_n_h0]);
    A.plotLine(ax2b,v2f,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.80 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.75 fig_n_wbrain fig_n_h0]);
    % --- env pair ---
    vfe = A.env(vf);
    v2fe = A.env(v2f);
    ax3 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.70 fig_n_w fig_n_h0]);
    A.plotLine(ax3,vfe,i);
    ax3b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.65 fig_n_w fig_n_h0]);
    A.plotLine(ax3b,v2fe,i);
    hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.70 fig_n_wbrain fig_n_h0]);
    hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.65 fig_n_wbrain fig_n_h0]);
    r_s = corr(v', v2', 'Type','Spearman');
    r_senv = corr(vfe', v2fe', 'Type','Spearman');
    A.saveFig(h,sprintf('./figures/voltages/filtering_%s_chans%i-%i_samp-%i_w-%i_rs_%i_rsenv_%i_gam',...
        sid,bchan1,b2chan1,samp,w,round(1000*r_s),round(1000*r_senv)))
    close(h);

end