close all;

fprintf('--- figure T1 ---\n')

if (~exist('AT','var'))
    load('xsub2/xsub_coh_all_v6.mat');
end

if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix    
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        resultsDir = '/media/klab/KLAB101/results/coh_w10';
        h5Dir = '/media/klab/KLAB101/h5_notch20';
        dir_art = '/media/klab/internal/data/h5_notch20/art';
        %resultsDir = '/media/klab/44/data/results';
        %h5Dir = '/media/klab/44/h5';
    elseif strcmp(strip(hname),'ubuntu_1604')
        resultsDir = '/nas_share/RawData/data/results';
        h5Dir = '/nas_share/RawData/scripts/synth/out';
    else
        resultsDir = '/mnt/cuenap2/data/results';
        h5Dir = '/mnt/cuenap2/scripts/synth/out';
    end
end

clear A;
if (~exist('A','var'))
    A = Analysis(resultsDir,h5Dir);
end

system('mkdir figures/figure_T1');
                                                                 
%                                                                          
% number of seconds to display
w = 60; % number of seconds used for coherence
w_show = 10; % number of seconds to plot
w_shift = 3600; % number of seconds for time shift
pick_fresh = false; % whether to compute coherence fresh
cohBro_threshold = 0.3; % only plot interactions with coh above this thresh
plot_diagnostic = false;
plot_t3 = false;

% ROI pair (according to xsub/xsub_coh_all.mat)
%
% AT.P.AtlNames: 
%    'Destrieux-2009'    'Desikan-Killiany'    'Brodmann'    'Yeo-2011-7'    'Yeo-2011-17'    'HCP-MMP1'    'M132'    'FVE91'    'FVEall'
%    'LVE00'    'PHT00'    'Brodmann'    'BoninBailey'    'FerryEtAl'    'UD86'    'PGR91'    'LANDMARK'    'LyonKaas'    'BRL87'
%    'SP78'
atl_idx = 2;

% load ROI labels
roi_labels = AT.P.AtlROIs{atl_idx}.LH.struct_names;
% Apply ROI clustering
roi_labels = roi_labels(cluster_i);
% Remove extra ROIs
roi_labels = roi_labels(com_i);
% 24    - DK:parsorbitalis
% 9     - DK:inferiorparietal
roi_i = 24;
roi_j = 9;

fprintf('Atlas:\t%s\n',AT.P.AtlNames{atl_idx})
fprintf('ROI 1:\t%s\n',roi_labels{roi_i})
fprintf('ROI 2:\t%s\n',roi_labels{roi_j})
roi_one = roi_labels{roi_i};
roi_two = roi_labels{roi_j};

Subjects = adjct_sub{roi_i,roi_j};
Subjects = reshape(Subjects,6,[])';
bchan = adjct_bchan{roi_i,roi_j};
ct = adjct_xsub{roi_i,roi_j};
mag = adjct_xsubM{roi_i,roi_j};

% get number of bipolar pairs responsible for roi pair
[~,n_bpair] = size(bchan);
I2_shuff = randperm(n_bpair);
for i2_i = 1:n_bpair

    
    i2 = I2_shuff(i2_i);
    sid = Subjects(i2,:);
    i = find(strcmp(A.h5eeg.subjects,sid),1);
    Fs = round(A.h5eeg.fs{i});
    if (isempty(i))
        fprintf(2,'E> Could not find subject %s\n',sid);
        break;
    end
    %sid = A.h5eeg.subjects{i};
    dmat = A.getDistanceMatrix(i);
    
    % Channel selection
%     chanBip = randi([1 A.h5eeg.n_bchan{i}]);
%     chan2Bip = randi([1 A.h5eeg.n_bchan{i}]);
    chanBip = bchan(1,i2);
    chan2Bip = bchan(2,i2);

    % Sample selection, remove artifacts
    if (pick_fresh)
        arts = h5read(A.h5eeg.filenames{i},'/h5eeg/artifacts');
        arts_w = h5readatt(A.h5eeg.filenames{i},'/h5eeg/artifacts','width');
        
        % select artifacts that are above coherence threshold
        isart = true;
        isbelowt = true;
        while ((isart ) || (isbelowt))
            samp = randi([1 (A.h5eeg.n_samples{i}-w*round(A.h5eeg.fs{i}))]);
            isart = (arts(1,ceil(samp/arts_w)) ~= 0);
            
            bchan1 = A.h5eeg.bip{i}(chanBip,1);
            bchan2 = A.h5eeg.bip{i}(chanBip,2);
            vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan1 samp],[1 w*round(A.h5eeg.fs{i})]);
            vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan2 samp],[1 w*round(A.h5eeg.fs{i})]);
            v = vb1 - vb2;
            
            b2chan1 = A.h5eeg.bip{i}(chan2Bip,1);
            b2chan2 = A.h5eeg.bip{i}(chan2Bip,2);
            vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan1 samp],[1 w*round(A.h5eeg.fs{i})]);
            vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan2 samp],[1 w*round(A.h5eeg.fs{i})]);
            v2 = vb1 - vb2;
            
            [ R ] = coherence( v, v2, Fs );
            %isbelowt = R(6) < cohBro_threshold;
            isbelowt = sum(R < cohBro_threshold) >= 1;
            fprintf('> isArt: %i\tcohBro: %.6f / %.6f\n',isart,R(6),cohBro_threshold)
        end
        
    else
        gfname = sprintf('%s/%s_graph-pcBroadband.h5',resultsDir,sid);
        h5i = h5info(gfname);
        ds = (h5i.Datasets.Dataspace);
        Ri = nan;
        Ric = 1;
        for ii = 1:(A.h5eeg.n_bchan{i}-1)
            for jj = (ii+1):A.h5eeg.n_bchan{i}
                if (((ii == chanBip) && (jj == chan2Bip)) || ((jj == chanBip) && (ii == chan2Bip)) )
                    Ri = Ric;
                    break
                end
                Ric = Ric + 1;
            end
            if (~ isnan(Ri))
                break
            end
        end
        R1 = h5read(gfname,'/R',[Ri 1],[1 ds.Size(2)]);
        arts = h5read(A.h5eeg.filenames{i},'/h5eeg/artifacts');
        arts_w = h5readatt(A.h5eeg.filenames{i},'/h5eeg/artifacts','width');
        R1arts = false(size(R1));
        for iii = 1:length(R1)
            arts_start = (iii-1) * w + 1;
            arts_end = arts_start + w - 1;
            if (arts_end > length(arts(1,:)))
                arts_end = length(arts(1,:));
            end
            art_ratio = sum(arts(1,arts_start:arts_end) ~= 0)/(arts_end-arts_start+1);
            R1arts(iii) = (art_ratio > 0);
        end
        % testing
        if (plot_diagnostic)
            system('mkdir figures/figure_T1/test');
            
            for i4 = 1:length(R1)

                samp = (i4 - 1)*w*round(A.h5eeg.fs{i}) + 1;
                
                bchan1 = A.h5eeg.bip{i}(chanBip,1);
                bchan2 = A.h5eeg.bip{i}(chanBip,2);
                vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan1 samp],[1 w*round(A.h5eeg.fs{i})]);
                vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan2 samp],[1 w*round(A.h5eeg.fs{i})]);
                v = vb1 - vb2;

                b2chan1 = A.h5eeg.bip{i}(chan2Bip,1);
                b2chan2 = A.h5eeg.bip{i}(chan2Bip,2);
                vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan1 samp],[1 w*round(A.h5eeg.fs{i})]);
                vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan2 samp],[1 w*round(A.h5eeg.fs{i})]);
                v2 = vb1 - vb2;
        
                arts_start = (i4-1) * w + 1;
                arts_end = arts_start + w - 1;
                if (arts_end > length(arts(1,:)))
                    arts_end = length(arts(1,:));
                end
                V_isart = (arts(1,arts_start:arts_end) ~= 0);
                v_isart = zeros(size(v));
                for i5 = 1:length(v_isart)
                    v_isart(i5) = V_isart(ceil(i5/A.h5eeg.fs{i}));
                end
                
                h = figure;
                %set(h,'PaperPositionMode','auto');         
                set(h,'PaperOrientation','landscape');
                set(h,'PaperUnits','normalized');
                set(h,'PaperPosition', [0 0 1 1]);
                %set(h,'Position',[0 0 1920 1200]);
                subplot(3,1,1)
                R1b = R1;
                plot(R1,'.','color',0.7*[1 1 1],'MarkerSize',1); hold on; 
                R1b(R1arts) = nan; 
                plot(R1b,'black.','MarkerSize',2); hold on;
                plot(i4,R1(i4),'redo','MarkerSize',1);
                axis([1 length(R1) 0 1]);
                ylabel('Coherence');xlabel('Minutes');
                title(sprintf('%s: %s %s',sid,A.h5eeg.bip_labels{i}{chanBip},A.h5eeg.bip_labels{i}{chan2Bip}));
                
                subplot(3,1,2)
                t = linspace(0,w,w*round(A.h5eeg.fs{i}));
                plot(t,v,'black'); hold on;
                vb = v;
                vb(~v_isart) = nan;
                plot(t,vb,'color',0.7*[1 1 1]);
                axis tight;
                xlabel('Seconds')
                ylabel(sprintf('%s IFP (uV)',A.h5eeg.bip_labels{i}{chanBip}))
                subplot(3,1,3)
                t = linspace(0,w,w*round(A.h5eeg.fs{i}));
                plot(t,v2,'black'); hold on;
                vb = v2;
                vb(~v_isart) = nan;
                plot(t,vb,'color',0.7*[1 1 1]);
                axis tight;
                xlabel('Seconds')
                ylabel(sprintf('%s IFP (uV)',A.h5eeg.bip_labels{i}{chan2Bip}))
                print(h,sprintf('figures/figure_T1/test/%s_%s_%s_coh-%i_min-%i',...
                    sid,A.h5eeg.bip_labels{i}{chanBip},A.h5eeg.bip_labels{i}{chan2Bip},round(1000*R1(i4)),i4),'-dpdf');
                close(h);
                %return
            end
        end

        R1rm = R1;
        R1rm(R1arts) = nan;
        samp_min = find(R1rm > cohBro_threshold,1);
        samp =  (samp_min - 1) * w * round(A.h5eeg.fs{i}) + 1; %randi([1 (A.h5eeg.n_samples{i}-w*round(A.h5eeg.fs{i}))]);
    end
    fprintf('Time:\t%i (%.2f hrs)\n',samp,samp/(A.h5eeg.fs{i}*3600))

    if (isempty(samp_min))
        fprintf(2,'W> No coherence values found above threshold: %.6f\n',cohBro_threshold);
    else
        % Read
        if (~ pick_fresh)
            bchan1 = A.h5eeg.bip{i}(chanBip,1);
            bchan2 = A.h5eeg.bip{i}(chanBip,2);
            vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan1 samp],[1 w*round(A.h5eeg.fs{i})]);
            vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[bchan2 samp],[1 w*round(A.h5eeg.fs{i})]);
            v = vb1 - vb2;

            b2chan1 = A.h5eeg.bip{i}(chan2Bip,1);
            b2chan2 = A.h5eeg.bip{i}(chan2Bip,2);
            vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan1 samp],[1 w*round(A.h5eeg.fs{i})]);
            vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan2 samp],[1 w*round(A.h5eeg.fs{i})]);
            v2 = vb1 - vb2;
        end

        vb1 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan1 samp+w_shift*round(A.h5eeg.fs{i})],[1 w*round(A.h5eeg.fs{i})]);
        vb2 = h5read(A.h5eeg.filenames{i},'/h5eeg/eeg',[b2chan2 samp+w_shift*round(A.h5eeg.fs{i})],[1 w*round(A.h5eeg.fs{i})]);
        v2_0 = vb1 - vb2;

        fprintf('\n%s\n',A.h5eeg.subjects{i})
        fprintf('bchan1: %i\n',bchan1)
        fprintf('bchan2: %i\n',bchan2)
        fprintf('dist bchan1,bchan2: %.4f\n',A.h5eeg.bip{i}(chanBip,3))
        fprintf('b2chan1: %i\n',b2chan1)
        fprintf('b2chan2: %i\n',b2chan2)
        fprintf('dist b2chan1,b2chan2: %.4f\n',A.h5eeg.bip{i}(chan2Bip,3))

        % --- Plot ------------------------------------------------------------ 
        h = figure;
        fig_w = 8.5;
        fig_h = 11.0;
        fig_n_w0 = 0.05; % width starting position
        fig_n_w = 0.6; % width
        fig_n_ws = 0.001; % space between brain
        fig_n_wbrain = 0.2; % width brain
        fig_n_h0 = 0.048; % height per axis
        atlas = '';
        set(h,'Position',[0 0 fig_w*100 fig_h*100])
        set(h, 'PaperUnits', 'Inches')
        set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
        set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
        % --- raw pair ---
        ax1 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.90 fig_n_w fig_n_h0]);
        title(sprintf('%s',sid));
        A.plotLine(ax1,v(1:(round(w_show*A.h5eeg.fs{i}))),i);    
        ax1b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.80 fig_n_w fig_n_h0]);
        A.plotLine(ax1b,v2(1:(round(w_show*A.h5eeg.fs{i}))),i);
        hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
        set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.90 fig_n_wbrain fig_n_h0]);
        hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
        set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.80 fig_n_wbrain fig_n_h0]);

        % --- time shifted pair ---
        ax1 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.60 fig_n_w fig_n_h0]);
        title(sprintf('%s',sid));
        A.plotLine(ax1,v(1:(round(w_show*A.h5eeg.fs{i}))),i);    
        ax1b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.50 fig_n_w fig_n_h0]);
        A.plotLine(ax1b,v2_0(1:(round(w_show*A.h5eeg.fs{i}))),i);
        hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
        set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.60 fig_n_wbrain fig_n_h0]);
        hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
        set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.50 fig_n_wbrain fig_n_h0]);


    %     
    %     % --- filt pair ---
    %     vf = A.del(v,i);
    %     v2f = A.del(v2,i);
    %     ax2 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.80 fig_n_w fig_n_h0]);
    %     A.plotLine(ax2,vf,i);
    %     ax2b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.75 fig_n_w fig_n_h0]);
    %     A.plotLine(ax2b,v2f,i);
    %     hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    %     set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.80 fig_n_wbrain fig_n_h0]);
    %     hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    %     set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.75 fig_n_wbrain fig_n_h0]);
    %     % --- env pair ---
    %     vfe = A.env(vf);
    %     v2fe = A.env(v2f);
    %     ax3 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.70 fig_n_w fig_n_h0]);
    %     A.plotLine(ax3,vfe,i);
    %     ax3b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.65 fig_n_w fig_n_h0]);
    %     A.plotLine(ax3b,v2fe,i);
    %     hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
    %     set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.70 fig_n_wbrain fig_n_h0]);
    %     hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
    %     set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.65 fig_n_wbrain fig_n_h0]);
    %     r_s = corr(v', v2', 'Type','Spearman');
    %     r_senv = corr(vfe', v2fe', 'Type','Spearman');


        [ R ] = coherence( v, v2, Fs );
        coh = round(1000*R);
        [ R0 ] = coherence( v, v2_0, Fs );
        coh0 = round(1000*R0);

        % Text

        A.saveFig(h,sprintf('./figures/figure_T1/%s_%s_%s_%i_%s_%s_%i_samp-%i_w-%i_coh_%i_%i_%i_%i_%i_%i_coh0_%i_%i_%i_%i_%i_%i',...
            sid,roi_one,A.h5eeg.bip_labels{i}{chanBip},chanBip,roi_two,A.h5eeg.bip_labels{i}{chan2Bip},chan2Bip,...
            samp,w,coh(1),coh(2),coh(3),coh(4),coh(5),coh(6),...
            coh0(1),coh0(2),coh0(3),coh0(4),coh0(5),coh0(6)))
        close(h);
        
        % --- Plot how to compute coherence -------------------------------
        if (plot_t3)
            % Constants
            PLI_S = 56;      % Power line interference frequency                          
            PLI_E = 64;                                                                   
            PL2_S = (Fs-180)-2;     % Power line interference frequency second band       
            PL2_E = (Fs-180)+2;                                                           
            PL3_S = 117;      % Power line interference frequency third band              
            PL3_E = 123; 
            %
            DEL_S = 0.5;     % Delta wave                                                 
            DEL_E = 3;                                                                    
            THE_S = 3;       % Theta wave                                                 
            THE_E = 8;                                                                    
            ALP_S = 8;       % Alpha wave                                                 
            ALP_E = 12;                                                                   
            BET_S = 12;      % Beta wave                                                  
            BET_E = 25;                                                                   
            GAM_S = 25;      % Gamma wave                                                 
            GAM_E = 100;                                                                  
            BRO_S = 0.5;     % Broadband                                                  
            BRO_E = 125; 

            [cxx,f] = mscohere(v', v2',hamming(2*Fs),[],[],Fs);
            cxx = sqrt(cxx);

            system('mkdir figures/figure_T3');

            % --- DELTA -------------------------------------------------------
            h = figure;
            set(h,'PaperOrientation','landscape');
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);

            xlabel('Frequency (Hz)')
            ylabel('Magnitude of Coherence')
            axis([0 max(f) 0 1])
            rectangle('Position',[0 0 DEL_S 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            rectangle('Position',[DEL_E 0 max(f) 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            plot(f, cxx,'black'); hold on;

            print(h,sprintf('figures/figure_T3/del_%s_%s_%s_samp-%i',...
                sid,A.h5eeg.bip_labels{i}{chanBip},A.h5eeg.bip_labels{i}{chan2Bip},samp),'-dpdf');
            close(h);

            % --- THETA -------------------------------------------------------
            h = figure;
            set(h,'PaperOrientation','landscape');
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);

            xlabel('Frequency (Hz)')
            ylabel('Magnitude of Coherence')
            axis([0 max(f) 0 1])
            rectangle('Position',[0 0 THE_S 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            rectangle('Position',[THE_E 0 max(f) 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            plot(f, cxx,'black'); hold on;

            print(h,sprintf('figures/figure_T3/the_%s_%s_%s_samp-%i',...
                sid,A.h5eeg.bip_labels{i}{chanBip},A.h5eeg.bip_labels{i}{chan2Bip},samp),'-dpdf');
            close(h);

            % --- ALPHA -------------------------------------------------------
            h = figure;
            set(h,'PaperOrientation','landscape');
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);

            xlabel('Frequency (Hz)')
            ylabel('Magnitude of Coherence')
            axis([0 max(f) 0 1])
            rectangle('Position',[0 0 ALP_S 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            rectangle('Position',[ALP_E 0 max(f) 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            plot(f, cxx,'black'); hold on;

            print(h,sprintf('figures/figure_T3/alp_%s_%s_%s_samp-%i',...
                sid,A.h5eeg.bip_labels{i}{chanBip},A.h5eeg.bip_labels{i}{chan2Bip},samp),'-dpdf');
            close(h);

            % --- BETA --------------------------------------------------------
            h = figure;
            set(h,'PaperOrientation','landscape');
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);

            xlabel('Frequency (Hz)')
            ylabel('Magnitude of Coherence')
            axis([0 max(f) 0 1])
            rectangle('Position',[0 0 BET_S 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            rectangle('Position',[BET_E 0 max(f) 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            plot(f, cxx,'black'); hold on;

            print(h,sprintf('figures/figure_T3/bet_%s_%s_%s_samp-%i',...
                sid,A.h5eeg.bip_labels{i}{chanBip},A.h5eeg.bip_labels{i}{chan2Bip},samp),'-dpdf');
            close(h);

            % --- GAMMA --------------------------------------------------------
            h = figure;
            set(h,'PaperOrientation','landscape');
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);

            xlabel('Frequency (Hz)')
            ylabel('Magnitude of Coherence')
            axis([0 max(f) 0 1])
            rectangle('Position',[0 0 GAM_S 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            rectangle('Position',[GAM_E 0 max(f) 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            rectangle('Position',[PLI_S 0 (PLI_E-PLI_S) 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            rectangle('Position',[PL2_S 0 (PL2_E-PL2_S) 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            plot(f, cxx,'black'); hold on;

            print(h,sprintf('figures/figure_T3/gam_%s_%s_%s_samp-%i',...
                sid,A.h5eeg.bip_labels{i}{chanBip},A.h5eeg.bip_labels{i}{chan2Bip},samp),'-dpdf');
            close(h);

            % --- BROADBAND ---------------------------------------------------
            h = figure;
            set(h,'PaperOrientation','landscape');
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);

            xlabel('Frequency (Hz)')
            ylabel('Magnitude of Coherence')
            axis([0 max(f) 0 1])
            rectangle('Position',[0 0 BRO_S 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            rectangle('Position',[BRO_E 0 max(f) 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            rectangle('Position',[PLI_S 0 (PLI_E-PLI_S) 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            rectangle('Position',[PL2_S 0 (PL2_E-PL2_S) 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            rectangle('Position',[PL3_S 0 (PL3_E-PL3_S) 1],'FaceColor',0.7*[1 1 1],'LineStyle','none'); hold on;
            plot(f, cxx,'black'); hold on;

            print(h,sprintf('figures/figure_T3/bro_%s_%s_%s_samp-%i',...
                sid,A.h5eeg.bip_labels{i}{chanBip},A.h5eeg.bip_labels{i}{chan2Bip},samp),'-dpdf');
            close(h);
        end
        
        %return
    end
    
    %break;
end

% 
% % Plot Theta
% h = figure;
% fig_w = 8.5;
% fig_h = 11.0;
% fig_n_w0 = 0.05; % width starting position
% fig_n_w = 0.4; % width
% fig_n_ws = 0.001; % space between brain
% fig_n_wbrain = 0.2; % width brain
% fig_n_h0 = 0.048; % height per axis
% %atlas = 'd';
% set(h,'Position',[0 0 fig_w*100 fig_h*100])
% set(h, 'PaperUnits', 'Inches')
% set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
% set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
% % --- raw pair ---
% ax1 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.90 fig_n_w fig_n_h0]);
% A.plotLine(ax1,v,i);    
% ax1b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.85 fig_n_w fig_n_h0]);
% A.plotLine(ax1b,v2,i);
% hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.90 fig_n_wbrain fig_n_h0]);
% hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.85 fig_n_wbrain fig_n_h0]);
% % --- filt pair ---
% vf = A.the(v,i);
% v2f = A.the(v2,i);
% ax2 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.80 fig_n_w fig_n_h0]);
% A.plotLine(ax2,vf,i);
% ax2b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.75 fig_n_w fig_n_h0]);
% A.plotLine(ax2b,v2f,i);
% hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.80 fig_n_wbrain fig_n_h0]);
% hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.75 fig_n_wbrain fig_n_h0]);
% % --- env pair ---
% vfe = A.env(vf);
% v2fe = A.env(v2f);
% ax3 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.70 fig_n_w fig_n_h0]);
% A.plotLine(ax3,vfe,i);
% ax3b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.65 fig_n_w fig_n_h0]);
% A.plotLine(ax3b,v2fe,i);
% hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.70 fig_n_wbrain fig_n_h0]);
% hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.65 fig_n_wbrain fig_n_h0]);
% r_s = corr(v', v2', 'Type','Spearman');
% r_senv = corr(vfe', v2fe', 'Type','Spearman');
% A.saveFig(h,sprintf('./figures/voltages/filtering_%s_chans%i-%i_samp-%i_w-%i_rs_%i_rsenv_%i_the',...
%     sid,bchan1,b2chan1,samp,w,round(1000*r_s),round(1000*r_senv)))
% close(h);
% 
% % Plot Alpha
% h = figure;
% fig_w = 8.5;
% fig_h = 11.0;
% fig_n_w0 = 0.05; % width starting position
% fig_n_w = 0.4; % width
% fig_n_ws = 0.001; % space between brain
% fig_n_wbrain = 0.2; % width brain
% fig_n_h0 = 0.048; % height per axis
% %atlas = 'd';
% set(h,'Position',[0 0 fig_w*100 fig_h*100])
% set(h, 'PaperUnits', 'Inches')
% set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
% set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
% % --- raw pair ---
% ax1 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.90 fig_n_w fig_n_h0]);
% A.plotLine(ax1,v,i);    
% ax1b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.85 fig_n_w fig_n_h0]);
% A.plotLine(ax1b,v2,i);
% hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.90 fig_n_wbrain fig_n_h0]);
% hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.85 fig_n_wbrain fig_n_h0]);
% % --- filt pair ---
% vf = A.alp(v,i);
% v2f = A.alp(v2,i);
% ax2 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.80 fig_n_w fig_n_h0]);
% A.plotLine(ax2,vf,i);
% ax2b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.75 fig_n_w fig_n_h0]);
% A.plotLine(ax2b,v2f,i);
% hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.80 fig_n_wbrain fig_n_h0]);
% hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.75 fig_n_wbrain fig_n_h0]);
% % --- env pair ---
% vfe = A.env(vf);
% v2fe = A.env(v2f);
% ax3 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.70 fig_n_w fig_n_h0]);
% A.plotLine(ax3,vfe,i);
% ax3b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.65 fig_n_w fig_n_h0]);
% A.plotLine(ax3b,v2fe,i);
% hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.70 fig_n_wbrain fig_n_h0]);
% hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.65 fig_n_wbrain fig_n_h0]);
% r_s = corr(v', v2', 'Type','Spearman');
% r_senv = corr(vfe', v2fe', 'Type','Spearman');
% A.saveFig(h,sprintf('./figures/voltages/filtering_%s_chans%i-%i_samp-%i_w-%i_rs_%i_rsenv_%i_alp',...
%     sid,bchan1,b2chan1,samp,w,round(1000*r_s),round(1000*r_senv)))
% close(h);
% 
% % Plot Beta
% h = figure;
% fig_w = 8.5;
% fig_h = 11.0;
% fig_n_w0 = 0.05; % width starting position
% fig_n_w = 0.4; % width
% fig_n_ws = 0.001; % space between brain
% fig_n_wbrain = 0.2; % width brain
% fig_n_h0 = 0.048; % height per axis
% %atlas = 'd';
% set(h,'Position',[0 0 fig_w*100 fig_h*100])
% set(h, 'PaperUnits', 'Inches')
% set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
% set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
% % --- raw pair ---
% ax1 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.90 fig_n_w fig_n_h0]);
% A.plotLine(ax1,v,i);    
% ax1b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.85 fig_n_w fig_n_h0]);
% A.plotLine(ax1b,v2,i);
% hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.90 fig_n_wbrain fig_n_h0]);
% hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.85 fig_n_wbrain fig_n_h0]);
% % --- filt pair ---
% vf = A.bet(v,i);
% v2f = A.bet(v2,i);
% ax2 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.80 fig_n_w fig_n_h0]);
% A.plotLine(ax2,vf,i);
% ax2b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.75 fig_n_w fig_n_h0]);
% A.plotLine(ax2b,v2f,i);
% hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.80 fig_n_wbrain fig_n_h0]);
% hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.75 fig_n_wbrain fig_n_h0]);
% % --- env pair ---
% vfe = A.env(vf);
% v2fe = A.env(v2f);
% ax3 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.70 fig_n_w fig_n_h0]);
% A.plotLine(ax3,vfe,i);
% ax3b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.65 fig_n_w fig_n_h0]);
% A.plotLine(ax3b,v2fe,i);
% hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.70 fig_n_wbrain fig_n_h0]);
% hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.65 fig_n_wbrain fig_n_h0]);
% r_s = corr(v', v2', 'Type','Spearman');
% r_senv = corr(vfe', v2fe', 'Type','Spearman');
% A.saveFig(h,sprintf('./figures/voltages/filtering_%s_chans%i-%i_samp-%i_w-%i_rs_%i_rsenv_%i_bet',...
%     sid,bchan1,b2chan1,samp,w,round(1000*r_s),round(1000*r_senv)))
% close(h);
% 
% % Plot Gamma
% h = figure;
% fig_w = 8.5;
% fig_h = 11.0;
% fig_n_w0 = 0.05; % width starting position
% fig_n_w = 0.4; % width
% fig_n_ws = 0.001; % space between brain
% fig_n_wbrain = 0.2; % width brain
% fig_n_h0 = 0.048; % height per axis
% %atlas = 'd';
% set(h,'Position',[0 0 fig_w*100 fig_h*100])
% set(h, 'PaperUnits', 'Inches')
% set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
% set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
% % --- raw pair ---
% ax1 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.90 fig_n_w fig_n_h0]);
% A.plotLine(ax1,v,i);    
% ax1b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.85 fig_n_w fig_n_h0]);
% A.plotLine(ax1b,v2,i);
% hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.90 fig_n_wbrain fig_n_h0]);
% hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.85 fig_n_wbrain fig_n_h0]);
% % --- filt pair ---
% vf = A.gam(v,i);
% v2f = A.gam(v2,i);
% ax2 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.80 fig_n_w fig_n_h0]);
% A.plotLine(ax2,vf,i);
% ax2b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.75 fig_n_w fig_n_h0]);
% A.plotLine(ax2b,v2f,i);
% hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.80 fig_n_wbrain fig_n_h0]);
% hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.75 fig_n_wbrain fig_n_h0]);
% % --- env pair ---
% vfe = A.env(vf);
% v2fe = A.env(v2f);
% ax3 = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.70 fig_n_w fig_n_h0]);
% A.plotLine(ax3,vfe,i);
% ax3b = axes('Parent',h,'Units','Normalized','Position',[fig_n_w0 0.65 fig_n_w fig_n_h0]);
% A.plotLine(ax3b,v2fe,i);
% hb = brainplot_one(A.h5eeg.subjects{i},[bchan1,bchan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.70 fig_n_wbrain fig_n_h0]);
% hb = brainplot_one(A.h5eeg.subjects{i},[b2chan1,b2chan2],atlas); copyobj(hb.Children,h); close(hb);
% set(h.Children(1),'Position',[fig_n_w0+fig_n_ws+fig_n_w 0.65 fig_n_wbrain fig_n_h0]);
% r_s = corr(v', v2', 'Type','Spearman');
% r_senv = corr(vfe', v2fe', 'Type','Spearman');
% A.saveFig(h,sprintf('./figures/voltages/filtering_%s_chans%i-%i_samp-%i_w-%i_rs_%i_rsenv_%i_gam',...
%     sid,bchan1,b2chan1,samp,w,round(1000*r_s),round(1000*r_senv)))
% close(h);
