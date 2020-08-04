close all;
clear;

metric = 'pcBroadband';
res_dir = '/media/klab/D0BA8713BA86F56E/data/results/coh_w10_2';
h5_dir = '/media/klab/KLAB101/h5_notch20'; %/mnt/cuenap/data/h5_notch20
art_dir = '/media/klab/D0BA8713BA86F56E/data/h5_art_50uv/art_idx';
coreg_dir = '/mnt/cuenap_ssd/coregistration';

dir_times = './figure_times';

system(sprintf('mkdir %s',dir_times));

% Get list of subjects with complete computations
%   S           List of subject names
%   S_art       List of artifact file names (.mat)
%   S_graph     List of graph file names (.h5)
%   S_dist      List of dist file names (.mat)
%   S_h5        List of h5 file names (.h5)
load_flist;

% Loop through subjects
for j = 1:19
for i = 1:length(S)
    tic;
    
    % Init
    fprintf('Processing: %s\n',S{i})
    ecog = H5eeg(S_h5{i});
    bip = h5readatt(S_h5{i},'/h5eeg/eeg','bip');
    [n_bchan,~] = size(bip);
    n_comb = nchoosek(n_bchan,2);
    
    % Read artifact matrix
    A = load(S_art{i});
    [n_comb2, n_w] = size(A.art_idx);
    
    % Read graph
    G.R = h5read(S_graph{i},'/R');
    G.chan1 = h5read(S_graph{i},'/chan1') + 1;
    G.chan2 = h5read(S_graph{i},'/chan2') + 1;
    G.w = h5read(S_graph{i},'/w');
    [n_comb3, n_w2] = size(G.R);
    
    % Read dist
    D = load(S_dist{i});
    n_comb4 = length(D.d);
    
    % Check number of combinations
    if (length(unique([n_comb,n_comb2,n_comb3,n_comb4])) ~= 1)
        fprintf(2,'E> Mismatch in number of combinations.\n')
        return;
    end
    
    % Check number of metrics
    if (abs(n_w2  - n_w) > 1)
        fprintf(2,'W> %s number of graph points (%i) >> artifacts (%i).\n',S{i},n_w2,n_w);
    else
        % Show graph
        h = figure;
        set(h,'Units','Normalized','Position',[0 0 1 1]);
        set(h,'PaperUnits','inches');
        set(h,'PaperPosition',[0 0 20 20]);
        comb_i = randi([1 n_comb]);
        g = G.R(comb_i,1:(end-1));
        t = linspace(0,n_w*(10/3600),n_w);
        art_i = A.art_idx(comb_i,:);
        
        % Get timeseries of maximum coherence
        g_clean = g(~art_i);
        [max_g,max_i] = max(g_clean);
        if (isempty(max_i))
            max_i = 1;
            max_g = g(max_i);
        else
            max_i = find(g == g_clean(max_i),1);
        end
        start_idx = double((round(ecog.fs)*G.w) * (max_i-1) + 1);
        
        % Plot coherence over time
        subplot(4,1,1)
        plot(t(art_i),g(art_i),'.','color',1*[1 0 0],'MarkerSize',5);
        hold on;
        plot(t(~art_i),g(~art_i),'black.','MarkerSize',5);
        hold on;
        plot(t(max_i),max_g,'o','color',[0 1 0],'MarkerSize',5);
        b1l = A.bchan_labels{G.chan1(comb_i)};
        b2l = A.bchan_labels{G.chan2(comb_i)};
        title(sprintf('%s %s:%s',S{i},b1l,b2l))
        axis([min(t) max(t) 0 1])
        xlabel('Hours')
        ylabel(sprintf('Coherence (%s)',metric))
        
        % Plot null
        subplot(4,1,2)
        g0 = random(D.d{comb_i},size(g));
        plot(t,g0,'black.','MarkerSize',5);
        axis([min(t) max(t) 0 1])
        title('Null')
        xlabel('Hours')
        ylabel(sprintf('Coherence (%s)',metric))
        
        
        b1c1 = bip(G.chan1(comb_i),1);
        b1c2 = bip(G.chan1(comb_i),2);
        b2c1 = bip(G.chan2(comb_i),1);
        b2c2 = bip(G.chan2(comb_i),2);
        w_samples = double(round(ecog.fs)*G.w);
        t_w = linspace(0,double(G.w),w_samples);
        v1c1 = h5read(S_h5{i},'/h5eeg/eeg',[b1c1 start_idx],[1 w_samples]);
        v1c2 = h5read(S_h5{i},'/h5eeg/eeg',[b1c2 start_idx],[1 w_samples]);
        v2c1 = h5read(S_h5{i},'/h5eeg/eeg',[b2c1 start_idx],[1 w_samples]);
        v2c2 = h5read(S_h5{i},'/h5eeg/eeg',[b2c2 start_idx],[1 w_samples]);
        v1 = v1c1 - v1c2;
        v2 = v2c1 - v2c2;
        
        subplot(4,1,3);
        plot(t_w,v1,'black-','LineWidth',1)
        xlabel('Seconds')
        ylabel(sprintf('%s-%s (uV)',A.chan_labels{b1c1},A.chan_labels{b1c2}))
        axis tight;
        title(sprintf('%s: %.4f',metric,max_g))
        
        subplot(4,1,4);
        plot(t_w,v2,'black-','LineWidth',1)
        xlabel('Seconds')
        ylabel(sprintf('%s-%s (uV)',A.chan_labels{b2c1},A.chan_labels{b2c2}))
        axis tight;
        
        print(h,sprintf('%s/%s_%s_combi-%i_maxi-%i_maxg-%i',dir_times,S{i},metric,comb_i,max_i,round(max_g*1000)),'-dpng');
        close(h)
        %return
    end
    %n_w = min([n_w,n_w2]);
    
    
    
    % Time
    t_sing = toc;
    fprintf('\tETA: %.2f min\n',(t_sing/60)*(length(S)-i));
end
end

fprintf('[!] All Done.\n')