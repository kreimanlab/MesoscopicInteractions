clear;
close all;

% Patients
Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'}; 
dir_results = '../results/coh_w10';
dir_art = '../h5_notch20/art_nosz';

atl = 2;
iM = 1;
Mags = {};
Pvals = {};
Rhos = {};
Dists = {};
for iSub = [3] %1:length(Subjects)
    sid = Subjects{iSub};
    
    % Load coherences
    fn_graph = sprintf('%s/%s_graph-%s.h5',dir_results,sid,metrics{iM});
    R = h5read(fn_graph,'/R');
    chan1 = h5read(fn_graph,'/chan1') + 1; % convert python to matlab
    chan2 = h5read(fn_graph,'/chan2') + 1;
    r_cols = h5read(fn_graph,'/r_cols');
    r_rows = h5read(fn_graph,'/r_rows');
    w = h5read(fn_graph,'/w');
    [n_comb,n_graph] = size(R);
    
    % Load artifacts
    fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
    art_idx = h5read(fn_art,'/art_idx');
    [~,n_art] = size(art_idx);
    if (n_art < n_graph)
        R = R(:,(1:n_art));
        [~,n_graph] = size(R);
    end
    
    % Load coherence thresholds
    fn_cache = sprintf('cache/xsub_out_%s_%i_atl%i',sid,iM,atl);
    C = load(fn_cache);
    
    % Binarize R
    R_bin = R;
    R_bin(art_idx == 1) = NaN;
    R_bin = (R_bin > C.coh_thresh);
    
    % Get significant magnitudes
    mag = C.mag;
    mag(C.mag < C.coh_thresh) = 0;
    mag(C.Dmats < C.dist_thresh) = NaN;
    
    % Compute usage correlations
    rhos = nan(1,n_comb);
    pvals = nan(1,n_comb);
    parfor i_comb = 1:n_comb
        if (mag(i_comb) > 0 )
            fprintf('[%s] %i of %i\n',sid,i_comb,n_comb);
            idx_comb = true(1,n_comb);
            idx_comb(i_comb) = false;
            usage = sum(R_bin(idx_comb,:))/n_comb;
            [rho,pval] = corr(usage',R_bin(i_comb,:)','Type','Spearman');
            rhos(i_comb) = rho;
            pvals(i_comb) = pval;
            
        end
    end
    
    Mags{iSub} = mag;
    Pvals{iSub} = pvals;
    Rhos{iSub} = rhos;
    Dists{iSub} = C.Dmats;
    
    %%
    % plot
    close all;
    h = figure('Position',[0,0,300,260]);
    X = mag(((mag>0)&(pvals'<0.05)));
    Y = rhos((mag>0)&(pvals'<0.05));
    plot(X,Y,'black.','MarkerSize',3);
    xlabel('Coherence');
    ylabel('Network Similarity');
    set(gca,'TickDir','out');
    box off;
    [r,p] = corr(X,Y','type','Spearman');
    print(h,'figures/constraint_coh_sim','-dpng','-r900');
    
    
    % Network Saturation
    h2 = figure('Position',[0,0,600,260]);
    i_comb = 1;
    idx_comb = true(1,n_comb);
    idx_comb(i_comb) = false;
    usage = sum(R_bin(idx_comb,:))/n_comb;
    %[rho,pval] = corr(usage',R_bin(i_comb,:)','Type','Spearman');
    T = ((1:length(usage))-1) .* (double(w)/(24*3600));
    plot(T,usage,'black-');
    xlabel('Time (days)');
    ylabel('Network Saturation');
    set(gca,'TickDir','out');
    box off;
    axis tight;
    print(h2,'figures/constraint_usage','-dpng','-r900');
    
    
    % Bipolar pair significance
    h3 = figure('Position',[0,0,600,160]);
    i_comb = 3;
    idx_comb = true(1,n_comb);
    idx_comb(i_comb) = false;
    usage = sum(R_bin(idx_comb,:))/n_comb;
    [rho,pval] = corr(usage',R_bin(i_comb,:)','Type','Spearman');
    
    fprintf('[!] corr: %.3f, p=%.3d, n=%i\n',rho,pval,length(usage))
    
    plot(T,R_bin(i_comb,:),'black.');
    %imagesc(R_bin(i_comb,:));
    xlabel('Time (days)');
    ylabel('Significance');
    set(gca,'TickDir','out');
    box off;
    axis tight;
    print(h3,'figures/constraint_rbin','-dpng','-r900');
    
    
    
    
    
    
    % Bipolar pair significance
    h3 = figure('Position',[0,0,600,160]);
    i_comb = 3;
    idx_comb = true(1,n_comb);
    idx_comb(i_comb) = false;
    usage = sum(R_bin(idx_comb,:))/n_comb;
    [rho,pval] = corr(usage',R_bin(i_comb,:)','Type','Spearman');
    
    fprintf('[!] corr: %.3f, p=%.3d, n=%i\n',rho,pval,length(usage))
    
    plot(T,R(i_comb,:),'black.');
    %imagesc(R_bin(i_comb,:));
    xlabel('Time (days)');
    ylabel('Coherence');
    set(gca,'TickDir','out');
    box off;
    axis tight;
    print(h3,'figures/constraint_r','-dpng','-r900');
    return
end