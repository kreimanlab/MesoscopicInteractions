close all;

if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        %resultsDir = '/media/klab/D0BA8713BA86F56E/data/results';
        resultsDir = '/mnt/cuenap2/data/results/coh_w30';
        h5Dir = '/media/klab/44/h5';
    elseif strcmp(strip(hname),'ubuntu_1604')
        resultsDir = '/nas_share/RawData/data/results/coh_w30';
        h5Dir = '/nas_share/RawData/scripts/synth/out';
    else
        resultsDir = '/mnt/cuenap2/data/results';
        h5Dir = '/mnt/cuenap2/scripts/synth/out';
    end
end

if (~exist('A','var'))
    A = Analysis(resultsDir,h5Dir);
end

system('mkdir figures/check_fit');

for fidx = 1:A.r.n_f

    % index into results
    %fidx = 26;

    graphfname = sprintf('%s/%s',resultsDir,A.r.graph_f{fidx});
    permfname = sprintf('%s/%s',resultsDir,A.r.perm_f{fidx});
    distfname = sprintf('%s/%s',resultsDir,A.r.dist_f{fidx});
    metric = A.r.metrics{fidx};
    
    if (startsWith(metric,'pc'))
        sid = A.r.subjects{fidx};
        sub_idx = find(strcmp(A.h5eeg.subjects,sid),1);
        n_chan = A.h5eeg.n_bchan{sub_idx};

        % Read permutations
        R = h5read(permfname,'/R');

        % Read fitted distribution
        D = load(distfname);

        % Read graph
        Rg = h5read(graphfname,'/R');
        
        % Threshold artifacts from graph
        Rg(Rg >= 0.8) = NaN;

        % Plot distributions
        [n_comb,n_graph] = size(Rg);

        % Threshold by distance
        dmat = A.getDistanceMatrix(sub_idx);
        is_dthresh = false(1,n_comb);
        dvec = nan(1,n_comb);
        ic = 1;
        for i1 = 1:(n_chan - 1)
            for i2 = (i1 + 1):n_chan
                if (dmat(i1,i2) > A.const.DIST_THRESH_MM)
                    is_dthresh(ic) = true;
                    dvec(ic) = dmat(i1,i2);
                end
                ic = ic + 1;
            end
        end

        h = figure;
        set(h,'Position',[0 0 1440 900])
        x = linspace(0,0.8,400);
        n_plot = n_comb;
        cc = gray(n_plot);
        cc = cc(randperm(n_plot),:);
        Iperm = 1:n_comb;
        Iperm = Iperm(randperm(n_comb));
        for i = 1:n_plot
            if (is_dthresh(Iperm(i)))
                [f_null,~] = hist(R(Iperm(i),:),x);
                [f_data,~] = hist(Rg(Iperm(i),:),x);
                
                col_i = round((dvec(Iperm(i)) - min(dvec) )/( max(dvec) - min(dvec) )*(length(cc)-1)+1); 
                
                plot(x,f_null/trapz(x,f_null),'-','color',0.4*[1 0.2 0.2],'LineWidth',1); hold on;
                plot(x,pdf(D.d{Iperm(i)},x),'-','color',0.4*[1 0.2 0.2],'LineWidth',1); hold on;
                plot(x,f_data/trapz(x,f_data),'-','color',1*cc(col_i,:),'LineWidth',1); hold on;
            end
        end
        set(gca,'TickDir','out')
        xlabel(sprintf('Metric: %s',metric))
        ylabel('Probability Density Function')
        title(sprintf('Subject: %s',sid))
        axis tight
        print(h,sprintf('figures/check_fit/%s_%s',sid,metric),'-depsc');
        %return;
        close(h);

    end

end