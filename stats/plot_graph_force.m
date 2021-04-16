close all;
clear;

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

for iM = 1:length(metrics)
    metric = metrics{iM};
    fprintf('[*] Metric: %s\n',metric);
    %Ca = load(sprintf('cache/fig_cluster2_sworld_circle_subjects_bip-%i.mat',iM));
    Ca = load(sprintf('cache/fig_cluster2_sworld_circle-%i.mat',iM));
    
    % Adjacency matrix
    A = Ca.Agr;
    for i = 1:length(A)
        A(i,i) = 0;
    end
    % define source and sink nodes
    idx_src = [8 11 19];
    idx_snk = [20];
    
    for k = 1:5
        Agr = A;
        
        % k-means cluster
        if (k>0)
            idx_mid = true(1,length(Agr));
            idx_mid(idx_src) = false;
            idx_mid(idx_snk) = false;
            
            numK = 1:length(Agr);
            numK = numK(idx_mid);
            
            AgrK = Agr(idx_mid,idx_mid);
            
            % Kmeans
            stream = RandStream('mlfg6331_64');  % Random number stream
            options = statset('UseParallel',1,'UseSubstreams',1,...
                'Streams',stream);
            kIdx = kmeans(AgrK,k,'Options',options,'MaxIter',10000,'Replicates',12*200);
            
            
            % shuffle order for plot
            [kIdx2,sIdx] = sort(kIdx);
            %numK2 = 1:length(Agr);
            %numK2(idx_mid) = sIdx;
            
            idx_src2 = 1:length(idx_src); % first come the sources
            idx_snk2 = (1:length(idx_snk)) + length(idx_src); % next come the sink

            % save k-means cluster index
            Clus = cell(1,k);
            for j = 1:k
                Clus{j} = find(kIdx2==j)+length(idx_src2)+length(idx_snk2);
            end
            
            idx_fstein = [idx_src,idx_snk,numK(sIdx)];
            Agr = Agr(idx_fstein,idx_fstein);
            %return
        end
        
        
        % Plot graph
        H = figure('visible','off','Position',[0 0 900 800]);
        G = graph(Agr);

        % strcmp(Ca.rois2_ab,['LIN'])
        h2 = plot(G,'Layout','layered','Source',idx_src2,'Sink',idx_snk2,...
            'Direction','right','NodeColor',[1 1 1]*0.1,'EdgeColor',[1 1 1]*0.1,'LineWidth',2);

        % Color weights
        % Determine node locations:
        cmap = inferno(1000);
        colormap(cmap);
        h2.EdgeCData = G.Edges.Weight;

        roi = Ca.rois2; %_ab
        roi = roi(idx_fstein);
        labelnode(h2,1:length(roi),roi);
        

        % Highlight source and sink nodes
        highlight(h2,idx_src2,'NodeColor',[1 1 1]*0.4); % 0.9*[13, 191, 126]/255
        highlight(h2,idx_snk2,'NodeColor',[1 1 1]*0.05); % 0.9*[25, 113, 255]/255
        node_cmap = turbo(k+2);
        node_cmap = node_cmap(2:(end-1),:);
        %node_cmap = turbo(k);
        for j = 1:k
            highlight(h2,Clus{j},'NodeColor',node_cmap(j,:));
        end
        
        
        col_var_min = min(G.Edges.Weight);
        col_var_max = max(G.Edges.Weight);
        cb = colorbar;
        caxis([col_var_min col_var_max])
        set(cb,'YTick',[col_var_min,col_var_max]);
        set(cb,'YTickLabel',{sprintf('%.2f',col_var_min),sprintf('%.2f',col_var_max)});
        set(cb,'TickLength',0);
        %set(cb,'Position',[0.8 0.1 0.01 0.4]);
        set(cb,'FontSize',10);
        if (iM == 1)
            cb.Label.String = sprintf('Coherence');
        else
            cb.Label.String = sprintf('%s Coherence',metric(3:end));
        end
        %set(gca,'Color',[1 1 1]*0.9)
        %colorbarlabel('Coherence');
        set(gca, 'DefaultFigureRenderer', 'painters');
        print(H,sprintf('./figures/plot_graph_force-%s_k-%i',metric(3:end),k),'-dpng','-r600'); % ,'-r900'
        close(H);
    end
    
end