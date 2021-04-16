close all;
clear;
rng('shuffle');

%/media/klab/internal/data/coreg/fsaverage_sym/label/all_surf_ielvis_m00037.label
dir_artLp = '/media/klab/internal/data/h5_notch20/art';
dir_corLp = '/media/klab/internal/data/coreg';
%dir_resLp = '/home/jerry/data/results/coh_w10';
setenv('SUBJECTS_DIR',dir_corLp);
dir_cacheLp = './cache';
dir_h5Lp = '/media/klab/internal/data/h5_notch20';
mkdir('figures/T17');

metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

% Colors
BRAIN_MESH_ALPHA = 1;
COLOR_MEDIAL_WALL = 0.1*[1 1 1];
COLOR_NOSIG = 0.5*[1 1 1];
COLOR_BLACK = 0.3*[1 1 1];
SURFACE_TYPE = 'pial';
paper_width = 8.5;
paper_height = 6.5;

% Options
trig_show_roi_label = false;
trig_blot_nosig = false;

% --- register t14 adjacency matrix to 150 parcellation -------------------
cluster_i = load('cache/fig_cluster3_cluster_i.mat');
cluster_i = cluster_i.cluster_i;
CaT14 = load(sprintf('./cache/figure_t14_%i_150',1));
cluster_i2 = zeros(size(cluster_i));
for i = 1:length(cluster_i)
    idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(i))));
    cluster_i2(i) = idx;
end
% cluster_i = cluster_i2;
% cluster_i2: converts Es indices to final adjacency matrix indices
% -------------------------------------------------------------------------


% load average surface
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_corLp,'fsaverage_sym','r',SURFACE_TYPE));

% load atlas definition
CaC1 = load('cache/fig_cluster2_reduce_new_annot_vertex_color2.mat');
CaC2 = load('cache/fig_cluster2_reduce_annot');


for i = 1 %1:length(metricsp)
    metric = metricsp{i};
    
    % Load adjacency matrix
    CaT14 = load(sprintf('%s/figure_t14_%i_150.mat',dir_cacheLp,i));
    CaR = load(sprintf('%s/fig_cluster2_reduce_%i_new.mat',dir_cacheLp,i));
    
    % Define clustering variables
    %cl_adj = CaT14.Adj_plt2;
    
    % marginalized distances clustering
    cl_adj = CaT14.Adj_plt2_cl(CaT14.cluster_i,CaT14.cluster_i);
    cl_adj_nosig = all(isnan(cl_adj) | (cl_adj == 0));
    cl_adj(isnan(cl_adj)) = -1;
    
    % Fill missing data
    %cl_adj = fillmissing(cl_adj,'nearest');
    %cl_adj(isnan(cl_adj)) = 0;
    
    % nosig
    %cl_adj_nosig = all(isnan(cl_adj) | (cl_adj == 0));
    
    
    %======================================================================
    % k-means clustering
    fn_kmeans = sprintf('cache/figure_t17_%i_kmedoids',i);
    
    if (exist([fn_kmeans,'.mat'],'file'))
        load([fn_kmeans,'.mat']);
    else
        K = 2:17;
        Ssumd = zeros(length(K),1);
        Res = cell(length(K),4);
        A2km = cl_adj;
        %A2km(isnan(A2km)) = 0; % set diagonals to zeros;
        for ki = 1:length(K)
            k = K(ki);
            stream = RandStream('mlfg6331_64');  % Random number stream
            options = statset('UseParallel',1,'UseSubstreams',1,'Streams',stream);
            %[idx,C,sumd,D] = kmeans(A2km,k,'Options',options,'MaxIter',10000,'Display','final','Replicates',18000);
            
            [idx,C,sumd,D] = kmedoids(A2km,k,'distance',@nanDist4clustering_neg,'Options',options,'OnlinePhase','on','Replicates',4*4*4*4*4);
            Res{ki,1} = idx;
            Res{ki,2} = C;
            Res{ki,3} = sumd;
            Res{ki,4} = D;

            % cross-val metric
            if (length(K) > 1)
                cv = sum(sumd); 
                Ssumd(ki) = cv; %sum(sumd);
            end
            
            fprintf('[*] k = %i, %i of %i computed.\n',K(ki),ki,length(K))
        end
        save(fn_kmeans,'K','Ssumd','Res','A2km');
        fprintf('[*] %s\n',fn_kmeans);
    end
    %======================================================================
    
    
    
    
    % 3d coordinates on fsaverage_sym pial surface
    cl_coord = CaR.E(:,3:5);
    %cl_coord = cl_coord(cluster_i2,:);
    
    for k = 1:length(K)
        
        % kmeans region colors
        col_offset = 1;
        n_col = K(k) + 1 + col_offset;
        cc_k = getPyPlot_cMap('nipy_spectral', n_col, false, 'python3'); % 'gist_rainbow'
        cc_k = cc_k((1+col_offset):(K(k)+col_offset),:);
        [n_cck,~] = size(cc_k);
        for icol = 1:n_cck
            co = cc_k(icol,:);
            if (all(co==0))
                cc_k(icol,:) = COLOR_BLACK;
            end
        end
        
        h = figure('visible','off');
        set(h,'PaperUnits','inches');
        set(h,'PaperPosition',[0 0 paper_width paper_height]);
        axis off;
        hold all;
        
        %subplot(2,1,1);
        [ha, pos] = tight_subplot(2,2,[0.01 0.01],8*[0.01 0.01],4*[0.01 0.01]);
        axes(ha(1));
        
        p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
            'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
        %set(p,'FaceVertexCData',0.8*[1 1 1]);
        
        % region highlight by number
        vcolor = 0.8*ones(length(s_vert),3);
        roi_cl = Res{k,1};
        region_const = 74; % <- region to highlight
        for j = 1:length(s_vert)
            r = CaC2.region(j);
            idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(r))));
%             if (idx == region_const)
%                 vcolor(j,:) = [1 0.2 0.2];
%             end

            if ((trig_blot_nosig) && (cl_adj_nosig(idx)))
                col_kmeans = COLOR_NOSIG;
            else
                col_kmeans = cc_k(roi_cl(idx),:);
            end

            vcolor(j,:) = col_kmeans;
            if (CaC1.region_medial(j) == 1)
                vcolor(j,:) = COLOR_MEDIAL_WALL;
            end
        end
        set(p,'FaceVertexCData',vcolor);
        brainlight;
        %plot3(cl_coord(:,1),cl_coord(:,2),cl_coord(:,3),'black.')

        cl_roi_labels = zeros(length(cl_adj),1);
        for j = 1:length(cl_adj)
            idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(j))));
            if (trig_show_roi_label)
                text(cl_coord(j,1),cl_coord(j,2),cl_coord(j,3),sprintf('%i',idx));
            end
            cl_roi_labels(j) = idx;
        end
        view(90,0);
        shading(gca,'interp');
        axis off;
        
        % Plot legend
        if (trig_blot_nosig)
            n_legend = (K(k)) + 1;
        else
            n_legend = (K(k));
        end
        for il = 1:n_legend
            box_width = 0.016;
            box_vert_spacer = 0.16 * box_width;
            horiz_offset = 0;
            horiz_text_offset = -1 * (box_width) - 0.001;
            horiz_pos = pos{1}(1) + pos{1}(3) + horiz_offset;
            vert_offset =  - il * (box_width + box_vert_spacer);
            vert_pos = pos{1}(2) + pos{1}(4) + vert_offset;
            if (il == ((K(k)) + 1))
                legendStr = 'N.S.';
                legendCol = COLOR_NOSIG;
            else
                legendStr = sprintf('%i',il);
                legendCol = cc_k(il,:);
            end
            annotation('rectangle',[horiz_pos + horiz_text_offset, vert_pos, box_width, 0.9*box_width*(paper_width/paper_height)],...
                'FaceColor',legendCol,'EdgeColor','black','LineWidth',1);
            annotation('textbox',[horiz_pos,vert_pos-0.002, box_width, box_width*(paper_width/paper_height)],...
                'String',legendStr,'EdgeColor','none',...
                'FontSize',10,'FontName','Helvetica','VerticalAlignment','middle')
            %text(horiz_pos + horiz_text_offset,vert_pos,sprintf('%i',il));
        end
        
        %return
        %subplot(2,1,2);
        axes(ha(2));
        hold all;
        p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
            'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
        set(p,'FaceVertexCData',vcolor);
        brainlight;
        view(-90,0);
        shading(gca,'interp');
        axis off;
        
        axes(ha(3));
        hold all;
        p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
            'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
        set(p,'FaceVertexCData',vcolor);
        brainlight;
        view(90,-90);
        shading(gca,'interp');
        axis off;
        
        
        
        % ------------ tsne -----------------------------------------------
        axes(ha(4));
        dimout = 2;
        options = statset('MaxIter',10000,'TolFun',1e-14);
        % 'distance',@nanDist4clustering
        %return
        Y = tsne(A2km,'distance',@nanDist4clustering_neg_inf,'Options',options,'Algorithm','exact','Standardize',true,'NumDimensions',dimout,'Verbose',1);
        %Y = tsne(A2km,'Options',options,'Algorithm','exact','Standardize',true,'NumDimensions',dimout,'Verbose',1);
        %h = figure('visible','off');
        omit_tsne_subplot = true;
        if (~omit_tsne_subplot)
            if (dimout == 2)
                %cck = jet(k);
                %cc_p = corrcmap(K(k)+2);
                %cck = cc_p(2:(end-1),:);
                %cck = corrcmap(k);
                for i2 = 1:length(Y(:,1))
                    if ((trig_blot_nosig) && (cl_adj_nosig(i2)))
                        col_tsne = COLOR_NOSIG;
                    else
                        col_tsne = cc_k(Res{k,1}(i2),:);
                    end
                    plot(Y(i2,1),Y(i2,2),'.','MarkerSize',15,'Color',col_tsne);
                    hold on;
                end
            else
                plot3(Y(:,1),Y(:,2),Y(:,3),'black.');
            end
            xlabel('t-SNE first dimension');
            ylabel('t-SNE second dimension');
            set(gca,'Box','off')
            set(gca,'TickDir','out')
            ax = gca;
            XLim2 = ax.XLim;
            YLim2 = ax.YLim;
            XLim2 = sign(XLim2).*(ceil(abs(XLim2))); %round(ax.XLim);
            YLim2 = sign(YLim2).*(ceil(abs(YLim2))); %round(ax.YLim);
            try
                set(gca,'XTick',unique([XLim2(1),0,XLim2(2)]))
                set(gca,'YTick',unique([YLim2(1),0,YLim2(2)]))
                set(gca,'XLim',XLim2);
                set(gca,'YLim',YLim2);
            catch
            end

            set(gca,'FontName','Helvetica');
            set(gca,'FontSize',10);
        else
            axis off;
        end



        print(h,sprintf('figures/T17/kmeans-%i_%s_%s_nosig-%i',K(k),metric(3:end),SURFACE_TYPE,trig_blot_nosig),'-depsc');
        print(h,sprintf('figures/T17/kmeans-%i_%s_%s_nosig-%i',K(k),metric(3:end),SURFACE_TYPE,trig_blot_nosig),'-dpng','-r400');
        close(h);
        
        
        
        
        % Separately plot tsne eps
        h = figure('visible','off');
        set(h,'PaperUnits','inches');
        set(h,'PaperPosition',[0 0 paper_width/2 paper_height/2]);
        hold all;
        if (dimout == 2)
            %cck = jet(k);
            %cc_p = corrcmap(K(k)+2);
            %cck = cc_p(2:(end-1),:);
            %cck = corrcmap(k);
            for i2 = 1:length(Y(:,1))
                if ((trig_blot_nosig) && (cl_adj_nosig(i2)))
                    col_tsne = COLOR_NOSIG;
                else
                    col_tsne = cc_k(Res{k,1}(i2),:);
                end
                plot(Y(i2,1),Y(i2,2),'.','MarkerSize',15,'Color',col_tsne);
                hold on;
            end
        else
            plot3(Y(:,1),Y(:,2),Y(:,3),'black.');
        end
        xlabel('t-SNE first dimension');
        ylabel('t-SNE second dimension');
        set(gca,'Box','off')
        set(gca,'TickDir','out')
        ax = gca;
        XLim2 = ax.XLim;
        YLim2 = ax.YLim;
        XLim2 = sign(XLim2).*(ceil(abs(XLim2))); %round(ax.XLim);
        YLim2 = sign(YLim2).*(ceil(abs(YLim2))); %round(ax.YLim);
        try
            set(gca,'XTick',unique([XLim2(1),0,XLim2(2)]))
            set(gca,'YTick',unique([YLim2(1),0,YLim2(2)]))
            set(gca,'XLim',XLim2);
            set(gca,'YLim',YLim2);
        catch
        end
        set(gca,'FontName','Helvetica');
        set(gca,'FontSize',10);
        
        print(h,sprintf('figures/T17/kmeans-%i_%s_%s_nosig-%i_tsne',K(k),metric(3:end),SURFACE_TYPE,trig_blot_nosig),'-depsc','-r400');
        close(h);
    end
    
    %[~,sIdx] = sort(cl_roi_labels);
    
end