close all;
clear;
rng('shuffle');

%/media/klab/internal/data/coreg/fsaverage_sym/label/all_surf_ielvis_sub18.label
dir_artLp = '/media/klab/internal/data/h5_notch20/art';
dir_corLp = '/media/klab/internal/data/coreg';
dir_resLp = '/media/klab/KLAB101/results/coh_w10';
setenv('SUBJECTS_DIR',dir_corLp);
dir_cacheLp = './cache';
dir_h5Lp = '/media/klab/KLAB101/h5_notch20';


%metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
metricsp = {'pcBroadband'};

tube_radius = 0.15;
tube_n_edges = 3;

% Corresponding variable (11: CT, 12: Mag)
variable_idx = 12;
        
%SurfT = {'pial','inflated','sphere','smoothwm'};
SurfT = {'pial'};

for ist = 1:length(SurfT)

    % Load fsaverage_sym surface
    trig_plot_all = true;
    surf_type = SurfT{ist}; %'pial';
    [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_corLp,'fsaverage_sym','r',surf_type));

    trig_plot_mag_dist = false;
    trig_plot_arc = true;
    trig_kmeans = false;

    for iM = 1:length(metricsp)
        fn_cache = [dir_cacheLp,'/xsub_out_all_',num2str(iM)];
        Ca = load(fn_cache);


        % build master coherence table
        D = [];
        Dpath = {};
        E = [];
        
        CaR = load(sprintf('%s/fig_cluster2_reduce.mat',dir_cacheLp));
        %CaR = load(sprintf('%s/fig_cluster2_reduce_clean.mat',dir_cacheLp));
        E = CaR.E;
        [n_E,~] = size(E);
        
        Cb = load('brainexport/red_6all_fsaverage.mat');
        
        count = 1;
        for iE = 1:(n_E-1)
            for jE = (iE+1):n_E
                mag = CaR.A(iE,jE);
                dis = CaR.Ad(iE,jE);
                coord1 = E(iE,3:5);
                coord2 = E(jE,3:5);
                if ((~isnan(mag)) && (mag ~= 0))
                    D = [D; [E(iE,1), E(iE,2), E(jE,2), dis, coord1, coord2, mag, mag]];
                    Dpath = [Dpath; {Cb.Paths{count}}];
                end
                count = count + 1;
            end
        end
        
        % Compute spatial variances of nodes
        Es = CaR.Es;
        [n_Es,~] = size(CaR.Es);
        Es_var = zeros(n_Es,4);
        for iv = 1:n_Es
            Es_var(iv,1:3) = var(Es{iv}(:,3:5));
            Es_var(iv,4) = (det(cov(Es{iv}(:,3:5))))^(1/3);
            % hist(sqrt(Es_var(:,4)));
        end
        
        
%         for iS = 1:length(Ca.Subjects)
%             % read subject cache
%             sid = Ca.Subjects{iS};
%             sidint = str2num(sid(2:end));
%             fprintf('[%s]\n',sid)
%             Cas = load(sprintf('cache/xsub_out_%s_%i.mat',sid,iM));
%             La = read_label('fsaverage_sym',sprintf('ielvis_%s',sid));
% 
% 
%             if isempty(La)
%                 fprintf(2,'[*] Skip %s, no fsaverage_sym label\n',sid)
%             else
% 
%                 % Get all electrode coordinates
%                 for j = 1:Cas.ecog.n_bchan
%                     b1c1 = Cas.ecog.bip(j,1);
%                     % get fsaverage_sym coordinates of electrodes
%                     %e_coords = s_vert(La(:,1) + 1,:);
%                     lidx = La(:,end) == b1c1;
%                     %coord1 = La(lidx,2:4);
%                     coord1 = s_vert(La(lidx,1) + 1,:);
%                     E = [E; [sidint, j, coord1]];
%                 end
% 
%                 % Get significant interaction coordinates
%                 for j = 1:length(Cas.ct)
%                     % Consider bip pairs in subject that are significant and
%                     % far enough apart
%                     if ((Cas.ct(j) > Cas.ct_thresh) && (Cas.Dmats(j) > Cas.dist_thresh))
%                         b1 = Cas.chan1(j) + 1;
%                         b2 = Cas.chan2(j) + 1;
%                         % get electrode numbers of first in bipolar pair
%                         b1c1 = Cas.ecog.bip(b1,1);
%                         b2c1 = Cas.ecog.bip(b2,1);
%                         % get fsaverage_sym coordinates of electrodes
%                         %e_coords = s_vert(La(:,1) + 1,:);
%                         lidx = La(:,end) == b1c1;
%                         %coord1 = La(lidx,2:4);
%                         coord1 = s_vert(La(lidx,1) + 1,:);
% 
%                         lidx = La(:,end) == b2c1;
%                         %coord2 = La(lidx,2:4);
%                         coord2 = s_vert(La(lidx,1) + 1,:);
% 
% 
% 
%                         if (isempty(coord1) || isempty(coord2))
%                             fprintf(2,'[!] coordinates empty for %s.\n',sid)
%                         else
%                             % --------------------------------------------------------------------------------
%                             D = [D; [sidint, b1, b2, Cas.Dmats(j), coord1, coord2, Cas.ct(j), Cas.mag(j)]];
%                         end
% 
%     %                     if (iS == 6)
%     %                         n_chan = max(max(Cas.ecog.bip(:,1:2)));
%     %                         [n_chan_l,~] = size(La);
%     %                         fprintf('n_chan: %i\nn_chan in label: %i \n',n_chan,n_chan_l)
%     %                         plot3(La(:,2),La(:,3),La(:,4),'blacko')
%     %                         return;
%     %                     end
%                     end
%                 end
%             end
% 
%         end
        
        sidintu = unique(D(:,1));
        n_sidintu = length(unique(D(:,1)));

    %     for i1 = 1:(Ca.n_rois-1)
    %         for i2 = (i1+1):Ca.n_rois
    %             AA = Ca.AdjAtl{i1,i2};
    %             AAs = Ca.AdjAtl_sid{i1,i2};
    %             for j = 1:length(AA)
    %                 % only consider significant interactions
    %                 if (AA(j) > 0)
    %                 end
    %             end
    %         end
    %     end



        if (trig_plot_mag_dist)
            h2 = figure;

            cc_sub = corrcmap(n_sidintu);
            
            
            subplot(2,1,1)
            % Plot magnitude distance for individual subjects
            %cc_sub = corrcmap(n_sidintu);
            for i = 1:n_sidintu
                Ds = D(D(:,1) == sidintu(i),:);
                dists = Ds(:,4);
                vars = Ds(:,variable_idx);
                plot(dists,vars,'.','color',cc_sub(i,:));
                hold on;
                % fit trendline
                p = polyfit(dists,vars,1);
                x = [Cas.dist_thresh max(dists)];
                plot(x,p(2) + x.*p(1),'-','color',cc_sub(i,:));
                hold on;
            end
            lz = cell(n_sidintu,1);
            for i = 1:length(lz)
                lz{i} = sprintf('m%i',sidintu(i));
            end
            %legend(lz)
            axis([Cas.dist_thresh 120 0 1]);
            set(gca,'TickDir','Out');
            set(gca,'Box','Off');
            xlabel('Distance (mm)');
            ylabel('Coherence');


            subplot(2,1,2)
            variable_idx = 11;
            % Plot magnitude distance for individual subjects
            %cc_sub = jet(n_sidintu);
            for i = 1:n_sidintu
                Ds = D(D(:,1) == sidintu(i),:);
                dists = Ds(:,4);
                vars = Ds(:,variable_idx);
                plot(dists,vars,'.','color',cc_sub(i,:));
                hold on;
                % fit trendline
                p = polyfit(dists,vars,1);
                x = [Cas.dist_thresh max(dists)];
                plot(x,p(2) + x.*p(1),'-','color',cc_sub(i,:));
                hold on;
            end
            lz = cell(n_sidintu,1);
            for i = 1:length(lz)
                lz{i} = sprintf('m%i',sidintu(i));
            end
            %legend(lz)
            axis([Cas.dist_thresh 120 0 1]);
            set(gca,'TickDir','Out');
            set(gca,'Box','Off');
            xlabel('Distance (mm)');
            ylabel('Consistency across time');
            
            print(h2,sprintf('figures/cluster1/mag_dist_nsub-%i',n_sidintu),'-depsc');
            close(h2);
        end






%         % K-means
%         if (trig_kmeans)
%             [n_E,~] = size(E);
%             [n_D,~] = size(D);
% 
%             K = [4:2:33];
%             n_inits = 20;
%             Iidx = cell(1,length(K));
%             Iidx_mindist = cell(1,length(K));
%             Sum_mindist = zeros(1,length(K));
% 
%             parfor i = 1:length(K)
%                 k = K(i);
% 
%                 var_init_prev = Inf;
%                 var_init_final = NaN;
%                 iidx_final = [];
%                 iidx_mindist_final = [];
%                 % init
%                 for j = 1:n_inits
%                     iidx = zeros(n_E,1);
%                     iidx(1:k) = 1:k;
%                     iidx = iidx(randperm(n_E));
%                     inits = find(iidx);
% 
%                     % Extras
%                     iidx_mindist = zeros(n_E,1);
% 
%                     % Coords of each init electrode
%                     coord_init = nan(k,3);
%                     for l = 1:k
%                         e = E(inits(l),:);close all
%                         e_sidint = e(1);
%                         e_b = e(2);
%                         % Find all interactions for this subject-bipolar set
%     %                     ee = find((D(:,1) == e_sidint) & ((D(:,2)==e_b) | (D(:,3)==e_b)));
%     %                     coords_all = D(ee,1:7);
%                         ee1 = find((D(:,1) == e_sidint) & ((D(:,3)==e_b)));
%                         ee2 = find((D(:,1) == e_sidint) & ((D(:,2)==e_b)));
%                         coords_all = [D(ee1,5:7); D(ee2,8:10)];
%                         % collapse all electrode pairs to single average
%                         if (isempty(coords_all))
%                             coord_init(l,:) = nan(1,3);
%                         else
%                             coord_init(l,:) = mean(coords_all,1);
%                         end
%                         %disp(coords_all)
%                     end
%                     %disp(coord_init)
%                     %return
% 
% 
%                     % first cluster loop
%                     % assign all points to one of k clusters
%                     for k2 = 1:n_E
%                         if (iidx(k2) == 0) % skip clustering inits
%                             % calculate distance to each init cluster
%                             cdist = nan(k,1);
% 
%                             coord_q = nan(1,3);
% 
%                             e = E(k2,:);
%                             e_sidint = e(1);
%                             e_b = e(2);
%                             % Find all interactions for this subject-bipolar set
%         %                     ee = find((D(:,1) == e_sidint) & ((D(:,2)==e_b) | (D(:,3)==e_b)));
%         %                     coords_all = D(ee,1:7);
%                             ee1 = find((D(:,1) == e_sidint) & ((D(:,3)==e_b)));
%                             ee2 = find((D(:,1) == e_sidint) & ((D(:,2)==e_b)));
%                             coords_all = [D(ee1,5:7); D(ee2,8:10)];
%                             % collapse all electrode pairs to single average
%                             if (isempty(coords_all))
%                                 coord_q(1,:) = nan(1,3);
%                             else
%                                 coord_q(1,:) = mean(coords_all,1);
%                             end
% 
%                            % cdist(l) = cluster_dist;
% 
%                             % compute distances
%                             for l = 1:k
%                                 if (isnan(coord_q(1)) && isnan(coord_init(l,1)))
%                                     % Both cluster and query have no interactions
%                                     cdist(l) = 0;
%                                 elseif (~isnan(coord_q(1)) && ~isnan(coord_init(l,1)))
%                                     % Both cluster and query have interactions
%                                     % return euclidean distances between
%                                     % geometric centers of interacting
%                                     % electrode set
%                                     cdist(l) = sqrt(sum((coord_q - coord_init(l,:)).^2));
%                                 else
%                                     % Only one has no interactions
%                                     cdist(l) = Inf;
%                                 end
%                             end
% 
%                             % Assign minimum distance cluster
%                             [mindist,cluster_id] = min(cdist);
%                             iidx(k2) = cluster_id;
%                             iidx_mindist(k2) = mindist;
% 
% 
%                         end
%                     end
% 
%                     % Calculate goodness of init
%                     var_init = 0;
%                     for k3 = 1:k
%                         var_init = var_init + mean(iidx_mindist(iidx==k3));
%                     end
%                     var_init = var_init / k;
%                     %var_init = sum(iidx_mindist);
% 
% 
%                     if (var_init < var_init_prev)
%                         var_init_final = var_init;
%                         var_init_prev = var_init;
% 
%                         iidx_final = iidx;
%                         iidx_mindist_final = iidx_mindist;
% 
%                         fprintf('[*] k = %i init %i of %i, log sum mindist: %.6f\n',k,j,n_inits,log(var_init_final));
%                     end
% 
%                 end % init loop
% 
%                 Iidx{i} = iidx_final;
%                 Iidx_mindist{i} = iidx_mindist_final;
%                 Sum_mindist(i) = var_init_final;
% 
%             end % large K loop
% 
%             % Plot variance against k
%             h = figure;
%             plot(K,Sum_mindist,'black-');
%             xlabel('k')
%             ylabel('Average cluster variance (mm)')
% 
%             [~,mIdx] = min(Sum_mindist);
% 
%             iidx = Iidx{mIdx};
%             iidx_mindist = Iidx_mindist{mIdx};
%             k = K(mIdx);
% 
%         end



        if (trig_plot_all)
            hfig = figure('Visible','off');
            set(hfig,'Position',1*[0 0 1920 1080])
            set(hfig,'PaperSize',[7.111 4]);
            set(hfig,'PaperUnits','inches');

            %hax1 = subplot(2,2,1,'Position',[0 0.5 0.5 0.5]);
            hax1 = axes;
            
            % plot surface
            SMOOTH_FACE = true;
            BRAIN_MESH_ALPHA = 0.1;
            p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
            set(p,'FaceVertexCdata',0.7*[1 1 1]);
            ax = gca;
            ax.Clipping = 'off';
            hold all;

            % colormap params
            n_cmap = 100;
            %cc = jet(n_cmap);
            cc = corrcmap(n_cmap);

            col_var_min = min(D(:,variable_idx)); %0.1;
            col_var_max = max(D(:,variable_idx)); %1;
            % Plot lines
            for i = 1:length(D)
                % Plot first electrode
                %plot3(D(i,5),D(i,6),D(i,7),'black.')
                % plot second electrode
                %plot3(D(i,8),D(i,9),D(i,10),'black.')
                % plot interaction
                % - magnitude -------------------------------------------------
                col_var = D(i,variable_idx);
                %col_line = cc(round((col_var/col_var_max)*n_cmap),:);
                col_line = cc(round(( (col_var - col_var_min)/(col_var_max - col_var_min) ) *(n_cmap-1)+1),:);


                path = Dpath{i};
                
                if (trig_kmeans)
                    % find E index for this D
                    %cc_kmeans = jet(k);
                    k = 12;
                    cc_kmeans = corrcmap(k);
                    

                    d_sidint = D(i,1);
                    d_e = D(i,2);
                    e_idx = find( (E(:,1) == d_sidint) & (E(:,2) == d_e) );

                    % find cluster id for this line
                    col_tube = cc_kmeans(iidx(e_idx),:);
                    %plot3([D(i,5),D(i,8)],[D(i,6),D(i,9)],[D(i,7),D(i,10)],'Color',cc_kmeans(iidx(e_idx),:));
                else
                    col_tube = col_line;
                    %plot3([D(i,5),D(i,8)],[D(i,6),D(i,9)],[D(i,7),D(i,10)],'Color',col_line);
                end
                
                if (trig_plot_arc)
                    path = downsample(path,1);
                    [n_path,~] = size(path);
                    [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
                    surf(X,Y,Z,'FaceColor',col_tube,'EdgeColor','none'); hold on;
                else
                    plot3([D(i,5),D(i,8)],[D(i,6),D(i,9)],[D(i,7),D(i,10)],'Color',col_tube);
                end
               
            end

            %plot electrodes
            for i = 1:length(E)
                if (trig_kmeans)
                    plot3(E(i,3),E(i,4),E(i,5),'.','Color',cc_kmeans(iidx(i),:));
                else
                    plot3(E(i,3),E(i,4),E(i,5),'black.');
                end
            end

            axis tight;
            daspect([1 1 1]);
            view(90,0);
            if (SMOOTH_FACE)
                p.FaceLighting = 'gouraud';
            end
            cb = colorbar;
            if (trig_kmeans)
                colormap(cc_kmeans)
                caxis([1 k])
            else
                colormap(cc)
                caxis([col_var_min col_var_max])
                set(cb,'YTick',[col_var_min,col_var_max]);
                set(cb,'TickLength',0);
            end
            
            
            % Copy ax before scale
            set(hax1,'Visible','off');
            hax_o = copyobj(hax1,hfig);
%             ax_clone = copyobj(hax1, hfig);
%             ax_clone2 = copyobj(hax1, hfig);
%             ax_clone3 = copyobj(hax1, hfig);
%             ax_clone4 = copyobj(hax1, hfig);
%             ax_clone5 = copyobj(hax1, hfig);
%             ax_clone6 = copyobj(hax1, hfig);
            
            % plot scale cube
            c = min(E(:,3:5)) + 10*[0 1 1];
            w = 10;
            cx = [c(1),c(1)  ,c(1)  ,c(1)  ,c(1),c(1),c(1),c(1)];
            cy = [c(2),c(2)+w,c(2)+w,c(2)  ,c(2),c(2),c(2),c(2)];
            cz = [c(3),c(3)  ,c(3)+w,c(3)+w,c(3),c(3),c(3),c(3)];
            text(c(1),c(2)+0.5*w,c(3)+0.5*w,sprintf('%i mm',w),'HorizontalAlignment','center','FontSize',6)
            plot3(cx,cy,cz,'black-')
            
            
            %margin
            m = 0.05;
            
            % Lateral view
            set(hax1,'Position',[0+m/2, 0.5+m/2, 0.5-m, 0.5-m]);
            
            
            % Medial view
            %ax_clone = copyobj(hax_o, hfig);
            set(hax_o,'Position',[0+m/2, 0+m/2, 0.5-m, 0.5-m]);
            set(hax_o,'View',[-90,0]);
            %axes(hax_o);
            %view(-90,0);
            
            % Superior view
            ax_clone2 = copyobj(hax_o, hfig);
            set(ax_clone2,'Position',[0.5+m/2, 0.25+m/2, 0.5-m, 0.25-m]);
            set(ax_clone2,'View',[-90,90]);
            %axes(ax_clone2);
            %view(-90,90);

            % Inferior view
            ax_clone3 = copyobj(hax_o, hfig);
            %hax4 = subplot(2,2,4,ax_clone3,'Position',[0.5 0 0.5 0.25]);
            set(ax_clone3,'Position',[0.5+m/2, 0+m/2, 0.5-m, 0.25-m]);
            set(ax_clone3,'View',[-90,-90]);
            %axes(ax_clone3);
            %view(-90,-90);
            
            % Posterior view
            ax_clone4 = copyobj(hax_o, hfig);
            %hax2 = subplot(2,2,2,ax_clone4,'Position',[0.5 0.5 0.25 0.5]);
            set(ax_clone4,'Position',[0.70+m/2, 0.5+m/2, 0.25-m, 0.5-m]);
            set(ax_clone4,'View',[0,0]);
            %axes(ax_clone4);
            %view(0,0);
            
            
            % Anterior view
            ax_clone5 = copyobj(hax_o, hfig);
            %hax2 = subplot(2,2,2,ax_clone4,'Position',[0.5 0.5 0.25 0.5]);
            set(ax_clone5,'Position',[0.55+m/2, 0.5+m/2, 0.25-m, 0.5-m]);
            set(ax_clone5,'View',[180,0]);
            %axes(ax_clone5);
            %view(180,0);
            

            % return
            if (trig_kmeans)
                surf_type = sprintf('%s_k-%i',surf_type,k);
            end
            print(hfig,sprintf('figures/cluster1/%s_nbip-%i_nedge-%i_arc',surf_type,length(E),length(D)),'-dpng','-r300');
            close(hfig);
        end



        % Show numbers
        fprintf('Number of electrodes: %i\n',length(E))
        fprintf('Number of significant interactions: %i\n',length(D))

    end

end