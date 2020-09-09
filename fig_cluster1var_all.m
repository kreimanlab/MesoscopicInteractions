close all;
clear;
rng('shuffle');

%/media/klab/internal/data/coreg/fsaverage_sym/label/all_surf_ielvis_sub18.label
dir_artLp = '/media/klab/internal/data/h5_notch20/art';
dir_corLp = '/media/klab/internal/data/coreg';
%dir_resLp = '/media/klab/KLAB101/results/coh_w10';
setenv('SUBJECTS_DIR',dir_corLp);
dir_cacheLp = './cache';
dir_h5Lp = '/media/klab/KLAB101/h5_notch20';

%metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
metricsp = {'pcBroadband'};

% Corresponding variable (11: CT, 12: Mag)
% variable_idx = 12;

%SurfT = {'pial','inflated','sphere','smoothwm'};
SurfT = {'pial'};

%Caa = load('cache/fig_cluster2_reduce_annot.mat');
Caa2 = load('cache/fig_cluster2_reduce_new_annot_vertex_color2.mat');
%vcolor = Caa2.vertex_color2_bw;
vcolor = Caa2.vertex_color2;
vcolor = vcolor * 0.6 + 0.6*[1 1 1]*0.4;
Caa = load('cache/fig_cluster2_reduce_annot.mat');
cc_cluster = Caa.cc;

% bypass face color
%vcolor = 0.7*ones(size(vcolor));
trig_color_mwall = true; % parfor


Caci = load('cache/fig_cluster3_cluster_i');
cluster_i = Caci.cluster_i; 

% Read medial wall
hemi = 'r';
l_mwall = read_label('fsaverage_sym',sprintf('%sh.Medial_wall',hemi));
vert_mwall = l_mwall(:,1) + 1;



% strip color
% col_const = 0.6*[1 1 1];
% trig_strip_col = true;
% if (trig_strip_col)
%     [n_v,~] = size(vcolor);
%     for i = 1:n_v
%         col = vcolor(i,:);
%         if (diff(diff(col)) ~= 0)
%             vcolor(i,:) = col_const;
%         end
%     end
% end


col_coeff = 0.2;
BRAIN_MESH_ALPHA = 1; % 0.3
vcolor = (1*[1 1 1] * col_coeff + vcolor * (1-col_coeff));


trig_plot_scale_cube = false;

for ist = 1:length(SurfT)

    % Load fsaverage_sym surface
    trig_plot_all = true;
    surf_type = SurfT{ist}; %'pial';
    [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_corLp,'fsaverage_sym','r',surf_type));
    
    % load spherical surface
    [s_vert_sph, faces_sph] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_corLp,'fsaverage_sym','r','sphere'));

    
    % Color medial wall
    if (trig_color_mwall)
        col_mwall = [1 1 1]*0.2;
        [n_vert,~] = size(s_vert);
        parfor i = 1:n_vert
            in_mwall = sum(vert_mwall == i) > 0;
            if (in_mwall)
                vcolor(i,:) = col_mwall;
            end
        end
    end
    
    
    % compute medial wall spherical center
    RADIUS_SPHERE = 100; % mm
    mwall_center = mean(s_vert_sph(vert_mwall,:));
    mwall_center = mwall_center / sqrt(sum(mwall_center.^2));
    mwall_center = RADIUS_SPHERE * mwall_center;
    
    trig_plot_mag_dist = false;

    trig_kmeans = false;

    for iM = 1:length(metricsp)
        fn_cache = [dir_cacheLp,'/xsub_out_all_',num2str(iM)];
        Ca = load(fn_cache);
        
        % Build fsaverage electrode list
        fprintf('[*] Building subject on fsaverage bipolar coord list..\n')
        Coord = [];
        for iI = 1:length(Ca.Subjects)
            sid = Ca.Subjects{iI};
            sid_int = str2double(sid(2:end));
            Cs = load(sprintf('%s/xsub_out_%s_%i',dir_cacheLp,sid,iM));
            l = read_label(sprintf('fsaverage_sym'),sprintf('ielvis_%s',sid));
            l(:,1) = l(:,1) + 1;
            for iJ = 1:Cs.ecog.n_bchan
                b1c1 = Cs.ecog.bip(iJ,1);
                lL = l((l(:,end) == b1c1),:);
                xyz = s_vert(lL(1),:);
                Coord = [Coord; [sid_int,iJ,xyz]];
            end
%             trisurf(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3),'EdgeColor','none','FaceColor',[1 1 1]*0.5); hold on;
%             plot3(Coord(:,3),Coord(:,4),Coord(:,5),'black.'); daspect([1 1 1]);
%             return
        end
        


        % build master coherence table
        D = [];
        E = [];
        
        CaR = load(sprintf('%s/fig_cluster2_reduce_%i_new.mat',dir_cacheLp,iM));
        E = CaR.E;
        [n_E,~] = size(E);
        
%         for iE = 1:(n_E-1)
%             for jE = (iE+1):n_E
%                 mag = CaR.A(iE,jE);
%                 dis = CaR.Ad(iE,jE);
%                 coord1 = E(iE,3:5);
%                 coord2 = E(jE,3:5);
%                 if ((~isnan(mag)) && (mag ~= 0))
%                     D = [D; [E(iE,1), E(iE,2), E(jE,2), dis, coord1, coord2, mag, mag]];
%                 end
%             end
%         end
        
        % Compute spatial variances of nodes
        Es = CaR.Es;
        [n_Es,~] = size(Es);
        Es_var = zeros(n_Es,4);
        Es_var_max = zeros(n_Es,4);
        for iv = 1:n_Es
            Es_var(iv,1:3) = var(Es{iv}(:,3:5));
            Es_var_max(iv,1:3) = var(Es{iv}(:,3:5));
            %Es_var(iv,4) = (det(cov(Es{iv}(:,3:5))))^(1/3);
            x = Es{iv}(:,3:5);
            [n_x,~] = size(x);
            xA = []; %zeros(1,nchoosek(n_x,2));
            xAc = 1;
            for ix = 1:(n_x-1)
                for ix2 = (ix+1):n_x
                    xd = sqrt(sum((x(ix,:)-x(ix2,:)).^2));
                    xA(xAc) = xd;
                    xAc = xAc + 1;
                end
            end
            if (isempty(xA))
                Es_var(iv,4) = 0;
                Es_var_max(iv,4) = 0;
            else
                Es_var(iv,4) = nanmean(xA(:));
                Es_var_max(iv,4) = nanmax(xA(:));
            end

            % hist(sqrt(Es_var(:,4)));
        end
        Es_stdev = Es_var(:,4); %sqrt(real(Es_var(:,4)));
        Es_stdev_max = Es_var_max(:,4);
        
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
        
        sidintu = unique(E(:,1));
        n_sidintu = length(unique(E(:,1)));

    %     for i1 = 1:(Ca.n_rois-1)
    %         for i2 = (i1+1):Ca.n_rois
    %             AA = Ca.AdjAtl{i1,i2};
    %             AAs = Ca.AdjAtl_sid{i1,i2};
    %             for j = 1:length(AA)
    %                 % only consider significant interactions
    %                 if (AA(j) > 0)
    %                 endDisplay image as texture-mapped surface
    %             end
    %         end
    %     end






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
            set(hfig,'PaperSize',3*[7.111 4]);
            set(hfig,'PaperUnits','inches');

            %hax1 = subplot(2,2,1,'Position',[0 0.5 0.5 0.5]);
            hax1 = axes;
            
            % plot surface
            SMOOTH_FACE = true;
            
            p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
                'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
            
            % Colorful regional surface
            set(p,'FaceVertexCData',vcolor);
            shading('interp');
            
            % Simple surface
            %set(p,'FaceVertexCdata',0.7*[1 1 1]);
            
            ax = gca;
            ax.Clipping = 'off';
            hold all;

            % colormap params
            n_cmap = 100;
            %cc = jet(n_cmap);
            %cc = corrcmap(n_cmap);
            %cc = viridis(n_cmap);
            cc = 0*ones(n_cmap,3);
%             col_var_min = min(D(:,variable_idx)); %0.1;
%             col_var_max = max(D(:,variable_idx)); %1;
%             % Plot lines
%             for i = 1:length(D)
%                 % Plot first electrode
%                 %plot3(D(i,5),D(i,6),D(i,7),'black.')
%                 % plot second electrode
%                 %plot3(D(i,8),D(i,9),D(i,10),'black.')
%                 % plot interaction
%                 % - magnitude -------------------------------------------------
%                 col_var = D(i,variable_idx);
%                 %col_line = cc(round((col_var/col_var_max)*n_cmap),:);
%                 col_line = cc(round(( (col_var - col_var_min)/(col_var_max - col_var_min) ) *(n_cmap-1)+1),:);
% 
% 
%                 if (trig_kmeans)
%                     % find E index for this D
%                     cc_kmeans = jet(k);
% 
%                     d_sidint = D(i,1);
%                     d_e = D(i,2);
%                     e_idx = find( (E(:,1) == d_sidint) & (E(:,2) == d_e) );
% 
%                     % find cluster id for this line
%                     plot3([D(i,5),D(i,8)],[D(i,6),D(i,9)],[D(i,7),D(i,10)],'Color',cc_kmeans(iidx(e_idx),:));
%                 else
%                     plot3([D(i,5),D(i,8)],[D(i,6),D(i,9)],[D(i,7),D(i,10)],'Color',col_line);
%                 end
%         %         if (i == 8000)
%         %             break
%         %         end
%             end

            %plot electrodes
            % Es_stdev - 1 standard deviation
            elec_size = 30;
            col_var_min = min(Es_stdev);
            col_var_max = max(Es_stdev);
            fprintf('Mean distance (mm) between electrodes in region mean: %.4f\n',mean(Es_stdev));
            fprintf('Mean distance (mm) between electrodes in region median: %.4f\n',median(Es_stdev));
            fprintf('Mean distance (mm) between electrodes in region std: %.4f\n',std(Es_stdev));
            fprintf('Mean distance (mm) between electrodes in region min: %.4f\n',min(Es_stdev));
            fprintf('Mean distance (mm) between electrodes in region max: %.4f\n',max(Es_stdev));
            fprintf('Max distance (mm) between electrodes in region mean: %.4f\n',mean(Es_stdev_max));
            fprintf('Max distance (mm) between electrodes in region median: %.4f\n',median(Es_stdev_max));
            fprintf('Max distance (mm) between electrodes in region std: %.4f\n',std(Es_stdev_max));
            fprintf('Max distance (mm) between electrodes in region min: %.4f\n',min(Es_stdev_max));
            fprintf('Max distance (mm) between electrodes in region max: %.4f\n',max(Es_stdev_max));
            
            % electrode mesh
            r = 0.5; n = 8;   
            [x,y,z] = sphere(n);
            hold all;
            % Main representative electrode color (of 150)
            %col_elec_loc = [0.8 0 0]; 
            col_elec_loc = [1 1 1]*0;
            
            
            for i = 1:length(E)
                
                fprintf('[*] Start %i of %i, region %i\n',i,length(E),cluster_i(i));
                
                % =====================================================
                % Medial wall remap for main electrode
                ei_coord = E(i,3:5);
                d2 = sum((s_vert - ei_coord).^2,2);
                [~,ei_idx] = min(d2);
                % Check if vertex in medial wall
                in_mwall = sum(vert_mwall == ei_idx) > 0;
                if (in_mwall)
                    % center of sphere at mwall_center
                    start = s_vert_sph(ei_idx,:);
                    direction = start - mwall_center;
                    direction = direction / sqrt(sum(direction.^2));
                    step = 1;
                    sat = true;
                    curr = start;
                    c = 1;
                    while (sat)
                        curr = start + direction*(step*c);
                        d2 = sum((s_vert_sph - curr).^2,2);
                        [~,finish_i] = min(d2);
                        in_mwall = sum(vert_mwall == finish_i) > 0;
                        sat = in_mwall; % run as long as finish point is in medial wall
                        c = c + 1;
                    end
                    Ei_new = s_vert(finish_i,:);
                    fprintf('\t[!] medial wall remap: [%.1f,%.1f,%.1f] -> [%.1f,%.1f,%.1f], %.1f mm\n',...
                        E(i,3),E(i,4),E(i,5),Ei_new(1),Ei_new(2),Ei_new(3),c-1);
                    E(i,3:5) = Ei_new;
                end
                % =====================================================

                col_var = Es_stdev(i);
                col_elec = cc(round(( (col_var - col_var_min)/(col_var_max - col_var_min) ) *(n_cmap-1)+1),:);
                %plot3(E(i,3),E(i,4),E(i,5),'o','Color',col_elec,'MarkerSize',elec_size*0.2); hold on;

                %fac = 2;
                fac = 1;
                X = x*r*fac + E(i,3);
                Y = y*r*fac + E(i,4);
                Z = z*r*fac + E(i,5);
                
                % Color main electrode
                %col_elec_loc = cc_cluster(i,:);
                
                q = surf(X,Y,Z,'facecolor',col_elec_loc,'LineStyle','none');
                q.SpecularStrength = 0;
                q.SpecularExponent = 15;
                q.SpecularColorReflectance = 0;

                % plot member bipolar electrodes
                blist = CaR.Es{i};
                [n_blist,~] = size(blist);
                for ib = 1:n_blist
                    
                    
                    


                    sid_int = blist(ib,1);
                    b1 = blist(ib,2);
                    sIdx = (Coord(:,1) == sid_int) & (Coord(:,2) == b1);
                    xyz = Coord(sIdx,3:5);
                    
                    
                    
                    
                    
                    
                    % =====================================================
                    % Medial wall remap for member electrodes
                    ei_coord = xyz; %blist(ib,3:5);
                    d2 = sum((s_vert - ei_coord).^2,2);
                    [~,ei_idx] = min(d2);
                    % Check if vertex in medial wall
                    in_mwall = sum(vert_mwall == ei_idx) > 0;
                    if (in_mwall)
                        % center of sphere at mwall_center
                        start = s_vert_sph(ei_idx,:);
                        direction = start - mwall_center;
                        direction = direction / sqrt(sum(direction.^2));
                        step = 1;
                        sat = true;
                        curr = start;
                        c = 1;
                        while (sat)
                            curr = start + direction*(step*c);
                            d2 = sum((s_vert_sph - curr).^2,2);
                            [~,finish_i] = min(d2);
                            in_mwall = sum(vert_mwall == finish_i) > 0;
                            sat = in_mwall; % run as long as finish point is in medial wall
                            c = c + 1;
                        end
                        Ei_new = s_vert(finish_i,:);
                        %fprintf('\t\t(*) child remap: [%.1f,%.1f,%.1f] -> [%.1f,%.1f,%.1f], %.1f mm\n',...
                        %    blist(ib,3),blist(ib,4),blist(ib,5),Ei_new(1),Ei_new(2),Ei_new(3),c-1);
                        
                        fprintf('\t\t(*) child remap: [%.1f,%.1f,%.1f] -> [%.1f,%.1f,%.1f], %.1f mm\n',...
                            xyz(1),xyz(2),blist(3),Ei_new(1),Ei_new(2),Ei_new(3),c-1);
                        
                        xyz = Ei_new;
                        %blist(i,3:5) = Ei_new;
                    end
                    % =====================================================
                    
                    
                    
                    
                    X = x*r + xyz(1);
                    Y = y*r + xyz(2);
                    Z = z*r + xyz(3);
                    q = surf(X,Y,Z,'facecolor',col_elec,'LineStyle','none');
                    q.SpecularStrength = 0;
                    q.SpecularExponent = 15;
                    q.SpecularColorReflectance = 0;
                    %plot3(xyz(1),xyz(2),xyz(3),'.','Color',col_elec,'MarkerSize',0.4*elec_size); hold on;
                end

                

                
                
            end

            axis tight;
            daspect([1 1 1]);
            view(90,0);
            if (SMOOTH_FACE)
                p.FaceLighting = 'gouraud';
            end
            
            % lighting
            p.AmbientStrength = 0.5 ; %0.3
            p.DiffuseStrength = 0.4 ;
            p.SpecularStrength = 0;
            p.SpecularExponent = 1;
            p.BackFaceLighting = 'lit';
            p.FaceLighting = 'gouraud';
            cam_elev = 0;
            camlight(-135,cam_elev);
            camlight(45,cam_elev);
            camlight(-225,cam_elev);
            camlight(-45,cam_elev);
            
            
            % colorbar
            cb = colorbar;
            if (trig_kmeans)
                colormap(cc_kmeans)
                caxis([1 k])
            else
                colormap(cc)
                caxis([col_var_min col_var_max])
                set(cb,'YTick',[col_var_min,col_var_max]);
                set(cb,'YTickLabel',{sprintf('%.1f',col_var_min),sprintf('%.1f',col_var_max)});
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
            if (trig_plot_scale_cube)
                c = min(E(:,3:5)) + 10*[0 1 1];
                w = 10;
                cx = [c(1),c(1)  ,c(1)  ,c(1)  ,c(1),c(1),c(1),c(1)];
                cy = [c(2),c(2)+w,c(2)+w,c(2)  ,c(2),c(2),c(2),c(2)];
                cz = [c(3),c(3)  ,c(3)+w,c(3)+w,c(3),c(3),c(3),c(3)];
                text(c(1),c(2)+0.5*w,c(3)+0.5*w,sprintf('%i mm',w),'HorizontalAlignment','center','FontSize',6)
                plot3(cx,cy,cz,'black-')
            end
            
            
            %margin
            m = 0.05;
            
            % Lateral view
            fprintf('[*] Making lateral view..\n');
            set(hax1,'Position',[0+m/2, 0.5+m/2, 0.5-m, 0.5-m]);
            
            
            % Medial view
            %ax_clone = copyobj(hax_o, hfig);
            fprintf('[*] Making medial view..\n');
            set(hax_o,'Position',[0+m/2, 0+m/2, 0.5-m, 0.5-m]);
            set(hax_o,'View',[-90,0]);
            %axes(hax_o);
            %view(-90,0);
            
%             % Superior view
%             fprintf('[*] Making superior view..\n');
%             ax_clone2 = copyobj(hax_o, hfig);
%             set(ax_clone2,'Position',[0.5+m/2, 0.25+m/2, 0.5-m, 0.25-m]);
%             set(ax_clone2,'View',[-90,90]);
%             %axes(ax_clone2);
%             %view(-90,90);

            % Inferior view
            fprintf('[*] Making inferior view..\n');
            ax_clone3 = copyobj(hax_o, hfig);
            %hax4 = subplot(2,2,4,ax_clone3,'Position',[0.5 0 0.5 0.25]);
            set(ax_clone3,'Position',[0.5+m/2, 0+m/2, 0.5-m, 0.25-m]);
            set(ax_clone3,'View',[-90,-90]);
            %axes(ax_clone3);
            %view(-90,-90);
            
%             % Posterior view
%             fprintf('[*] Making posterior view..\n');
%             ax_clone4 = copyobj(hax_o, hfig);
%             %hax2 = subplot(2,2,2,ax_clone4,'Position',[0.5 0.5 0.25 0.5]);
%             set(ax_clone4,'Position',[0.70+m/2, 0.5+m/2, 0.25-m, 0.5-m]);
%             set(ax_clone4,'View',[0,0]);
%             %axes(ax_clone4);
%             %view(0,0);
            
            
%             % Anterior view
%             fprintf('[*] Making anterior view..\n');
%             ax_clone5 = copyobj(hax_o, hfig);
%             %hax2 = subplot(2,2,2,ax_clone4,'Position',[0.5 0.5 0.25 0.5]);
%             set(ax_clone5,'Position',[0.55+m/2, 0.5+m/2, 0.25-m, 0.5-m]);
%             set(ax_clone5,'View',[180,0]);
%             %axes(ax_clone5);
%             %view(180,0);
            

            fprintf('[!] Saving..\n');
           % return
            print(hfig,sprintf('figures/cluster1/%s_nbip-%i_nedge-%i_1stdev',surf_type,length(E),length(D)),'-dpng','-r400');
            close(hfig);
        end



        % Show numbers
        fprintf('Number of electrodes: %i\n',length(E))
        fprintf('Number of significant interactions: %i\n',length(D))

    end

end

% Mean distance (mm) between electrodes in region mean: 7.7043
% Mean distance (mm) between electrodes in region median: 7.4710
% Mean distance (mm) between electrodes in region std: 1.8626
% Mean distance (mm) between electrodes in region min: 3.4138
% Mean distance (mm) between electrodes in region max: 12.9894
% Max distance (mm) between electrodes in region mean: 17.2329
% Max distance (mm) between electrodes in region median: 16.9921
% Max distance (mm) between electrodes in region std: 4.3361
% Max distance (mm) between electrodes in region min: 8.8188
% Max distance (mm) between electrodes in region max: 27.6423