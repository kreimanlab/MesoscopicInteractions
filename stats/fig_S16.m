close all; clear; rng('shuffle');

dir_artLp = '/media/klab/internal/data/h5_notch20/art';
dir_corLp = '/media/klab/internal/data/coreg';
setenv('SUBJECTS_DIR',dir_corLp);
dir_cacheLp = './cache';
dir_h5Lp = '/media/klab/KLAB101/h5_notch20';

%metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
metricsp = {'pcBroadband'};

SurfT = {'pial'};

Caa2 = load('cache/fig_cluster2_reduce_new_annot_vertex_color2.mat');
%vcolor = Caa2.vertex_color2_bw;
vcolor = Caa2.vertex_color2;
Caa = load('cache/fig_cluster2_reduce_annot.mat');
cc_cluster = Caa.cc;


% Make figure and subplots
h = figure('Position',round(0.8*[0, 0, 8.5, 11]*100));
set(h,'PaperUnits','inches');
set(h,'PaperPosition',[0, 0, 8.5, 11]);
[ha, pos] = tight_subplot(2,2,[.01 .01],[.01 .01],[.01 .01]);

trig_plot_marker_not_surf = false;
msize_elec = 12;
msize_elec_loc = 12;


axes(ha(1));
% --- First subplot: just electrodes --------------------------------------

%trig_all_elec_same = true;

% bypass face color
vcolor = 0.7*ones(size(vcolor));
trig_color_mwall = true; % parfor

Caci = load('cache/fig_cluster3_cluster_i');
cluster_i = Caci.cluster_i;

% Read medial wall
hemi = 'r';
l_mwall = read_label('fsaverage_sym',sprintf('%sh.Medial_wall',hemi));
vert_mwall = l_mwall(:,1) + 1;

col_coeff = 0.2;
BRAIN_MESH_ALPHA = 1; % 0.3
vcolor = (1*[1 1 1] * col_coeff + vcolor * (1-col_coeff));

trig_plot_scale_cube = false;

ist = 1;


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

    end

    % build master coherence table
    D = [];
    E = [];

    CaR = load(sprintf('%s/fig_cluster2_reduce_%i_new.mat',dir_cacheLp,iM));
    E = CaR.E;
    [n_E,~] = size(E);

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
    end
    Es_stdev = Es_var(:,4); %sqrt(real(Es_var(:,4)));
    Es_stdev_max = Es_var_max(:,4);

    sidintu = unique(E(:,1));
    n_sidintu = length(unique(E(:,1)));

    if (trig_plot_all)
        %hfig = figure('Visible','off');
        %set(hfig,'Position',1*[0 0 1920 1080])
        %set(hfig,'PaperSize',3*[7.111 4]);
        %set(hfig,'PaperUnits','inches');

        %hax1 = subplot(2,2,1,'Position',[0 0.5 0.5 0.5]);
        %hax1 = axes;

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
        cc = 0*ones(n_cmap,3); % for electrodes

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
        % Electrode marker radius
        r = 0.5;
        % Electrode marker quality
        n = 8;   
        [x,y,z] = sphere(n);
        [x4,y4,z4] = sphere(4);
        hold all;
        % Main representative electrode color (of 150)
        col_elec_loc = [0.8 0 0]; %[1 1 1]*0;


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

            
            
            
            % Color main electrode
            col_elec_loc = cc_cluster(i,:);
            col_elec_loc = col_elec;
            if (trig_plot_marker_not_surf)
                plot3(E(i,3),E(i,4),E(i,5),'.','Color',col_elec_loc,'MarkerSize',msize_elec_loc);
            else
                fac = 1;
                % Use low quality polygons
                X = x4*r*fac + E(i,3);
                Y = y4*r*fac + E(i,4);
                Z = z4*r*fac + E(i,5);
                q = surf(X,Y,Z,'facecolor',col_elec_loc,'LineStyle','none');
                q.SpecularStrength = 0;
                q.SpecularExponent = 15;
                q.SpecularColorReflectance = 0;
            end

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

                if (trig_plot_marker_not_surf)
                    plot3(xyz(1),xyz(2),xyz(3),'.','Color',col_elec,'MarkerSize',msize_elec);
                else
                    X = x4*r + xyz(1);
                    Y = y4*r + xyz(2);
                    Z = z4*r + xyz(3);
                    q = surf(X,Y,Z,'facecolor',col_elec,'LineStyle','none');
                    q.SpecularStrength = 0;
                    q.SpecularExponent = 15;
                    q.SpecularColorReflectance = 0;
                end
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
        
        % Cache plotted axis
        set(gca,'Visible','off');
        ax_1_1 = gca;
        save('cache/fig_S15_ax_1_1','ax_1_1');
    end
end
% -------------------------------------------------------------------------


print(h,sprintf('figures/figure_S16_%s_nbip-%i_nedge-%i_1stdev',surf_type,length(E),length(D)),'-dpng','-r400');
print(h,sprintf('figures/figure_S16_%s_nbip-%i_nedge-%i_1stdev',surf_type,length(E),length(D)),'-depsc','-r400');
close(h);