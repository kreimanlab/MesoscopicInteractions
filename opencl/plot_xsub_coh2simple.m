close all;
clear;

BRAIN_MESH_ALPHA = 0.4; % Surface transparancy between 0,1
SMOOTH_FACE = true;
BLUR_ROI = true;
ELEC_RADIUS = 0.7;
COLOR_METRIC_MIN = 0.15;
COLOR_METRIC_MAX = 0.45;
M_LINE_WIDTH = 3;

THETA_TINT_COLOR = [0 1 0]; % lab
GAMMA_TINT_COLOR = [0 0 1]; %lab
N_TINT = 2;

%DELTA_TINT_COLOR = [30 50 0]; % lab
%THETA_TINT_COLOR = [30 50 45]; % lab
%ALPHA_TINT_COLOR = [30 50 90]; % lab
%BETA_TINT_COLOR = [30 50 135]; % lab
%GAMMA_TINT_COLOR = [30 50 180]; %lab
%N_TINT = 5;


% Path
[~,hname] = system('hostname');
if strcmp(strip(hname),'ubuntu_1604')
    resultsDir = '/nas_share/RawData/data/results';
    h5Dir = '/nas_share/RawData/scripts/synth/out';
else
    h5Dir = '/media/klab/44/h5';
end


%fsaverage_sym_dir = 'home/klab/data/h5eeg/artifact_removed/test/opencl/coreg/fsaverage_sym';
fsaverage_sym_dir = 'fsaverage_sym';

% Colormap
n_colormap = 100;
cc = corrcmap(n_colormap);
%cc = 1 - gray(n_colormap);

% Read
load('coreg/fsaverage_sym_flat');
load('xsub/xsub_coh_all','Br','De','Th','Al','Be','Ga');
%[vtx_infl, fac_infl] = read_surf('coreg/fsaverage_sym/surf/lh.inflated');
[vtx_pial, fac_pial] = read_surf('coreg/fsaverage_sym/surf/lh.inflated');

% %fv = patch('Faces',fac_pial,'Vertices',vtx_pial);
% fv.vertices = vtx_pial;
% fv.faces = fac_pial+1;
% N = patchnormals(fv);
% patch(fv,'FaceColor','red','EdgeColor','none');
% daspect([1,1,1])
% axis tight
% camlight
% camlight(-80,-10)
% lighting gouraud

%N = isonormals(D,FV.vertices);

% Color
C = ones(size(vtx_pial));

% ROI border
% fname = sprintf('coreg/fsaverage_sym/label/lh.aparc.annot');
% [vertices,label,colortable] = read_annotation(fname);
% vertices = vertices + 1;
% for i = 1:length(fac_pial)
%     vB = fac_pial(i,:);
%     tv = roi_id(vB);
%     if (length(unique(tv)) > 1)
%         plot(mean(pa.x(vB)),mean(pa.y(vB)),'black.')
%     end
% end

% Plot surface
p = trisurf(fac_pial + 1,vtx_pial(:,1),vtx_pial(:,2),vtx_pial(:,3),'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
set(p,'FaceVertexCdata',C);
ax = gca;
ax.Clipping = 'off';
hold on;
% Lighting
axis tight;
daspect([1 1 1]);
view(90,0);
if (SMOOTH_FACE)
    p.FaceLighting = 'gouraud';
end
% Intensity is reduced for 2 hemisphere plots
p.AmbientStrength = 0.3;
p.DiffuseStrength = 0.4;
p.SpecularStrength = 0;
p.SpecularExponent = 1;
p.BackFaceLighting = 'lit';
cam_elev = 0;
camlight(-135,cam_elev);
camlight(45,cam_elev);
camlight(-225,cam_elev);
camlight(-45,cam_elev);
axis vis3d off;
%lighting phong; % smooth lighting
%material dull;
if (BLUR_ROI)
    shading interp;
else
    shading flat;
end

n_rois = length(Br.bchan);
n_tot_rois = nchoosek(n_rois,2);
n_i = 1;
for i = 1:(n_rois-1)
    for j = (i+1):n_rois
        tic;
        % If is an interaction
        if (Br.isI(i,j))
            n_pairs = length(Br.xsub{i,j});
            
            % Save bip
            prev_sid = Br.sub{i,j}(1:6);
            h5fn = sprintf('%s/%s.h5',h5Dir,prev_sid);
            bip = h5readatt(h5fn,'/h5eeg/eeg','bip');
            l = read_label(fsaverage_sym_dir,sprintf('ielvis_%s',prev_sid)); %all_surf_
            
            B1c = nan(n_pairs,3);
            B2c = nan(n_pairs,3);
            for k = 1:n_pairs
                sub_start = (k-1)*6+1;
                sid = Br.sub{i,j}(sub_start:(sub_start+5));
                % only read if different
                if (~strcmp(sid,prev_sid))
                    h5fn = sprintf('%s/%s.h5',h5Dir,sid);
                    bip = h5readatt(h5fn,'/h5eeg/eeg','bip');
                    l = read_label(fsaverage_sym_dir,sprintf('ielvis_%s',sid)); %all_surf_
                    prev_sid = sid;
                end
                b1 = Br.bchan{i,j}(1,k);
                b1c1 = bip(b1,1);
                b1c2 = bip(b1,2);
                b2 = Br.bchan{i,j}(2,k);
                b2c1 = bip(b2,1);
                b2c2 = bip(b2,2);
                
                % Find electrodes
                if (~isempty(l))
                    b1c1_i = find(l(:,end) == b1c1,1);
                    b1c2_i = find(l(:,end) == b1c2,1);
                    b2c1_i = find(l(:,end) == b2c1,1);
                    b2c2_i = find(l(:,end) == b2c2,1);
                    % ==== TEMPORARY - THIS IS WRONG ========================== 
    %                 b1c1_i = b1c1;
    %                 b1c2_i = b1c2;
    %                 b2c1_i = b2c1;
    %                 b2c2_i = b2c2;
    %                 

                    % Electrode coordinates
                    b1c1_c = vtx_pial(l(b1c1_i,1),:);
                    b1c2_c = vtx_pial(l(b1c2_i,1),:);
                    b2c1_c = vtx_pial(l(b2c1_i,1),:);
                    b2c2_c = vtx_pial(l(b2c2_i,1),:);

                    % Average coordinates to get bip electrode
                    b1_c = (b1c1_c + b1c2_c)/2;
                    b2_c = (b2c1_c + b2c2_c)/2;
                    B1c(k,:) = b1_c;
                    B2c(k,:) = b2_c;
                end
                
                %fprintf('%s:\t%i-%i\t\t%i-%i\n',sid,b1c1,b1c2,b2c1,b2c2)
            end
            
            b1_c = nanmedian(B1c);
            b2_c = nanmedian(B2c);
            %COLOR_METRIC_MIN = 0.2;
            %COLOR_METRIC_MAX = 0.4;
            % Plot interactivity line
            m_norm = (nanmean(Br.xsubM{i,j}) - COLOR_METRIC_MIN)/(COLOR_METRIC_MAX - COLOR_METRIC_MIN);
            m_norm_th = (nanmean(Th.xsubM{i,j}) - COLOR_METRIC_MIN)/(COLOR_METRIC_MAX - COLOR_METRIC_MIN);
            m_norm_ga = (nanmean(Ga.xsubM{i,j}) - COLOR_METRIC_MIN)/(COLOR_METRIC_MAX - COLOR_METRIC_MIN);

            cc_i = ceil(m_norm*n_colormap);
            if (~isnan(cc_i) && (cc_i ~= 0))
                icolor = cc(cc_i,:);
                % Band tint
                icolor = icolor + (THETA_TINT_COLOR - icolor) * m_norm_th*(1/(N_TINT));
                icolor = icolor + (GAMMA_TINT_COLOR - icolor) * m_norm_ga*(1/(N_TINT));

                plot3([b1_c(1) b2_c(1)],[b1_c(2) b2_c(2)],[b1_c(3) b2_c(3)],'Color',icolor,'LineWidth',M_LINE_WIDTH)

                % Plot electrode 1 as sphere
                e_x = b1_c(1);
                e_y = b1_c(2);
                e_z = b1_c(3);
                n = 21;
                s_color = zeros(n+1,n+1,3);
                [x,y,z] = sphere(n);
                radius = ELEC_RADIUS;
                s = surf(radius*x+e_x,radius*y+e_y,radius*z+e_z,s_color);
                hold on;
                s.FaceLighting = 'none';
                s.EdgeColor = 'none';

                % Plot electrode 2 as sphere
                e_x = b2_c(1);
                e_y = b2_c(2);
                e_z = b2_c(3);
                n = 21;
                s_color = zeros(n+1,n+1,3);
                [x,y,z] = sphere(n);
                radius = ELEC_RADIUS;
                s = surf(radius*x+e_x,radius*y+e_y,radius*z+e_z,s_color);
                hold on;
                s.FaceLighting = 'none';
                s.EdgeColor = 'none';

            end
            %return
        end
        t_sing = toc;
        fprintf('%i of %i (ETA: %.2f min)\n',n_i,n_tot_rois,t_sing*(n_tot_rois-n_i)/60);
        n_i = n_i + 1;
    end
    %return
end
