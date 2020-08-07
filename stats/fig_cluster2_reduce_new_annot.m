close all;
clear;

% load ROI numbers
cluster_i = load('cache/fig_cluster3_cluster_i.mat');
cluster_i = cluster_i.cluster_i;

% create annot file for location areas
subjects_dir = '/media/jerry/internal/data/coreg';
fn_in = 'cache/fig_cluster2_reduce_1_new.mat';
fn_out = 'cache/fig_cluster2_reduce_annot';

% for arclength
sphere_radius = 100;

% Read sphere
hemi = 'r';
sid_const = 'fsaverage_sym';
[s_vert_sph, faces_sph] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,'sphere'));
faces_sph = faces_sph + 1 * (1 - min(faces_sph(:)));

% Read pial
hemi = 'r';
sid_const = 'fsaverage_sym';
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,'pial'));
faces = faces + 1 * (1 - min(faces(:)));

% Read medial wall
l_mwall = read_label('fsaverage_sym',sprintf('%sh.Medial_wall',hemi));
vert_mwall = l_mwall(:,1) + 1;

% Read locations
Ca = load(fn_in);

% generate colors
n_regions = length(Ca.A);
%cc0 = corrcmap(n_regions);
cc0 = getPyPlot_cMap('gist_rainbow', n_regions, false, 'python3');
%cc0 = 0.75*ones(n_regions,3);
% for i = 1:n_regions
%     f = (i-1)/n_regions;
%     a = 0.2 + 0.6*(cos(f*(pi)) + 1)/2;
%     b = 0 + 1*(sin(f*(pi)) + 0)/1;
%     c = 0 + 1*(cos(f*(pi)) + 0)/1;
%     %cc0(i,:) = cc0(i,:)*(1 - f) + [0.5 0 0]*(f);
%     col = lab2rgb([a*100 (0.5-b)*2*100 (0.5-c)*2*100],'OutputType','uint8');
%     cc0(i,:) = col / 256;
% end

%cc = cc(randperm(n_regions),:);

% organize colors
n_perm = 6*1200; %*2*2*2;
Arcl = zeros(1,n_perm);
Cc = cell(1,n_perm);
E_sphere0 = Ca.E_sphere;
fprintf('[*] Randomly shuffling regions until region colors no longer clash..\n')
n_nn = 6;
found_old = false;
if (exist('cache/fig_cluster2_reduce_annot.mat','file'))
    CaOld = load('cache/fig_cluster2_reduce_annot.mat');
    rIdx_old = CaOld.rIdx;
    keep_ratio = 1; %0.95;
    fprintf('\t* Found pre-existing shuffled file, keeping: %.3f.\n',keep_ratio)
    found_old = true;
    n_r = length(rIdx_old);
end
parfor k = 1:n_perm
    if (found_old)
        r2keep = false(1,n_r);
        r2keep(1:round(keep_ratio*n_r)) = true;
        r2keep = r2keep(randperm(n_r));
        rIdx(r2keep) = rIdx_old(r2keep);
        rIdxOldScr = rIdx_old(~r2keep);
        rIdx(~r2keep) = rIdxOldScr(randperm(length(rIdxOldScr)));
    else
        rIdx = randperm(n_regions);
    end
    cc = cc0(rIdx,:);
    E_sphere = E_sphere0(rIdx,:);
    %cc_arcl = nan(n_regions,n_regions);
    cc_arcl = nan(n_regions,1);
    for i = 1:(n_regions)
        u = E_sphere(i,3:5);
        
        iIdx = (i-n_nn):(i+n_nn);
        t_arcl = zeros(1,length(iIdx));
        for j = 1:length(iIdx)
            a = iIdx(j);
            if (a < 1)
                j2 = a + n_regions;
            elseif (a > n_regions)
                j2 = a - n_regions;
            else
                j2 = a;
            end
            
            v = E_sphere(j2,3:5);

            % Arclength
            theta = atan2d(norm(cross(u,v)),dot(u,v));
            theta = theta * (2*pi/360);
            arcl = theta*sphere_radius;

            t_arcl(j) = arcl;
            %cc_arcl(i,j) = arcl;
            %cc_arcl(j,i) = arcl;
        end
        cc_arcl(i,:) = nansum(t_arcl);
    end
    mean_arcl = nanmean(cc_arcl(:));
    Arcl(k) = mean_arcl;
    Cc{k} = rIdx;
end

[best,mIdx] = max(Arcl);
rIdx = Cc{mIdx};
cc = cc0(rIdx,:);

fprintf('Best arcl: %.6f mm\n',best/(2*n_nn+1));
%return

% Tint color by RAS
Rcol = [1 1 1]*1; %[1 0 0];
Acol = [1 1 1]*0; %[1 1 0];
Scol = [1 1 1]*0.5; %[0 0 1];
[n_E,~] = size(Ca.E);
alpha = 0; %0.9;
for i = 1:n_E
    r = Ca.E(i,3);
    a = Ca.E(i,4);
    s = Ca.E(i,5);
    col = cc(i,:);
    % mix
    r_f = alpha*( (max(Ca.E(:,3)) - r)/(max(Ca.E(:,3))-min(Ca.E(:,3))) );
    col = (1 - r_f)*col + r_f*Rcol;
    a_f = alpha*( (max(Ca.E(:,4)) - a)/(max(Ca.E(:,4))-min(Ca.E(:,4))) );
    col = (1 - a_f)*col + a_f*Acol;
    s_f = alpha*( (max(Ca.E(:,5)) - s)/(max(Ca.E(:,5))-min(Ca.E(:,5))) );
    col = (1 - s_f)*col + s_f*Scol;
    cc(i,:) = col;
end

% Vertex-wise location assignment
fprintf('[*] vertex-wise location assignment..\n')
V = [];
L = [];
for i = 1:length(Ca.Es_sphere)
    es = Ca.Es_sphere{i}(:,3:5);
    [n_es,~] = size(es);
    V = [V; es];
    L = [L; ones(n_es,1)*i];
end

[n_vert,~] = size(s_vert_sph);
vertex = 1:n_vert;
region = zeros(n_vert,1);
vertex_color = zeros(n_vert,3);
fprintf('[*] assigning..\n')

for i = 1:n_vert
    
%     % Arclength
%     theta = atan2d(norm(cross(u,v)),dot(u,v));
%     theta = theta * (2*pi/360);
%     arcl = theta*sphere_radius;
    
    coord = s_vert_sph(i,:);
    distance = (sum((coord - V).^2,2)); % sqrt
    [~,mIdx] = min(distance);
    region(i) = L(mIdx);
    vertex_color(i,:) = cc(L(mIdx),:);
end

fprintf('[*] saving ..\n')
save(fn_out,'vertex_color','region','s_vert','faces','cc','rIdx','cluster_i');
%return

%%


% -------------------------------------------------------------------------
% Redraw border -- takes a while
% -------------------------------------------------------------------------

trig_draw_border = true;
if (trig_draw_border)
    % region borders
    col_border = [1 1 1]*0.2;
    col_mwall = col_border; %[1 1 1]*0.1;
    region_medial = zeros(n_vert,1);
    region_border = zeros(n_vert,1);
    vertex_color2 = (vertex_color);
    vertex_color2_bw = 0.7*ones(size(vertex_color));
    [n_vert,~] = size(vertex_color);
    fprintf('[*] fetching borders..\n')
    parfor i = 1:n_vert
        j = (i == faces(:,1)) | (i == faces(:,2)) | (i == faces(:,3));
        v_neigh = unique(faces(j,:));
        n_uregion = length(unique(region(v_neigh)));
        
        % check if is in medial wall
        in_mwall = sum(vert_mwall == i) > 0;
        
        if (in_mwall)
            vertex_color2(i,:) = [col_mwall];
            vertex_color2_bw(i,:) = [col_mwall];
            region_medial(i) = 1;
        end
        
        if (n_uregion > 1)
            vertex_color2(i,:) = [col_border];
            vertex_color2_bw(i,:) = [col_border];
            region_border(i) = 1;
        end
        
    end
    save('cache/fig_cluster2_reduce_new_annot_vertex_color2','vertex_color2','vertex_color2_bw','region_border','region_medial');
else
    %load('cache/fig_cluster2_reduce_new_annot_vertex_color2.mat');
    vertex_color2 = vertex_color;
end


h = figure('PaperUnits','inches','PaperPosition',[0 0 4 4],'Position',[0 0 800 800]);
fprintf('[*] plotting..\n')
p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),'EdgeColor','none','FaceVertexCData',vertex_color2);
daspect([1 1 1]);
hold all;
shading('interp');
p.AmbientStrength = 0.3 ;
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
view(90,0);
axis off;

% color legend
ox = min(s_vert(:,1));
oy = min(s_vert(:,2));
oz = min(s_vert(:,3));
offset = 10;
for i = 1:length(cluster_i)
    c = Ca.E(i,3:5); %[ox,oy,oz + i*offset]; %
    b = text(c(1),c(2),c(3),sprintf('%i',cluster_i(i))); %,'Color',cc(i,:)); 
end
colormap(cc0);
colorbar;


print(h,'figures/fig_cluster2_reduce_new_annot','-depsc','-r400');

% Read annotation
% atlname = 'HCP-MMP1';
% hemi = 'r';
% sid_const = 'fsaverage_sym';
% fn_annot = sprintf('%s/%s/label/%sh.%s.annot',subjects_dir,sid_const,hemi,atlname);
% [vertices,label,ctab] = read_annotation(fn_annot);
% xIdx = ~ (all(ctab.table==0,2));
% roi_table = ctab.table(xIdx,:);
% roi_names = ctab.struct_names(xIdx,:);

fprintf('[!] Done.\n')
