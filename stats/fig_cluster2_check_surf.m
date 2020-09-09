close all;
clear;


cmap = corrcmap(100);
alpha_pial = 0.3;
col_pial = 0.7*[1 1 1];
SURFACE_TYPE = 'pial';
hemi = 'r';
sid_const = 'fsaverage_sym';
subjects_dir = '/media/jerry/internal/data/coreg';
CaR = load('cache/fig_cluster2_reduce.mat');

[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,SURFACE_TYPE));
p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
    'EdgeColor','none','FaceColor',col_pial,'FaceAlpha',alpha_pial);
hold all;
daspect([1 1 1]);
p.AmbientStrength = 0.3 ;
p.DiffuseStrength = 0.4 ;
p.SpecularStrength = 0;
p.SpecularExponent = 1;
p.BackFaceLighting = 'lit';
cam_elev = 0;
camlight(-135,cam_elev);
camlight(45,cam_elev);
camlight(-225,cam_elev);
camlight(-45,cam_elev);

view(90,0);
axis off;
set(gcf,'Position',[0 0 1920 1080]);



% Plot electrodes
[n_E,~] = size(CaR.E);
for i = 1:n_E
    e = CaR.E(i,:);
    plot3(e(3),e(4),e(5),'black.','MarkerSize',40); hold on;
end