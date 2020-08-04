close all;
clear;

% dir def
dir_cache = './cache';
dir_cor = '/media/klab/internal/data/coreg';

% Find overlap with these areas
dk_roi1 = 'superiortemporal';
dk_roi2 = 'parsopercularis';

iM = 1;
CaC = load('cache/fig_cluster2_reduce_annot');
CaC_border = load('cache/fig_cluster2_reduce_new_annot_vertex_color2.mat');
Ca2 = load(sprintf('%s/fig_cluster2_reduce_%i_new.mat',dir_cache,iM));

% Read fsaverage_sym
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_cor,'fsaverage_sym','r','pial'));
faces = faces + 1;

p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),'EdgeColor','None');
set(p,'FaceVertexCData',CaC.vertex_color);
brainlight;