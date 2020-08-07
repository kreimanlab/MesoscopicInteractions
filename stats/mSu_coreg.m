close all;
clear;


% Register to macaque labels
D_anat = load('../coreg/AdjMKFV.mat');
D = load('../coreg/all_surf_ielvis_macaque.mat');
load('SuMAP.mat');
crop_pix = 50; % vertical crop
crop_pix_h = 20; % horizontal crop
Su.I = Su.I(crop_pix:(end-crop_pix), crop_pix_h:(end-crop_pix_h),:);
%dist_thresh = 20;
atl = 1;
atl_name = D.AtlasName{atl};
atl_labels = D.AtlasLabels{atl};