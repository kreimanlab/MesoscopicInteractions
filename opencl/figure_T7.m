close all;

fprintf('--- figure T3 ---\n')

if (~exist('AT','var'))
    load('xsub2/xsub_coh_all_v6.mat');
end

if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix    
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        resultsDir = '/media/klab/44/data/results';
        h5Dir = '/media/klab/44/h5';
    elseif strcmp(strip(hname),'ubuntu_1604')
        resultsDir = '/nas_share/RawData/data/results';
        h5Dir = '/nas_share/RawData/scripts/synth/out';
    else
        resultsDir = '/mnt/cuenap2/data/results';
        h5Dir = '/mnt/cuenap2/scripts/synth/out';
    end
end

clear A;
if (~exist('A','var'))
    A = Analysis(resultsDir,h5Dir);
end

system('mkdir figures/figure_T3');


% ROI pair (according to xsub/xsub_coh_all.mat)
%
% AT.P.AtlNames: 
%    'Destrieux-2009'    'Desikan-Killiany'    'Brodmann'    'Yeo-2011-7'    'Yeo-2011-17'    'HCP-MMP1'    'M132'    'FVE91'    'FVEall'
%    'LVE00'    'PHT00'    'Brodmann'    'BoninBailey'    'FerryEtAl'    'UD86'    'PGR91'    'LANDMARK'    'LyonKaas'    'BRL87'
%    'SP78'
atl_idx = 2;

% load ROI labels
roi_labels = AT.P.AtlROIs{atl_idx}.LH.struct_names;
% Apply ROI clustering
roi_labels = roi_labels(cluster_i);
% Remove extra ROIs
roi_labels = roi_labels(com_i);
% 24    - DK:parsorbitalis
% 9     - DK:inferiorparietal
roi_i = 24;
roi_j = 9;