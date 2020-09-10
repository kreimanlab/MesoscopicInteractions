%analysis1

close all; clear;

% Load results
[~,hname] = system('hostname');
if strcmp(strip(hname),'hopper')
    %resultsDir = '/media/klab/44/data/results';
    resultsDir = '/mnt/cuenap2/data/results/coh_w30';
    h5Dir = '/media/klab/44/h5';
    artDir = '/mnt/cuenap2/scripts/synth/out_art';
elseif strcmp(strip(hname),'ubuntu_1604')
    resultsDir = '/nas_share/RawData/data/results/coh_w30';
    h5Dir = '/nas_share/RawData/scripts/synth/out';
    artDir = '/nas_share/RawData/scripts/synth/out_art';
end
A = Analysis(resultsDir,h5Dir);


