close all;
clear;
clc;
rng shuffle;

% Get atlas names
CaAtl = load(sprintf('./cache/xsub_out_all_%i',1));
AtlNames = CaAtl.C.AtlNames;

% Read figure f14 results
iM = 1;
atl = 2;
Ca14 = load(sprintf('./cache/figure_t14_%i_atl%i_%s',iM,atl,AtlNames{atl}));


% Init paths
% dir_art = '/media/klab/KLAB101/h5_notch20/art_nosz'; %'/media/klab/internal/data/h5_notch20/art';
% dir_res = '/media/klab/KLAB101/results/coh_w10';
% dir_cor = '/media/klab/internal/data/coreg'; %'/mnt/cuenap_ssd/coregistration';
% dir_h5 = '/media/klab/KLAB101/h5_notch20';
% setenv('SUBJECTS_DIR',dir_cor);
% Patients
Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};


for i = 1:length(Subjects)
    sid = Subjects{i};
    fprintf('[%s]\n',sid);
    
    % Read xsub_out_stats
    ca_fname = sprintf('cache/xsub_out_%s_%i_atl%i',sid,iM,atl);
    Ca = load(ca_fname);
    
    ecog = H5eeg(sprintf('%s/%s.h5',Ca.dir_h5,sid));
    
    % subfile indices
    ecog_files = h5readatt(ecog.filename,'/h5eeg','files');
    ecog_n_samples = h5readatt(ecog.filename,'/h5eeg','n_samples');
    
    % start and stops
    idx_starts = cumsum([1; ecog_n_samples]);
    idx_starts = idx_starts(1:(end-1));
    idx_stops = cumsum([1; ecog_n_samples])-1;
    idx_stops = idx_stops(2:end);
    
    % Find times and breaks
    for j = 1:length(ecog_n_samples)
        try
            time_start = ecog.getTimeFromSample('/media/klab/KLAB101/stamps',idx_starts(j));
            time_stop = ecog.getTimeFromSample('/media/klab/KLAB101/stamps',idx_stops(j));
            fprintf('\t(%i/%i) %s - %s\n',j,length(ecog_n_samples),time_start,time_stop);
        catch
            return
            fprintf(2,'\t(%i/%i) %i - %i\n',j,length(ecog_n_samples),idx_starts(j),idx_stops(j));
        end
    end
    %ecog.getSampleFromTime('/media/klab/KLAB101/stamps','04/04/2007 19:10:55')
    
    %return
end