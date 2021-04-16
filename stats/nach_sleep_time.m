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
Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};


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
            %return
            fprintf(2,'\t(%i/%i) %i - %i\n',j,length(ecog_n_samples),idx_starts(j),idx_stops(j));
        end
    end
    %ecog.getSampleFromTime('/media/klab/KLAB101/stamps','04/04/2007 19:10:55')
    
    %return
end