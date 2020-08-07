close all;
clear;

dir_h5 = '/media/klab/internal/data/h5_notch20'; %'/media/klab/KLAB101/h5_notch20';
dir_art = sprintf('%s/art_nosz',dir_h5); %dir_artLp = '/media/klab/KLAB101/h5_notch20/art_nosz';
dir_res = '/home/jerry/data/results/coh_w10'; %dir_resLp = '/media/klab/KLAB101/results/coh_w10'; 
dir_cor = '/media/klab/internal/data/coreg';
dir_cache = './cache';
dir_stamp = '/media/klab/internal/data/stamps';
dir_video = '/media/klab/internal/data/videos';

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

% Patients
Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124',...
    'mSu'};

% Exclude monkey
Subjects = Subjects(1:(end-1));

for iM = 1:5
    coh_thresh = [];
    coh_thresh_5 = [];
    for i = 1:length(Subjects)
        sid = Subjects{i};
        %fprintf('[%s-%s]\n',metrics{iM},sid);
        Ca = load(sprintf('cache/xsub_out_%s_%i',sid,iM));
        Ca5 = load(sprintf('cache/archive/xsub_out_%s_%i',sid,iM));
        sat_dist = Ca.Dmats > Ca.dist_thresh;
        coh_thresh = [coh_thresh; Ca.coh_thresh(sat_dist)];
        coh_thresh_5 = [coh_thresh_5; Ca5.coh_thresh(sat_dist)];
    end
    fprintf('[%s] n=%i \n',metrics{iM},length(coh_thresh));
    fprintf('[*] res 0.5hz\n');
    fprintf('\tmean: %.6f\n',mean(coh_thresh));
    fprintf('\tstd: %.6f\n',std(coh_thresh));
    fprintf('\tmedian: %.6f\n',median(coh_thresh));
    fprintf('[*] res 5hz\n');
    fprintf('\tmean: %.6f\n',mean(coh_thresh_5));
    fprintf('\tstd: %.6f\n',std(coh_thresh_5));
    fprintf('\tmedian: %.6f\n',median(coh_thresh_5));
    fprintf('[*] difference: %.6f\n',mean(coh_thresh_5) - mean(coh_thresh));
    %return
end