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
Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48',...
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