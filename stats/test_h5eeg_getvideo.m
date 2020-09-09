close all;
clear all;

fn_h5 = '/media/klab/KLAB101/h5_notch20/sub14.h5';
dir_stamp = '/media/klab/internal/data/stamps';
dir_video = '/media/klab/internal/data/videos';
% fn_stamp_1 = '/media/klab/KLAB101/stamps/sub1_2710fd11_stamps.txt';
% f = fopen(fn_stamp_1,'r');
% fclose(f);
ecog = H5eeg(fn_h5);

sample = randi([1 ecog.n_samples]);
%sample = 45727847;
time = ecog.getTimeFromSample(dir_stamp,sample);
vids = ecog.getVideo(dir_video,sample,time);

fprintf('Sample: %i\tTime: %s\n',sample,time)
fprintf('Video: %s\n',vids)
