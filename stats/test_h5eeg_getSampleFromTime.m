close all; clear;

fn_h5 = '/media/klab/KLAB101/h5_notch20/sub1.h5';
dir_stamp = '/media/klab/internal/data/stamps';
dir_video = '/media/klab/internal/data/videos';
ecog = H5eeg(fn_h5);
ecog.getSampleFromTime(dir_stamp, '07-Apr-2007 01:55:09')