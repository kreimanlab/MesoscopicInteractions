close all;
clear;

D = dir('./h5');
Dname = {D.name};
Disdir = [D.isdir];
h5fs = Dname(~Disdir);

for i = 1:length(h5fs)
    h5fn = sprintf('./h5/%s',h5fs{i});

    fs = h5readatt(h5fn,'/h5eeg/eeg','rate');
    n_samples = h5readatt(h5fn,'/h5eeg/eeg','n_samples');
    days = round(n_samples/fs)/(3600*24);
    fprintf('%s\t%.2f\n',h5fn,days);
end