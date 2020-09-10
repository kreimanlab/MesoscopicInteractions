close all;
clear;

h5_dir = './h5';

% Get list of .h5 files
d = dir(h5_dir);
dname = {d.name};
disdir = [d.isdir];
d = dname(~disdir);
H5 = {};
j = 1;
for i = 1:length(d)
    if (endsWith(d{i},'.h5'))
        H5{j} = d{i};
        j = j + 1;
    end
end

% Calculate mu and sigma for each .h5 file
for i = 1:length(H5)
    infn = [h5_dir,'/',H5{i}];
    ecog = H5eeg(infn);
    size_mb = ecog.n_samples*32*(1/8)*(1e-6);
    mu = zeros(1,ecog.n_chan);
    sigma = zeros(1,ecog.n_chan);
    for j = 1:ecog.n_chan
        fprintf('Procesing chan %i of %i in %s (%.2f MB)\n',j,ecog.n_chan,infn,size_mb)
        tic;
        %v = h5read(infn,'/h5eeg/eeg',[j 1],[1 ecog.n_samples]);
        
        t_sing = toc;
        fprintf('\t%.2f MB/s\n',size_mb/t_sing)
    end
    return
end