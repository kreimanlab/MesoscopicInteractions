close all;
clear;

h5_dir = './h5';
hdf5_dir = '/mnt/cuenap2/data/h5eeg';

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

for i = 1:length(H5)
    
    sid = strsplit(H5{i},'.h5');
    sid = sid{1};
    % Get list of .hdf5 files
    d = dir(sprintf('%s/%s',hdf5_dir,sid));
    dname = {d.name};
    disdir = [d.isdir];
    d = dname(~disdir);
    Hd5 = {};
    j = 1;
    for i2 = 1:length(d)
        if (endsWith(d{i2},'.hdf5'))
            Hd5{j} = d{i2};
            j = j + 1;
        end
    end

    for i3 = 1:length(Hd5)
        infn = sprintf('%s/%s/%s',hdf5_dir,sid,Hd5{i3});
        aux_labels = h5readatt(infn,'/h5eeg/aux','labels');
        ekg_i = contains(lower(aux_labels),'ekg') | contains(lower(aux_labels),'ecg');
        n_ekg = sum(ekg_i);
        fprintf('%s\t%i\n',infn,n_ekg)
        %disp(aux_labels)
    end
end
