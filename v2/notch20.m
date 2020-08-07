close all;
clear;

[~,host] = system('hostname');
if contains(host,'ubuntu_1604')
    h5dir_in = '/nas_share/cuenap/data/h5';
    h5dir_out = '/nas_share/cuenap/data/h5_notch20';
elseif contains(host,'o2.rc.hms.harvard.edu')
    h5dir_in = '/n/scratch2/jw324/data/h5';
    h5dir_out = '/n/scratch2/jw324/data/h5_notch20';
else
    h5dir_in = '/media/klab/44/h5';
    %h5dir_out = '/media/klab/44/h5_notch20';
    h5dir_out = '/media/klab/D0BA8713BA86F56E/data/h5_notch20';
end
toggle_rsync = false;

system(sprintf('mkdir %s', h5dir_out));

fprintf('Apply 20 Hz notch filter\n')
d = dir(h5dir_in);
disdir = [d.isdir];
dname = {d.name};
Subs = dname(~disdir);

%=============================================================================
Subs = {'m00043.h5'};

Sub = cell(size(Subs));
Subo = cell(size(Subs));
for i = 1:length(Sub)
    Sub{i} = sprintf('%s/%s',h5dir_in,Subs{i});
    Subo{i} = sprintf('%s/%s',h5dir_out,Subs{i});
    cpycmd = sprintf('rsync -avPI %s %s',Sub{i},Subo{i});
    fprintf('\t%s\n',Sub{i});
    
    if (toggle_rsync)
        fprintf('\t%s\n',cpycmd);
        system(cpycmd);
    end
    
    % Get chunksize
    I = h5info(Sub{i});
    Ics = {I.Groups.Datasets.ChunkSize};
    In = {I.Groups.Datasets.Name};
    Ics = Ics{4};
    In = In{4};
    fprintf('\tchunksize %s: %i, %i\n',In,Ics(1),Ics(2));
    
    % Build start indices
    ecog = H5eeg(Sub{i});
    Fs = round(ecog.fs);
    Starts = 1:Ics(2):ecog.n_samples;
    Ends = Starts + Ics(2) - 1;
    Ends(Ends > ecog.n_samples) = ecog.n_samples;
    for j = 1:length(Starts)
        
        % Read
        read_chunk = [Ics(1) (Ends(j) - Starts(j) + 1)];
        sizemb = prod(read_chunk)*32*(1/8)*(1e-6);
        tic;
        V = h5read(Sub{i},'/h5eeg/eeg',[1 Starts(j)],read_chunk);
        t_read = toc;
        
        % Filter
        Vf = single(zeros(size(V)));
        for k = 1:ecog.n_chan
            Vf(k,:) = single(FilterData(double(V(k,:)),Fs,'notch',20));
        end
        
        % Write
        tic;
        h5write(Subo{i},'/h5eeg/eeg',Vf,[1 Starts(j)],read_chunk);
        t_write = toc;
        
        fprintf('\t\t[%.0f%%]\tRead %s at %.2f MB/s, Wrote at %.2f MB/s\n',100*Starts(j)/Starts(end),Subs{i},sizemb/t_read,sizemb/t_write)
        
    end
    %return
    
    
    
end
