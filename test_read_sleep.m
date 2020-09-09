close all;
clear;

% Directory definitions
dir_annot = '/home/jerry/data/behavior/combined';
dir_h5 = '/home/jerry/data/h5_notch20';
dir_stamps = '/media/jerry/internal/data/stamps';
fn_annot = 'sub21_4ce6814b';



% Read in .sleep file
f_sleep = fopen(sprintf('%s/%s.sleep',dir_annot,fn_annot),'r');
sl = textscan(f_sleep,'%s','Delimiter','\n');
sl = sl{1};
fclose(f_sleep);

ss = strsplit(fn_annot,'_');
sid = ss{1};
subfile = ss{2};
ecog = H5eeg(sprintf('%s/%s.h5',dir_h5,sid));
subfiles = h5readatt(ecog.filename,'/h5eeg','files');
n_samples = h5readatt(ecog.filename,'/h5eeg','n_samples');
ds_factors = h5readatt(ecog.filename,'/h5eeg','ds_factors');
sIdx = find(strcmp(subfiles,subfile),1);
subf_start = [0;cumsum(n_samples)]+1;
subf_start = subf_start(1:(end-1));
if (~isempty(sIdx))
    start_idx = subf_start(sIdx);
    start_eidx = h5read(ecog.filename,'/h5eeg/aux',[1 start_idx],[1 1]);
    start_etime_absolute = ecog.getTimeFromSample(dir_stamps, start_idx);
    start_etime_absolute_dt = datetime(start_etime_absolute,'format','MM/dd/uuuu HH:mm:ss');
    % Time drift error occurs only once for calculating starting etime
    start_etime_seconds = start_eidx / (ecog.fs * ds_factors(sIdx));
else
    fprintf(2,'[!] subject %s, subfile %s not found.\n',sid,subfile)
end

Behavior = {};
Time = {};
Index = [];
count = 0;
for i = 1:length(sl)
%     try
        line = strsplit(sl{i},',');
        etimestr = line{1};
        behavior = line{2};
        ee = strsplit(etimestr,'.');
        etimesec = str2double(ee{2}) / 1000; % Milliseconds
        ee2 = strsplit(ee{1},':');
        etimesec = etimesec + str2double(ee2{3}); % Seconds
        etimesec = etimesec + 60*str2double(ee2{2}); % Minutes
        etimesec = etimesec + 60*60*str2double(ee2{1}); % Hours
        
        etime_absolute = datestr(start_etime_absolute_dt + seconds(etimesec-start_etime_seconds));
        etime_index = ecog.getSampleFromTime(dir_stamps,etime_absolute);
        
        Behavior = [Behavior; behavior];
        Time = [Time; etime_absolute];
        Index = [Index; etime_index];
        
        fprintf('[%i of %i] %s\t%i\t%s\n',i,length(sl),etime_absolute,etime_index,behavior)
%     catch
%         fprintf('[!] Skipped\n')
%     end
end

save(sprintf('%s_h5index',fn_annot))
