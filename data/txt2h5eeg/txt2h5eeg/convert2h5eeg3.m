clear all;

% Display header
fprintf('o>-----------------------------<o>-----------------------------<o\n');
fprintf('|                ECOG PUTATIVE INTERACTION MAPPER               |\n');
fprintf('|                              ~=~                              |\n');
fprintf('|                   .mat to H5EEG Conversion                    |\n');
fprintf('o>-----------------------------<o>-----------------------------<o\n');

% Get directory information
listing = dir;
cwd_names = {listing.name};
cwd_isdir = [listing.isdir];
Dirs = cwd_names(cwd_isdir);
Dirs = Dirs(3:end); % remove '.' and '..' from directory list

% Select directories to use
Dirs = Dirs(1);
fprintf('[!] Using directories:\n')
disp(Dirs)

% read electrode labels
channel_label = textscan(fopen('channel_labels.txt','r'),'%s','delimiter','\n');
channel_label = channel_label{1};
aux_label = textscan(fopen('aux_labels.txt','r'),'%s','delimiter','\n');
aux_label = aux_label{1};
isaux = textscan(fopen('isaux.txt','r'),'%f','delimiter','\n');
isaux = logical(isaux{1});
if ((length(aux_label)+length(channel_label)) ~= length(isaux))
    fprintf(2,'[!] Warning: channel and aux label files do not match length in isaux.txt\n');
end

% assemble a list of all available files
Flist = {};
Flistnaked = {};
N_chan = [];
N_samples = [];
Fs_all = [];
k = 1;
for i = 1:length(Dirs)
    % read directory information from peter
    flist = textscan(fopen([Dirs{i} '/peter_filelist'],'r'),'%s');
    fsumm = textscan(fopen([Dirs{i} '/peter_summary'],'r'),'%s');
    flist = flist{1};
    fsumm = fsumm{1};
    
    % accumulate file list
    for j = 1:length(flist)
        Flist{k,1} = [Dirs{i}, '/', flist{j}, '.mat'];
        ft = strsplit(flist{j},'_');
        Flistnaked{k,1} = [Dirs{i}, '/',ft{1},'_',ft{2},'_',ft{3},'_',ft{4}];
        %Flistnaked{k,1} = [Dirs{i}, '/',ft{1},'_',ft{2},'_',ft{3},'_',ft{4},'_',ft{5}];
        N_samples(k) = str2double(fsumm{2});
        N_chan(k) = str2double(fsumm{4});
        Fs_all(k) = str2double(fsumm{6});
        k = k + 1;
    end
end
Fs = Fs_all(1);
N_CHAN = N_chan(1)-1;

% Warn against channel label file
if (length(isaux) ~= N_CHAN)
    fprintf(2,'[!] Warning: number of channel labels is not equal to N_CHAN.\n');
end

% Check for data consistency
if ~(length(unique(N_chan)) == 1)
    fprintf(2, '[!] Error: number of channels for input files are different.\n');
    return
end

% Check for sampling frequency consistency
if ~(length(unique(Fs_all)) == 1)
    fprintf(2, '[!] Error: sampling frequencies for input files are different.\n');
    return
end

% Display information
cdir = strsplit(pwd,'/');
patient_name = cdir{end};
fprintf('[#] Patient: %s\n',cdir{end})
fprintf('[#] Current working directory: %s\n',pwd)

% calculate start time by extrapolating sample change
fprintf('[#] Calculating start time\n')
file = 1;
%disp([Flistnaked{file},'_100000','.txt'])
%firsts = textscan(fopen([Flistnaked{file},'_100000','.txt'],'r'),'%s','delimiter', '\n');
firsts = textscan(fopen([Flistnaked{file},'.txt'],'r'),'%s','delimiter', '\n');
firsts = firsts{1};
firstt = strsplit(firsts{1,1});
subtime = 0;
for i = 1:Fs
    fir = strsplit(firsts{i,1});
    if (strcmp(firstt{2},fir{2}))
        subtime = subtime + 1;
    end
end
msec = (Fs-subtime)/Fs;
start_datenum = datenum([firstt{1},' ',firstt{2}],'mm/dd/yyyy HH:MM:SS');
start_time = datestr(start_datenum,'yyyy-mm-ddTHH:MM:SS');
start_hdate = firstt{1};
start_htime = firstt{2};
start_etime = eval(firstt{3});

% read annotations
fprintf('[*] Reading in annotations.txt\n')
annotations = textscan(fopen('annotations.txt','r'),'%s','delimiter','\n');
annotations = annotations{1};
adate = 0;
rem = 0;
Annot = {};
j = 1;
eTimeBypass = false;
for i = 1:length(annotations)
    
    % Account for versions of XLTEK that include day number
    dayNum = 0;
    if (startsWith(annotations{i},'d'))
        lineL = strsplit(annotations{i});
        dayNum = eval(lineL{1}(2:end));
        n2skip = length(lineL{1});
        linet = annotations{i};
        annotations{i} = linet((n2skip+2):end);
    end
    
    lineL = strsplit(annotations{i});
    
    try
        %fprintf('%s\n',lineL{1});
        if (sum(isspace(lineL{1})==0)>0)
            if (~isempty(strfind(lineL{1},'.')))
                % When annot file contains elapsed time
                eTimeBypass = true;
                %fprintf('[!] Annot file contains elapsed time. This format is not implemented yet.\n')
                ett = strsplit(lineL{1},':');
                etime = round((60*60*eval(ett{1}) + 60*eval(ett{2}) + eval(ett{3}))*Fs);
                etdiff = (etime-annot_start_etime)/(Fs*60*60*24);
                Annot{j,1} = datestr(datenum([datestr(annot_start_datenum),' ',annot_start_htime])+etdiff,'mm/dd/yyyy');
                Annot{j,2} = datestr(datenum([datestr(annot_start_datenum),' ',annot_start_htime])+etdiff,'HH:MM:SS');
                Annot{j,3} = etime;
                Annot{j,4} = strjoin({lineL{2:end}});
                j = j + 1;
            else
                % When annot file contains absolute time
                atime = datenum(lineL{1},'HH:MM:SS');
                if (atime < rem)
                    adate = adate + 1;
                end
                if (dayNum == 0)
                    Annot1 = datestr(annot_start_datenum+adate,'mm/dd/yyyy');
                else
                    % Override date in new version of XLTEK
                    Annot1 = datestr(annot_start_datenum+(dayNum-1),'mm/dd/yyyy');
                end
                Annot2 = lineL{1};
                Annot{j,1} = Annot1;
                Annot{j,2} = Annot2;

                Annot{j,3} = datenum([Annot1,' ',Annot2],'mm/dd/yyyy HH:MM:SS');
                Annot{j,4} = strjoin({lineL{2:end}});
                Atime = strsplit(lineL{1},':');

                rem = atime;
                j = j + 1;
            end
        end
        %fprintf(['    ',num2str(dtime),'\n']);
    catch
        annotemp = strsplit(annotations{i},':');
        if (strcmp(annotemp{1},'Creation Date'))
            annotempspace = strsplit(annotations{i},'\t');
            annotempspace = annotempspace{2};
            fprintf('[*] %s found: %s\n',annotemp{1},annotempspace);
            annot_start_datenum = datenum(annotempspace(9:end),'mmm dd, yyyy');
            annot_start_htime = annotempspace(1:8);
            annot_start_hdate = datestr(annot_start_datenum,'mm/dd/yyyy');
            annot_start_etime = eval(firstt{3});
        else
            fprintf('-!- ignoring line: %s\n',annotations{i});
        end
    end
end

% skip through annotations to first sample in file
if (eTimeBypass)
    % outputs annot_marker, the index of first legit annotation
    fprintf('[!] Elapsed time found in annotation file.\n')
    n_annot = j - 1;
    edget = annot_start_etime + sum(N_samples);
    for annot_marker = 1:n_annot
        if (Annot{annot_marker,3} > edget)
            fprintf(2,'[!] Warn - annot etime: %i > last sample etime: %i\n',Annot{annot_marker,3},edget);
            fprintf(2,'\tTrimming annotations after this point.\n');
            % actual trimming
            Annot(annot_marker:end,:) = [];
            break;
        end
    end
    % first annotation follows creation date in file
    annot_marker = 1;
else
    fprintf('[!] Absolute time found in annotation file.\n')
    annot_marker = 1;
    annot_datenum = Annot{annot_marker,3};
    while (start_datenum > annot_datenum)
        annot_marker = annot_marker + 1;
        try
            annot_datenum = Annot{annot_marker,3};
        catch
            fprintf(2,'[!] Warning: annotation file event times not found in data.\n');
            annot_marker = annot_marker - 1;
            break
        end
        if (annot_marker == 1e6)
            fprintf(2,'[!] Warning: over 1e6 annotations were skipped, breaking..\n');
            break;
        end
    end
end
fprintf('[#] Skipped to line %i in annotation file (of %i)\n',annot_marker,j-1);

% ------------------ Begin conversion ------------------------------------%
h5fname = [patient_name,'.hdf5'];
system(['rm ',h5fname]);
% Add aux channels
AUX_CHAN = sum(isaux);

% save electrophysiological data
% Main EEG data
% * *`h5eeg/`* -- *Group* containing the data for this file.
%h5create(h5fname,'/h5eeg/eeg',[N_CHAN Inf],'ChunkSize',[1 Fs]);
%h5create(h5fname,'/h5eeg/eeg',[N_CHAN-AUX_CHAN Inf],'ChunkSize',[1 Fs]);

% --------------------------- CHANGED NOV28,2017 ---------------------------
% Increased chunksize, changed to single type
h5create(h5fname,'/h5eeg/eeg',[N_CHAN-AUX_CHAN Inf],...
    'ChunkSize',[(N_CHAN-AUX_CHAN),round(Fs*60)],'Datatype','single');

% Create hard link /h5eeg/raw to /h5eeg/eeg
fid = H5F.open(h5fname,'H5F_ACC_RDWR','H5P_DEFAULT');
gid = H5G.open(fid,'/h5eeg');
H5L.create_hard(gid,'eeg',gid,'raw','H5P_DEFAULT','H5P_DEFAULT');
H5G.close(gid);
H5F.close(fid);

% Aux
% * **`aux`** -- **Dataset** (Samples by Channels) containing non-electrophysiologial data used to mark/give meaning to the **`eeg`** dataset.
%h5create(h5fname,'/h5eeg/aux',[1 Inf],'ChunkSize',[1 Fs]);
% --------------------------- CHANGED NOV28,2017 ---------------------------
%% Increased chunksize, changed to single type
%h5create(h5fname,'/h5eeg/aux',[(1+AUX_CHAN) Inf],'ChunkSize',[1 Fs]);
h5create(h5fname,'/h5eeg/aux',[(1+AUX_CHAN) Inf],...
    'ChunkSize',[(1+AUX_CHAN),Fs*60],'Datatype','single');

place = 1;
TOTAL_SAMPLES = 0;
k = 0;
conversion_ratio = Fs*60*60*24;
for file = 1:length(Flist)
    % Development shortcut
    %if (file >= 3)
    %    break
    %end
    
    fname = Flist{file};
    fprintf('\t* saving file %s ..', fname);
    load(fname);
    
    % Write to /h5eeg/eeg
    [N_SAMPLES, ~] = size(data);
    start = [1 place];
    count = [(N_CHAN-AUX_CHAN) N_SAMPLES];
    %count = [N_CHAN N_SAMPLES];
    
    % Write ECoG electrodes
    %h5write(h5fname,'/h5eeg/eeg',data(:,[false,~isaux'])',start,count); %2:end
    % --------------------------- CHANGED NOV28,2017 ---------------------------
    % changed from double to single type
    h5write(h5fname,'/h5eeg/eeg',single(data(:,[false,~isaux'])'),start,count); %2:end
    
    % Write elapsed time and aux electrodes to /h5eeg/aux
    %h5write(h5fname,'/h5eeg/aux',data(:,[true,isaux'])',start,[(1+AUX_CHAN) N_SAMPLES]);
    h5write(h5fname,'/h5eeg/aux',single(data(:,[true,isaux'])'),start,[(1+AUX_CHAN) N_SAMPLES]);
    
    % h5write(h5fname,'/h5eeg/eeg',data(:,2:end)',start,count);
    %  Write to /h5eeg/aux
    % h5write(h5fname,'/h5eeg/aux',data(:,1)',start,[1 N_SAMPLES]);

    % Hold on to events for /h5eeg/events
    try
        current_datenum = Annot{annot_marker,3};
        if (eTimeBypass)
            data_datenum = data(:,1);
        else
            data_datenum = start_datenum + (data(:,1) - start_etime)/(conversion_ratio);
        end
        for i = 1:N_SAMPLES
            if (data_datenum(i) >= current_datenum)
                k = k + 1;
                annot_Names{k} = Annot{annot_marker,4};
                if (eTimeBypass)
                    annot_Idx(k) = Annot{annot_marker,3} - start_etime + 1;
                    if (annot_Idx(k) < 1)
                        fprintf(2,'[!] Warn: annotation time leads to negative index. Check if annotation file matches data file.\n');
                    end
                else
                    annot_Idx(k) = place + i - 1;
                end
                annot_Times{k} = strjoin(Annot(annot_marker,1:2));
                annot_marker = annot_marker + 1;
                try
                    current_datenum = Annot{annot_marker,3};
                catch
                    %fprintf(2,'[!] Warning: Annot inaccessible 2\n');
                    break
                end
            end
        end
    catch
        %fprintf(2,'[!] Warning: Annot inaccessible 1.\n');
    end
    % update place
    place = place + N_SAMPLES;
    clear data;
    fprintf(' Done.\n');
    TOTAL_SAMPLES = TOTAL_SAMPLES + N_SAMPLES;
end

% * Attribute: `subject` -- Subject identifier (optional)
h5writeatt(h5fname,'/h5eeg/','subject',patient_name);
% * Attribute: `experiment` -- Experiment identifier (optional)
h5writeatt(h5fname,'/h5eeg/','experiment','unknown');
% * Attribute: `datetime` -- Human readable date/time of start of recording (optional)
h5writeatt(h5fname,'/h5eeg/','datetime',start_time);
% * Attribute: `timestamp` -- POSIX formatted timestamp of start of recording (optional)
h5writeatt(h5fname,'/h5eeg/','timestamp',(start_datenum-datenum(1970,1,1))*86400+msec);

% * Attribute: `labels` -- One dimensional array of strings for channel labels
labels = [];
%for i = 1:N_CHAN
for i = 1:(N_CHAN-AUX_CHAN)
    newlabel = channel_label{i};
    while (length(newlabel) < 8)
        newlabel = [newlabel,' '];
    end
    labels = [labels;newlabel];
end
h5attput(h5fname,'/h5eeg/eeg','labels',channel_label);

% * Attribute: `offsets` -- One dimensional array of offsets in digital A/D units
h5writeatt(h5fname,'/h5eeg/eeg','offsets',zeros(1,N_CHAN-AUX_CHAN));
% * Attribute: `gains` -- One dimensional array of gains to convert channels to muV
h5writeatt(h5fname,'/h5eeg/eeg','gains',ones(1,N_CHAN-AUX_CHAN));
% * Attribute: `rate` -- Sampling Rate
h5writeatt(h5fname,'/h5eeg/eeg','rate',Fs);
% Write dimensions
h5writeatt(h5fname,'/h5eeg/eeg','n_chan',N_CHAN-AUX_CHAN);
h5writeatt(h5fname,'/h5eeg/eeg','n_samples',TOTAL_SAMPLES);

% % Aux
% % * **`aux`** -- **Dataset** (Samples by Channels) containing non-electrophysiologial data used to mark/give meaning to the **`eeg`** dataset.
% * Attribute: `labels` -- One dimensional array of strings for auxiliary channel labels.
%labels = {'Elapsed time'};
labels = ['Elapsed time', aux_label'];
h5attput(h5fname,'/h5eeg/aux','labels',labels);
% * Attribute: `rate` -- Sampling Rate
h5attput(h5fname,'/h5eeg/aux','rate',Fs);
% Write dimensions
h5writeatt(h5fname,'/h5eeg/aux','n_chan',AUX_CHAN+1);
h5writeatt(h5fname,'/h5eeg/aux','n_samples',TOTAL_SAMPLES);

% Events
% * **`events`** -- **Dataset** (Compound Datatype) containing indices of events in the **`eeg`** dataset.
annot = struct('name','FILE_BEGIN','start_idx',1,'duration',0);
annot.name = {annot.name};

% -- LOW LEVEL COMPOUND DATA TYPE MAGIC ----------------------------------%
%
%Create the required data types
%
strType = H5T.copy('H5T_C_S1');
H5T.set_size(strType, 'H5T_VARIABLE');
sz(1) = H5T.get_size(strType);
doubleType=H5T.copy('H5T_NATIVE_DOUBLE');
sz(2) = H5T.get_size(doubleType);
sz(3) = H5T.get_size(doubleType);

%
% Computer the offsets to each field. The first offset is always zero.
%
offset(1)=0;
offset(2:3)=cumsum(sz(1:2));

%
% Create the compound datatype for memory.
%
memtype = H5T.create ('H5T_COMPOUND', sum(sz));
H5T.insert (memtype, 'name',offset(1),strType);
H5T.insert (memtype, 'start_idx',offset(2), doubleType);
H5T.insert (memtype, 'duration',offset(3), doubleType);

%
% Create the compound datatype for the file.  Because the standard
% types we are using for the file may have different sizes than
% the corresponding native types, we must manually calculate the
% offset of each member.
%
filetype = H5T.create ('H5T_COMPOUND', sum(sz));
H5T.insert (filetype, 'name', offset(1),strType);
H5T.insert (filetype, 'start_idx', offset(2), doubleType);
H5T.insert (filetype, 'duration',offset(3), doubleType);

%
% Create dataspace.
%
dims = 1;
max_dims = {'H5S_UNLIMITED'};
chunk_dims = 1;
space = H5S.create_simple(1, dims, max_dims);
prop = H5P.create('H5P_DATASET_CREATE');
H5P.set_chunk(prop, chunk_dims);

%
% Create the dataset and write the compound data to it.
%
DATASET = '/h5eeg/events';
file = H5F.open (h5fname, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
dset = H5D.create (file, DATASET, filetype, space, 'H5P_DEFAULT',prop,'H5P_DEFAULT');
H5D.write (dset, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', annot);

%
% Extend the dataset
%
for i = 1:k
    % modify annotations to store
    
    annot.name = annot_Names(i);
    annot.start_idx = annot_Idx(i);
    annot.duration = 0;
    
    % write
    newind = i;
    H5D.set_extent (dset, newind);
    filespace = H5D.get_space (dset);
    H5S.select_hyperslab (filespace, 'H5S_SELECT_SET', i-1, 1, 1, 1);
    memspace = H5S.create_simple (1, 1, 1);
    H5D.write (dset, memtype, memspace, filespace, 'H5P_DEFAULT', annot);
end

%
% Clean up
%
H5D.close (dset);
H5S.close (space);
H5T.close (filetype);
H5F.close (file);
fprintf('[!] All Done.\n')

% -- LOW LEVEL COMPOUND DATA TYPE MAGIC ----------------------------------%

% Display final file
%h5disp(h5fname) % Warning takes a long time
%z = h5read('sub21.hdf5','/h5eeg/events')
%z.name

exit()
