clear;
close all;

% --- TESTS BEFORE RUN ---
%       getCoreg.sh
%       bip_coverage.m
%       find $SUBJECTS_DIR | grep mmp
%       find $SUBJECTS_DIR | grep macaque.mat

% Stitch files all into one hdf5 file

% History of successfully ran patients
% Subs = {'m00006','m00019','m00023','m00024','m00026','m00030','m00037','m00038','m00043','m00060','m00068','m00083'};

% Queue
%Subs = {'m00032','m00039'}; % Skipped m00001 no artifacts found
%h5dir = '/share2/CACHEDEV1_DATA/RawData/h5eeg';
%h5dirArt = '/mnt/cuenap/data/h5eeg/artifact_removed';

%Subs = {'m00047','m00048','m00049','m00052','m00053','m00055'}; % skip 45 cant find n_chan. 55 didnt work. 47, 48 need redo

% Feb 28, 2018
%Subs = {'m00001','m00035','m00039','m00044','m00045','m00047','m00048','m00056','m00058','m00059','m00061','m00073','m00075','m00079','m00084','m00095','m00096','m00100'};

% Mar 6, 2018
%Subs = {'m00071'};

% Mar 8, 2018
%Subs = {'m00097','m00107'};

% Mar 9, 2018
%Subs = {'m00124'};

% Mar 12, 2018
%Subs = {'m00003'};

% Mar 16, 2018
%Subs = {'m00004','m00005','m00017','m00018','m00021','m00022','m00025','m00027','m00028','m00033','m00122'};

% Mar 22, 2018
%Subs = {'m00055'};

% Mar 23, 2018
%Subs = {'mSu'};

% Aug 21, 2018
%Subs = {'m00043'};

% May 9, 2019 for Adam, 2000Hz
%Subs = {'m00055'};

% June 11, 2019 Annabelle
%Subs = {'mChibi'};

% July 27, 2019
%Subs = {'m00068'};

% July 31, 2019
Subs = {'mChibi'};

%h5dir = '/share2/CACHEDEV1_DATA/RawData/data/h5eeg';
%h5dirArt = '/share2/CACHEDEV1_DATA/RawData/data/h5eeg/artifact_removed';
h5dir = '/Volumes/cuenap/data/h5eeg';
h5dirArt = '/Volumes/cuenap/data/h5eeg/artifact_removed';

FINAL_FS = 250; % target final sampling frequency (Hz)
DRY_RUN = false;
%DRY_RUN = true;
SAVE_FILE_LIST = true;

% ----- Find and check for hdf5 files ------------------------------------------
Files = cell(length(Subs),2);
for subI = 1:length(Subs)
    sid = Subs{subI};


    % inputs
    pid = sid;
    inDir = h5dir;

    % build input file list "infiles"
    pDir = [inDir,'/',pid,'/'];
    fprintf('Looking for clean h5eeg files in: %s\n',pDir)
    pd = dir(pDir);
    plist = {pd.name};
    ifi = 1;
    infiles = {};
    Ecogs = {};
    start_datenum = [];
    end_datenum = [];
    start_dstr = {};
    end_dstr = {};

    for i = 1:length(plist)
        if (endsWith(plist{i},'.hdf5'))
            infiles{ifi} = [pDir,plist{i}];
            ecog = H5eeg(infiles{ifi});

            % Save ecogs
            Ecogs{ifi} = ecog;

            % Calculate start and end dates
            start_datenum = [start_datenum; ecog.start_time];
            start_dstr{ifi} = datestr((ecog.start_time/ecog.fs)/86400);
            end_datenum = [end_datenum; ecog.end_time];
            end_dstr{ifi} = datestr((ecog.end_time/ecog.fs)/86400);
            ifi = ifi + 1;
        end
    end

    % sort infiles by starting time
    [~,sortIdx] = sort(start_datenum);
    infiles = infiles(sortIdx);
    Ecogs = Ecogs(sortIdx);
    start_dstr = start_dstr(sortIdx);
    start_datenum = start_datenum(sortIdx);
    end_dstr = end_dstr(sortIdx);
    end_datenum = end_datenum(sortIdx);
    n_chan = zeros(length(Ecogs),1);
    fs = zeros(length(Ecogs),1);
    total_time = zeros(length(Ecogs),1);
    for i = 1:length(Ecogs)
        n_chan(i) = Ecogs{i}.n_chan;
        fs(i) = round(Ecogs{i}.fs);
        total_time(i) = Ecogs{i}.n_samples / Ecogs{i}.fs;
    end

    % choose fs and chan for the whole subject based on the longest file
    [~,mIdx] = max(total_time);
    fs = fs(mIdx);
    n_chan = n_chan(mIdx);
    fprintf('[*     ] %s: %.4f Days\n',sid,sum(total_time)/(24*3600));


    if (DRY_RUN)
        % Check bip
        bipChan = bip(sid);

        chkFiles = {sprintf('%s/label/all_parcellation.mat',sid),...
                sprintf('%s/label/all_surf_ielvis.label',sid),...
                sprintf('%s/label/all_surf_ielvis_mmp.mat',sid),...
                sprintf('%s/label/all_surf_ielvis_macaque.mat',sid)};
        for i2 = 1:length(chkFiles)
            chkFile = chkFiles{i2};
            if (~exist(chkFile))
                fprintf(2,'[ERROR ] Could not find necessary file: %s\n',chkFile);
            end
        end
    end


    isOK = true;
    offset = 0;
    for iT = 1:length(infiles)
        i = iT + offset;
        if (i > 1)
            % Check to make sure times are ok
            isOK = (start_datenum(i) > end_datenum(i-1));

            % Check number of channels
            ecog = Ecogs{i};
            isOK = (isOK) & (ecog.n_chan == n_chan);
            if (~ (ecog.n_chan == n_chan))
                fprintf('[  !!  ] ecog.n_chan: %i, subject n_chan: %i\n',ecog.n_chan,n_chan);
            end

            % Check sampling frequency - try to salvage, don't drop
            %isOK = (isOK) & (round(ecog.fs) == fs);
            if (~ (round(ecog.fs) == fs))
                fprintf('[  !!  ] ecog.fs: %i, subject fs: %i\n',round(ecog.fs),fs);
            end
        end
        if (isOK)
            fprintf('[  OK  ] %s - %s (%.4f Hz) %s \n',start_dstr{i},end_dstr{i},Ecogs{i}.fs,infiles{i});
        else
            % if times overlap, pop out the smaller file
            popIdx = i - (Ecogs{i-1}.n_samples < Ecogs{i}.n_samples);
            fprintf(2,'[  !!  ] %s - %s %s \n',start_dstr{i},end_dstr{i},infiles{i});
            fprintf(2,'[ !FIX ] Dropping: %s\n',infiles{popIdx});
            infiles(popIdx) = [];
            Ecogs(popIdx) = [];
            offset = offset - 1;
        end
    end

    % find companion artifact files
    fprintf('Looking for companion artifact removed .hdf5 files\n')
    artfiles = {};
    offset = 0;
    for iT = 1:length(infiles)
        i = iT + offset;
        suffix = strsplit(Ecogs{i}.filename,'_');
        suffix = suffix{end};
        h5artname = sprintf('%s/%s/%s_%s',h5dirArt,sid,sid,suffix);
        if (exist(h5artname,'file'))
            artfiles{i} = h5artname;
        else
            fprintf(2,'[ !FIX ] Dropping companion to: %s\n',h5artname);
            % pop out filename
            infiles(i) = [];
            Ecogs(i) = [];
            offset = offset - 1;
        end
    end

    % gather variables: Ecogs, artfiles
    Files{subI,1} = Ecogs;
    Files{subI,2} = artfiles;
end

fprintf('[  OK  ] Filenames checked and ready.\n')



if (SAVE_FILE_LIST)
    % save mat
    save('synth_Files','Files','-V7.3');
    % save csv
    osf = fopen('synth_Files.csv','w');
    for i = 1:length(Subs)
        Ecogs = Files{i,1};
        artfname = Files{i,2};
        for j = 1:length(Ecogs)
            e = Ecogs{j};
            fprintf(osf,'%s,%s,%s,%i,%i,%.8f,%.8f,%.8f\n',...
            Subs{i},artfname{j},e.filename,e.n_chan,e.n_samples,e.start_time,e.end_time,e.fs);
        end
    end
end

if (DRY_RUN)
    return
end

% DEBUG
%return

% ----- Combine files ----------------------------------------------------------
[n_Files,~] = size(Files);
for fIdx = 1:n_Files

    % Pre-run checks
    skipSubject = false;
    % Check if subject is actionable
    if (length(Ecogs{1}) == 0)
        fprintf(2,'[  !!  ] Subject %s has no available files to convert.\n',sid);
        skipSubject = true;
    end

    if (~skipSubject)

        % parameters
        sid = Subs{fIdx};
        Ecogs = Files{fIdx,1};
        artfname = Files{fIdx,2};
        outfilename = sprintf('out2/%s.h5',sid);
        chunk_seconds = 60*60; % seconds
        datatype = 'single';
        max_mem = 1000; % maximum memory usage in MB
        bytes_double = 8;
        %n_chan_aux = h5readatt(Ecogs{1}.filename,'/h5eeg/aux','n_chan');
        n_chan_aux = 1;

        % Clear existing file
        % Check files and make parameters
        h5_n_chan = Ecogs{1}.n_chan;
        h5_n_samples = zeros(1,length(Ecogs));
        infiles = cell(1,length(Ecogs));
        suffixes = cell(1,length(Ecogs));
        start_timestrings = cell(1,length(Ecogs));
        timestamps = zeros(1,length(Ecogs));
        ds_factors = zeros(1,length(Ecogs));
        h5_fs = zeros(1,length(Ecogs));
        for i = 1:length(Ecogs)
            ds_factors(i) = round(Ecogs{i}.fs/FINAL_FS);
            h5_n_samples(i) = ceil(Ecogs{i}.n_samples / ds_factors(i));
            infiles{i} = Ecogs{i}.filename;
            h5_fs(i) = Ecogs{i}.fs / ds_factors(i);

            % Build extra atts
            suffix = strsplit(Ecogs{i}.filename,'_');
            suffix = strsplit(suffix{end},'.hdf5');
            suffixes{i} = suffix{1};
            att = Ecogs{i}.readAttributes();
            start_timestrings{i} = att.datetime;
            timestamps(i) = att.timestamp;
        
            fprintf('\tPlan downsample %.4f Hz --> %.4f Hz: %s\n',...
                Ecogs{i}.fs,Ecogs{i}.fs/ds_factors(i),Ecogs{i}.filename);
        end
        %h5_fs = mode(h5_fs);
        h5_fs = h5_fs(1);
        chunk_height = h5_n_chan;
        att = Ecogs{1}.readAttributes();
        chunk_len = chunk_seconds * round(h5_fs);

        % Build Bipolar ref
        bipChan = bip_chibi(sid);

        fprintf('[!] Clearing existing outfile: %s\n',outfilename)
        system(['rm ',outfilename]);
        fprintf('[#] Creating combined h5 file\n')
        fprintf('[*]\toutfilename: %s\n',outfilename)
        fprintf('[*]\tchunk height: %i\n',chunk_height)
        fprintf('[*]\tchunk width: %i\n',chunk_len)
        fprintf('[*]\tchunk width (sec): %i\n',chunk_seconds)
        fprintf('[*]\ttotal time: %.2f days\n',sum(h5_n_samples)/(h5_fs*3600*24))
        fprintf('[*]\tsum(h5_n_samples): %i\n',sum(h5_n_samples))
        fprintf('[*]\th5_n_chan: %i\n',h5_n_chan)
        fprintf('[*]\tn_chan_aux: %i\n',n_chan_aux)
        fprintf('[*]\th5_fs: %.4f\n',h5_fs)
        
        % create h5 file
        h5create(outfilename,'/h5eeg/eeg',[h5_n_chan,sum(h5_n_samples)],'ChunkSize',[chunk_height,chunk_len],'Datatype',datatype);
        h5create(outfilename,'/h5eeg/aux',[n_chan_aux,sum(h5_n_samples)],'ChunkSize',[1,chunk_len],'Datatype',datatype);
        % Common average referencing
        h5create(outfilename,'/h5eeg/cavg',[1,sum(h5_n_samples)],'ChunkSize',[1,chunk_len],'Datatype',datatype);
        h5writeatt(outfilename,'/h5eeg','subject',att.subject);
        h5writeatt(outfilename,'/h5eeg','experiment',att.experiment);
        h5attput(outfilename,'/h5eeg','datetime',start_timestrings);
        h5attput(outfilename,'/h5eeg','timestamp',timestamps);
        h5attput(outfilename,'/h5eeg','files',suffixes);
        h5attput(outfilename,'/h5eeg','n_samples',h5_n_samples);
        h5attput(outfilename,'/h5eeg','ds_factors',ds_factors);
        % att.events = h5read(self.filename{fIdx},'/h5eeg/events'); % copy events further down
        h5writeatt(outfilename,'/h5eeg/eeg','n_chan',h5_n_chan);
        h5writeatt(outfilename,'/h5eeg/eeg','n_samples',sum(h5_n_samples));
        h5writeatt(outfilename,'/h5eeg/eeg','rate',h5_fs);
        h5writeatt(outfilename,'/h5eeg/aux','n_chan',n_chan_aux);
        h5writeatt(outfilename,'/h5eeg/aux','n_samples',sum(h5_n_samples));
        h5writeatt(outfilename,'/h5eeg/aux','rate',h5_fs);
        
        % copy label cells
        atlasF = load(sprintf('%s/label/all_surf_ielvis_mmp',sid));
        atlasFmac = load(sprintf('%s/label/all_surf_ielvis_macaque',sid));
        atlasAll = load(sprintf('%s/label/all_parcellation',sid));

        % clean empty labels
        for i = 1:length(atlasF.AtlasLabels{1})
            if isempty(atlasF.AtlasLabels{1}{i})
                atlasF.AtlasLabels{1}{i} = 'UNKNOWN';
            end
        end
        for i = 1:length(atlasFmac.AtlasLabels{1})
            if isempty(atlasFmac.AtlasLabels{1}{i})
                atlasFmac.AtlasLabels{1}{i} = 'UNKNOWN';
            end
        end
        for i = 1:length(atlasFmac.AtlasLabels{2})
            if isempty(atlasFmac.AtlasLabels{2}{i})
                atlasFmac.AtlasLabels{2}{i} = 'UNKNOWN';
            end
        end
        for ii = 1:20
            for i = 1:length(atlasAll.AtlLabels{ii})
                if isempty(atlasAll.AtlLabels{ii}{i})
                    atlasAll.AtlLabels{ii}{i} = 'UNKNOWN';
                end
            end
        end

        channel_label = h5readatt(Ecogs{1}.filename,'/h5eeg/eeg','labels');
        h5attput(outfilename,'/h5eeg/eeg','labels',channel_label);
        h5attput(outfilename,'/h5eeg/eeg','bip',bipChan);
        h5attput(outfilename,'/h5eeg/eeg','labels_mmp',atlasF.AtlasLabels{1});
        h5attput(outfilename,'/h5eeg/eeg','labels_m132',atlasFmac.AtlasLabels{1});
        h5attput(outfilename,'/h5eeg/eeg','labels_fve',atlasFmac.AtlasLabels{2});
        for ii = 1:20
            h5attput(outfilename,'/h5eeg/eeg',sprintf('labels_%i',ii),atlasAll.AtlLabels{ii});
            h5attput(outfilename,'/h5eeg/eeg',sprintf('labels_%i_atlname',ii),atlasAll.AtlNames{ii});
            h5attput(outfilename,'/h5eeg/eeg',sprintf('labels_%i_roinames_l',ii),atlasAll.AtlROIs{ii}.LH.struct_names);
            h5attput(outfilename,'/h5eeg/eeg',sprintf('labels_%i_roinames_r',ii),atlasAll.AtlROIs{ii}.RH.struct_names);
        end
        h5attput(outfilename,'/h5eeg/eeg','hemi',atlasAll.EleHemi);
        aux_label = h5readatt(Ecogs{1}.filename,'/h5eeg/aux','labels');
        h5attput(outfilename,'/h5eeg/aux','labels',aux_label(1));


        % Copy events
        write_offset = [0,cumsum(h5_n_samples)];
        write_offset = write_offset(1:length(h5_n_samples));
        annot = struct('name','FILE_BEGIN','start_idx',1,'duration',0);
        annot.name = {annot.name};
        annotCount = 1;
        for subfIdx = 1:length(Ecogs)
            infilename = Ecogs{subfIdx}.filename;
            annotAll = h5read(infilename,'/h5eeg/events');
            
            for aIdx = 1:length(annotAll.start_idx)
                annot.name{annotCount} = annotAll.name{aIdx};

                % First downsample to file-specific index
                ev_start_idx = ceil(annotAll.start_idx(aIdx)/ds_factors(subfIdx));

                % Adjust for file offset
                ev_start_idx_offset = ev_start_idx + write_offset(subfIdx);

                annot.start_idx(annotCount) = ev_start_idx_offset;
                annot.duration(annotCount) = annotAll.duration(aIdx);
                %fprintf('[DEBUG ] %s annot start_idx: %i --> %i\n',suffixes{subfIdx},annotAll.start_idx(aIdx),annot.start_idx(annotCount))
                % TODO: fix sampling rate drifted indices with elapsed time files
                annotCount = annotCount + 1;
            end
        end
            
        % -- LOW LEVEL COMPOUND DATA TYPE MAGIC ----------------------------------%
        %Create the required data types
        strType = H5T.copy('H5T_C_S1');
        H5T.set_size(strType, 'H5T_VARIABLE');
        sz(1) = H5T.get_size(strType);
        doubleType=H5T.copy('H5T_NATIVE_DOUBLE');
        sz(2) = H5T.get_size(doubleType);
        sz(3) = H5T.get_size(doubleType);
        % Compute the offsets to each field. The first offset is always zero.
        offset(1)=0;
        offset(2:3)=cumsum(sz(1:2));
        % Create the compound datatype for memory.
        memtype = H5T.create ('H5T_COMPOUND', sum(sz));
        H5T.insert (memtype, 'name',offset(1),strType);
        H5T.insert (memtype, 'start_idx',offset(2), doubleType);
        H5T.insert (memtype, 'duration',offset(3), doubleType);
        % Write events
        filetype = H5T.create ('H5T_COMPOUND', sum(sz));
        H5T.insert (filetype, 'name', offset(1),strType);
        H5T.insert (filetype, 'start_idx', offset(2), doubleType);
        H5T.insert (filetype, 'duration',offset(3), doubleType);
        % Create dataspace
        dims = 1;
        max_dims = {'H5S_UNLIMITED'};
        chunk_dims = 1;
        space = H5S.create_simple(1, dims, max_dims);
        prop = H5P.create('H5P_DATASET_CREATE');
        H5P.set_chunk(prop, chunk_dims);
        % Create the dataset and write the compound data to it.
        DATASET = '/h5eeg/events';
        file = H5F.open (outfilename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
        dset = H5D.create (file, DATASET, filetype, space, 'H5P_DEFAULT',prop,'H5P_DEFAULT');
        H5D.write (dset, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', annot);
        % Extend the dataset
        for i = 1:length(annotAll.name)
            % modify annotations to store
            annot.name = annotAll.name(i); % must be cell, not string
            annot.start_idx = annotAll.start_idx(i);
            annot.duration = annotAll.duration(i);
            % write
            newind = i;
            H5D.set_extent (dset, newind);
            filespace = H5D.get_space (dset);
            H5S.select_hyperslab (filespace, 'H5S_SELECT_SET', i-1, 1, 1, 1);
            memspace = H5S.create_simple (1, 1, 1);
            H5D.write (dset, memtype, memspace, filespace, 'H5P_DEFAULT', annot);
        end
        % Clean up
        H5D.close (dset);
        H5S.close (space);
        H5T.close (filetype);
        H5F.close (file);
        fprintf('[+] Finished writing events.\n')
        % -- LOW LEVEL COMPOUND DATA TYPE MAGIC ----------------------------------%

        write_offset = [0,cumsum(h5_n_samples)];
        write_offset = write_offset(1:length(h5_n_samples));
        widths = zeros(1,length(Ecogs));
        all_arts = [];
        all_artifact_starts = [];
        for subfIdx = 1:length(Ecogs)
            infilename = Ecogs{subfIdx}.filename;
            infilenameArt = artfname{subfIdx};
            
            % copy artifacts
            tic;
            art = h5read(infilenameArt,'/h5eeg/artifacts'); %,[1 start_idx],[1 end_idx-start_idx+1]
            width = h5readatt(infilenameArt,'/h5eeg/artifacts','width');
            artifact_starts = ceil((1:width:(Ecogs{subfIdx}.n_samples))/ds_factors(subfIdx)) + write_offset(subfIdx);
            widths(subfIdx) = round(width/ds_factors(subfIdx));
            all_arts = [all_arts, art];
            all_artifact_starts = [all_artifact_starts, artifact_starts];
        end

        % Create artifact
        if (size(all_arts) ~= size(all_artifact_starts))
            fprintf(2,'[!] ERR: existing artifact annotations may be missing bins.\n')
        end
        
        h5create(outfilename,'/h5eeg/artifacts',[2,length(all_arts)],'Datatype',datatype);
        h5write(outfilename,'/h5eeg/artifacts',single([all_arts;all_artifact_starts])); %,[1 start_idx],[1 end_idx-start_idx+1]
        h5writeatt(outfilename,'/h5eeg/artifacts','width',widths(1));
        t_art = toc;
        fprintf('[+] Artifacts copied (%.2f sec).\n',t_art)
        
        % Create unplugged
        h5create(outfilename,'/h5eeg/unplugged',[2,length(all_arts)],'Datatype',datatype);
        % Write placeholder unplugged
        h5write(outfilename,'/h5eeg/unplugged',single([zeros(size(all_arts));all_artifact_starts])); %,[1 start_idx],[1 end_idx-start_idx+1]
        h5writeatt(outfilename,'/h5eeg/unplugged','width',widths(1));
        fprintf('[+] Wrote zeros to unplugged.\n')
        %return

        % Calculate pre-copy params
        write_offset = [0,cumsum(h5_n_samples)];
        write_offset = write_offset(1:length(h5_n_samples));
        read_bytes = 0;
        for iTemp = 1:length(Ecogs)
            read_bytes = read_bytes + Ecogs{iTemp}.n_samples * (Ecogs{iTemp}.n_chan+2) * bytes_double;
        end
        write_bytes = sum(h5_n_samples)*(2+h5_n_chan)*(bytes_double/2);
        fprintf('[#] Start EEG copy\n');
        fprintf('[*]\tApprox. read bytes: %.3f GB\n',read_bytes * 1e-9);
        fprintf('[*]\tApprox. write bytes: %.3f GB\n',write_bytes * 1e-9);

        % Estimate run time
        tic;
        testBlock = [Ecogs{1}.n_chan,round(Ecogs{1}.fs)];
        n_test = 10;
        for iT2 = 1:n_test
            testDat = h5read(Ecogs{1}.filename,'/h5eeg/eeg',[1 1],testBlock);
        end
        tRead = toc;
        testSpeed = (n_test*prod(testBlock)*(bytes_double/2))/tRead;
        fprintf('[*]\tTested read speed: %.3f MB/s\n',testSpeed*1e-6);
        fprintf('[*]\tTheoretical read time: %.1f hrs\n',(read_bytes/testSpeed)/(3600));
        fprintf('[*]\tTheoretical write time: %.1f hrs\n',(write_bytes/testSpeed)/(3600));
        fprintf('[*]\tTheoretical total I/O time: %.1f hrs\n',((read_bytes+write_bytes)/testSpeed)/(3600));
        
        %return 

        % Copy
        for subfIdx = 1:length(Ecogs)
            infilename = Ecogs{subfIdx}.filename;
            infilenameArt = artfname{subfIdx};
            suffix = strsplit(Ecogs{subfIdx}.filename,'_');
            suffix = suffix{end};
                
            % Round max_mem_samples to the nearest integer multiple of downsample factor
            max_mem_samples = round((max_mem*1e6/bytes_double)/h5_n_chan);

            %max_mem_samples = 100*chunk_len;
            max_mem_samples = floor(max_mem_samples/ds_factors(subfIdx))*ds_factors(subfIdx);

            % Start indices
            cIdx = 1:max_mem_samples:Ecogs{subfIdx}.n_samples;
            n_cIdx = length(cIdx);
            for i = 1:n_cIdx
                start_idx = cIdx(i);
                end_idx = start_idx + max_mem_samples - 1;
                if (end_idx > Ecogs{subfIdx}.n_samples)
                    end_idx = Ecogs{subfIdx}.n_samples;
                end
                fprintf('\t(%s_%s) start: %i\tend: %i (%.2f %%)\n',sid,suffix,start_idx,end_idx,100*i/n_cIdx) % dev

                % copy cavg
                tic;
                cavg = h5read(infilenameArt,'/h5eeg/cavg',[1 start_idx],[1 end_idx-start_idx+1]);
                % downsample
                cavg_ds = downsample(cavg,ds_factors(subfIdx));
                start_idx_ds = ceil(start_idx/ds_factors(subfIdx)) + write_offset(subfIdx);
                h5write(outfilename,'/h5eeg/cavg',single(cavg_ds),[1 start_idx_ds],size(cavg_ds))
                t_cavg = toc;
                %fprintf('\t\tcavg copied: %.2f s\n',t_cavg);
                
                % copy aux
                for j = 1:n_chan_aux
                    tic;
                    aux = h5read(infilename,'/h5eeg/aux',[j start_idx],[1 end_idx-start_idx+1]);
                    t_aux_read = toc;
                    aux_ds = downsample(aux,ds_factors(subfIdx));
                    h5write(outfilename,'/h5eeg/aux',single(aux_ds),[j start_idx_ds],size(aux_ds));
                    t_aux = toc;
                    fprintf('\t\tAux read: %.2f MB/s\n',((end_idx-start_idx+1)*bytes_double*1e-6)/t_aux_read)
                    %fprintf('\t\tAux copied: %.2f s, %i of %i, ETA: %.2f min\n',t_aux,j,n_chan_aux,(n_chan_aux-j)*t_aux/60);
                end

                % copy eeg

                % Read
                tic;
                eeg = h5read(infilename,'/h5eeg/eeg',[1 start_idx],[h5_n_chan end_idx-start_idx+1]);
                t_read = toc;
                fprintf('\t\tEEG read: %.2f MB/s\n',h5_n_chan*(end_idx-start_idx+1)*bytes_double*1e-6/t_read)

                % Downsample
                eeg_ds = downsample(eeg',ds_factors(subfIdx))';

                % Write
                tic;
                h5write(outfilename,'/h5eeg/eeg',single(eeg_ds),[1 start_idx_ds],size(eeg_ds));
                t_write = toc;
                fprintf('\t\tEEG write: %.2f MB/s\n',numel(eeg_ds)*(bytes_double/2)*1e-6/t_write)
                fprintf('\t(*) ETA: %.2f hrs.\n',(n_cIdx-i)*((t_read+t_write+t_cavg)/3600))
            end
        end



    end

end

%end
fprintf('[!] All Done.\n')
