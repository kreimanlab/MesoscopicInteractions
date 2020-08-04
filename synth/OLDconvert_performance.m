% convert_performance
clear; clc; close all;

% inputs
pid = 'm00083';
inDir = '../../';
%inDir = '/mnt/cuenap/data/h5eeg/artifact_removed/';

% build input file list
pDir = [inDir,pid,'/'];
pd = dir(pDir);
plist = {pd.name};
ifi = 1;
infiles = {};
outfiles = {};
for i = 1:length(plist)
    if (endsWith(plist{i},'.hdf5'))
        infiles{ifi} = [pDir,plist{i}];
        outfiles{ifi} = strrep(infiles{ifi},'.hdf5','.h5');
        ifi = ifi + 1;
    end
end

for fIdx = 1:length(infiles)

    infilename = infiles{fIdx};
    outfilename = outfiles{fIdx};
    %infilename = '../../m00023/m00023_0c36dc94.hdf5';
    %outfilename = '../../m00023/m00023_0c36dc94.h5';
    chunk_seconds = 60; % seconds
    
    datatype = 'single';
    max_mem = 50; % maximum memory usage in MB
    bytes_double = 8;
    n_chan_aux = h5readatt(infilename,'/h5eeg/aux','n_chan');

    % clear existing file
    fprintf('[!] infile: %s\n[!] outfile: %s\n',infilename,outfilename)
    system(['rm ',outfilename]);

    % read
    ecog = H5eeg2({infilename});
    chunk_height = ecog.n_chan; % --- Chunk height ------------------------
    att = ecog.readAttributes(1);
    chunk_len = chunk_seconds * round(ecog.fs);
    if (ecog.n_samples <= 0)
        fprintf(2,'[WTF?] input n_samples not a positive number.\n');
    end

    h5create(outfilename,'/h5eeg/eeg',[ecog.n_chan,ecog.n_samples],'ChunkSize',[chunk_height,chunk_len],'Datatype',datatype);
    h5create(outfilename,'/h5eeg/aux',[n_chan_aux,ecog.n_samples],'ChunkSize',[1,chunk_len],'Datatype',datatype);
    h5create(outfilename,'/h5eeg/cavg',[1,ecog.n_samples],'ChunkSize',[1,chunk_len],'Datatype',datatype);
    h5writeatt(outfilename,'/h5eeg','subject',att.subject);
    h5writeatt(outfilename,'/h5eeg','experiment',att.experiment);
    h5writeatt(outfilename,'/h5eeg','datetime',att.datetime);
    h5writeatt(outfilename,'/h5eeg','timestamp',att.timestamp);
    % att.events = h5read(self.filename{fIdx},'/h5eeg/events'); % TODO: copy events
    h5writeatt(outfilename,'/h5eeg/eeg','n_chan',ecog.n_chan);
    h5writeatt(outfilename,'/h5eeg/eeg','n_samples',ecog.n_samples);
    h5writeatt(outfilename,'/h5eeg/eeg','rate',ecog.fs);
    h5writeatt(outfilename,'/h5eeg/aux','n_chan',h5readatt(infilename,'/h5eeg/aux','n_chan'));
    h5writeatt(outfilename,'/h5eeg/aux','n_samples',h5readatt(infilename,'/h5eeg/aux','n_samples'));
    h5writeatt(outfilename,'/h5eeg/aux','rate',h5readatt(infilename,'/h5eeg/aux','rate'));
    
    % copy label cells
    h5attput(outfilename,'/h5eeg/eeg','labels',h5readatt(infilename,'/h5eeg/eeg','labels'));
    h5attput(outfilename,'/h5eeg/aux','labels',h5readatt(infilename,'/h5eeg/aux','labels'));
   % h5writeatt(outfilename,'/h5eeg/aux','rate',h5readatt(infilename,'/h5eeg/aux','labels'));

    % copy artifacts
    tic;
    art = h5read(infilename,'/h5eeg/artifacts'); %,[1 start_idx],[1 end_idx-start_idx+1]
    width = h5readatt(infilename,'/h5eeg/artifacts','width');
    h5create(outfilename,'/h5eeg/artifacts',size(art),'Datatype',datatype);
    h5write(outfilename,'/h5eeg/artifacts',single(art)); %,[1 start_idx],[1 end_idx-start_idx+1]
    h5writeatt(outfilename,'/h5eeg/artifacts','width',width);
    t_art = toc;
    fprintf('\tArtifacts copied: %.2f s\n',t_art)

    % Copy
    max_mem_samples = round(max_mem*1e6/bytes_double);
    cIdx = 1:max_mem_samples:ecog.n_samples; % copy start indices
    n_cIdx = length(cIdx);
    for i = 1:n_cIdx
        start_idx = cIdx(i);
        end_idx = start_idx + max_mem_samples - 1;
        if (end_idx > ecog.n_samples)
            end_idx = ecog.n_samples;
        end
        fprintf('\t(>) start: %i\tend: %i (%.2f %%)\n',start_idx,end_idx,100*i/n_cIdx) % dev

        % copy cavg
        tic;
        cavg = h5read(infilename,'/h5eeg/cavg',[1 start_idx],[1 end_idx-start_idx+1]);
        h5write(outfilename,'/h5eeg/cavg',single(cavg),[1 start_idx],[1 end_idx-start_idx+1])
        t_cavg = toc;
        fprintf('\t\t cavg copied: %.2f s\n',t_cavg);
        
        % copy aux
        for j = 1:n_chan_aux
            tic;
            aux = h5read(infilename,'/h5eeg/aux',[j start_idx],[1 end_idx-start_idx+1]);
            h5write(outfilename,'/h5eeg/aux',single(aux),[j start_idx],[1 end_idx-start_idx+1]);
            t_aux = toc;
            fprintf('\t\t Aux copied: %.2f s, %i of %i, ETA: %.2f min\n',t_aux,j,n_chan_aux,(n_chan_aux-j)*t_aux/60);
        end

        % copy eeg
        for chan = 1:ecog.n_chan
            % Read
            tic;
            eeg = h5read(infilename,'/h5eeg/eeg',[chan start_idx],[1 end_idx-start_idx+1]);
            t_read = toc;
            fprintf('\t\tRead: %.2f MB/s\t',(end_idx-start_idx+1)*bytes_double*1e-6/t_read)

            % Write
            tic;
            h5write(outfilename,'/h5eeg/eeg',single(eeg),[chan start_idx],[1 end_idx-start_idx+1]);
            t_write = toc;
            fprintf('Wrote: %.2f MB/s\tETA: %.2f min\n',...
                (end_idx-start_idx+1)*bytes_double*1e-6/t_write,(t_read+t_write)*(ecog.n_chan-chan)/60)
        end
    end
    
    % Copy events
    annotAll = h5read(infilename,'/h5eeg/events');
    annot = struct('name','FILE_BEGIN','start_idx',1,'duration',0);
    annot.name = {annot.name};
    
    % -- LOW LEVEL COMPOUND DATA TYPE MAGIC ----------------------------------%
    %Create the required data types
    strType = H5T.copy('H5T_C_S1');
    H5T.set_size(strType, 'H5T_VARIABLE');
    sz(1) = H5T.get_size(strType);
    doubleType=H5T.copy('H5T_NATIVE_DOUBLE');
    sz(2) = H5T.get_size(doubleType);
    sz(3) = H5T.get_size(doubleType);
    % Computer the offsets to each field. The first offset is always zero.
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
    fprintf('[!] All Done.\n')
    % -- LOW LEVEL COMPOUND DATA TYPE MAGIC ----------------------------------%

end

fprintf('[!] Done.\n')