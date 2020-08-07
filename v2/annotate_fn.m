function [ ] = annotate_fn( sid )


vers = 2;

[~,host] = system('hostname');
if contains(host,'ubuntu_1604')
    h5dir = '/nas_share/RawData/scripts/synth/out';
    outdir = '/nas_share/RawData/scripts/synth/out_art';
elseif contains(host,'hopper')
    h5dir = '/mnt/cuenap2/scripts/synth/out';
    outdir = '/mnt/cuenap2/scripts/synth/out_art';
elseif contains(host,'o2.rc.hms.harvard.edu')
    %h5dir = '/n/groups/kreiman/jerry/data/h5';
    h5dir = '/n/scratch2/jw324/data/h5_notch20';
    %outdir = '/n/groups/kreiman/jerry/data/h5_art';
    outdir = '/n/groups/kreiman/jerry/data/h5_art_50uv_2';
end

system(['mkdir ',outdir]);

% Get list of .h5 files
if (true)
    %H5 = {'m00083.h5'};
    H5 = {[sid,'.h5']};
else

    d = dir(h5dir);
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

end

for i = 1:length(H5)
    sid = strsplit(H5{i},'.h5');
    sid = sid{1};
    h5fn = sprintf('%s/%s',h5dir,H5{i});
    fprintf('%s\n',h5fn);
    ecog = H5eeg(h5fn);
    bip = h5readatt(h5fn,'/h5eeg/eeg','bip');
    [n_bchan,~] = size(bip);
    Fs = round(ecog.fs);
    outfn = sprintf('%s/%s_art.h5',outdir,sid);

    % Optimize chunk
    fid = H5F.open(h5fn);
    dset_id = H5D.open(fid,'/h5eeg/eeg');
    dcpl = H5D.get_create_plist(dset_id);
    [~,chunk_dims] = H5P.get_chunk(dcpl);
    H5D.close(dset_id);
    H5F.close(fid);
    
    % Check that annotation width falls neatly within chunk size
    if (mod(chunk_dims(1),Fs) ~= 0)
        fprintf(2,'E> Chunk size divided by width is not an integer\n');
    end
    % Read chunk for cache
    %Cache = zeros(chunk_dims);
    
    % Make starting and ending indices
    Starts = 1:Fs:ecog.n_samples;
    Ends = Starts + Fs - 1;
    Ends(end) = ecog.n_samples;
    Starts = Starts(1:(end-1));
    Ends = Ends(1:(end-1));
    
    % Make h5 file (delete existing)
    if (exist(outfn,'file'))
        system(['rm -f',outfn]);
    end
    h5create(outfn,'/artifacts',[n_bchan length(Starts)],'Datatype','uint8')
    h5writeatt(outfn,'/artifacts','width',Fs);
    h5writeatt(outfn,'/artifacts','n_bchan',n_bchan);
    h5writeatt(outfn,'/artifacts','n_arts',length(Starts));
    h5writeatt(outfn,'/artifacts','version',vers);
    
    % Annotate
    % Read cache on first go
%     cache_start = 1;
%     cache_end = cache_start + chunk_dims(1) - 1;
%     tic;
%     eeg = ecog.readEEG({cache_start cache_end});
%     t_read = toc;
%     Cache = eeg.data;
%     fprintf('Read: %.2f MB/s\n',(prod(chunk_dims)*4*(1e-6))/t_read);
    %return
    for j = 1:length(Starts)
        % If index is not in cache, read the appropriate chunk
        % NOTE: if Starts(j) to Ends(j) spans 2 chunks, this will break
%         if ~((Starts(j) >= cache_start) && (Ends(j) <= cache_end))
%             cache_start = (ceil(Starts(j)/chunk_dims(1)) - 1)*chunk_dims(1) + 1;
%             cache_end = cache_start + chunk_dims(1) - 1;
%             tic;
%             eeg = ecog.readEEG({cache_start cache_end});
%             t_read = toc;
%             Cache = eeg.data;
%             fprintf('Read: %.2f MB/s\n',(prod(chunk_dims)*4*(1e-6))/t_read);
%         else
%             tic;
%         end

        tic;
        eeg = ecog.readEEG({Starts(j) Ends(j)});

%         start_i = Starts(j) - cache_start + 1;
%         end_i = start_i + (Ends(j)-Starts(j)+1) - 1;
%         eegdata = Cache(start_i:end_i,:);
%         tic;
%         eeg = ecog.readEEG({Starts(j) Ends(j)});
%         t_read = toc;
        
        A = uint8(zeros(n_bchan,1));
        for k = 1:n_bchan
            V = eeg.data(:,bip(k,1)) - eeg.data(:,bip(k,2));
            V = V - mean(V);
            X = features_v2(V, Fs);
            Y = classify_v2(X);
            A(k) = uint8(Y);
        end
        
        % Write
        h5write(outfn,'/artifacts',A,[1 j],[n_bchan 1]);
        
        % For testing
        %h5read(outfn,'/artifacts',[1 1],[n_bchan 10])
        
        % Print
        t_sing = toc;
        eta = (length(Starts) - j) * t_sing;
        size_mb = ecog.n_chan * (Ends(j) - Starts(j)) * 32 * (1/8) * (1e-6);
        fprintf('%s: %.6f %% ETA: %.2f hr\n',...
            sid,100*j/length(Starts),eta/3600);
        
        
%         if (j == 100)
%             return
%         end
    end
end

end
