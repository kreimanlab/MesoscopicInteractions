function [ ] = annotate_opti_fn( sid )

vers = 2;

[~,host] = system('hostname');
if contains(host,'ubuntu_1604')
    h5dir = '/nas_share/RawData/scripts/synth/out';
    outdir = '/nas_share/RawData/scripts/synth/out_art';
elseif contains(host,'hopper')
    %h5dir = '/mnt/cuenap2/scripts/synth/out';
    %outdir = '/mnt/cuenap2/scripts/synth/out_art';
    h5dir = '/media/klab/KLAB101/h5_notch20';
    outdir = '/media/klab/KLAB101/h5_art_opti';
elseif contains(host,'o2.rc.hms.harvard.edu')
    %h5dir = '/n/groups/kreiman/jerry/data/h5';
    h5dir = '/n/scratch2/jw324/data/h5_notch20';
    %outdir = '/n/groups/kreiman/jerry/data/h5_art';
    outdir = '/n/groups/kreiman/jerry/data/h5_art_opti';
end

system(['mkdir ',outdir]);

% Get list of .h5 files
% if (true)
    %H5 = {'sub40.h5'};
H5 = {[sid,'.h5']};
% else
% 
%     d = dir(h5dir);
%     dname = {d.name};
%     disdir = [d.isdir];
%     d = dname(~disdir);
%     H5 = {};
%     j = 1;
%     for i = 1:length(d)
%         if (endsWith(d{i},'.h5'))
%             H5{j} = d{i};
%             j = j + 1;
%         end
%     end
% 
% end

%for i = 1:length(H5)
i = 1;
sid = strsplit(H5{i},'.h5');
sid = sid{1};
h5fn = sprintf('%s/%s',h5dir,H5{i});
fprintf('%s\n',h5fn);
ecog = H5eeg(h5fn);
bip = h5readatt(h5fn,'/h5eeg/eeg','bip');
[n_bchan,~] = size(bip);
Fs = round(ecog.fs);
outfn = sprintf('%s/%s_art.h5',outdir,sid);

% Multiple parameters
X1_t1 = [10,25,50,75,100];
X1_t2 = [1000,1250,1500,1750,2000,2500,3000];
X2_t = [200,250,300,350,400,500]/4;
n_opti = length(X1_t1)*length(X1_t2)*length(X2_t);
%x1_t1 = 25; % uV
%x1_t2 = 1000; % uV
%x2_t = 300/4; % uV / msec

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
    system(['rm -f ',outfn]);
end
h5create(outfn,'/artifacts',[n_bchan*n_opti length(Starts)],'Datatype','uint8')
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
    A = uint8(zeros(n_bchan*n_opti,1));
   
    for k = 1:n_bchan
        V = eeg.data(:,bip(k,1)) - eeg.data(:,bip(k,2));
        V = V - mean(V);
        X = features_v2_l(V, Fs);

        c_i = 1;
        for i1 = 1:length(X1_t1)
            for i2 = 1:length(X1_t2)
                for i3 = 1:length(X2_t)

                    Y = classify_v2_l(X,X1_t1(i1),X1_t2(i2),X2_t(i3));
                    A(k + (c_i - 1)*n_bchan) = uint8(Y);
                    c_i = c_i + 1;
                    
                end
            end
        end
        
    end

    
    %return
    h5write(outfn,'/artifacts',A,[1 j],[n_bchan*n_opti 1]);
    
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

%end


function [ X ] = features_v2_l( V, Fs )
% V - 1 by n vector

X(1) = max(V) - min(V); % uV
X(2) = max(abs(diff(V))) / (1e3/Fs); % uV / msec

end

function [ Y ] = classify_v2_l( X, x1_t1, x1_t2, x2_t )

%x1_t1 = 25; % uV
%x1_t2 = 1000; % uV
%x2_t = 300/4; % uV / msec

%x3_t = 0.75; % ratio [0,1]
%x3_t = 0.6;
%x4_t = 0.4;

if (X(1) < x1_t1)
    Y = 1;
elseif (X(1) > x1_t2)
    Y = 2;
elseif (X(2) > x2_t)
    Y = 3;
% elseif (X(3) > x3_t)
%     Y = 4;
% elseif (X(4) > x4_t)
%     Y = 5;
else
    Y = 0;
end

end


