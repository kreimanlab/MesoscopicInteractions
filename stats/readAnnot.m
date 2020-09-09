[~, txt, ~] = xlsread("annot_6_11.xlsx");
h5dir = '/mnt/cuenap/data/h5_notch20'; %'/Volumes/cuenap/data/h5_notch20';
stampsDir = '/media/jerry/internal/data/stamps'; %'/Users/jerrywang/Desktop/Annabelle/summer/stamps';

Subjects = unique(txt(:, 1));
subjId = 0;
annot = {};
fmt = 'MM/dd/yy HH:mm:ss';

% DEBUGGING!
subjId = 1;
for i = [19, 21]
% for i = 1:length(txt)

    ctr = 3;
    start = [];
    stop = [];
    fn = [h5dir '/' txt{i, 1} '.h5'];
    ecog = H5eeg([h5dir '/' txt{i, 1} '.h5']);
    while true
        subfile = txt{i, 2};
        if size(txt, 2) > ctr && ~isnat(datetime(txt(i, ctr), 'InputFormat', fmt))
            startDate = datetime(txt(i, ctr), 'InputFormat', fmt);
            startSamp = ecog.getSampleFromTime(stampsDir, startDate);
            stopDate = datetime(txt(i, ctr+1), 'InputFormat', fmt);
            stopSamp = ecog.getSampleFromTime(stampsDir, stopDate);
            start = [start; startSamp];
            stop = [stop; stopSamp];
            ctr = ctr + 2;
        else
            break
        end
    end
    
    if i == 1 || ~strcmp(txt{i, 1}, txt{i-1, 1})
        subjId = subjId + 1;
        annot{subjId, 1} = txt{i, 1};
        annot{subjId, 2} = subfile;
        annot{subjId, 3} = start;
        annot{subjId, 4} = stop;
    else
        annot{subjId, 2} = [annot{subjId, 2}; subfile];
        annot{subjId, 3} = [annot{subjId, 3}; start];
        annot{subjId, 4} = [annot{subjId, 4}; stop];
    end
end

save('annot_temp.mat', 'annot')
for i = 1:length(annot)

    fn = [h5dir '/' annot{i, 1} '.h5'];
    ecog = H5eeg([h5dir '/' txt{i, 1} '.h5']);
    subfile = annot{i, 2};

    freq = h5readatt(fn, '/h5eeg/eeg', 'rate');
    file_start = h5readatt(fn, '/h5eeg', 'datetime');
    len = cumsum(h5readatt(fn, '/h5eeg', 'n_samples'));
    files = h5readatt(fn, '/h5eeg', 'files');
    
    file_ind = zeros(size(files, 1), 1);
    for j = 1:size(files, 1)
        file_ind(j) = sum(strcmp(subfile, files(j)));
    end
    total_len = sum(h5readatt(fn, '/h5eeg', 'n_samples'));
    inds = zeros(total_len, 1);

    for j = 1:length(files)
        if ~file_ind(j)
           if j == 1
               start = 1;
           else
               start = len(j-1)+1;
           end
           inds(start:len(j)) = 2;       
        end
    end
    starts = annot{i, 3};
    stops = annot{i, 4};
    for j = 1:length(starts)
        inds(starts(j):stops(j)) = 1;
    end
    annot{i, 5} = inds;
end

% save('annot_6_11.mat', 'annot', 'txt', '-v7.3')

for i = 1:length(annot)
   a = annot{i, 5};
   plot(a)
   
   b=1;
end

