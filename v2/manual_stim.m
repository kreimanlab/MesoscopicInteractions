% remove seizures by adding it to art_idx in artifact .h5 files

close all;
clear;

load('annotation_paths2_stim.mat');
ARTIFACT_CODE_SEIZURE = 7;

art_idx_dir = 'art';
art_idx_dir_out = 'art_stim';
trig_write = true;

outL = carveF(art_idx_dir, '.h5');

%Buffer = [0.0625,0.125,0.25,0.5,1,2,4,8,16,32,64,128]; % minutes
Buffer = [45];
marked_a = {};
marked_ai = {};
for i_buf = 1:length(Buffer)

    i_c = 1;
    pct_marked_already_a = [];
    pct_marked_already_ai = [];
    for i = 1:length(outL)
        % build input and output paths
        infn = outL{i};
        fn_art = strsplit(infn,'/');
        fn_art = fn_art{end};
        sid = strsplit(fn_art,'_art.h5');
        sid = sid{1};
        outfn = sprintf('%s/%s',art_idx_dir_out,fn_art);
        fprintf('%s:\n',sid)

        % search for patient in seizure annotations
        sid_idx = find(strcmp(f_sid,sid));
        %if (~ isempty(sid_idx))
        % read existing artifacts
        artifacts = h5read(infn,'/artifacts');
        art_idx = h5read(infn,'/art_idx');
        width = round(h5readatt(infn,'/artifacts','width'));
        w = round(h5readatt(infn,'/art_idx','w'));
        n_arts = round(h5readatt(infn,'/artifacts','n_arts'));
        n_ai = floor(n_arts/w);
        %end

        is_sz_a = false(size(artifacts));
        is_sz_ai = false(size(art_idx));
        % mark annotated seizures
        for j = 1:length(sid_idx)
            j_f = sid_idx(j);
            %fprintf('\t[%s] %s %s\n',f_prefix{j_f},datestr(f_annot_dt{j_f}),f_annot_str{j_f});

            % convert from index to artifacts index
            %fprintf('Buffer: %.2f\n',Buffer(i_buf));
            idx_a = ceil(f_annot_idx(j_f) / width);
            offset_a = floor(Buffer(i_buf)*60/(2)); % floor: odd buffer lengths
            idx_ai = ceil(f_annot_idx(j_f) / (width*w));
            offset_ai = floor(Buffer(i_buf)*60/(w*2)); % floor: odd buffer lengths

            % start and stop indices
            a_start = idx_a - round(0.1*offset_a);
            a_end = idx_a + offset_a;
            ai_start = idx_ai - round(0.1*offset_ai);
            ai_end = idx_ai + offset_ai;

            if (a_start < 1)
                a_start = 1;
            end
            if (ai_start < 1)
                ai_start = 1;
            end

            % get chunk
            a_existing = artifacts(:,a_start:a_end);
            ai_existing = art_idx(:,ai_start:ai_end);
            is_sz_a(:,a_start:a_end) = true;
            is_sz_ai(:,ai_start:ai_end) = true;

            % calculate marked seizures
            pct_marked_already_a(i_c) = sum(a_existing(:) ~= 0)/numel(a_existing);
            pct_marked_already_ai(i_c) = sum(ai_existing(:) ~= 0)/numel(ai_existing);

            i_c = i_c + 1;
        end

        frac_sz_a = sum(is_sz_a(:))/numel(is_sz_a);
        frac_sz_ai = sum(is_sz_ai(:))/numel(is_sz_ai);
        fprintf('\tFraction marked as seizure: %.12f\n',frac_sz_a);
        fprintf('\tFraction marked as seizure (pairwise): %.12f\n',frac_sz_ai);

        if (trig_write)
            artifacts(is_sz_a) = uint8(ARTIFACT_CODE_SEIZURE);
            art_idx(is_sz_a) = uint8(1);
            h5write(outfn,'/artifacts',artifacts);
            h5writeatt(outfn,'/artifacts','frac_sz',frac_sz_a);
            h5write(outfn,'/art_idx',art_idx);
            h5writeatt(outfn,'/art_idx','frac_sz',frac_sz_ai);
        end

    end


    marked_a{i_buf} = pct_marked_already_a;
    marked_ai{i_buf} = pct_marked_already_ai;
end

marked_a_mean = zeros(1,length(marked_a));
marked_ai_mean = zeros(1,length(marked_ai));
marked_a_std = zeros(1,length(marked_a));
marked_ai_std = zeros(1,length(marked_ai));
for i = 1:length(marked_a)
    marked_a_mean(i) = mean(marked_a{i});
    marked_ai_mean(i) = mean(marked_ai{i});
    marked_a_std(i) = std(marked_a{i});
    marked_ai_std(i) = std(marked_ai{i});
end

