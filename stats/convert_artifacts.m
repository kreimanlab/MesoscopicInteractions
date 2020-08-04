close all;
clear;

% Convert art_idx

dir_h5_in = '/media/klab/KLAB101/h5_notch20/art_nosz';
dir_h5_out = '/media/klab/KLAB101/h5_notch20/art_nosz2';
dir_h5 = '/media/klab/KLAB101/h5_notch20';

h5_in_filenames = carveF(dir_h5_in,'.h5');
h5_out_filenames = carveF(dir_h5_out,'.h5');

for i = 1:length(h5_in_filenames)
    fn_in = h5_in_filenames{i};
    fn_out = h5_out_filenames{i};
    sid = fn_in((end-12):(end-7));
    fprintf('\tIn:  %s\n',fn_in);
    
    % Load
    artifacts = h5read(fn_in,'/artifacts');
    artifacts_w = h5readatt(fn_in,'/artifacts','width');
    n_arts = h5readatt(fn_in,'/artifacts','n_arts');
    art_idx = h5read(fn_in,'/art_idx');
    art_idx_w = h5readatt(fn_in,'/art_idx','w');
    ecog = H5eeg(sprintf('%s/%s.h5',dir_h5,sid));
    
    % First, figure out what fraction are sz
    art_idx_notsz = zeros(size(art_idx));
    art_idx_new = zeros(size(art_idx));
    art_idx_new2 = zeros(size(art_idx));
    i3 = 1;
    for i1 = 1:(ecog.n_bchan)
        for i2 = (i1 + 1):ecog.n_bchan
            a_raw_b1 = artifacts(i1,:);
            a_raw_b2 = artifacts(i2,:);
            ai_raw = art_idx(i3,:);
            ai_notsz = zeros(size(ai_raw));
            ai_new = zeros(size(ai_raw));
            ai_new2 = zeros(size(ai_raw));
            starts = 1:art_idx_w:(n_arts-art_idx_w);
            parfor j = 1:length(starts)
                istart = starts(j);
                istop = istart + art_idx_w - 1;
                % For each 20-second segment
                a1 = a_raw_b1(istart:istop) > 0;
                a2 = a_raw_b2(istart:istop) > 0;
                is_ai = any(a1) || any(a2);
                is_ai_new = sum([a1, a2]) > 1;
                is_ai_new2 = sum([a1, a2]) > 2;
                is_sz = (ai_raw(j) && (~ is_ai ));
                ai_notsz(j) = ~ is_sz;
                ai_new(j) = is_ai_new || is_sz;
                ai_new2(j) = is_ai_new2 || is_sz;
            end
            %return
            art_idx_notsz(i3,:) = ai_notsz;
            art_idx_new(i3,:) = ai_new;
            art_idx_new2(i3,:) = ai_new2;
            i3 = i3 + 1;
        end
    end
    return
    
    fprintf('\tOut: %s\n',fn_out);
end