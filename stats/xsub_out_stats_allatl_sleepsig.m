clear;
close all;
rng shuffle;

[~,host] = system('hostname');
if contains(host,'hopperu')
    % Path definitions
    dir_art = '/mnt/cuenap/data/scripts/v2/art_szs';%'/media/klab/KLAB101/h5_notch20/art_nosz';%'/mnt/cuenap/data/h5_notch20/art_nosz';
    dir_res = '/mnt/cuenap/data/results/coh_w20';
    dir_cor = '/mnt/cuenap_ssd/coregistration';
    dir_h5 = '/mnt/cuenap/data/h5_notch20';
    %dir_out = '/';
elseif contains(host,'o2.rc.hms.harvard.edu')
    % Path definitions
    dir_art = '/n/groups/kreiman/jerry/data/h5/art';
    dir_res = '/n/scratch2/jw324/data/results/coh_w10';
    dir_cor = '/n/scratch2/jw324/data/coreg';
    dir_h5 = '/n/groups/kreiman/jerry/data/h5';
elseif contains(host,'kraken')
    dir_art = '/media/klab/untitled/h5_notch20/art_nosz'; %'/media/klab/internal/data/h5_notch20/art';
    dir_res = '/media/klab/untitled/results/coh_w10';
    dir_cor = '/media/klab/internal/data/coreg'; %'/mnt/cuenap_ssd/coregistration';
    dir_h5 = '/media/klab/untitled/h5_notch20';
elseif contains(host,'ubuntu_1604')
    dir_art = '/nas_share/RawData/data/h5_notch20/art_nosz';
    dir_res = '/nas_share/RawData/scripts/opencl/results_res5hz';
    dir_cor = '/mnt/cuenap_ssd/coregistration';
    %dir_h5 = '/nas_share/RawData/data/h5_notch20';
    dir_h5 = '/mnt/cuenap/data/h5_notch20';
end

setenv('SUBJECTS_DIR',dir_cor);

% Patients
Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

SubjectsSleep = {'m00019','m00023','m00024','m00026','m00030','m00032','m00035',...
    'm00039','m00043','m00044','m00049','m00079','m00083','m00084'};

%Subjects = {'mSu'};

% Randomly rearrange subjects
%Subjects = Subjects(randperm(length(Subjects)));

% Exclude monkey
% Subjects = Subjects(1:(end-1));
% Subjects = Subjects(1:2);
% Subjects{3} = 'mSu';

% Constants
metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'}; % 'pcDelta'

% Randomly rearrange metrics
%metrics = metrics(randperm(length(metrics)));

% atlas index: 1 - Destrieux, 2 - Desikan-Killiany, 6 - HCP-MMP1
% Update: 1 for 150 area parcellation

fprintf('[*] Loading ../sleep/annot_all.mat\n');
CaSleep = load('../sleep/annot_all.mat');
fprintf('\tDone.\n');

%%
trig_permutation = false;

atl = 2;

n_perm = 10000;
for metrici = [1 5]
    metric = metrics{metrici};

    n_samples_all = nan(length(SubjectsSleep),1);
    for i = 1:length(SubjectsSleep)
        sid = SubjectsSleep{i};
        sid2 = CaSleep.annot{i,1};
        sid3 = CaSleep.Subjs{i};

        % check sid match
        if (~strcmp(sid,sid2))
            fprintf(2,'[!] Warn. mismatched subject strings: %s, %s\n',sid,sid2);
        end
        if (~strcmp(sid,sid3))
            fprintf(2,'[!] Warn. mismatched subject strings: %s, %s\n',sid,sid3);
        end

        % load H5eeg
        fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
        ecog = H5eeg(fn_h5);
        n_samples_all(i) = ecog.n_samples;

        % Chunk behaviors
        label_sleep_raw = CaSleep.annot{i,5};
        w_seconds = 10; % seconds
        w = round(ecog.fs)*w_seconds;
        starts = (1:w:ecog.n_samples);
        starts = starts(1:(end-1));
        label_sleep = nan(size(starts));
        c = 1;
        for j = 1:length(starts)
            idx_start = starts(j);
            idx_end = idx_start + w - 1;
            seg_sleep = label_sleep_raw(idx_start:idx_end);

            % check that the chunk is valid
            isvalid = (sum(mode(seg_sleep) == seg_sleep) == length(seg_sleep));
            if (isvalid)
                label_sleep(c) = mode(seg_sleep);
            end
            c = c + 1;
        end
        % trim last
        label_sleep = label_sleep(1:(end-1));


        % Load
        % --- Begin copypasta ---------------------------------------------
        fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
        fn_dist = sprintf('%s/%s_dists-%s-%i.mat',dir_res,sid,metric,n_perm);
        fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
        fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
        fn_coreg = sprintf('%s/%s/label/all_parcellation.mat',dir_cor,sid);
        
        % --- Load graph and artifacts ---
        if (~strcmp(sid,'mSu'))
            C = load(fn_coreg);
            atl_labels = C.AtlLabels{atl};
            atl_labels_hemi = C.EleHemi;
        else
            % special case for mSu
            C = load(sprintf('%s/%s/label/all_parcellation.mat',dir_cor,Subjects{1}));
            atl_labels = C.AtlLabels{atl};
            atl_labels_hemi = C.EleHemi;
        end
        art_idx = h5read(fn_art,'/art_idx');
        art_idxB = (art_idx == 1);
        frac_art = h5readatt(fn_art,'/art_idx','frac_art');
        w_art = h5readatt(fn_art,'/art_idx','w');
        % trim last time sample of graph
        R = h5read(fn_graph,'/R',[1 1],size(art_idx));
        
        if (isa(R,'uint16'))
            R = single(R) ./ (2^16 - 1);
            if ( (max(R(:)) > 1) || (min(R(:)) < 0) )
                fprintf(2,'WARN: graph file values not in range 0 - 1.\n');
            end
        end
        
        w = double(h5read(fn_graph,'/w'));
        [~,n_graph] = size(R);
        chan1 = h5read(fn_graph,'/chan1');
        chan2 = h5read(fn_graph,'/chan2');
        
        % Filter out high coherences
        high_coh_thresh = 1; % or not
        frac_high_coh = sum(R(:) > high_coh_thresh) / numel(R);
        R(R > high_coh_thresh) = -1;

        D = load(fn_dist);
        % --- end copypasta ---------------------------------------------
        
        
        % load xsub_out_all
        Ca = load(sprintf('./cache/xsub_out_%s_%i_atl%i.mat',sid,metrici,atl));
        
        % p-values
        % sleep labels: label_sleep 1x30240
        % coherence: R 1953x30239
        [n_pairs,~] = size(R);
        pvals = nan(n_pairs,1);
        pvals_perm = nan(n_pairs,1);
        ns_sleep = nan(n_pairs,1);
        ns_wake = nan(n_pairs,1);
        means_sleep = nan(n_pairs,1);
        stds_sleep = nan(n_pairs,1);
        means_wake = nan(n_pairs,1);
        stds_wake = nan(n_pairs,1);
        ranksums = nan(n_pairs,1);
        obsdiff = nan(n_pairs,1);
        effsize = nan(n_pairs,1);
        for k = 1:n_pairs
            R1 = R(k,:);
            
            % remove artifacts
            R1(art_idxB(k,:)) = NaN;
            
            % remove sub-threshold coherences
            C1 = Ca.coh_thresh(k);
            L0 = label_sleep;
%             if (length(C1) > length(R1))
%                 C1 = C1(1:(length(R1)));
%             elseif (length(C1) < length(R1))
%                 R1 = R1(1:(length(C1)));
%                 L0 = L0(1:(length(C1)));
%             end
            R1(R1 < C1) = NaN;
            
            Rsleep = R1(L0 == 1);
            Rwake = R1(L0 == 0);
            ns_sleep(k) = sum(~isnan(Rsleep));
            ns_wake(k) = sum(~isnan(Rwake));
            means_sleep(k) = nanmean(Rsleep);
            stds_sleep(k) = nanstd(Rsleep);
            means_wake(k) = nanmean(Rwake);
            stds_wake(k) = nanstd(Rwake);
            
            % ranksum test
            X1 = Rsleep(~isnan(Rsleep));
            X2 = Rwake(~isnan(Rwake));
            try
                [p,h,stats] = ranksum(double(X1),double(X2));
                ranksums(k) = stats.ranksum;
                pvals(k) = p;
            catch
                ranksums(k) = nan;
                pvals(k) = nan;
            end
            
            fprintf('[%i/%i %s] %i of %i\n',i,length(SubjectsSleep),sid,k,n_pairs);
            
            if (trig_permutation)
                [p, observeddifference, effectsize] = permutationTest(X1, X2, 1e7);
                pvals_perm(k) = p;
                obsdiff(k) = observeddifference;
                effsize(k) = effectsize;
            end
            %return
        end
        save(sprintf('./cache/%s_sleepsig_metric-%i_atl-%i',sid,metrici,atl),'-v6','-regexp','^(?!(CaSleep|R|R_null|R4m|R4mN|D|d_i|art_idx|art_idxB)$).');
        %return
    end % for sleep subjects
end % for frequency band


fprintf('All Done.\n')


function save_tmp(sid, R_null, atl, metrici)
    save(sprintf('cache/%s_R_null_atl%i_%i',sid,atl,metrici),'R_null','-v6');
end

function R_null = load_tmp(R_null_fn)
    R_null = load(R_null_fn);
end
