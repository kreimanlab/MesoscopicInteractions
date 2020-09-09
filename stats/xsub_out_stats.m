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
    dir_art = '/media/klab/KLAB101/h5_notch20/art_nosz'; %'/media/klab/internal/data/h5_notch20/art';
    dir_res = '/media/klab/KLAB101/results/coh_w10';
    dir_cor = '/media/klab/internal/data/coreg'; %'/mnt/cuenap_ssd/coregistration';
    dir_h5 = '/media/klab/KLAB101/h5_notch20';
elseif contains(host,'ubuntu_1604')
    dir_art = '/nas_share/RawData/data/h5_notch20/art_nosz';
    dir_res = '/nas_share/RawData/scripts/opencl/results_res5hz';
    dir_cor = '/mnt/cuenap_ssd/coregistration';
    dir_h5 = '/nas_share/RawData/data/h5_notch20';
end

setenv('SUBJECTS_DIR',dir_cor);

% Patients
Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

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

atl = 2; % atlas index: 1 - Destrieux, 2 - Desikan-Killiany, 6 - HCP-MMP1
n_perm = 10000;
dist_thresh = 17; % mm

% Thresholds
% coh_thresh = 0.3;
% ct_thresh = 0.1;
% cp_thresh = 0.1;

% Thresholds drawn using statistics
p_val = 0.01;
p_val_ct = 0.01;
coh_thresh = NaN;
%ct_thresh = NaN;
ct_thresh = 0.05;
%cp_thresh = NaN;
cp_thresh = 0.1;
art_thresh = 0.25;

% coh_thresh_o = [0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3];
% ct_thresh_o = [0.05 0.1 0.2 0.3];
% cp_thresh_o = [0.05 0.1 0.2 0.3];
% n_params = length(coh_thresh_o)*length(ct_thresh_o)*length(cp_thresh_o);
% all_o = zeros(n_params,3);
% all_i = 1;
% for ii = 1:length(coh_thresh_o)
%     for iii = 1:length(ct_thresh_o)
%         for iv = 1:length(cp_thresh_o)
%             all_o(all_i,1) = coh_thresh_o(ii);
%             all_o(all_i,2) = ct_thresh_o(iii);
%             all_o(all_i,3) = cp_thresh_o(iv);
%             all_i = all_i + 1;
%         end
%     end
% end

% Cache files
cfiles = cell(1,length(Subjects));
for i = 1:length(cfiles)
    cfiles{i} = sprintf('cache/%s_R_null.mat',Subjects{i});
end

% parfor
trig_cache = false;
trig_xsub = true;
curr_d = pwd;
% delete(gcp('nocreate'));
% poolobj = parpool;
% if (trig_cache)
%     addAttachedFiles(poolobj,cfiles);
% end

for metrici = length(metrics) %1:length(metrics)
    
    tic;

    metric = metrics{metrici};
    % coh_thresh = 0.3;
    % ct_thresh = 0.1;
    % cp_thresh = 0.5;
%     coh_thresh = coh_thresh_o(ii);
%     ct_thresh = ct_thresh_o(iii);
%     cp_thresh = cp_thresh_o(iv);
%     coh_thresh = all_o(v,1);
%     ct_thresh = all_o(v,2);
%     cp_thresh = all_o(v,3);

    % Load atlas
    C1 = load(sprintf('%s/%s/label/all_parcellation.mat',dir_cor,Subjects{1}));
    atl_name = C1.AtlNames{atl};
    n_rois = C1.AtlROIs{atl}.LH.numEntries;
    rois = C1.AtlROIs{atl}.LH.struct_names;
    n_comb_atl = nchoosek(n_rois,2);
    AdjAtl = cell(n_rois,n_rois);
    AdjAtlN = cell(n_rois,n_rois);
    AdjAtl_sid = cell(n_rois,n_rois);
%     AdjAtl_xyz = cell(n_rois,n_rois);
    adjct_dist = cell(n_rois,n_rois);
    %AdjAtlN_sid = cell(n_rois,n_rois);

    % Print header
    %fprintf('pt\tmetric\tcoh_thresh\tct_thresh\tisx\tisxN\n')

    % Main loop
    for i = 1:length(Subjects)
        
        %tic;
        subject_i = i;
        sid = Subjects{i};
        fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
        fn_dist = sprintf('%s/%s_dists-%s-%i.mat',dir_res,sid,metric,n_perm);
        fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
        fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
        fn_coreg = sprintf('%s/%s/label/all_parcellation.mat',dir_cor,sid);

        % Check if files exist
        ckf = {fn_art,fn_dist,fn_graph,fn_h5,fn_coreg};
        for j = 1:length(ckf)
            if (~exist(ckf{j},'file'))
                fprintf(2,'E> File not found: %s\n',ckf{j});
                %return
            end
        end

        % Init H5eeg
        ecog = H5eeg(fn_h5);
        try
            Dmat = ecog.getDist2();
        catch
            fprintf(2,'[!] Warning: older .h5 format found. Continuing with no ground-truth neighbor detection.\n')
            Dmat = ecog.getDist();
        end
        
        Dmats = zeros(ecog.n_comb,1);
        count = 1;
        for i1 = 1:(ecog.n_bchan-1)
            for i2 = (i1+1):ecog.n_bchan
                Dmats(count) = Dmat(i1,i2);
                count = count + 1;
            end
        end

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
        %t_load = toc;
        %fprintf('[%s] Loading took: %.3f seconds\n',sid,t_load)
        
        % Filter out high coherences
        high_coh_thresh = 1; % or not
        frac_high_coh = sum(R(:) > high_coh_thresh) / numel(R);
        R(R > high_coh_thresh) = -1;

        D = load(fn_dist);
        
        % --- Prepare null data ---
        if (~ trig_cache)
            
            R_null = zeros(size(R));
            %D = load(fn_dist);
            for j = 1:ecog.n_comb
                R_null(j,:) = random(D.d{j},[1 n_graph]);
            end
            
        else
            
            fn_tmp = sprintf('%s/cache/%s_R_null.mat',curr_d,sid);
            if (~exist(fn_tmp,'file'))
                R_null = zeros(size(R));
                %D = load(fn_dist);
                for j = 1:ecog.n_comb
                    R_null(j,:) = random(D.d{j},[1 n_graph]);
                end
                save_tmp(sid,single(R_null));
            else
                try
                    R_null = load_tmp(fn_tmp);
                    R_null = R_null.R_null;
                catch
                    % Re-sample and save if file is corrupt
                    R_null = zeros(size(R));
                    for j = 1:ecog.n_comb
                        R_null(j,:) = random(D.d{j},[1 n_graph]);
                    end
                    save_tmp(sid,single(R_null));
                end
            end
        
        end
        %t_rand = toc;
        %fprintf('[%s] Random sampling took: %.3f seconds\n',sid,t_rand-t_load)

        % ----------------------------------------------------------------------
        % Compute coherence thresholds
        coh_thresh = zeros(ecog.n_comb,1);
        for i3 = 1:ecog.n_comb
            % Multiple Hypothesis Test Correction: Bonferroni every hour
            bonf_seconds = 1*60;%3600/2;
            coh_thresh(i3) = icdf(D.d{i3},1 - p_val/(bonf_seconds/w));
            
            if (coh_thresh(i3) > 1)
                coh_thresh(i3) = 1;
            end
        end
        
        % Compute consistency across time thresholds
        %ct_thresh = binoinv(1-p_val_ct,n_graph,p_val/(bonf_seconds/w))/n_graph;
        %ct_thresh = binoinv(1-(p_val_ct/ecog.n_comb),n_graph,p_val/(bonf_seconds/w))/n_graph;
        
        
        
        %%
        % SLEEP
        %
        % R - coherences over time. Script ignores NaNs
        %
        %%
        
        
        % --- Build CT ---
        % distance threshold
        d_i = Dmats > dist_thresh;
        
        n_yes = sum((R >= coh_thresh) & (~art_idxB),2); % broadcast thresholding, future releases of matlab may not support
        n_no = sum((R < coh_thresh) & (~art_idxB),2);
        ct = n_yes./(n_yes + n_no);
        is_x = (ct >= ct_thresh);
        frac_x = sum(is_x(d_i))/numel(is_x(d_i));

        % artifact fraction threshold
        not_art = ((n_no + n_yes) ./ n_graph);
        not_art_threshed = ((n_no + n_yes) ./ n_graph) > art_thresh;

        
        
        %R4m = double(R);
        R4m = R;
        R4m(R < coh_thresh) = NaN; % --- Comment this line to include all points in average magnitude ---
        R4m(art_idxB) = NaN;
        mag = double(nanmean(R4m,2));
        mag_std = double(nanstd(R4m')');

        n_yes_null = sum((R_null >= coh_thresh) & (~art_idxB),2);
        n_no_null = sum((R_null < coh_thresh) & (~art_idxB),2);
        ct_null = n_yes_null./(n_yes_null + n_no_null);
        is_x_null = (ct_null >= ct_thresh);
        frac_x_null = sum(is_x_null(d_i))/numel(is_x_null(d_i));
        
        
        
        %R4mN = double(R_null);
        R4mN = R_null;
        R4mN(R_null < coh_thresh) = NaN; % --- Comment this line to include all points in average magnitude ---
        R4mN(art_idxB) = NaN;
        mag_null = double(nanmean(R4mN,2));
        mag_null_std = double(nanstd(R4m')');
        %fprintf('%i\t%s\t%s\t%.8f\t%.8f\t%.8f\t%.8f\n',v,sid,metric,coh_thresh,ct_thresh,frac_x,frac_x_null);
        %fprintf('%s\t%s\t%.8f\t%.8f\t%.8f\t%.8f\n',sid,metric,coh_thresh,ct_thresh,frac_x,frac_x_null);
        fprintf('%s\t%s\t%.8f\t%.8f\t%.8f\t%.8f\n',sid,metric,median(coh_thresh),ct_thresh,frac_x,frac_x_null);

        
        % --- DEBUG ---
%         i_rnd = randi([1 ecog.n_comb]);
%         subplot(2,1,1)
%         Rplt = R(i_rnd,:);
%         Rplt(art_idxB(i_rnd,:)) = NaN;
%         R0plt = R_null(i_rnd,:);
%         R0plt(art_idxB(i_rnd,:)) = NaN;
%         plot(Rplt,'black.'); hold on; plot([0 n_graph],[1 1]*coh_thresh(i_rnd),'--'); axis([0 n_graph 0 1]);
%         title(sprintf('%.2f mm',Dmats(i_rnd)));
%         %hist(ct(:),50);
%         subplot(2,1,2)
%         plot(R0plt,'black.'); hold on; plot([0 n_graph],[1 1]*coh_thresh(i_rnd),'--'); axis([0 n_graph 0 1]);
%         %hist(ct_null(:),50);
%         return
        % --- END DEBUG ---

        if ((~ strcmp(sid,'mSu')) && trig_xsub)
            % Map ct onto atlas
            b1c1 = zeros(ecog.n_comb,1);
            b1c2 = zeros(ecog.n_comb,1);
            b2c1 = zeros(ecog.n_comb,1);
            b2c2 = zeros(ecog.n_comb,1);
            for k = 1:ecog.n_comb
                if (d_i(k)) % only map if distance threshold is met
                    b1c1(k) = ecog.bip(chan1(k)+1,1);
                    b1c2(k) = ecog.bip(chan1(k)+1,2);
                    b2c1(k) = ecog.bip(chan2(k)+1,1);
                    b2c2(k) = ecog.bip(chan2(k)+1,2);

                    % Get ROI names
                    roi_b1c1 = atl_labels{b1c1(k)};
                    roi_b1c2 = atl_labels{b1c2(k)};
                    roi_b2c1 = atl_labels{b2c1(k)};
                    roi_b2c2 = atl_labels{b2c2(k)};
                    
                    % Get ROI hemispheres
                    roi_b1c1_hemi = atl_labels_hemi{b1c1(k)};
                    roi_b1c2_hemi = atl_labels_hemi{b1c2(k)};
                    roi_b2c1_hemi = atl_labels_hemi{b2c1(k)};
                    roi_b2c2_hemi = atl_labels_hemi{b2c2(k)};
                    
                    % Check that every roi is on same hemisphere
                    hemi_list = {roi_b1c1_hemi,roi_b1c2_hemi,roi_b2c1_hemi,roi_b2c2_hemi};
                    isSameHemi = (length(unique(hemi_list)) == 1);
                    
                    if (~isSameHemi)
                        fprintf(2,'[*] skipped xhemi roi pair: %s bip pair %i,%i\n',sid,chan1(k)+1,chan2(k)+1)
                    else

                        % Get ROI indices
                        i_b1c1 = find(strcmp(rois,roi_b1c1),1);
                        i_b1c2 = find(strcmp(rois,roi_b1c2),1);
                        i_b2c1 = find(strcmp(rois,roi_b2c1),1);
                        i_b2c2 = find(strcmp(rois,roi_b2c2),1);
                        i_b1 = [i_b1c1,i_b1c2];
                        i_b2 = [i_b2c1,i_b2c2];

                        % Get value to write to AdjAtl
                        if (is_x(k))
                            el = mag(k);
                        else
                            el = zeros(1,1);
                        end

                        if (is_x_null(k))
                            el_null = mag_null(k);
                        else
                            el_null = zeros(1,1);
                        end

                        % Map onto atlas adjacency matrix
                        try
                            for k1 = 1:1%length(i_b1)
                                for k2 = 1:1%length(i_b2)
                                    ki1 = i_b1(k1);
                                    ki2 = i_b2(k2);

                                    AdjAtl{ki1,ki2} = [AdjAtl{ki1,ki2}; el];
                                    AdjAtl{ki2,ki1} = [AdjAtl{ki2,ki1}; el];

                                    AdjAtlN{ki1,ki2} = [AdjAtlN{ki1,ki2}; el_null];
                                    AdjAtlN{ki2,ki1} = [AdjAtlN{ki2,ki1}; el_null];

                                    AdjAtl_sid{ki1,ki2} = [AdjAtl_sid{ki1,ki2},subject_i];
                                    AdjAtl_sid{ki2,ki1} = [AdjAtl_sid{ki2,ki1},subject_i];

                                    adjct_dist{ki1,ki2} = [adjct_dist{ki1,ki2}; Dmats(k)];
                                    adjct_dist{ki2,ki1} = [adjct_dist{ki2,ki1}; Dmats(k)];

    %                                 AdjAtl_xyz{ki1,ki2} = [AdjAtl_xyz{ki1,ki2}; ecog.bip()];
    %                                 AdjAtl_xyz{ki2,ki1} = [AdjAtl_xyz{ki2,ki1}; ecog.bip()];
                                end
                            end

                        catch
                            %fprintf('ct-atl map skip.\n');
                        end
                    end % check hemi
                end
            end
        end
        
        % build adjacency matrices
        AdjMagRaw = nan(ecog.n_bchan,ecog.n_bchan);
        AdjCTRaw = nan(ecog.n_bchan,ecog.n_bchan);
        AdjIsSig = nan(ecog.n_bchan,ecog.n_bchan);
        AdjIsDistOk = nan(ecog.n_bchan,ecog.n_bchan);
        AdjArtFraC = nan(ecog.n_bchan,ecog.n_bchan);
        k0 = 1;
        for i0 = 1:(ecog.n_bchan-1)
            for j0 = (i0 + 1):ecog.n_bchan
                AdjMagRaw(i0,j0) = mag(k0);
                AdjMagRaw(j0,i0) = mag(k0);
                AdjCTRaw(i0,j0) = ct(k0);
                AdjCTRaw(j0,i0) = ct(k0);
                AdjIsSig(i0,j0) = is_x(k0);
                AdjIsSig(j0,i0) = is_x(k0);
                AdjIsDistOk(i0,j0) = d_i(k0);
                AdjIsDistOk(j0,i0) = d_i(k0);
                AdjArtFraC(i0,j0) = 1 - not_art(k0);
                AdjArtFraC(j0,i0) = 1 - not_art(k0);
                k0 = k0 + 1;
            end
        end
        
        AdjCT = AdjCTRaw;
        AdjCT(AdjIsSig~=1) = 0;
        AdjCT(AdjIsDistOk~=1) = NaN;
        AdjMag = AdjMagRaw;
        AdjMag(AdjIsSig~=1) = 0;
        AdjMag(AdjIsDistOk~=1) = NaN;
        %return


        % Predict ETA
        %t_sing = toc;
        %fprintf('%s\tETA: %.2f hrs\n',sid,t_sing*(length(Subjects) - i)/3600)
        
        % Save individual subjects
        save(sprintf('cache/xsub_out_%s_%i',sid,metrici),'-v6','-regexp','^(?!(AdjAtl|AdjAtlN|AdjAtl_sid|R|R_null|R4m|R4mN|D|d_i|art_idx|art_idxB)$).');

    end
    
    %if (~strcmp(sid,'mSu'))

    % xsub
    cp = zeros(n_comb_atl,1);
    cp_null = zeros(n_comb_atl,1);
    count = 1;
    for i = 1:(n_rois-1)
        for j = (i+1):n_rois
            l = AdjAtl{i,j};
            l_null = AdjAtlN{i,j};
            cp(count) = sum(l ~= 0)/numel(l);
            cp_null(count) = sum(l_null ~= 0)/numel(l_null);
            count = count + 1;
        end
    end
    
    for cp_thresh = [0.001 0.01 0.03 0.05 0.07 0.1 0.15 0.2 0.25 0.3 0.35 0.4]
    
        frac_y = sum(cp > cp_thresh) / sum(~isnan(cp));
        frac_y_null = sum(cp_null > cp_thresh) / sum(~isnan(cp_null));
        fprintf('%s\t%s\t%.8f\t%.8f\t%.8f\n','xsub',metric,cp_thresh,frac_y,frac_y_null);

    end

    % Save final result
    save(sprintf('cache/xsub_out_all_%i',metrici),'-v6','-regexp','^(?!(R|R_null|R4m|R4mN|D|d_i|art_idx|art_idxB)$).');

    %end
    %return;
    
    % Clear loop indices
    %clear i;
    %clear j;
    
    %t_sing = toc;
    %fprintf('Est. total time: %.2f hrs\n',t_sing * n_params * (1/3600));
end

% Print finish message
fprintf('Done.\n')


function save_tmp(sid, R_null)
    save(sprintf('cache/%s_R_null',sid),'R_null','-v6');
end

function R_null = load_tmp(R_null_fn)
    R_null = load(R_null_fn);
end
