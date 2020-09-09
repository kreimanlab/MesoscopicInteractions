close all;
clear;
rng shuffle;


[~,host] = system('hostname');
if contains(host,'hopper')
    % Path definitions
    dir_art = '/media/klab/internal/data/h5_notch20/art';
    dir_res = '/media/klab/internal/data/results/coh_w10';
    dir_cor = '/media/klab/internal/data/coreg';
    dir_h5 = '/media/klab/KLAB101/h5_notch20';
elseif contains(host,'o2.rc.hms.harvard.edu')
    % Path definitions
    dir_art = '/n/groups/kreiman/jerry/data/h5/art';
    dir_res = '/n/scratch2/jw324/data/results/coh_w10';
    dir_cor = '/n/scratch2/jw324/data/coreg';
    dir_h5 = '/n/groups/kreiman/jerry/data/h5';
end

% Patients
Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48',...
    'mSu'};

% Exclude monkey
Subjects = Subjects(1:(end-1));

% Constants
metrics = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};
metric = metrics{1};
atl = 2; % atlas index: 1 - Destrieux, 2 - Desikan-Killiany, 6 - HCP-MMP1
n_perm = 10000;
dist_thresh = 20; % mm

% Thresholds
coh_thresh_o = [0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3];
ct_thresh_o = [0.05 0.1 0.2 0.3];
cp_thresh_o = [0.05 0.1 0.2 0.3];
n_params = length(coh_thresh_o)*length(ct_thresh_o)*length(cp_thresh_o);
all_o = zeros(n_params,3);
all_i = 1;
for ii = 1:length(coh_thresh_o)
    for iii = 1:length(ct_thresh_o)
        for iv = 1:length(cp_thresh_o)
            all_o(all_i,1) = coh_thresh_o(ii);
            all_o(all_i,2) = ct_thresh_o(iii);
            all_o(all_i,3) = cp_thresh_o(iv);
            all_i = all_i + 1;
        end
    end
end

% Cache files
cfiles = cell(1,length(Subjects));
for i = 1:length(cfiles)
    cfiles{i} = sprintf('cache/%s_R_null.mat',Subjects{i});
end

% parfor
trig_cache = true;
curr_d = pwd;
% delete(gcp('nocreate'));
% poolobj = parpool;
% if (trig_cache)
%     addAttachedFiles(poolobj,cfiles);
% end
for v = 1:n_params
    tic;

    % coh_thresh = 0.3;
    % ct_thresh = 0.1;
    % cp_thresh = 0.5;
%     coh_thresh = coh_thresh_o(ii);
%     ct_thresh = ct_thresh_o(iii);
%     cp_thresh = cp_thresh_o(iv);
    coh_thresh = all_o(v,1);
    ct_thresh = all_o(v,2);
    cp_thresh = all_o(v,3);

    % Load atlas
    C1 = load(sprintf('%s/%s/label/all_parcellation.mat',dir_cor,Subjects{1}));
    atl_name = C1.AtlNames{atl};
    n_rois = C1.AtlROIs{atl}.LH.numEntries;
    rois = C1.AtlROIs{atl}.LH.struct_names;
    n_comb_atl = nchoosek(n_rois,2);
    AdjAtl = cell(n_rois,n_rois);
    AdjAtlN = cell(n_rois,n_rois);

    % Print header
    %fprintf('pt\tmetric\tcoh_thresh\tct_thresh\tisx\tisxN\n')

    % Main loop
    for i = 1:length(Subjects)
        %tic;
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
        Dmat = ecog.getDist();
        Dmats = zeros(ecog.n_comb,1);
        count = 1;
        for i1 = 1:(ecog.n_bchan-1)
            for i2 = (i1+1):ecog.n_bchan
                Dmats(count) = Dmat(i1,i2);
                count = count + 1;
            end
        end

        % --- Load graph and artifacts ---
        C = load(fn_coreg);
        atl_labels = C.AtlLabels{atl};
        art_idx = h5read(fn_art,'/art_idx');
        art_idxB = (art_idx == 1);
        frac_art = h5readatt(fn_art,'/art_idx','frac_art');
        w_art = h5readatt(fn_art,'/art_idx','w');
        % trim last time sample of graph
        R = h5read(fn_graph,'/R',[1 1],size(art_idx));
        [~,n_graph] = size(R);
        chan1 = h5read(fn_graph,'/chan1');
        chan2 = h5read(fn_graph,'/chan2');
        t_load = toc;
        %fprintf('[%s] Loading took: %.3f seconds\n',sid,t_load)

        % --- Prepare null data ---
        if (~ trig_cache)
            R_null = zeros(size(R));
            D = load(fn_dist);
            for j = 1:ecog.n_comb
                R_null(j,:) = random(D.d{j},[1 n_graph]);
            end
        else
            
            fn_tmp = sprintf('%s/cache/%s_R_null.mat',curr_d,sid);
            if (~exist(fn_tmp,'file'))
                R_null = zeros(size(R));
                D = load(fn_dist);
                for j = 1:ecog.n_comb
                    R_null(j,:) = random(D.d{j},[1 n_graph]);
                end
                save_tmp(sid,single(R_null));
            else
                R_null = load_tmp(fn_tmp);
            end
        
        end
        
        %return
        
        %t_rand = toc;
        %fprintf('[%s] Random sampling took: %.3f seconds\n',sid,t_rand-t_load)

        % --- Build CT ---
        d_i = Dmats > dist_thresh;

        n_yes = sum((R >= coh_thresh) & (~art_idxB),2);
        n_no = sum((R < coh_thresh) & (~art_idxB),2);
        ct = n_yes./(n_yes + n_no);
        is_x = (ct >= ct_thresh);
        frac_x = sum(is_x(d_i))/numel(is_x(d_i));
        %R4m = double(R);
        R4m = R;
        R4m(R < coh_thresh) = NaN; % --- Comment this line to include all points in average magnitude ---
        R4m(art_idxB) = NaN;
        mag = double(nanmean(R4m,2));

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
        fprintf('%i\t%s\t%s\t%.8f\t%.8f\t%.8f\t%.8f\n',v,sid,metric,coh_thresh,ct_thresh,frac_x,frac_x_null);


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
                    for k1 = 1:length(i_b1)
                        for k2 = 1:length(i_b2)
                            ki1 = i_b1(k1);
                            ki2 = i_b2(k2);
                            AdjAtl{ki1,ki2} = [AdjAtl{ki1,ki2}; el];
                            AdjAtl{ki2,ki1} = [AdjAtl{ki2,ki1}; el];

                            AdjAtlN{ki1,ki2} = [AdjAtlN{ki1,ki2}; el_null];
                            AdjAtlN{ki2,ki1} = [AdjAtlN{ki2,ki1}; el_null];
                        end
                    end

                catch
                    %fprintf('ct-atl map skip.\n');
                end
            end
        end

        %return

        % Predict ETA
        %t_sing = toc;
        %fprintf('%s\tETA: %.2f hrs\n',sid,t_sing*(length(Subjects) - i)/3600)

        %return

    end

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
    frac_y = sum(cp > cp_thresh) / sum(~isnan(cp));
    frac_y_null = sum(cp_null > cp_thresh) / sum(~isnan(cp_null));
    fprintf('%i\t%s\t%s\t%.8f\t%.8f\t%.8f\n',v,'xsub',metric,cp_thresh,frac_y,frac_y_null);
    %return;
    
    % Clear loop indices
    %clear i;
    %clear j;
    
    t_sing = toc;
    fprintf('Est. total time: %.2f hrs\n',t_sing * n_params * (1/3600));
end

% Print finish message
fprintf('Done.\n')

function save_tmp(sid, R_null)
    save(sprintf('cache/%s_R_null',sid),'R_null','-v6');
end

function R_null = load_tmp(R_null_fn)
    R_null = load(R_null_fn);
end
