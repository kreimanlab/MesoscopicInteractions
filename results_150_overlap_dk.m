close all;
clear;

% load surf
dir_cacheLp = './cache';
dir_corLp = '/media/klab/internal/data/coreg';
metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
metrics_suffix = {'0.5-125 Hz','3-8 Hz','8-12 Hz','12-30 Hz','30-100 Hz'};
SURFACE_TYPE = 'pial';
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_corLp,'fsaverage_sym','r',SURFACE_TYPE));
iM = 1;

% original subject
sid = 'sub3';
sid_i = str2double(sid(2:end));
b1 = 35;
b2 = 84;
%CaR = load(sprintf('%s/fig_cluster2_reduce_%i_new',dir_cache,iM));
CaA = load(sprintf('%s/xsub_out_all_%i',dir_cacheLp,iM));

% convert to integer subject number
sid_int = find(strcmp(CaA.Subjects,sid));
Ca = load(sprintf('%s/xsub_out_%s_%i',dir_cacheLp,sid,iM));

% Find ROIs bipolar electrodes map onto
rois = Ca.C.AtlLabels{2};
b1c1 = Ca.ecog.bip(b1,:);
b2c1 = Ca.ecog.bip(b2,:);
b1_roi = rois{b1c1};
b2_roi = rois{b2c1};
%--------------------------------------------------------------------------

% Loop through all rois
dk_rois = Ca.C.AtlROIs{2}.LH.struct_names;

for idk = 1:length(dk_rois)
    bx_roi = dk_rois{idk};


    % Begin copy pasta from figure_t18.m
    
    % Load adjacency matrix
    CaT14 = load(sprintf('%s/figure_t14_%i_150.mat',dir_cacheLp,iM));

    % region highlight by number

    % Color DK region
    const_col_surface = 0.9 * ones(1,3);
    [v_a, l, ct] = read_annotation(sprintf('%s/%s/label/%sh.aparc.annot',dir_corLp,'fsaverage_sym','r'));
    [n_vert,~] = size(s_vert);
    const_col_surface2 = repmat(const_col_surface,[n_vert 1]);
    vtx_r1 = false(length(l),1);
    vtx_r2 = false(length(l),1);
    vtx_rx = false(length(l),1);
    for ihi = 1:ct.numEntries
        roi_id = ct.table(ihi,end);
        cond_rx = strcmp(ct.struct_names{ihi},bx_roi);
        cond_r1 = strcmp(ct.struct_names{ihi},b1_roi);
        cond_r2 = strcmp(ct.struct_names{ihi},b2_roi);
        if (cond_r1)
            vtx_r1(l==roi_id) = true;
        end
        if (cond_r2)
            vtx_r2(l==roi_id) = true;
        end
        if (cond_rx)
            vtx_rx(l==roi_id) = true;
        end
        cond_to_color = cond_r1 | cond_r2;
        if (cond_to_color)
            roi_col = ct.table(ihi,1:3) / 255;
            n_vtx = sum(l==roi_id);
            const_col_surface2(l==roi_id,:) = repmat(roi_col,[n_vtx 1]);
        end
    end
    col_r1 = ct.table(strcmp(ct.struct_names,b1_roi),1:3) / 255;
    col_r2 = ct.table(strcmp(ct.struct_names,b2_roi),1:3) / 255;
    col_rx = ct.table(strcmp(ct.struct_names,bx_roi),1:3) / 255;
    vcolor = const_col_surface2;
    %vcolor = 0.8*ones(length(s_vert),3);
    %roi_cl = Res{k,1};

    % deprecated
    cluster_i = load('cache/fig_cluster3_cluster_i.mat');
    cluster_i = cluster_i.cluster_i;
    % end deprecated
    %     CaT14 = load(sprintf('./cache/figure_t14_%i_150',iM));
    %     cluster_i2 = zeros(size(cluster_i));
    %     for i = 1:length(cluster_i)
    %         idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(i))));
    %         cluster_i2(i) = idx;
    %     end
    %     cluster_i = cluster_i2; 


    CaC2 = load('cache/fig_cluster2_reduce_annot');

    % 150 region - superior temporal overlap
    tint1_factor = 0.1;
    tint1_col = [0 0.4 0];
    tint1_border_col = [0 0.2 0];
    % 150 region - pars op overlap
    tint2_factor = 0.1;
    tint2_col = [0.4 0 0];
    tint2_border_col = [0.2 0 0];

    % First loop to find overlap areas
    %   label each vertex in s_vert with one of 150 areas
    vtx_region = nan(length(s_vert),1);
    for j = 1:length(s_vert)
        r = CaC2.region(j);
        idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(r))));
        vtx_region(j) = idx;
    end
    %fprintf('[*] First loop complete.\n')

    r1_overlap_150 = zeros(length(CaT14.rois_plt),1);
    for r1 = 1:length(r1_overlap_150)
        %f_overlap = mean(vtx_r1(vtx_r1) & (vtx_region(vtx_r1) == r1));
        f_overlap = mean(vtx_r1(vtx_region == r1));
        r1_overlap_150(r1) = f_overlap;
    end
    %plot(r1_overlap_150,'black.')

    r2_overlap_150 = zeros(length(CaT14.rois_plt),1);
    for r2 = 1:length(r2_overlap_150)
        f_overlap = mean(vtx_r2(vtx_region == r2));
        r2_overlap_150(r2) = f_overlap;
    end

    rx_overlap_150 = zeros(length(CaT14.rois_plt),1);
    for rx = 1:length(rx_overlap_150)
        f_overlap = mean(vtx_rx(vtx_region == rx));
        rx_overlap_150(rx) = f_overlap;
    end

    thresh_overlap = 0.5; % figure_t18.m: 0.3
    region_const1 = find(r1_overlap_150 > thresh_overlap);
    region_const2 = find(r2_overlap_150 > thresh_overlap);
    region_constx = find(rx_overlap_150 > thresh_overlap);
    
    outstr = '';
    for ix = 1:length(region_constx)
        if (ix == 1)
            outstr = sprintf('%s%i',outstr,region_constx(ix));
        else
            outstr = sprintf('%s,%i',outstr,region_constx(ix));
        end
    end
    fprintf('[%s] 150 areas: %s\n',bx_roi,outstr);
end

% April 8, 2020

% overlap threshold: 0.5

% [unknown] 150 areas: 25,30,83,86,88,93,144
% [bankssts] 150 areas: 85
% [caudalanteriorcingulate] 150 areas: 
% [caudalmiddlefrontal] 150 areas: 57,142,150
% [corpuscallosum] 150 areas: 
% [cuneus] 150 areas: 55,84
% [entorhinal] 150 areas: 10,106
% [fusiform] 150 areas: 4,5,6,75,105,113,123
% [inferiorparietal] 150 areas: 2,15,21,59,73,87,89,95,104,132,141
% [inferiortemporal] 150 areas: 13,16,17,23,68,98,100,101,102,107,108,109
% [isthmuscingulate] 150 areas: 69
% [lateraloccipital] 150 areas: 22,42,49,79,80,112
% [lateralorbitofrontal] 150 areas: 36,114,120
% [lingual] 150 areas: 44,47,71
% [medialorbitofrontal] 150 areas: 34,50
% [middletemporal] 150 areas: 7,14,18,19,20,78,128,129
% [parahippocampal] 150 areas: 
% [paracentral] 150 areas: 24
% [parsopercularis] 150 areas: 31,125
% [parsorbitalis] 150 areas: 
% [parstriangularis] 150 areas: 32,39,130
% [pericalcarine] 150 areas: 82
% [postcentral] 150 areas: 29,48,67,131,136,139
% [posteriorcingulate] 150 areas: 27
% [precentral] 150 areas: 119,121,122,126,127,133,135,137,143
% [precuneus] 150 areas: 41,43,66
% [rostralanteriorcingulate] 150 areas: 
% [rostralmiddlefrontal] 150 areas: 28,35,51,53,54,70,90,91,103,116,118,140
% [superiorfrontal] 150 areas: 58,60,62,63,72,74,76,94,115
% [superiorparietal] 150 areas: 3,45,64,92,117
% [superiortemporal] 150 areas: 9,46,77,97,110,111,145
% [supramarginal] 150 areas: 8,33,40,65,96,138,149
% [frontalpole] 150 areas: 
% [temporalpole] 150 areas: 147
% [transversetemporal] 150 areas: 
% [insula] 150 areas: 