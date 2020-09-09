close all;
clear;


% Atlas atl1-20 file definition
D = load('./cache/xsub_out_all_1_atl1');
atl_names = D.C.AtlNames;
atl_fn = { ...
    'aparc.a2009s.annot'; ... %Destrieux-2009
    'aparc.annot'; ... %Desikan-Killiany
    'PALS_B12_Brodmann.annot'; ... % Brodmann
    'Yeo2011_7Networks_N1000.annot'; ... % Yeo-2011-7
    'Yeo2011_17Networks_N1000.annot'; ... % Yeo-2011-17
    'HCP-MMP1.annot'; ... % HCP-MMP1
    'MACAQUE_M132.annot'; ... % M132
    'MACAQUE_FVE91.annot'; ... % FVE91
    'MACAQUE_FVEall.annot'; ... % FVEall
    'MACAQUE_LVE00.annot'; ... % LVE00
    'MACAQUE_PHT00.annot'; ... % PHT00
    'MACAQUE_Brodmann.annot'; ... % Brodmann
    'MACAQUE_BoninBailey.annot'; ... % BoninBailey
    'MACAQUE_FerryEtAl.annot'; ... % FerryEtAl
    'MACAQUE_UD86.annot'; ... % UD86
    'MACAQUE_PGR91.annot'; ... % PGR91
    'MACAQUE_LANDMARK.annot'; ... % LANDMARK
    'MACAQUE_LyonKaas.annot'; ... % LyonKaas
    'MACAQUE_BRL87.annot'; ... % BRL87
    'MACAQUE_SP78.annot' % SP78
};

metrici = 1;
% Load 150 locations
D = load(sprintf('brainexport/red_6all_fsaverage_%i.mat',metrici));

% Colors
BRAIN_MESH_ALPHA = 1;
COLOR_MEDIAL_WALL = 0.1*[1 1 1];
COLOR_NOSIG = 0.5*[1 1 1];
COLOR_BLACK = 0.3*[1 1 1];
SURFACE_TYPE = 'pial';
paper_width = 8.5;
paper_height = 6.5;


% load surf
dir_cacheLp = './cache';
dir_corLp = '/media/klab/internal/data/coreg';
metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
metrics_suffix = {'0.5-125 Hz','3-8 Hz','8-12 Hz','12-30 Hz','30-100 Hz'};
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_corLp,'fsaverage_sym','r',SURFACE_TYPE));



%i = 1;
%iM = i;
iM = 1;

metric = metricsp{iM};


% plot 150 boundaries on DK areas

%h = figure('visible','off');
%set(h,'PaperUnits','inches');
%set(h,'PaperPosition',[0 0 paper_width paper_height]);
%axis off;
%hold all;

%subplot(2,1,1);
%[ha, pos] = tight_subplot(2,2,[0.01 0.01],8*[0.01 0.01],4*[0.01 0.01]);
%axes(ha(1));
%p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
%    'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
%hold all;



%--------------------------------------------------------------------------
%
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



% Load
CaT14 = load(sprintf('%s/figure_t14_%i_150.mat',dir_cacheLp,iM));
cluster_i = load('cache/fig_cluster3_cluster_i.mat');
cluster_i = cluster_i.cluster_i;
CaC2 = load('cache/fig_cluster2_reduce_annot');

% region highlight by number
AChans = cell(1,length(atl_names));
AChansR = cell(1,length(atl_names));
for atl = 1:length(atl_names)
    atl_ifn = sprintf('%s/%s/label/%sh.%s',dir_corLp,'fsaverage_sym','r',atl_fn{atl});
    % Read atlas
    [v_a, l, ct] = read_annotation(atl_ifn);
    idx_clean = ~(all(ct.table == 0,2));
    ct.table = ct.table(idx_clean,:);
    ct.struct_names = ct.struct_names(idx_clean,:);
    ct.numEntries = sum(idx_clean);
    [n_reg,~] = size(ct.table);
    fprintf('[*] %s: %s\n',atl_names{atl},atl_ifn);
    fprintf('\t# regions: %i\n',n_reg)

    % Match atlas region to pial vertex
    vtx_reg = false(length(l),n_reg);
    for ihi = 1:ct.numEntries
        roi_id = ct.table(ihi,end);
        %cond_r1 = strcmp(ct.struct_names{ihi},b1_roi);
        vtx_reg(l==roi_id,ihi) = true;
    end

    % Check each vertex is mapped
    if (~all(sum(vtx_reg,2) == 1))
        fprintf(2,'[!] Warn: not all vertices labeled by atl.\n');
    end
    
    %   label each vertex in s_vert with one of 150 areas
    fprintf('\t(*) 150-vertex suchen..\n')
    vtx_region = nan(length(s_vert),1);
    for j = 1:length(s_vert)
        r = CaC2.region(j);
        idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(r))));
        vtx_region(j) = idx;
    end
    fprintf('\t(*) fertig.\n')
    
    % Assign each one of 150 areas to region in atl
    atl_chans = 1:length(cluster_i);
    overlap_thresh = 0.1; % Default 0.1
    for i = 1:length(atl_chans)
        idx_vr = (vtx_region == i);
        % Score overlap with regions, normalize
        vtx_reg_r = vtx_reg(idx_vr,:);
        reg_score = sum(vtx_reg_r)./sum(sum(vtx_reg_r));
        % reg_score(reg_score < overlap_thresh) = 0;
        [~,reg_max] = max(reg_score);
        atl_chans(i) = reg_max;
    end
    
    % Assign each one of the N areas in atl to region in 150-parcellation
    atl_chansR = cell(1,n_reg);
    for i = 1:length(atl_chansR)
        idx_vr = vtx_reg(:,i);
        vtx_reg_r = vtx_region(idx_vr,:);
        reg_scoreR = nan(1,length(cluster_i));
        for j = 1:length(reg_scoreR)
            reg_scoreR(j) = sum(vtx_reg_r == j);
        end
        reg_scoreR = reg_scoreR / sum(reg_scoreR);
        
        % Only map 1 of 150 regions to atlas region if overlap threshold
        atl_chansR{i} = find(reg_scoreR > overlap_thresh);
    end
    
    AChans{atl} = atl_chans;
    AChansR{atl} = atl_chansR;
    
end

save('./cache/custom_parcellation_atlas_overlap.mat','AChans','AChansR');
fprintf('[!] Done.\n');
return





% Color DK region
const_col_surface = 0.9 * ones(1,3);
[v_a, l, ct] = read_annotation(sprintf('%s/%s/label/%sh.aparc.annot',dir_corLp,'fsaverage_sym','r'));
[n_vert,~] = size(s_vert);
const_col_surface2 = repmat(const_col_surface,[n_vert 1]);
vtx_r1 = false(length(l),1);
vtx_r2 = false(length(l),1);
for ihi = 1:ct.numEntries
    roi_id = ct.table(ihi,end);
    cond_r1 = strcmp(ct.struct_names{ihi},b1_roi);
    cond_r2 = strcmp(ct.struct_names{ihi},b2_roi);
    if (cond_r1)
        vtx_r1(l==roi_id) = true;
    end
    if (cond_r2)
        vtx_r2(l==roi_id) = true;
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
vcolor = const_col_surface2;
%vcolor = 0.8*ones(length(s_vert),3);
%roi_cl = Res{k,1};

cluster_i = load('cache/fig_cluster3_cluster_i.mat');
cluster_i = cluster_i.cluster_i;
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
fprintf('[*] 150-vertex suchen..\n')
vtx_region = nan(length(s_vert),1);
for j = 1:length(s_vert)
    r = CaC2.region(j);
    idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(r))));
    vtx_region(j) = idx;
end
fprintf('[*] fertig.\n')

% Compute overlap of each 150 areas with atlas of interest

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



thresh_overlap = 0.3;
region_const1 = find(r1_overlap_150 > thresh_overlap);
region_const2 = find(r2_overlap_150 > thresh_overlap);
coord1 = cell(length(region_const1),1);
coord2 = cell(length(region_const2),1);

