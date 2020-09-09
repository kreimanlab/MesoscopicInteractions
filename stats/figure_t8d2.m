close all
clear
rng shuffle;


%metric = 'pcBroadband';
n_perm = 10000;
perm_alpha_sec = 12;
cp_thresh_override = 0.05; % 0.01, new 0.05
trig_zeros_to_nan = false;
trig_remove_diagonal = true;
n_pairs_thresh = 1; %round(20 * 1 * (35/90)); % at least this many electrode pairs to be considered 12 * 4 * (35/90)
n_subs_thresh = 1; % at least this many subjects to be considered
n_subs_ct_thresh = 1; % significant CTs in region pair must be from at least this many subjects
fig_fmt = '-depsc';
trig_eps = false;
trig_mag_no_cp = true;
n_colorbar_ticks = 5;
system('mkdir figures');
system('mkdir figures/T8d2');

% Fast i/o definitions
dir_artL = '/media/klab/internal/data/h5_notch20/art';
dir_resL = '/media/klab/internal/data/results/coh_w10';
dir_corL = '/media/klab/internal/data/coreg';
dir_cacheL = './cache'; % /old_coh-stats_ct-100_cp-100
subjects_dirL = '/mnt/cuenap_ssd/coregistration';

% Slow i/o definitions
dir_h5L = '/media/klab/KLAB101/h5_notch20';

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'}; % 'pcDelta'

banned_rois = {};%{'24a_M132'};

% Patients
%Subjects = {'mSu'};
%sid = Subjects{1};
%iM = 1;

for iM = 1:length(metrics)
    
metric = metrics{iM};



% fn_artL = sprintf('%s/%s_art.h5',dir_artL,sid);
% fn_distL = sprintf('%s/%s_dists-%s-%i.mat',dir_resL,sid,metric,n_perm);
% fn_graphL = sprintf('%s/%s_graph-%s.h5',dir_resL,sid,metric);
% fn_permL = sprintf('/mnt/cuenap2/data/results/coh_w10/%s_perm-%s-%i.h5',sid,metric,n_perm);
% fn_h5L = sprintf('%s/%s.h5',dir_h5L,sid);
% fn_coregL = sprintf('%s/%s/label/all_parcellation.mat',dir_corL,sid);
% fn_cacheL = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheL,sid,iM);
% 
% % Load mSu cache
% load(fn_cacheL);

% Load human cache
load(sprintf('%s/xsub_out_all_%i.mat',dir_cacheL,iM));

% Register to macaque labels
D_anat = load('../coreg/AdjMKFV.mat');
D = load('../coreg/all_surf_ielvis_macaque.mat');
load('SuMAP.mat');
crop_pix = 50; % vertical crop
crop_pix_h = 20; % horizontal crop
Su.I = Su.I(crop_pix:(end-crop_pix), crop_pix_h:(end-crop_pix_h),:);
%dist_thresh = 20;
atl = 1;
atl_name = D.AtlasName{atl};
atl_labels = D.AtlasLabels{atl};


% init
Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

rois = C1.AtlROIs{7}.LH.struct_names;
n_rois = length(rois);
AdjAtl_loc = cell(n_rois,n_rois);
AdjAtlN_loc = cell(n_rois,n_rois);
adjct_dist_loc = cell(n_rois,n_rois);
AdjAtl_sid_loc = cell(n_rois,n_rois);

for iS = 1:length(Subjects)
    
    sid = Subjects{iS};
    fn_cacheL = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheL,sid,iM);
    % 
    % % Load subject cache
    load(fn_cacheL);
    
    % Overwrite with macaque ROIs
    rois = C1.AtlROIs{7}.LH.struct_names;
    atl_labels = C.AtlLabels{7};
    
    % ban labels
    for ii = 1:length(atl_labels)
        if (sum(strcmp(banned_rois,atl_labels{ii})) ~= 0)
            atl_labels{ii} = 'BANNED';
        end
    end

    % Map onto macaque atlas
    %dist_thresh_mac = dist_thresh / 2; % macaque grids have 5mm spacing instead of 10mm
    d_i = Dmats > dist_thresh;
    %fprintf('Distance threshold fraction removed: %.4f\n',1-sum(d_i)/length(d_i))
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
                        AdjAtl_loc{ki1,ki2} = [AdjAtl_loc{ki1,ki2}; el];
                        AdjAtl_loc{ki2,ki1} = [AdjAtl_loc{ki2,ki1}; el];

                        AdjAtlN_loc{ki1,ki2} = [AdjAtlN_loc{ki1,ki2}; el_null];
                        AdjAtlN_loc{ki2,ki1} = [AdjAtlN_loc{ki2,ki1}; el_null];

                        adjct_dist_loc{ki1,ki2} = [adjct_dist_loc{ki1,ki2}; Dmats(k)];
                        adjct_dist_loc{ki2,ki1} = [adjct_dist_loc{ki2,ki1}; Dmats(k)];
    % 
                        AdjAtl_sid_loc{ki1,ki2} = [AdjAtl_sid_loc{ki1,ki2},subject_i];
                        AdjAtl_sid_loc{ki2,ki1} = [AdjAtl_sid_loc{ki2,ki1},subject_i];
                    end
                end

            catch e
                %fprintf('[W] ct-atl map skip.\n');
            end
        end
    end
    fprintf('%s\tfinished mapping %s to M132\n',metric,sid);
end
%adjct_dist_sav = adjct_dist_loc;

% M132 rois
rois = C1.AtlROIs{7}.LH.struct_names;
n_rois = length(rois);

% Calculate final functional interaction matrix
Adj = nan(n_rois,n_rois);
AdjMag = nan(n_rois,n_rois);
AdjCP = nan(n_rois,n_rois);
N_bchan = nan(n_rois,n_rois);
cp_thresh = cp_thresh_override;
Dmat = Inf(n_rois,n_rois); % set to inf to avoid removing from every instance

DistsAtl = [];
for i1 = 1:n_rois
    for i2 = 1:n_rois
        AA = AdjAtl_loc{i1,i2};
        AA_dist = adjct_dist_loc{i1,i2};
        
        AA_sub = AdjAtl_sid_loc{i1,i2};
        n_pairs = length(AA_sub);
        n_subs = length(unique(AA_sub));
        n_subs_ct = length(unique(AA_sub(AA ~= 0)));
        % ROI pair coverage condition (grey)
        if ( (~ isempty(AA))  && (n_pairs >= n_pairs_thresh) && (n_subs >= n_subs_thresh))
            N_bchan(i1,i2) = length(AA);
            frac_cp = sum(AA ~= 0)/length(AA);
            AdjCP(i1,i2) = frac_cp;
            % ROI pair significance condition (white)
            %if (frac_cp >= cp_thresh)
            if ( (frac_cp >= cp_thresh) && (n_subs_ct >= n_subs_ct_thresh) )
                Adj(i1,i2) = 1;
                %AdjMag(i1,i2) = mean(AA(frac_cp > cp_thresh));
                AdjMag(i1,i2) = mean(AA(AA ~= 0));
                DistsAtl = [DistsAtl; [mean(AA_dist(AA ~= 0)), mean(AA(AA ~= 0))]];
            else
                Adj(i1,i2) = 0;
                AdjMag(i1,i2) = 0;
            end
        end
    end
end
save('../coreg/AdjSu','Adj','AdjMag','AdjCP','Dmat');
%return

if (trig_mag_no_cp)
    AdjCP = AdjMag;
end

if (trig_zeros_to_nan)
    AdjCP(AdjCP == 0) = NaN;
end


frac_xsub_msu = sum(Adj(~isnan(Adj))) /numel((~isnan(Adj)));
fprintf('%s - human fraction of ROIs significant: %.4f\n',metric,frac_xsub_msu)


% compare to macaque
load('../coreg/AdjMKFV.mat')

%AdjMKneurons = AdjMKneurons + AdjMKneurons';
% make symmetric
a = AdjMKneurons;
b = zeros(size(a));
for i = 1:length(a)
    for j = i:length(a)
        val = nansum([a(i,j), a(j,i)]);
        b(i,j) = val;
        b(j,i) = val;
    end
end
AdjMKneurons = b;

%AdjMKflne = (AdjMKflne + AdjMKflne');
a = AdjMKflne;
b = zeros(size(a));
for i = 1:length(a)
    for j = i:length(a)
        val = nanmean([a(i,j), a(j,i)]);
        b(i,j) = val;
        b(j,i) = val;
    end
end
AdjMKflne = b;

% 
% if (trig_remove_diagonal)
%     for izz = 1:length(AdjCP)
%         AdjCP(izz,izz) = NaN;
%     end
%     for izz = 1:length(AdjMKflne)
%         AdjMKflne(izz,izz) = NaN;
%     end
% end
% 
% AdjMKneurons = AdjMKneurons + AdjMKneurons';
% 
% AdjMKneuronst = AdjMKneurons;
% AdjMKneuronst(AdjMKneuronst == 0) = NaN;
% AdjMK = log10(AdjMKneuronst);
% AdjMK_bin = ~isnan(AdjMK);
% AdjMK_bin = double(AdjMK_bin);


AdjMKflnet = AdjMKflne;
sav = AdjMKflnet;
AdjMKflnet(sav == 0) = NaN;
AdjMKflne = log10(AdjMKflnet);
AdjMKflneRaw = sav;
%AdjMKflne(isnan(AdjMKflne)) = 0;

AdjMKneuronst = AdjMKneurons;
AdjMKneuronst(AdjMKneuronst == 0) = NaN;
AdjMK = log10(AdjMKneuronst);
AdjMKRaw = AdjMKneuronst;

AdjMK_bin = ~isnan(AdjMK);
AdjMK_bin = double(AdjMK_bin);


% Align macaque map to human map
%MK_human_labels = AT.P.AtlROIs{atl_i}.LH.struct_names;
MK_human_labels = rois;
for i = 1:length(MK_human_labels)
    st = strsplit(MK_human_labels{i},'_M132');
    st = st{1};
    MK_human_labels{i} = st;
end
n_found = 0;
ind_hum2mac = zeros(1,length(AdjMK));
for i = 1:length(MK_labels)
    % Edits
    MK_labels{i} = replace(MK_labels{i},'/','_');
    if strcmp(MK_labels{i},'Core')
        MK_labels{i} = 'Aud-core';
    end
    if strcmp(MK_labels{i},'ENTO')
        MK_labels{i} = 'Entorhinal';
    end
    if strcmp(MK_labels{i},'INS')
        MK_labels{i} = 'Insula';
    end
    if strcmp(MK_labels{i},'PERI')
        MK_labels{i} = 'Perirhinal';
    end
    if strcmp(MK_labels{i},'PIR')
        MK_labels{i} = 'Piriform';
    end
    if strcmp(MK_labels{i},'POLE')
        MK_labels{i} = 'Temporal-pole';
    end
    if strcmp(MK_labels{i},'Pro.St.')
        MK_labels{i} = 'Prostriate';
    end
    if strcmp(MK_labels{i},'SII')
        MK_labels{i} = 'S2';
    end
    if strcmp(MK_labels{i},'SUB')
        MK_labels{i} = 'Subiculum';
    end
    if strcmp(MK_labels{i},'TEa_ma')
        MK_labels{i} = 'TEam-a';
    end
    if strcmp(MK_labels{i},'TEa_mp')
        MK_labels{i} = 'TEam-p';
    end
    
    la = find(strcmp(MK_human_labels,MK_labels{i}),1);
    if (~isempty(la))
        %AdjCS_mac(i,:) = AdjCS(la,:);
        ind_hum2mac(i) = la;
        n_found = n_found + 1;
    end
end

% Calculate correlations mSu ECoG with MK tractography
% Binary comparisons:
%   Adj(ind_hum2mac,ind_hum2mac), AdjMK_bin
AdjMK_bin_f = Adj(ind_hum2mac,ind_hum2mac);
%v1 = AdjMK_bin_f(~isnan(AdjMK_bin_f));
%v2 = AdjMK_bin(~isnan(AdjMK_bin_f));

% triangle conversion
v1t = [];
v2t = [];
for iz = 1:(length(AdjMK_bin_f)-1)
    for iiz = (iz + 1):length(AdjMK_bin_f)
        if ((~isnan(AdjMK_bin_f(iz,iiz))) && (~isnan(AdjMK_bin(iz,iiz))))
            v1t = [v1t; AdjMK_bin_f(iz,iiz)];
            v2t = [v2t; AdjMK_bin(iz,iiz)];
        end
    end
end

trig_plot_flne = true;
if (trig_plot_flne)
    AdjMKneuronsP = AdjMKflneRaw;
    AdjMKneuronsP((AdjMKflneRaw == 0)) = NaN;
    AdjMKneuronsP = log10(AdjMKneuronsP);
    AdjMKneuronsP((AdjMKflneRaw == 0)) = 0;
else
    AdjMKneuronsP = AdjMKneuronsRaw;
    AdjMKneuronsP((AdjMKneuronsRaw == 0)) = NaN;
    AdjMKneuronsP = log10(AdjMKneuronsP);
    AdjMKneuronsP((AdjMKneuronsRaw == 0)) = 0;
end


[r_bin,p_val_bin] = corr(v1t,v2t,'Type','Spearman');
fprintf('Binary r: %.4f, p = %d\n',r_bin,p_val_bin);
% Weight comparions:
%   AdjCP(ind_hum2mac,ind_hum2mac), AdjMKflne
%AdjMK_f = AdjCP(ind_hum2mac,ind_hum2mac);
AdjMK_f = AdjMag(ind_hum2mac,ind_hum2mac);
%v1 = AdjMK_f(~isnan(AdjMK_f));
%v2 = AdjMKneurons(~isnan(AdjMK_f));
%v2 = AdjMKflne(~isnan(AdjMK_f));
%v2(isnan(v2)) = 0;

% triangle conversion
v1t = [];
v2t = [];
for iz = 1:(length(AdjMK_f)-1)
    for iiz = (iz + 1):length(AdjMK_f)
        if ((~isnan(AdjMK_f(iz,iiz))) && (~isnan(AdjMKflne(iz,iiz))))
            v1t = [v1t; AdjMK_f(iz,iiz)];
            v2t = [v2t; AdjMKflne(iz,iiz)];
        end
    end
end

%[r,p_val] = corr(v1,v2);
[rS,p_valS] = corr(v1t,v2t,'type','Spearman');
[rS2,p_valS2] = corr(v1t,v2t,'type','Pearson');
fprintf('Mag vs flne Spearman: r = %.4f, p = %d\n',rS,p_valS);
fprintf('Mag vs flne Pearson: r = %.4f, p = %d\n',rS2,p_valS2);
%fprintf('CP vs Neurons Pearson r: %.4f, p = %d\n',rS,p_valS);




% Remove zeros
v1_no0 = AdjMK_f(~isnan(AdjMK_f));
%v2_no0 = AdjMKneurons(~isnan(AdjMK_f));
v2_no0 = AdjMKflne(~isnan(AdjMK_f));
v1_no0s = v1_no0;
v1_no0 = v1_no0(v1_no0s ~= 0);
v2_no0 = v2_no0(v1_no0s ~= 0);
%v2(isnan(v2)) = 0;
[rS_no0,p_valS_no0] = corr(v1_no0,v2_no0,'type','Spearman');
%[rS_no0,p_valS_no0] = corr(v1_no0,v2_no0,'type','Pearson');
fprintf('Mag vs flne, zeros removed: r = %.4f, p = %d\n',rS_no0,p_valS_no0);
%return


fprintf('[!] New correlations\n')

% triangle conversion
v1t = [];
v2t = [];
v3t = [];
v4t = [];
v5t = [];
for iz = 1:(length(AdjMK_f)-1)
    for iiz = (iz + 1):length(AdjMK_f)
        if ((~isnan(AdjMK_f(iz,iiz))) && (~isnan(AdjMKflne(iz,iiz))))
%             v1t = [v1t; AdjMK_f(iz,iiz)];
%             %v2t = [v2t; AdjMKflne(iz,iiz)];
%             v2t = [v2t; AdjMKneuronsP(iz,iiz)];
%             v3t = [v3t; AdjMKneurons(iz,iiz)];
%             v4t = [v4t; AdjMKflne(iz,iiz)];
            v1t = [v1t; AdjMK_f(iz,iiz)];
            v2t = [v2t; AdjMK(iz,iiz)];
            v3t = [v3t; AdjMKRaw(iz,iiz)];
            v4t = [v4t; AdjMKflne(iz,iiz)];
            v5t = [v5t; AdjMKflneRaw(iz,iiz)];
        end
    end
end

%[r,p_val] = corr(v1,v2);
[rS,p_valS] = corr(v1t,v2t,'type','Spearman');
[rS2,p_valS2] = corr(v1t,v2t,'type','Pearson');
[rSr,p_valSr] = corr(v1t,v3t,'type','Spearman');
[rS2r,p_valS2r] = corr(v1t,v3t,'type','Pearson');
fprintf('Mag vs neurons Spearman: r = %.4f, p = %d\n',rSr,p_valSr);
fprintf('Mag vs neurons Pearson: r = %.4f, p = %d\n',rS2r,p_valS2r);
fprintf('Mag vs log neurons Spearman: r = %.4f, p = %d\n',rS,p_valS);
fprintf('Mag vs log neurons Pearson: r = %.4f, p = %d\n',rS2,p_valS2);

[rSf,p_valSf] = corr(v1t,v4t,'type','Spearman');
[rS2f,p_valS2f] = corr(v1t,v4t,'type','Pearson');
[rSrf,p_valSrf] = corr(v1t,v5t,'type','Spearman');
[rS2rf,p_valS2rf] = corr(v1t,v5t,'type','Pearson');
fprintf('Mag vs FLNe Spearman: r = %.4f, p = %d\n',rSrf,p_valSrf);
fprintf('Mag vs FLNe Pearson: r = %.4f, p = %d\n',rS2rf,p_valS2rf);
fprintf('Mag vs log FLNe Spearman: r = %.4f, p = %d\n',rSf,p_valSf);
fprintf('Mag vs log FLNe Pearson: r = %.4f, p = %d\n',rS2f,p_valS2f);

return

% --- PLOT  --------------------------------------------------
colordef_adj;
color_nocov = COLOR_ADJ_NOCOV;
color_distt = color_nocov;
color_isnan = color_nocov;
color_not_sig = COLOR_ADJ_NOSIG;

%color_distt = 0.6*[1 1 1];
%color_isnan = 0.6*[1 1 1];
color_not_sig = 0.2*[1 1 1];
fontsz = 10;
fsize = fontsz;

% Calculate rows of nans
Adj2 = Adj(ind_hum2mac,ind_hum2mac);
ind_isnan = false(length(ind_hum2mac),1);
for i = 1:length(ind_isnan)
    ind_isnan(i) = (sum(isnan(Adj2(i,:))) == length(ind_isnan));
end

% Make figure
h = figure;
set(h,'Position',round(1.2*[0 0 1.4*1080 0.5*1080]))

subplot(1,2,1); % ---------------------------------------------------------------
Adj_plt = Adj(ind_hum2mac,ind_hum2mac);
Dmat_plt = Dmat(ind_hum2mac,ind_hum2mac);
rois_plt = rois(ind_hum2mac);
for i = 1:length(rois_plt)
    rois_plt{i} = replace(rois_plt{i}(1:(end-5)), '_','/');
end
[n_roi_p,~] = size(Adj_plt);
% generate image from matrix
%Adj_plt(Dmat_plt <= dist_thresh) = nanmean(Adj_plt(:));
Adj_plt(Dmat_plt <= dist_thresh) = nan;
v = Adj_plt;
[v_n, v_m] = size(v);
map = COLOR_ADJ_CMAP; %corrcmap(100);
minv = min(v(:));
maxv = max(v(:));
ncol = size(map,1);
s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
Im = ind2rgb(s,map);
for i = 1:n_roi_p
    for j = 1:n_roi_p
        if (isnan(Adj_plt(i,j)))
            Im(i,j,:) = color_isnan;
        end
        if (Dmat_plt(i,j) <= dist_thresh)
            Im(i,j,:) = color_distt;
        end
    end
end
Im = Im(~ind_isnan,~ind_isnan,:); % remove rows of nans
rois_plt = rois_plt(~ind_isnan);

% --- Clustering (output: cluster_i) --------------------------------------
%adjct_dist_loc = adjct_dist_sav; % overwrite adjct_dist from xsub cache
adjct_dist_cl = adjct_dist_loc(ind_hum2mac,ind_hum2mac);
adjct_dist_cl = adjct_dist_cl(~ind_isnan,~ind_isnan);
%adjct_dist = adjct_dist(~ind_isnan)
n_rois_cl = length(rois_plt);
Y = zeros(1,nchoosek(n_rois_cl,2));
yc = 1;
for m1 = 1:(n_rois_cl-1)
    for m2 = (m1+1):n_rois_cl
        ad = adjct_dist_cl{m2,m1};
        if (isempty(ad))
            Y(yc) = Inf;%Inf; 600;
        else
            Y(yc) = mean(ad);
        end
        yc = yc + 1;
    end
end
roi_dist = squareform(Y);
% average   2778.832643576550 mm
% centroid  2773.791496333409 mm
% complete  2808.894068525726 mm
% median    2767.246370358684 mm
% single    2836.305543514634 mm
% ward      2778.832643576550 mm
% weighted  2778.918494966365 mm

Z = linkage(Y); %,'centroid'
cluster_i = optimalleaforder(Z,Y,'transformation','inverse'); % ,'transformation','inverse'
roi_dist = roi_dist(cluster_i,cluster_i);
%roi_dist(isinf(roi_dist)) = ;
clash = nansum(nansum(triu(roi_dist,1) - triu(roi_dist,2)));
fprintf('cluster clash: %.12f mm\n',clash)
%cluster_i = optimalleaforder(Z,Y);

% -------------------------------------------------------------------------

% Apply clustering
Im = Im(cluster_i,cluster_i,:);
rois_plt = rois_plt(cluster_i);

imagesc(Im);
%imagesc(Im, 'Parent', ax);
colormap(map);
if (minv == maxv)
    cb = colorbar;
    set(cb,'TickLength',0);
else
    cb = colorbar('Xtick',linspace(minv,maxv,n_colorbar_ticks));
    set(cb,'TickLength',0);
end
% 
% % --- Axis label staggering---------------------
% bip_labels = rois_plt;
% stagger_labels = true;
% xtick_vec = 1:length(bip_labels);
% ytick_vec = 1:length(bip_labels);
% ax1_xind = true(1,length(xtick_vec));
% ax1_yind = true(1,length(ytick_vec));
% ax2_xind = true(1,length(xtick_vec));
% ax2_yind = true(1,length(ytick_vec));
% if (stagger_labels)
%     odd_i = downsample(1:length(ax1_xind),2,1);
%     even_i = downsample(1:length(ax1_xind),2,0);
%     ax1_xind(even_i) = false;
%     ax2_xind(odd_i) = false;
%     ax1_yind(odd_i) = false;
%     ax2_yind(even_i) = false;
% end
% 
% % Show axis labels
% %ax1 = ax;
% ax = gca;
% set(ax,'xtick',xtick_vec(ax1_xind),'Ticklength',[0.001, 0.001],...
%     'xticklabel',bip_labels(ax1_xind),'fontsize',fsize,'TickDir','out');
% xtickangle(90);
% set(ax,'ytick',ytick_vec(ax1_yind),'Ticklength',[0.001, 0.001],...
%     'yticklabel',bip_labels(ax1_yind),'fontsize',fsize,'TickDir','out');
% 
% ax2 = copyobj(ax,ax.Parent);
% set(ax2,'XAxisLocation','top');
% set(ax2,'YAxisLocation','right');
% set(ax2,'xtick',xtick_vec(ax2_xind),'xticklabel',bip_labels(ax2_xind));
% set(ax2,'ytick',ytick_vec(ax2_yind),'yticklabel',bip_labels(ax2_yind));
% return

% --- END Axis label staggering---------------------

 % colormap
% colormap(ax,map)
% colormap(ax2,map)
% if (minv == maxv)
%     caxis(ax,[(minv - 1e-6), (maxv + 1e-6)]); 
%     caxis(ax2,[(minv - 1e-6), (maxv + 1e-6)]); 
% else
%     caxis(ax,[minv maxv]);
%     caxis(ax2,[minv maxv]);
% end
%return
yticks(1:length(rois_plt));
yticklabels(rois_plt);
xticks(1:length(rois_plt));
xticklabels(rois_plt);
xtickangle(90);
set(gca,'tickdir','out');
set(gca,'fontsize',fontsz);
set(gca,'TickLength',[0.001, 0.001])
if (strcmp(metric,'pcBroadband'))
    title('Functional Interactions')
else
    title(sprintf('Functional Interactions - %s',metric(3:end)))
end

subplot(1,2,2) % ---------------------------------------------------------------
Adj_plt_mSu = Adj_plt;
Adj_plt = AdjMK_bin;
Dmat_plt = Dmat(ind_hum2mac,ind_hum2mac);
rois_plt = rois(ind_hum2mac);
for i = 1:length(rois_plt)
    %rois_plt{i} = rois_plt{i}(1:(end-5));
    rois_plt{i} = replace(rois_plt{i}(1:(end-5)), '_','/');
end
[n_roi_p,~] = size(Adj_plt);
% generate image from matrix
%Adj_plt(Dmat_plt <= dist_thresh) = nanmean(Adj_plt(:));
Adj_plt(Dmat_plt <= dist_thresh) = nan;
v = Adj_plt;
[v_n, v_m] = size(v);
map = COLOR_ADJ_CMAP; %corrcmap(100);
minv = min(v(:));
maxv = max(v(:));
ncol = size(map,1);
s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
Im = ind2rgb(s,map);
for i = 1:n_roi_p
    for j = 1:n_roi_p
        if (isnan(Adj_plt_mSu(i,j)))
            Im(i,j,:) = color_isnan;
        end
        if (Dmat_plt(i,j) <= dist_thresh)
            Im(i,j,:) = color_distt;
        end
    end
end
Im = Im(~ind_isnan,~ind_isnan,:); % remove rows of nans
rois_plt = rois_plt(~ind_isnan);

% Apply clustering
Im = Im(cluster_i,cluster_i,:);
rois_plt = rois_plt(cluster_i);

imagesc(Im);
colormap(map);
if (minv == maxv)
    cb = colorbar;
    set(cb,'TickLength',0);
else
    cb = colorbar('Xtick',linspace(minv,maxv,n_colorbar_ticks));
    set(cb,'TickLength',0);
end
yticks(1:length(rois_plt));
yticklabels(rois_plt);
xticks(1:length(rois_plt));
xticklabels(rois_plt);
xtickangle(90);
set(gca,'tickdir','out');
set(gca,'fontsize',fontsz);
set(gca,'TickLength',[0.001, 0.001])
title({'Anatomical Connections - Markov Kennedy';sprintf('r = %.2f (p = %.2d)',r_bin,p_val_bin)})

%return

print(h,sprintf('figures/T8d2/Corr_mSu_MK_bin_%s_r-%i',metric,round(1000*r_bin)),fig_fmt);
if (trig_eps)
    print(h,sprintf('figures/T8d2/Corr_mSu_MK_bin_%s_r-%i',metric,round(1000*r_bin)),'-depsc');
end
close(h);





% =======================================================================



%AdjCP(ind_hum2mac,ind_hum2mac), AdjMKflne
% Make figure
h = figure;
set(h,'Position',round(1.2*[0 0 1.4*1080 0.5*1080]))

subplot(1,2,1) % ---------------------------------------------------------------
Adj_plt = AdjCP(ind_hum2mac,ind_hum2mac);
Dmat_plt = Dmat(ind_hum2mac,ind_hum2mac);
rois_plt = rois(ind_hum2mac);
for i = 1:length(rois_plt)
    rois_plt{i} = replace(rois_plt{i}(1:(end-5)), '_','/');
    %rois_plt{i} = rois_plt{i}(1:(end-5));
end
[n_roi_p,~] = size(Adj_plt);
% generate image from matrix
%Adj_plt(Dmat_plt <= dist_thresh) = nanmean(Adj_plt(:));
Adj_plt(Dmat_plt <= dist_thresh) = nan;
v = Adj_plt;
[v_n, v_m] = size(v);
map = COLOR_ADJ_CMAP; %corrcmap(100);
minv = min(v(v~=0));
maxv = max(v(v~=0));
ncol = size(map,1);
s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
Im = ind2rgb(s,map);
for i = 1:n_roi_p
    for j = 1:n_roi_p
        if (isnan(Adj_plt(i,j)))
            Im(i,j,:) = color_isnan;
        end
        if (Dmat_plt(i,j) <= dist_thresh)
            Im(i,j,:) = color_distt;
        end
        if (Adj_plt(i,j) == 0)
            Im(i,j,:) = color_not_sig;
        end
    end
end
Im = Im(~ind_isnan,~ind_isnan,:); % remove rows of nans
rois_plt = rois_plt(~ind_isnan);

% Apply clustering
Im = Im(cluster_i,cluster_i,:);
rois_plt = rois_plt(cluster_i);

imagesc(Im);
colormap(map);
if (minv == maxv)
    cb = colorbar;
    set(cb,'TickLength',0);
else
    cb = colorbar('Xtick',linspace(minv,maxv,n_colorbar_ticks));
    set(cb,'TickLength',0);
end
%caxis([min(Adj_plt(:)) max(Adj_plt(:))])
caxis([minv maxv])
yticks(1:length(rois_plt));
yticklabels(rois_plt);
xticks(1:length(rois_plt));
xticklabels(rois_plt);
xtickangle(90);
set(gca,'tickdir','out');
set(gca,'fontsize',fontsz);
set(gca,'TickLength',[0.001, 0.001])
if (strcmp(metric,'pcBroadband'))
    title('Functional Interactions')
else
    title(sprintf('Functional Interactions - %s',metric(3:end)))
end

subplot(1,2,2) % ---------------------------------------------------------------
Adj_plt_mSu = Adj_plt;
AdjMKneuronsP = AdjMKneurons;
AdjMKneuronsP((AdjMKneurons == 0)) = NaN;
AdjMKneuronsP = log10(AdjMKneuronsP);
AdjMKneuronsP((AdjMKneurons == 0)) = 0;
Adj_plt = AdjMKneuronsP; %Adj_plt = AdjMKneurons;
Dmat_plt = Dmat(ind_hum2mac,ind_hum2mac);
rois_plt = rois(ind_hum2mac);
for i = 1:length(rois_plt)
    rois_plt{i} = replace(rois_plt{i}(1:(end-5)), '_','/');
    %rois_plt{i} = rois_plt{i}(1:(end-5));
end
[n_roi_p,~] = size(Adj_plt);
% generate image from matrix
%Adj_plt(Dmat_plt <= dist_thresh) = nanmean(Adj_plt(:));
Adj_plt(Dmat_plt <= dist_thresh) = nan;
v = Adj_plt;
[v_n, v_m] = size(v);
map = COLOR_ADJ_CMAP; %corrcmap(100);
minv = min(v(v~=0));
maxv = max(v(v~=0));
ncol = size(map,1);
s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
Im = ind2rgb(s,map);
for i = 1:n_roi_p
    for j = 1:n_roi_p
        if (Adj_plt(i,j) == 0)
            Im(i,j,:) = color_not_sig;
        end
        if (isnan(Adj_plt_mSu(i,j)))
            Im(i,j,:) = color_isnan;
        end
        if (Dmat_plt(i,j) <= dist_thresh)
            Im(i,j,:) = color_distt;
        end
        
    end
end
Im = Im(~ind_isnan,~ind_isnan,:); % remove rows of nans
rois_plt = rois_plt(~ind_isnan);

% Apply clustering
Im = Im(cluster_i,cluster_i,:);
rois_plt = rois_plt(cluster_i);

imagesc(Im);
colormap(map);
if (minv == maxv)
    cb = colorbar;
    set(cb,'TickLength',0);
else
    cb = colorbar('Xtick',linspace(minv,maxv,n_colorbar_ticks));
    set(cb,'TickLength',0);
end
%caxis([min(Adj_plt(:)) max(Adj_plt(:))])
caxis([minv maxv])
yticks(1:length(rois_plt));
yticklabels(rois_plt);
xticks(1:length(rois_plt));
xticklabels(rois_plt);
xtickangle(90);
set(gca,'tickdir','out');
set(gca,'fontsize',fontsz);
set(gca,'TickLength',[0.001, 0.001])
title({'Anatomical Connections - Markov Kennedy';sprintf('r = %.2f (p = %.2d)',rS,p_valS)})

%return

print(h,sprintf('figures/T8d2/Corr_mSu_MK_f_%s_r-%i',metric,round(1000*rS)),fig_fmt);
if (trig_eps)
    print(h,sprintf('figures/T8d2/Corr_mSu_MK_f_%s_r-%i',metric,round(1000*rS)),'-depsc');
end
close(h);



end