close all
clear
rng shuffle;


trig_plot_flne = true;


%metric = 'pcBroadband';
n_perm = 10000;
perm_alpha_sec = 12;
cp_thresh_override = 0.05; % 0.01, new 0.05
fig_fmt = '-dpng';
trig_eps = true;
trig_mag_no_cp = true;
n_colorbar_ticks = 5;
system('mkdir figures');
system('mkdir figures/T8d1');

% Fast i/o definitions
dir_artL = '/media/klab/internal/data/h5_notch20/art';
dir_resL = '/media/klab/internal/data/results/coh_w10';
dir_corL = '/media/klab/internal/data/coreg';
dir_cacheL = './cache';
subjects_dirL = '/mnt/cuenap_ssd/coregistration';

% Slow i/o definitions
dir_h5L = '/media/klab/KLAB101/h5_notch20';

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'}; %'pcDelta',

% Patients
Subjects = {'mSu'};
sid = Subjects{1};
%iM = 1;

for iM = 1:5 %1:length(metrics)
    
metric = metrics{iM};

fn_artL = sprintf('%s/%s_art.h5',dir_artL,sid);
fn_distL = sprintf('%s/%s_dists-%s-%i.mat',dir_resL,sid,metric,n_perm);
fn_graphL = sprintf('%s/%s_graph-%s.h5',dir_resL,sid,metric);
fn_permL = sprintf('/mnt/cuenap2/data/results/coh_w10/%s_perm-%s-%i.h5',sid,metric,n_perm);
fn_h5L = sprintf('%s/%s.h5',dir_h5L,sid);
fn_coregL = sprintf('%s/%s/label/all_parcellation.mat',dir_corL,sid);
fn_cacheL = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheL,sid,iM);

% Load mSu cache
load(fn_cacheL);

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


% Map onto atlas
%dist_thresh_mac = dist_thresh / 2; % macaque grids have 5mm spacing instead of 10mm
d_i = Dmats > dist_thresh;
%fprintf('Distance threshold fraction removed: %.4f\n',1-sum(d_i)/length(d_i))
b1c1 = zeros(ecog.n_comb,1);
b1c2 = zeros(ecog.n_comb,1);
b2c1 = zeros(ecog.n_comb,1);
b2c2 = zeros(ecog.n_comb,1);

rois = C1.AtlROIs{7}.LH.struct_names;
n_rois = length(rois);
AdjAtl = cell(n_rois,n_rois);
AdjAtlN = cell(n_rois,n_rois);
adjct_dist = cell(n_rois,n_rois);
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
                    
                    adjct_dist{ki1,ki2} = [adjct_dist{ki1,ki2}; Dmats(k)];
                    adjct_dist{ki2,ki1} = [adjct_dist{ki2,ki1}; Dmats(k)];
% 
%                     AdjAtl_sid{ki1,ki2} = [AdjAtl_sid{ki1,ki2},subject_i];
%                     AdjAtl_sid{ki2,ki1} = [AdjAtl_sid{ki2,ki1},subject_i];
                end
            end

        catch
            %fprintf('ct-atl map skip.\n');
        end
    end
end

% Calculate final functional interaction matrix
Adj = nan(n_rois,n_rois);
AdjMag = nan(n_rois,n_rois);
AdjCP = nan(n_rois,n_rois);
N_bchan = nan(n_rois,n_rois);
cp_thresh = cp_thresh_override;

for i1 = 1:n_rois
    for i2 = 1:n_rois
        AA = AdjAtl{i1,i2};
        if (~ isempty(AA))
            N_bchan(i1,i2) = length(AA);
            frac_cp = sum(AA ~= 0)/length(AA);
            AdjCP(i1,i2) = frac_cp;
            if (frac_cp >= cp_thresh)
                Adj(i1,i2) = 1;
                %AdjMag(i1,i2) = mean(AA(frac_cp > cp_thresh));
                AdjMag(i1,i2) = mean(AA(AA ~= 0));
            else
                Adj(i1,i2) = 0;
                AdjMag(i1,i2) = 0;
            end
        end
    end
end


if (trig_mag_no_cp)
    AdjCP = AdjMag;
end



frac_xsub_msu = sum(Adj(~isnan(Adj))) /numel((~isnan(Adj)));
fprintf('%s - mSu fraction of ROIs significant: %.4f\n',metric,frac_xsub_msu)


% compare to macaque
load('../coreg/AdjMKFV.mat')
AdjMKflne_orig = AdjMKflne;
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
v1 = AdjMK_bin_f(~isnan(AdjMK_bin_f));
v2 = AdjMK_bin(~isnan(AdjMK_bin_f));

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
%AdjMK_f = Adj(ind_hum2mac,ind_hum2mac);

AdjMK_f = AdjMag(ind_hum2mac,ind_hum2mac);
%v1 = AdjMK_f(~isnan(AdjMK_f));
%v2 = AdjMKneurons(~isnan(AdjMK_f));
%v2 = AdjMKflne(~isnan(AdjMK_f));
%v2(isnan(v2)) = 0;

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
fprintf('Mag vs FLNe Spearman: r = %.4f, p = %d, n=%i\n',rSrf,p_valSrf,length(v1t));
fprintf('Mag vs FLNe Pearson: r = %.4f, p = %d, n=%i\n',rS2rf,p_valS2rf,length(v1t));
fprintf('Mag vs log FLNe Spearman: r = %.4f, p = %d, n=%i\n',rSf,p_valSf,length(v1t));
fprintf('Mag vs log FLNe Pearson: r = %.4f, p = %d, n=%i\n',rS2f,p_valS2f,length(v1t));

% Calculate lognormal fit
X_ecog = v1t;
X_ecog = X_ecog(X_ecog>0);

%X_flne = v5t;
%v6t = [];
%for ii1 = 1:(AdjMKflne_orig;
%X_flne = AdjMKflneRaw(:);

cIdx = ~all(AdjMKflne_orig==0);
AdjMKflne_orig_comp = AdjMKflne_orig(cIdx,cIdx);
X_flne = AdjMKflne_orig(:);
X_flne = X_flne(X_flne>0);


% 
% % randomly shuffle
% n_perm = 1000;
% fprintf('[*] Shuffling correlation with n = %i..\n',n_perm);
% rng('shuffle');
% R = [];
% n_v1t = length(v1t);
% %R = zeros(4,n_perm);
% % P = zeros(4,n_perm);
% for i = 1:n_perm
%     sIdx = randperm(n_v1t);
%     sIdx2 = randperm(n_v1t);
%     sIdx3 = randperm(n_v1t);
%     v1t = v1t(sIdx);
%     v2t = v2t(sIdx2);
%     v3t = v3t(sIdx3);
%     
%     [rSt,p_valSt] = corr(v1t,v2t,'type','Spearman');
%     [rS2t,p_valS2t] = corr(v1t,v2t,'type','Pearson');
%     [rSrt,p_valSrt] = corr(v1t,v3t,'type','Spearman');
%     [rS2rt,p_valS2rt] = corr(v1t,v3t,'type','Pearson');
%     
%     % save
% %     R(1,i) = rSt;
% %     R(2,i) = rS2t;
% %     R(3,i) = rSrt;
% %     R(4,i) = rS2rt;
%     R = [R, [rSt; rS2t; rSrt; rS2rt] ];
%     
% %     P(1,i) = p_valS;
% %     P(2,i) = p_valS2;
% %     P(3,i) = p_valSr;
% %     P(4,i) = p_valS2r;
% end
% 
% p_S = sum(sort(R(1,:)) > rS)/n_perm;
% p_S2 = sum(sort(R(1,:)) > rS2)/n_perm;
% p_Sr = sum(sort(R(1,:)) > rSr)/n_perm;
% p_S2r = sum(sort(R(1,:)) > rS2r)/n_perm;
% 
% fprintf('Mag vs neurons Spearman: r = %.4f, p = %d\n',rSr,p_Sr);
% fprintf('Mag vs neurons Pearson: r = %.4f, p = %d\n',rS2r,p_S2r);
% fprintf('Mag vs log neurons Spearman: r = %.4f, p = %d\n',rS,p_S);
% fprintf('Mag vs log neurons Pearson: r = %.4f, p = %d\n',rS2,p_S2);


% --- PLOT  --------------------------------------------------
colordef_adj;
color_nocov = COLOR_ADJ_NOCOV;
color_distt = color_nocov;
color_isnan = color_nocov;
color_not_sig = COLOR_ADJ_NOSIG;
    
%color_distt = 0.6*[1 1 1];
%color_isnan = 0.6*[1 1 1];
%color_not_sig = 0.2*[1 1 1];
fontsz = 10;

% Calculate rows of nans
Adj2 = Adj(ind_hum2mac,ind_hum2mac);
ind_isnan = false(length(ind_hum2mac),1);
for i = 1:length(ind_isnan)
    ind_isnan(i) = (sum(isnan(Adj2(i,:))) == length(ind_isnan));
end

% Make figure
h = figure('visible', 'off');
set(h,'Position',round(1.2*[0 0 1.4*1080 0.5*1080]))

subplot(1,2,1) % ---------------------------------------------------------------
Adj_plt = Adj(ind_hum2mac,ind_hum2mac);
Adj_plt_clean = AdjMag(ind_hum2mac,ind_hum2mac);
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
map = COLOR_ADJ_CMAP; %inferno(100);
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
        if (Adj(i,j) == 0)
            Im(i,j,:) = color_not_sig;
        end
    end
end
Im = Im(~ind_isnan,~ind_isnan,:); % remove rows of nans
rois_plt = rois_plt(~ind_isnan);
Adj_plt_clean = Adj_plt_clean(~ind_isnan,~ind_isnan);

% --- Clustering (output: cluster_i) --------------------------------------
adjct_dist_cl = adjct_dist(ind_hum2mac,ind_hum2mac);
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

% Old linkage
%Z = linkage(Y); %,'centroid'
% New linkage, ignoring NaNs
Z = linkage(Adj_plt_clean,'ward',@nanDist4clustering);

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
map = COLOR_ADJ_CMAP; %inferno(100);
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
        if (Adj(i,j) == 0)
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

print(h,sprintf('figures/T8d1/Corr_mSu_MK_bin_%s_r-%i',metric,round(1000*r_bin)),fig_fmt);
if (trig_eps)
    print(h,sprintf('figures/T8d1/Corr_mSu_MK_bin_%s_r-%i',metric,round(1000*r_bin)),'-depsc');
end
close(h);





% =======================================================================



%AdjCP(ind_hum2mac,ind_hum2mac), AdjMKflne
% Make figure
h = figure('visible', 'off');
set(h,'Position',round(1.2*[0 0 1.4*1080 0.5*1080]))

subplot(1,2,1) % ---------------------------------------------------------------
%Adj_plt = AdjCP(ind_hum2mac,ind_hum2mac);
Adj_plt = AdjMag(ind_hum2mac,ind_hum2mac);
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
map = COLOR_ADJ_CMAP; %inferno(100);
%minv = min(v(:));
%maxv = max(v(:));
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
Adj_plt_mSu_red = Adj_plt(~ind_isnan,~ind_isnan);
rois_plt = rois_plt(~ind_isnan);

% Apply clustering
Im = Im(cluster_i,cluster_i,:);
rois_plt = rois_plt(cluster_i);
Adj_plt_mSu_red = Adj_plt_mSu_red(cluster_i,cluster_i);
Adj_plt_mSu_red_rois = rois_plt;

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




%============
% Plot log-normal fit
AdjV = nan(nchoosek(length(Adj_plt),2),1);
count = 1;
for i1 = 1:(length(Adj_plt)-1)
    for i2 = (i1+1):length(Adj_plt)
        AdjV(count) = Adj_plt(i1,i2);
        count = count + 1;
    end
end
AdjV = AdjV(AdjV>0);
hh = figure('visible','on');
set(hh,'Position',[0 0 150 150]);
%h = histogram(log10(AdjV),'Normalization','pdf','DisplayStyle','bar','FaceColor',0.5*[1 1 1]);
hm = histogram(AdjV,'Normalization','pdf','DisplayStyle','bar','FaceColor',0.5*[1 1 1]);
x_cen = linspace(hm.BinLimits(1),hm.BinLimits(2),100); %h.BinEdges + h.BinWidth/2;
%y_norm = normpdf(x_cen,(nanmean(log10(AdjV))),(nanstd(log10(AdjV))));
y_norm = lognpdf(x_cen,(nanmean(log(AdjV))),(nanstd(log(AdjV))));
xticks([x_cen(1), 0.5*(x_cen(1)+x_cen(end)), x_cen(end)]);
hold on;
plot(x_cen,y_norm,'black-');
xlabel('Coherence');
ylabel('pdf')
set(gca,'TickDir','out');
box off;

% --- KS test ---
Xt = (log10(AdjV));
[~,p,ksstat,cv] = kstest((Xt-mean(Xt))/std(Xt));
fprintf('[*] Log(C) kstest p=%.4d, n=%i\n',p,length(Xt));
fprintf('\tmin:%.8f, max:%.8f\n',min(Xt),max(Xt));

print(hh,sprintf('figures/T8d1/Corr_mSu_MK_f_%s_r-%i_hist',metric,round(1000*rSrf)),'-dsvg');
%saveas(gcf,sprintf('figures/T20_lognormal/metric%i_atl%i_hist',iM,atl),'epsc');
close(hh);
%return
%============






subplot(1,2,2) % ---------------------------------------------------------------
Adj_plt_mSu = Adj_plt;
% AdjMKneuronsP = AdjMKneurons;
% AdjMKneuronsP((AdjMKneurons == 0)) = NaN;
% AdjMKneuronsP = log10(AdjMKneuronsP);
% AdjMKneuronsP((AdjMKneurons == 0)) = 0;
% 
% AdjMKneuronsP = AdjMKflne;
% AdjMKneuronsP((AdjMKflne == 0)) = NaN;
% AdjMKneuronsP = log10(AdjMKneuronsP);
% AdjMKneuronsP((AdjMKflne == 0)) = 0;

%AdjMKflneRaw2 = AdjMKflneRaw;
%AdjMKflneRaw2((AdjMKflneRaw == 0)) = NaN;
Adj_plt = AdjMKneuronsP; % AdjMKflneRaw2

% Preserve raw connection weights for flne
Adj_plt_MKraw = AdjMKflne_orig;
Adj_plt_MKraw = (Adj_plt_MKraw + (Adj_plt_MKraw'))/2;
%Adj_plt_MKraw(all(Adj_plt_MKraw == 0),:) = NaN;
%Adj_plt_MKraw(:,all(Adj_plt_MKraw == 0)) = NaN;


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
map = COLOR_ADJ_CMAP; %inferno(100);
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
Adj_plt_MK = Adj_plt(~ind_isnan,~ind_isnan);
Adj_plt_MKraw = Adj_plt_MKraw(~ind_isnan,~ind_isnan);

% Apply clustering
Im = Im(cluster_i,cluster_i,:);
rois_plt = rois_plt(cluster_i);
Adj_plt_MK = Adj_plt_MK(cluster_i,cluster_i);
Adj_plt_MKraw = Adj_plt_MKraw(cluster_i,cluster_i);
Adj_plt_MK_rois = rois_plt;
Adj_plt_MKraw_rois = rois_plt;

imagesc(Im);
if (minv == maxv)
    cb = colorbar;
    set(cb,'TickLength',0);
else
    cb = colorbar('Xtick',linspace(minv,maxv,n_colorbar_ticks));
    set(cb,'TickLength',0);
end
caxis([minv maxv])
yticks(1:length(rois_plt));
yticklabels(rois_plt);
xticks(1:length(rois_plt));
xticklabels(rois_plt);
xtickangle(90);
set(gca,'tickdir','out');
set(gca,'fontsize',fontsz);
set(gca,'TickLength',[0.001, 0.001])
title({'Anatomical Connections - Markov Kennedy';sprintf('r = %.2f (p = %.2d)',rSrf,p_valSrf)})

%return

print(h,sprintf('figures/T8d1/Corr_mSu_MK_f_%s_r-%i',metric,round(1000*rSrf)),fig_fmt);
if (trig_eps)
    print(h,sprintf('figures/T8d1/Corr_mSu_MK_f_%s_r-%i',metric,round(1000*rSrf)),'-depsc');
end
close(h);

% Save dendrogram
h = figure('visible','off','Position',[0 0 100 400]);
hD = dendrogram(Z,0,'Reorder',fliplr(cluster_i),'Orientation','right','ColorThreshold',0);
for ihD = 1:length(hD)
    hD(ihD).Color = 0*[1 1 1];
end
%yticklabels(flipud(rois_plt));
axis tight;
axis off;
print(h,sprintf('figures/T8d1/Corr_mSu_MK_f_%s_r-%i_dendro',metric,round(1000*rSrf)),'-depsc');
print(h,sprintf('figures/T8d1/Corr_mSu_MK_f_%s_r-%i_dendro',metric,round(1000*rSrf)),'-dpng');
close(h);

% Save matrices
save(sprintf('figure_t8d1_adj_%i',iM),'-v7.3','-nocompression','Adj_plt_mSu_red','Adj_plt_mSu_red_rois','Adj_plt_MKraw','Adj_plt_MKraw_rois');


% Save roi names
ofn = sprintf(sprintf('figures/T8d1/Corr_mSu_MK_f_%s_r-%i.txt',metric,round(1000*rSrf)));
of = fopen(ofn,'w');
for irn = 1:length(rois_plt)
    fprintf(of,'%s\n',rois_plt{irn});
end
fclose(of);


end