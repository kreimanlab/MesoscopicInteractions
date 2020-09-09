close all;
clear;
rng shuffle;

metricpp = 'pcBroadband';
n_perm = 10000;
perm_alpha_sec = 12;
cp_thresh_override = 0.05;
p_val = 0.01;

% Fast i/o definitions
dir_art = '/media/jerry/KLAB101/h5_notch20/art_nosz';
dir_res = '/media/jerry/KLAB101/results/coh_w10';
dir_cor = '/media/jerry/internal/data/coreg';
dir_cache = './cache';
dir_stamp = '/media/klab/internal/data/stamps';
dir_video = '/media/klab/internal/data/videos';
%subjects_dirL = '/mnt/cuenap_ssd/coregistration';

% Slow i/o definitions
dir_h5 = '/media/klab/KLAB101/h5_notch20';
%dir_h5L = '/mnt/cuenap/data/h5_notch20';

Subjects = {'sub3'};
r_samp_const = [107865601, 72056321, 110476801, 25121281];
coh_const = [0.535, 0.364, 0.280, 0.093];
coh_annot = {'No interaction','Fig. 1','Strong interaction','Medium interaction'};
% High - 0.535, str996
%sub3__35_superiortemporal_PT7-PT8__84_parsopercularis_RP55-RP56__29mm_ct405_mag221_coh1t172_str996_Br535_Th727_Al286_Be154_Ga536_Samp107865601_BrN445_ThN467_AlN132_BeN574_GaN575_SampN107865601.eps
% Original - 0.364, str987
%sub3__35_superiortemporal_PT7-PT8__84_parsopercularis_RP55-RP56__29mm_ct405_mag221_coh1t172_str987_Br364_Th662_Al277_Be149_Ga394_Samp72056321_BrN337_ThN405_AlN162_BeN316_GaN440_SampN72056321
% Low - 0.280 - str920 
%sub3__35_superiortemporal_PT7-PT8__84_parsopercularis_RP55-RP56__29mm_ct405_mag221_coh1t172_str920_Br280_Th663_Al143_Be293_Ga211_Samp110476801_BrN247_ThN356_AlN218_BeN154_GaN278_SampN110476801
% Not sig - 0.093 - str100 nosig
%sub3__35_superiortemporal_PT7-PT8__84_parsopercularis_RP55-RP56__29mm_ct405_mag221_coh1t172_str100_Br93_Th237_Al77_Be98_Ga82_Samp25121281_BrN112_ThN136_AlN112_BeN70_GaN126_SampN25121281

bchan1_const = 35;
bchan2_const = 84;

% Sort time
[~,sIdx] = sort(r_samp_const);
r_samp_const = r_samp_const(sIdx);
coh_const = coh_const(sIdx);


system('mkdir figures/T1d2');

iM = 1;
metric = 1;
sid = Subjects{1};
fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
fn_dist = sprintf('%s/%s_dists-%s-%i.mat',dir_res,sid,metric,n_perm);
fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
fn_perm = sprintf('/mnt/cuenap2/data/results/coh_w10/%s_perm-%s-%i.h5',sid,metric,n_perm);
fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
fn_coreg = sprintf('%s/%s/label/all_parcellation.mat',dir_cor,sid);
fn_cache = sprintf('%s/xsub_out_%s_%i.mat',dir_cache,sid,iM);

ecog = H5eeg(fn_h5);
Ca = load(fn_cache);

% build raw voltages
n_seg = length(r_samp_const); % Number of 10-second segments to plot
n_plot_samp = round(ecog.fs)*Ca.w;
V1b = zeros(n_seg,n_plot_samp);
V2b = zeros(n_seg,n_plot_samp);
b1c1 = ecog.bip(bchan1_const,1);
b1c2 = ecog.bip(bchan1_const,2);
b2c1 = ecog.bip(bchan2_const,1);
b2c2 = ecog.bip(bchan2_const,2);
Timetxt = cell(1,n_seg);
t = linspace(0,Ca.w,n_plot_samp);
parfor i = 1:n_seg
    v_b1c1 = h5read(ecog.filename,'/h5eeg/eeg',[b1c1 r_samp_const(i)],[1 n_plot_samp]);
    v_b1c2 = h5read(ecog.filename,'/h5eeg/eeg',[b1c2 r_samp_const(i)],[1 n_plot_samp]);
    v_b2c1 = h5read(ecog.filename,'/h5eeg/eeg',[b2c1 r_samp_const(i)],[1 n_plot_samp]);
    v_b2c2 = h5read(ecog.filename,'/h5eeg/eeg',[b2c2 r_samp_const(i)],[1 n_plot_samp]);
    v_b1 = v_b1c1 - v_b1c2;
    v_b2 = v_b2c1 - v_b2c2;
    v_b1 = v_b1 - mean(v_b1);
    v_b2 = v_b2 - mean(v_b2);
    V1b(i,:) = v_b1;
    V2b(i,:) = v_b2;
    
    % Get raw times
    %fn_ca = sprintf('%s/%s',dir_cache,'figure_t1d2_one_Timetxt.mat');
    ttraw = ecog.getTimeFromSample(dir_stamp,r_samp_const(i));
    Timetxt{i} = ttraw(1:(end-1));
end


%%

% Set plot parameters
between_sec = 2;
text_yoffset = 10; % date text
text_xoffset = 0.2;
text_coh_yoffset = 140; % coherence value text
text_coh_xoffset = text_xoffset;
text_annot_yoffset = 250; % annotation text
text_annot_xoffset = text_xoffset;
ytxt = sprintf('IFP (\\muV)');
%xtxt = 'Time (seconds)';
xtxt = [];

% plot 
vmax = max([V1b(:), V2b(:)]);
vmin = min([V1b(:), V2b(:)]);
vnearest = 100;
vmax = ceil(vmax/vnearest)*vnearest;
vmin = floor(vmin/vnearest)*vnearest;

PLOT_W_PIX = 1080;
PLOT_H_PIX = 1080/2;
col_faint = 0.5*[1 1 1];
h = figure('visible','on','Position',[0  0  PLOT_W_PIX PLOT_H_PIX]);

subplot(2,1,1)
xtick_t = zeros(1,n_seg);
for i = 1:n_seg
    plot(t + (between_sec+Ca.w)*(i-1),V1b(i,:),'black'); hold on
    xtick_t(i) = (i-1)*((Ca.w)+between_sec);
    text(text_xoffset + (i-1)*((Ca.w)+between_sec),vmin(1)-text_yoffset,Timetxt{i},'HorizontalAlignment','left','VerticalAlignment','top');
end
box off;
set(gca,'TickDir','out');
ylabel(ytxt);
xlabel(xtxt);
arrow_x_w = 0.04;
arrow_y_w = 0.05;
% Arrows
axp = get(gca,'Position');
xs=axp(1);
xe=axp(1)+axp(3)+arrow_x_w;
ys=axp(2);
ye=axp(2)+axp(4)+arrow_y_w;
annotation('arrow', [xs xe],[ys ys]);
axis([0 n_seg*(Ca.w)+(n_seg-1)*between_sec vmin(1) vmax(1)])
%xticks([]);
xticks(xtick_t);
xticklabels('');

% plot coherence values
for i = 1:n_seg
    cohtxt = sprintf('\\rho=%.2f',coh_const(i));
    text(text_coh_xoffset + (i-1)*((Ca.w)+between_sec),vmin(1)-text_coh_yoffset,cohtxt,'HorizontalAlignment','left','VerticalAlignment','top');
    text(text_annot_xoffset+ (i-1)*((Ca.w)+between_sec),vmin(1)-text_annot_yoffset,coh_annot{i})
end


subplot(2,1,2)
for i = 1:n_seg
    plot(t + (between_sec+Ca.w)*(i-1),V2b(i,:),'black'); hold on
    text(text_xoffset + (i-1)*((Ca.w)+between_sec),vmin(2)-text_yoffset,Timetxt{i},'HorizontalAlignment','left','VerticalAlignment','top');
end
box off;
set(gca,'TickDir','out');
ylabel(ytxt);
xlabel(xtxt);
arrow_x_w = 0.04;
arrow_y_w = 0.05;
% Arrows
axp = get(gca,'Position');
xs=axp(1);
xe=axp(1)+axp(3)+arrow_x_w;
ys=axp(2);
ye=axp(2)+axp(4)+arrow_y_w;
annotation('arrow', [xs xe],[ys ys]);
axis([0 n_seg*(Ca.w)+(n_seg-1)*between_sec vmin(2) vmax(2)])
%xticks([]);
xticks(xtick_t);
xticklabels('');

ofstr = sprintf('%s_%i_%i_t1d2_one',sid,bchan1_const,bchan2_const);
print(h,sprintf('figures/T1d2/%s',ofstr),'-dpng');
print(h,sprintf('figures/T1d2/%s',ofstr),'-depsc');
close(h);

% Print finish message
fprintf('Done.\n')
