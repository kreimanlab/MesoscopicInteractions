

%figure_dist_mag
clear
close all;

system('mkdir figures/dist_mag');

dir_art = '/media/klab/internal/data/h5_notch20/art';
dir_cor = '/media/klab/internal/data/coreg';
dir_res = '/home/jerry/data/results/coh_w10';
setenv('SUBJECTS_DIR',dir_cor);
dir_cache = './cache';
dir_h5 = '/media/klab/internal/data/h5_notch20/';

metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
atl = 2; % 2 - DK
size_font = 14;

Subjects = {'mSu','sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};


for iM = [1 5]
    D = [];
    C = [];
    CT = [];
    mag_thresh = 0.3;%0.3581;
    S = [];
    R = [];
    
    
    for iSub = 1:length(Subjects)
        sid = Subjects{iSub};
        metric = metricsp{iM};
        Ca = load(sprintf('%s/xsub_out_%s_%i',dir_cache,sid,iM));
        iscov = (Ca.Dmats > Ca.dist_thresh) & (Ca.ct > Ca.ct_thresh);
        if (sum(iscov) ~= 0)
            [r,p] = corr(Ca.Dmats(iscov),Ca.mag(iscov),'type','Pearson');
            fprintf('[%i] Pearson r: %.3f, p: %.3d\n',iSub-1,r,p);
            if (iSub == 1)
                % mSu
                % Adjust units given interelectrode distance is 5mm
                Dsu = Ca.Dmats(iscov); % * (1/10) * (5/1.1);
                Csu = Ca.mag(iscov);
                CTsu = Ca.ct(iscov);
            else
                % humans
                D = [D; Ca.Dmats(iscov)];
                C = [C; Ca.mag(iscov)];
                CT = [CT; Ca.ct(iscov)];
                
                % store subject number
                S = [S; (iSub-1)*ones(sum(iscov),1)];
                
                % store region
                bchans = [double(Ca.chan1(iscov)+1),double(Ca.chan2(iscov)+1)];
                % convert bchan to region names
                rois = zeros(size(bchans));
                alabels = Ca.C.AtlLabels{atl};
                aRois = Ca.C.AtlROIs{atl}.LH.struct_names;
                for i = 1:sum(iscov)
                    b1c1 = Ca.ecog.bip(bchans(i,1),1);
                    b2c1 = Ca.ecog.bip(bchans(i,2),1);
                    
                    roi_b1c1 = alabels{b1c1};
                    roi_b2c1 = alabels{b2c1};
                    
                    
                    rois(i,1) = find(strcmp(roi_b1c1,aRois),1);
                    rois(i,2) = find(strcmp(roi_b1c1,aRois),1);
                end
                R = [R; rois];
            end
%             [r,p] = corr(Ca.Dmats(iscov),Ca.mag(iscov),'type','Spearman');
%             fprintf('[%i] Spearman r: %.3f, p: %.3d\n',iSub-1,r,p);
        end
        %return
    end
    
    h = figure;
    subplot(2,1,1)
    plot(D,C,'.','color',0*[1 1 1],'MarkerSize',1); hold on;
    set(gca,'tickdir','out');
    set(gca,'fontsize',size_font);
    box off;
    axis tight;
    xlabel('Distance (mm)')
    ylabel('Coherence')
    [r,p] = corr(D,C,'type','Pearson');
    fprintf('\t[dist-mag] Pearson r: %.3f, p: %.3d\n',r,p);
    [r,p] = corr(D,C,'type','Spearman');
    fprintf('\t[dist-mag] Spearman r: %.3f, p: %.3d\n',r,p);
    [r,p] = partialcorr(D,C,[S,R],'type','Spearman');
    fprintf('\t[dist-mag] Partial Spearman r | subject number, roi: %.3f, p: %.3d\n',r,p);
    ax = gca;
    YLim = ax.YLim;
    
    subplot(2,1,2)
    plot(Dsu,Csu,'.','color',0*[1 1 1],'MarkerSize',1);
    set(gca,'tickdir','out');
    set(gca,'fontsize',size_font);
    box off;
    axis([min(Dsu) max(Dsu) YLim(1) YLim(2)]);
    xlabel('Distance (mm)')
    ylabel('Coherence')
    [r,p] = corr(Dsu,Csu,'type','Pearson');
    fprintf('\t[mSu: dist-mag] Pearson r: %.3f, p: %.3d\n',r,p);
    [r,p] = corr(Dsu,Csu,'type','Spearman');
    fprintf('\t[mSu: dist-mag] Spearman r: %.3f, p: %.3d\n',r,p);
    
    % Find patterns in threshold
    fprintf('[!] Coherence threshold: %.2f\n',mag_thresh);
    pass_thresh = C > mag_thresh;
    n_sub_pass = length(unique(S(pass_thresh)));
    fprintf('\tNumber of subjects above threshold: %i of %i\n',n_sub_pass,length(Subjects)-1);
    n_roi_all = length(unique(reshape(R,[],1)));
    n_roi_pass = length(unique(reshape(R(pass_thresh,:),[],1)));
    fprintf('\tNumber of rois above threshold: %i of %i\n',n_roi_pass,n_roi_all);

    print(h,sprintf('figures/dist_mag/dist-mag_%s',metric(3:end)),'-depsc');
    
    
    h = figure;
    subplot(2,1,1)
    plot(CT,C,'.','color',0*[1 1 1],'MarkerSize',1);
    set(gca,'tickdir','out');
    set(gca,'fontsize',size_font);
    box off;
    axis([min(CT) max(CT) min([C; Csu]) max(C)]);
    ax = gca;
    YLim = ax.YLim;
    xlabel('Consistency across time')
    ylabel('Coherence')
    [r,p] = corr(CT,C,'type','Pearson');
    fprintf('\t[CT-mag] Pearson r: %.3f, p: %.3d\n',r,p);
    [r,p] = corr(CT,C,'type','Spearman');
    fprintf('\t[CT-mag] Spearman r: %.3f, p: %.3d\n',r,p);
    [r,p] = partialcorr(CT,C,[S,R],'type','Spearman');
    fprintf('\t[CT-mag] Partial Spearman r | subject number, roi: %.3f, p: %.3d\n',r,p);
    
    subplot(2,1,2)
    plot(CTsu,Csu,'.','color',0*[1 1 1],'MarkerSize',1);
    set(gca,'tickdir','out');
    set(gca,'fontsize',size_font);
    box off;
    %axis tight;
    axis([min(CTsu) max(CTsu) YLim(1) YLim(2)]);
    xlabel('Consistency across time')
    ylabel('Coherence')
    [r,p] = corr(CTsu,Csu,'type','Pearson');
    fprintf('\t[mSu: CT-mag] Pearson r: %.3f, p: %.3d\n',r,p);
    [r,p] = corr(CTsu,Csu,'type','Spearman');
    fprintf('\t[mSu: CT-mag] Spearman r: %.3f, p: %.3d\n',r,p);
    
    print(h,sprintf('figures/dist_mag/ct-mag_%s',metric(3:end)),'-depsc');
    
    %%
    % Show patient consistency
    Stats = [];
    for iSub = 1:length(Subjects)
        sIdx = (S == iSub);
        if (sum(sIdx) > 0)
            [r,p] = corr(C(sIdx),D(sIdx),'type','Pearson');
            [rS,pS] = corr(C(sIdx),D(sIdx),'type','Spearman');
            fprintf('%i,%s,%.3f,%.3d\n',iSub,Subjects{iSub},rS,pS);
            Stats = [Stats; [rS,pS]];
        end
    end
    
    return
end

% 
% % Pick coherence threshold
% 

h = figure;
[f,x] = hist(C,50);
f_norm = f/trapz(x,f);
plot(x,f_norm,'black-'); hold on;
plot([mag_thresh mag_thresh],[min(f_norm) max(f_norm)],'--','Color',0.5*[1 1 1]); hold on;
text([mag_thresh],[max(f_norm)*0.9],sprintf('%.2f',mag_thresh),'fontsize',size_font);
set(gca,'tickdir','out');
set(gca,'fontsize',size_font)
box off;
axis([min(x) max(x) min(f_norm) max(f_norm)]);
xlabel(sprintf('%s Coherence',metric(3:end)));
ylabel('Probability')
print(h,sprintf('figures/dist_mag/%s_mag_hist',metric(3:end)),'-depsc')



% 
% mkdir: cannot create directory ‘figures/dist_mag’: File exists
% [0] Pearson r: -0.351, p: 5.263e-59
% [1] Pearson r: -0.064, p: 2.817e-02
% [2] Pearson r: -0.280, p: 3.028e-18
% [3] Pearson r: -0.203, p: 1.405e-06
% [4] Pearson r: -0.029, p: 6.877e-01
% [5] Pearson r: 0.215, p: 1.122e-03
% [6] Pearson r: -0.267, p: 1.670e-02
% [7] Pearson r: -0.131, p: 4.080e-04
% [8] Pearson r: 0.458, p: 1.777e-03
% [9] Pearson r: -0.221, p: 1.049e-02
% [10] Pearson r: 0.054, p: 5.215e-01
% [11] Pearson r: 0.212, p: 8.216e-05
% [12] Pearson r: 0.169, p: 2.122e-01
% [13] Pearson r: 0.145, p: 9.645e-02
% [14] Pearson r: 0.116, p: 5.820e-01
% [15] Pearson r: 0.317, p: 7.400e-04
% [16] Pearson r: -0.172, p: 2.680e-14
% [17] Pearson r: -0.045, p: 5.513e-01
% [18] Pearson r: NaN, p: NaN
% [19] Pearson r: -0.118, p: 2.830e-01
% [21] Pearson r: 0.024, p: 8.968e-01
% [22] Pearson r: 0.788, p: 2.019e-02
% [23] Pearson r: 0.003, p: 9.785e-01
% [24] Pearson r: -0.350, p: 4.033e-04
% [25] Pearson r: 0.054, p: 4.488e-01
% [27] Pearson r: -0.089, p: 7.734e-01
% [28] Pearson r: 0.294, p: 5.219e-01
% [29] Pearson r: 0.274, p: 7.202e-02
% [30] Pearson r: -1.000, p: NaN
% [31] Pearson r: 0.293, p: 2.541e-01
% [32] Pearson r: -0.167, p: 2.458e-01
% [33] Pearson r: -0.300, p: 2.223e-02
% [34] Pearson r: 0.133, p: 1.753e-04
% [35] Pearson r: 0.052, p: 8.135e-01
% [37] Pearson r: 0.029, p: 7.844e-01
% [38] Pearson r: -0.485, p: 6.690e-02
% [39] Pearson r: 0.479, p: 9.745e-02
% [40] Pearson r: -0.201, p: 1.931e-15
% [41] Pearson r: 0.181, p: 5.710e-06
% [42] Pearson r: -0.144, p: 5.341e-03
% [43] Pearson r: 0.179, p: 7.339e-01
% [44] Pearson r: 0.200, p: 7.040e-01
% [45] Pearson r: -0.365, p: 1.392e-03
% [46] Pearson r: -0.516, p: 2.358e-01
% [47] Pearson r: -0.144, p: 1.634e-01
% [48] Pearson r: -0.144, p: 6.457e-03
% 	[dist-mag] Pearson r: 0.186, p: 9.662e-92
% 	[dist-mag] Spearman r: 0.173, p: 2.321e-79
% 	[dist-mag] Partial Spearman r | subject number, roi: 0.166, p: 1.262e-72
% 	[mSu: dist-mag] Pearson r: -0.351, p: 5.263e-59
% 	[mSu: dist-mag] Spearman r: -0.506, p: 2.473e-130
% [!] Coherence threshold: 0.30
% 	Number of subjects above threshold: 31 of 48
% 	Number of rois above threshold: 31 of 32
% 	[CT-mag] Pearson r: 0.078, p: 3.930e-17
% 	[CT-mag] Spearman r: 0.024, p: 9.273e-03
% 	[CT-mag] Partial Spearman r | subject number, roi: 0.025, p: 5.890e-03
% 	[mSu: CT-mag] Pearson r: 0.769, p: 000
% 	[mSu: CT-mag] Spearman r: 0.780, p: 000
% [0] Pearson r: -0.258, p: 1.864e-27
% [1] Pearson r: -0.069, p: 1.872e-02
% [2] Pearson r: -0.400, p: 1.920e-30
% [3] Pearson r: -0.189, p: 1.448e-03
% [4] Pearson r: -0.102, p: 1.823e-01
% [5] Pearson r: 0.207, p: 4.437e-03
% [6] Pearson r: -0.109, p: 5.821e-01
% [7] Pearson r: -0.133, p: 3.531e-02
% [8] Pearson r: 0.264, p: 1.377e-01
% [9] Pearson r: -0.151, p: 7.894e-02
% [10] Pearson r: 0.030, p: 7.634e-01
% [11] Pearson r: 0.267, p: 2.756e-07
% [12] Pearson r: 0.633, p: 1.508e-02
% [13] Pearson r: -0.049, p: 6.094e-01
% [14] Pearson r: -0.330, p: 3.223e-01
% [15] Pearson r: 0.318, p: 1.592e-02
% [16] Pearson r: -0.061, p: 5.721e-02
% [17] Pearson r: -0.176, p: 1.460e-02
% [18] Pearson r: 0.356, p: 6.441e-01
% [19] Pearson r: -0.057, p: 3.502e-01
% [21] Pearson r: 0.196, p: 3.588e-02
% [22] Pearson r: 0.695, p: 8.380e-03
% [23] Pearson r: 0.047, p: 5.479e-01
% [24] Pearson r: -0.290, p: 3.932e-04
% [25] Pearson r: 0.110, p: 5.439e-02
% [27] Pearson r: -0.101, p: 6.534e-01
% [28] Pearson r: 0.330, p: 1.959e-01
% [29] Pearson r: 0.123, p: 4.386e-01
% [30] Pearson r: -0.806, p: 4.034e-01
% [31] Pearson r: 0.021, p: 9.372e-01
% [32] Pearson r: -0.268, p: 6.822e-02
% [33] Pearson r: -0.060, p: 6.260e-01
% [34] Pearson r: 0.061, p: 1.114e-01
% [35] Pearson r: 0.113, p: 6.554e-01
% [37] Pearson r: -0.077, p: 3.966e-01
% [38] Pearson r: -0.239, p: 2.960e-01
% [39] Pearson r: 0.424, p: 1.938e-01
% [40] Pearson r: 0.138, p: 5.888e-07
% [41] Pearson r: 0.187, p: 2.562e-06
% [42] Pearson r: -0.101, p: 2.868e-02
% [43] Pearson r: -0.748, p: 1.279e-02
% [44] Pearson r: 0.500, p: 9.773e-02
% [45] Pearson r: -0.111, p: 3.184e-01
% [46] Pearson r: -0.070, p: 6.753e-01
% [47] Pearson r: -0.194, p: 1.372e-02
% [48] Pearson r: 1.000, p: 1.740e-02
% 	[dist-mag] Pearson r: 0.203, p: 4.095e-90
% 	[dist-mag] Spearman r: 0.258, p: 1.867e-145
% 	[dist-mag] Partial Spearman r | subject number, roi: 0.213, p: 1.121e-98
% 	[mSu: dist-mag] Pearson r: -0.258, p: 1.864e-27
% 	[mSu: dist-mag] Spearman r: -0.341, p: 6.494e-48
% [!] Coherence threshold: 0.30
% 	Number of subjects above threshold: 30 of 48
% 	Number of rois above threshold: 29 of 32
% 	[CT-mag] Pearson r: 0.257, p: 3.124e-144
% 	[CT-mag] Spearman r: 0.176, p: 1.010e-67
% 	[CT-mag] Partial Spearman r | subject number, roi: 0.149, p: 1.690e-48
% 	[mSu: CT-mag] Pearson r: 0.803, p: 000
% 	[mSu: CT-mag] Spearman r: 0.818, p: 000