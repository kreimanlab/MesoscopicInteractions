close all;
clear;
rng shuffle;

metricp = 'pcBroadband';
n_perm = 10000;
perm_alpha_sec = 12;

% Fast i/o definitions
dir_artLp = '/media/jerry/KLAB101/h5_notch20/art_nosz';
dir_resLp = '/media/jerry/KLAB101/results/coh_w10';
dir_corLp = '/media/jerry/internal/data/coreg';
dir_cacheLp = './cache';
%subjects_dirL = '/mnt/cuenap_ssd/coregistration';
subjects_dirL = '/media/jerry/internal/data/coreg';
setenv('SUBJECTS_DIR',subjects_dirL);

% Slow i/o definitions
dir_h5Lp = '/media/jerry/KLAB101/h5_notch20';
%dir_h5Lp = '/mnt/cuenap/data/h5_notch20';

metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'}; % ,'pcDelta'

% Patients
Subjectsp = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124',...
    'mSu'};

% Exclude monkey
Subjectsp = Subjectsp(1:(end-1));

% Triggers
system('mkdir figures');
system('mkdir figures/T7d2');

% only works for macaque, not human
trig_count_second_bip = false;

% Sorted Macaque mSu coverage:
% sort(D.AtlasLabels{1})

% Plot
% Roi1 = {{'V4'};            {'V1'}; {'V1'}; {'V1'}; {'V1'}... % {'V4'}
%         };
% Roi2 = {{'TEa_ma','TEav'}; {'V4'}; {'TEO'}; {'5'}; {'2'}... % {'8l','8m'} 'TEa_ma','TEa_mp'
%         };

% pos
% v4 - teo: 0.1614
% v1 - teo: 0.1688
% v1 - (v2,v4): 0.2394
% neg
% v4 - 9: 0.1713
% v4 - 7B: bad
% v4 - 5: bad
% v4 - 2: high
% v4 - 3: high
% v4 - ProM: hi
% v4 - PBc: hi
% v4 - F7: no human
% v4 - F2: no human
% v4 - F1: hi
% v4 - 1: hi
% v4 - v6a: no cov
% v4 - mip: hi
% v4 - 7op: hi
% v1 - 2: 0.1781
% v1 - 5: 0.2525
% v1 - MIP: bad
% v1 - V6A: bad
% v1 - 7B: 0.163
% v1 - 1: incomplete
% v1 - f1, f2: bad
% v1 - 3, ProM: high

%'TEa_ma','TE_mp','TEav','TEpv',
Roi1 = { {'V1'};  {'V1'};  {'V4'};          {'V1'};  {'V1'}; {'V4'};... % {'V4'}
        };
Roi2 = { {'V4'};  {'TEO'}; {'TEad','TEpd'}; {'2'};   {'5'};  {'F5'};... % {'8l','8m'} 'TEa_ma','TEa_mp'
        };
% Roi2 test
% V4
%   7B - 6/15 mon, 2/90 hum
%   F2 - error: no coverage in humans
%   F1 - 0/25 mon , 0/4 hum
%   F5 - 0/20 mon , 0/6 hum
%   F7 - error
%   10 - 0/5 mon, 15/316 hum
%   46d - 0/10 mon, 0/14 hum
%   24c - error
%   2 - 9/40 mon, 2/89 hum
%   5 - 1/15 mon, 0/40 hum
%   ProM - 1/10 mon, 0/42 hum
% V1
%   7B - 2/27 mon, 0/17 hum
%   F2 - error
%   F1 - 0/45 mon, 0/4 hum
%   F5 - 0/36 mon, 0/2 hum
%   F7 - error
%   10 - 0/9 mon, 5/174 hum
%   46d - error
%   24c - error
%   2 - 3/72 mon, 0/27 hum *
%   5 - 1/27 mon, 0/34 hum
%   ProM - 0/18 mon, 0/9 hum
    
% channels for first ROI2
chans2_1 = [112,113,114,125,128]; % [125, 128, 127,112,113]
% Annot = {'Zhou & Desimone 2016';...
%         'Bastos & Fries 2015';...
%         'Bastos & Fries 2015';...
%         'Markov & Kennedy 2014';...
%         'Markov & Kennedy 2014'
%         };
Annot = {'Bastos-Fries';...
        'Bastos-Fries';...
        'Zhou-Desimone';...
        'Markov-Kennedy';...
        'Markov-Kennedy';...
        'Markov-Kennedy'
        }; % 'Gregoriou-Desimone';...
is_pos = [true, true, true, false, false, false];

trig_avg_partial = false;
trig_subtract_coh_thresh = false;
trig_divide_avg_coh = false;

n_eg = length(Roi1);


% Load mSu labels
D = load('../coreg/all_surf_ielvis_macaque.mat'); % electrode labels
D_anat = load('../coreg/AdjMKFV.mat'); % Markov Kennedy Adj matrix
D_SuMAP = load('SuMAP.mat'); % brain map

% Display to check ROI selection:
roi2_1 = join(Roi2{1});
fprintf('Composite ROI: %s\n',roi2_1{1});
for i = 1:length(chans2_1)
    fprintf('Electrode %i\t%s\n',chans2_1(i),D.AtlasLabels{1}{chans2_1(i)});
end

% For each positive or negative example
iMp = 1;
atlM = 7;
Coh = cell(2,n_eg);
Coh_avg = zeros(2,n_eg);
Coh_std = zeros(2,n_eg);
Fracsig = zeros(2,n_eg);
Npairs = zeros(2,n_eg);
Npairs_sig = zeros(2,n_eg);
Elecs_mSu = cell(2,n_eg);
for i = 1:n_eg
    if (length(Roi2{i}) <= 1)
        roi2txt = Roi2{i}{1};
        fprintf('[*] Example %s, %s & %s\n',Annot{i},Roi1{i}{1},roi2txt)
    else
        roi2txt = '';
        for i2 = 1:length(Roi2{i})
            if (i2 > 1)
                roi2txt = [roi2txt, ', ', Roi2{i}{i2}];
            else
                roi2txt = [roi2txt, Roi2{i}{i2}];
            end
        end
        fprintf('[*] Example %s, %s & %s\n',Annot{i},Roi1{i}{1},roi2txt)
    end
    n_pairs = 0;
    coh_sig = [];
    coh_nosig = [];
    subs_eg = {};
    dists_mm = [];
    
    
    % Get average coherence for humans
    for j = 1:length(Subjectsp)
        sidp = Subjectsp{j};
        fn_artLp = sprintf('%s/%s_art.h5',dir_artLp,sidp);
        fn_distLp = sprintf('%s/%s_dists-%s-%i.mat',dir_resLp,sidp,metricp,n_perm);
        fn_graphLp = sprintf('%s/%s_graph-%s.h5',dir_resLp,sidp,metricp);
        fn_permLp = sprintf('/mnt/cuenap2/data/results/coh_w10/%s_perm-%s-%i.h5',sidp,metricp,n_perm);
        fn_h5Lp = sprintf('%s/%s.h5',dir_h5Lp,sidp);
        fn_coregLp = sprintf('%s/%s/label/all_parcellation.mat',dir_corLp,sidp);
        fn_cacheLp = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheLp,sidp,iMp);
        
        Ca = load(fn_cacheLp);
        
        
        
        ecog = Ca.ecog;
        %Ca.Dmats = Ca.Dmats / ((1/10) * (5/1.1));
        %ecog = H5eeg(fn_h5Lp);
        
        % Check for pairs that match description
        roi1 = Roi1{i};
        roi2 = Roi2{i};
        count = 1;
        for ii = 1:(ecog.n_bchan-1)
            for jj = (ii+1):ecog.n_bchan
                b1c1 = ecog.bip(ii,1);
                b1c2 = ecog.bip(ii,2);
                b2c1 = ecog.bip(jj,1);
                b2c2 = ecog.bip(jj,2);
                
                dist_count = Ca.Dmats(count);
                sat_dist = dist_count > Ca.dist_thresh;
                m132L = Ca.C.AtlLabels{atlM};
                roi1a = [roi1{1},'_M132'];
                roi2a = {};
                for i3 = 1:length(Roi2{i})
                    roi2a{i3} = [roi2{i3},'_M132'];
                end
%                 if (i == 1)
%                     % Manually pick electrodes for first positive example
%                     sat_fwd = any(strcmp(m132L{b1c1},roi1a)) && any(b2c1 == chans2_1); % && any(strcmp(m132L{b2c1},roi2a));
%                     sat_bak = any(strcmp(m132L{b2c1},roi1a)) && any(b1c1 == chans2_1); % && any(strcmp(m132L{b1c1},roi2a));
%                 else
                sat_fwd = any(strcmp(m132L{b1c1},roi1a)) && any(strcmp(m132L{b2c1},roi2a));
                sat_bak = any(strcmp(m132L{b2c1},roi1a)) && any(strcmp(m132L{b1c1},roi2a));
%                 end
                sat_notempty = (~isempty(m132L{b1c1})) && (~isempty(m132L{b2c1}));
%                 fname = sprintf('%s_%i_%s_%i_%s',sidp,ii,m132L{b1c1},jj,m132L{b2c1});
%                 fprintf('%s\n',fname);
                if (sat_dist && (sat_fwd || sat_bak) && sat_notempty)
                    
                    % Significant electrode pair in right areas
                    fname = sprintf('%s_%i_%s_%i_%s',sidp,ii,m132L{b1c1},jj,m132L{b2c1});
                    %fprintf('%s\n',fname);
                    n_pairs = n_pairs + 1;
                    subs_eg = {subs_eg{:}, sidp};
                    
                    % Get coherence, if pair is significant
                    sat_sig =  Ca.ct(count) > Ca.ct_thresh;
                    coh_mag = Ca.mag(count);
                    if (trig_subtract_coh_thresh)
                        coh_mag = coh_mag - Ca.coh_thresh(count);
                    end
                    if (sat_sig)
                        coh_sig = [coh_sig, coh_mag];
                    else
                        coh_sig = [coh_sig, NaN];
                        coh_nosig = [coh_nosig, coh_mag];
                    end
                    
                    dists_mm = [dists_mm, dist_count];
                    
                    %return
                end
                count = count + 1;
            end
        end
    end
    
    
    % Get average coherence for mSu
    sidp = 'mSu';
    fn_h5Lp = sprintf('%s/%s.h5',dir_h5Lp,sidp);
    fn_cacheLp = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheLp,sidp,iMp);
    fn_artLp = sprintf('%s/%s_art.h5',dir_artLp,sidp);
    fn_graphLp = sprintf('%s/%s_graph-%s.h5',dir_resLp,sidp,metricp);
    %ecog = H5eeg(fn_h5Lp);
    Ca = load(fn_cacheLp);
    ecog = Ca.ecog;
    % DEBUG adjust ct threshold
    Ca.ct_thresh = 0.05;
    
    if (i == 1)
        %return
    end
    
    n_pairs_mSu = 0;
    coh_sig_mSu = [];
    coh_nosig_mSu = [];
    dists_mm_mSu = [];
    % Copied from human section above
    % Check for pairs that match description
    roi1 = Roi1{i};
    roi2 = Roi2{i};
    count = 1;
    rois_idx_mSu = [];
    rois_idx_mSu2 = [];
    for ii = 1:(ecog.n_bchan-1)
        for jj = (ii+1):ecog.n_bchan
            b1c1 = ecog.bip(ii,1);
            b1c2 = ecog.bip(ii,2);
            b2c1 = ecog.bip(jj,1);
            b2c2 = ecog.bip(jj,2);

            dist_count = Ca.Dmats(count);
            sat_dist = dist_count > Ca.dist_thresh;
            %m132L = Ca.C.AtlLabels{atlM};
            m132L = D.AtlasLabels{1}; % replace with macaque labels
            
            manual_prune = find(strcmp(m132L,'V4_M132'));
            m132L{manual_prune(6)} = 'V3_M132';
            m132L{manual_prune(7)} = 'V3_M132';
            manual_prune = find(strcmp(m132L,'TEO_M132'));
            m132L{manual_prune(end)} = 'V4_M132';
            m132L{28} = '8l_M132';
            m132L{28-8} = '8l_M132';
            
            roi1a = [roi1{1},'_M132'];
            %roi2a = [roi2{1},'_M132'];
            roi2a = {};
            for i3 = 1:length(Roi2{i})
                roi2a{i3} = [roi2{i3},'_M132'];
            end
            if (i == 3)
                % Manually pick electrodes for first positive example
                sat_fwd = any(strcmp(m132L{b1c1},roi1a)) && (any(b2c1 == chans2_1)); % && any(strcmp(m132L{b2c1},roi2a));
                if (trig_count_second_bip)
                    sat_fwd = sat_fwd || (any(strcmp(m132L{b1c2},roi1a)) && (any(b2c1 == chans2_1)));
                    sat_fwd = sat_fwd || (any(strcmp(m132L{b1c1},roi1a)) && (any(b2c2 == chans2_1)));
                    sat_fwd = sat_fwd || (any(strcmp(m132L{b1c2},roi1a)) && (any(b2c2 == chans2_1)));
                end
                
                sat_bak = any(strcmp(m132L{b2c1},roi1a)) && (any(b1c1 == chans2_1)); % && any(strcmp(m132L{b1c1},roi2a));
                if (trig_count_second_bip)
                    sat_bak = sat_bak || (any(strcmp(m132L{b2c2},roi1a)) && (any(b1c1 == chans2_1))); 
                    sat_bak = sat_bak || (any(strcmp(m132L{b2c1},roi1a)) && (any(b1c2 == chans2_1))); 
                    sat_bak = sat_bak || (any(strcmp(m132L{b2c2},roi1a)) && (any(b1c2 == chans2_1))); 
                end
            else
                sat_fwd = any(strcmp(m132L{b1c1},roi1a)) && any(strcmp(m132L{b2c1},roi2a));
                
                if (trig_count_second_bip)
                    % consider second bipolar electrode
                    sat_fwd = sat_fwd || (any(strcmp(m132L{b1c2},roi1a)) && any(strcmp(m132L{b2c1},roi2a)));
                    sat_fwd = sat_fwd || (any(strcmp(m132L{b1c1},roi1a)) && any(strcmp(m132L{b2c2},roi2a)));
                    sat_fwd = sat_fwd || (any(strcmp(m132L{b1c2},roi1a)) && any(strcmp(m132L{b2c2},roi2a)));
                end
                
                sat_bak = any(strcmp(m132L{b2c1},roi1a)) && any(strcmp(m132L{b1c1},roi2a));
                
                if (trig_count_second_bip)
                    % consider second bipolar electrode
                    sat_bak = sat_bak || (any(strcmp(m132L{b2c2},roi1a)) && any(strcmp(m132L{b1c1},roi2a)));
                    sat_bak = sat_bak || (any(strcmp(m132L{b2c1},roi1a)) && any(strcmp(m132L{b1c2},roi2a)));
                    sat_bak = sat_bak || (any(strcmp(m132L{b2c2},roi1a)) && any(strcmp(m132L{b1c2},roi2a)));
                end
            end
                
            %sat_fwd = any(strcmp(m132L{b1c1},roi1a)) && any(strcmp(m132L{b2c1},roi2a));
            %sat_bak = any(strcmp(m132L{b2c1},roi1a)) && any(strcmp(m132L{b1c1},roi2a));
            sat_notempty = (~isempty(m132L{b1c1})) && (~isempty(m132L{b2c1}));
%                 fname = sprintf('%s_%i_%s_%i_%s',sidp,ii,m132L{b1c1},jj,m132L{b2c1});
%                 fprintf('%s\n',fname);
            

            if (sat_dist && (sat_fwd || sat_bak) && sat_notempty)
                % Significant electrode pair in right areas
                fname = sprintf('%s_%i_%s_%i_%s',sidp,ii,m132L{b1c1},jj,m132L{b2c1});
                %fprintf('%s\n',fname);
                n_pairs_mSu = n_pairs_mSu + 1;
                %subs_eg = {subs_eg{:}, sidp};
                
                coh_mag = Ca.mag(count);
                if (trig_subtract_coh_thresh)
                    coh_mag = coh_mag - Ca.coh_thresh(count);
                end
                    
                % Get coherence, if pair is significant
                sat_sig =  Ca.ct(count) > Ca.ct_thresh;
                if (sat_sig)
                    coh_sig_mSu = [coh_sig_mSu, coh_mag];
                else
                    coh_sig_mSu = [coh_sig_mSu, NaN];
                    coh_nosig_mSu = [coh_nosig_mSu, coh_mag];
                end
                
                dists_mm_mSu = [dists_mm_mSu, dist_count];
                
                if (sat_fwd)
                    c1_sav = b1c1;
                    c2_sav = b2c1;
                else
                    c1_sav = b2c1;
                    c2_sav = b1c1;
                end
                rois_idx_mSu = [rois_idx_mSu, c1_sav];
                rois_idx_mSu2 = [rois_idx_mSu2, c2_sav];
            end
            count = count + 1;
        end
    end
    Elecs_mSu{1,i} = sort(unique(rois_idx_mSu));
    Elecs_mSu{2,i} = sort(unique(rois_idx_mSu2));
    
    n_pairs_sig = sum(~isnan(coh_sig));
    avg_mag = nanmean(coh_sig);
    std_mag = nanstd(coh_sig);
    avg_mag_nosig = nanmean(coh_nosig);
    std_mag_nosig = nanstd(coh_nosig);
    
    n_pairs_sig_mSu = sum(~isnan(coh_sig_mSu));
    avg_mag_mSu = nanmean(coh_sig_mSu);
    std_mag_mSu = nanstd(coh_sig_mSu);
    avg_mag_nosig_mSu = nanmean(coh_nosig_mSu);
    std_mag_nosig_mSu = nanstd(coh_nosig_mSu);
    
    
    fprintf('\tn_pairs: %i\n',n_pairs);
    fprintf('\tn_unique_subs: %i\n',length(unique(subs_eg)));
    fprintf('\tavg mag: %.6f, std: %.6f\n',avg_mag,std_mag);
    fprintf('\tnosig avg mag: %.6f, std: %.6f\n',avg_mag_nosig, std_mag_nosig);
    fprintf('\tfrac sig: %.6f (%i of %i)\n',n_pairs_sig/n_pairs, n_pairs_sig, n_pairs);
    fprintf('\tdistance: %.2f mm, std: %.2f mm\n',mean(dists_mm),std(dists_mm))
    
    fprintf('\tn_pairs_mSu: %i\n',n_pairs_mSu);
    fprintf('\tavg mag mSu: %.6f, std: %.6f\n',avg_mag_mSu,std_mag_mSu);
    fprintf('\tnosig avg mag: %.6f, std: %.6f\n',avg_mag_nosig_mSu, std_mag_nosig_mSu);
    fprintf('\tfrac sig mSu: %.6f (%i of %i)\n',n_pairs_sig_mSu/n_pairs_mSu, n_pairs_sig_mSu, n_pairs_mSu);
    fprintf('\tdistance: %.2f mm?, std: %.2f mm?\n',mean(dists_mm_mSu),std(dists_mm_mSu))
    
%     if (is_pos(i))
%         % On positive examples
%         cond_human = ~isnan(avg_mag);
%         cond_macaque = ~isnan(avg_mag_mSu);
%     else
%         % On negative examples
%         cond_human = ~isnan(avg_mag);
%         cond_macaque = ~isnan(avg_mag_mSu);
%     end
    coh_all = [coh_nosig, coh_sig(~isnan(coh_sig))];
    coh_all_mSu = [coh_nosig_mSu, coh_sig_mSu(~isnan(coh_sig_mSu))];
    % Check coherence merge successful
    if (length(coh_all) == n_pairs)
        fprintf('length(coh_all) == n_pairs\n')
    end
    if (length(coh_all_mSu) == n_pairs_mSu)
        fprintf('length(coh_all_mSu) == n_pairs_mSu\n')
    end
    Coh{2,i} = coh_all;
    Coh{1,i} = coh_all_mSu;
    
    
    
    if (trig_avg_partial)
        if (is_pos(i))   
            Coh_avg(2,i) = avg_mag;
            Coh_std(2,i) = std_mag;
            Coh_avg(1,i) = avg_mag_mSu;
            Coh_std(1,i) = std_mag_mSu;
        else
            Coh_avg(2,i) = avg_mag_nosig;
            Coh_std(2,i) = std_mag_nosig;
            Coh_avg(1,i) = avg_mag_nosig_mSu;
            Coh_std(1,i) = std_mag_nosig_mSu;
        end
    else
        % Average regardless if significant or not
        
        
        Coh_avg(2,i) = nanmean(coh_all);
        Coh_std(2,i) = nanstd(coh_all);
        Coh_avg(1,i) = nanmean(coh_all_mSu);
        Coh_std(1,i) = nanstd(coh_all_mSu);
        
    end
    
    
    
    % Fraction significant
    Fracsig(2,i) = n_pairs_sig/n_pairs;
    Fracsig(1,i) = n_pairs_sig_mSu/n_pairs_mSu;
    Npairs(2,i) = n_pairs;
    Npairs_sig(2,i) = n_pairs_sig;
    Npairs(1,i) = n_pairs_mSu;
    Npairs_sig(1,i) = n_pairs_sig_mSu;
end

%%
col_bar = 0.85*[1 1 1];
col_bar_mSu = 0.5*[1 1 1];
col_elecs = {0.1*[1 1 1],0.8*[1 0 0]};
col_elec_border = [0 0 0];
size_elec = 12;
factor_elec_border = 0.7;
size_cap = 12; % error bar cap size
bar_center_adjust = 5.5;

PLOT_W_PIX = 1080;
PLOT_H_PIX = 1080*(3/4);
%h = figure('Position',[0 0 PLOT_W_PIX PLOT_H_PIX]);
h = figure('Position',[0 0 PLOT_W_PIX PLOT_H_PIX],'visible','on');

% Crop image
Im = D_SuMAP.Su.I;
Im = Im * 1.4;
%crop_pix = 0; % doesn't work yet
%crop_pix_h = 0;
%Im = Im(crop_pix_h:(end-crop_pix_h),crop_pix:(end-crop_pix),:);
for i = 1:n_eg
    subplot(3,n_eg,i)
    xticks([]);
    yticks([]);
    
    %Im(xs,ys,:) = [1 0 0];
    imagesc(Im); hold on;
    for ii = 1:2
        elecs = Elecs_mSu{ii,i};
        xs = D_SuMAP.Su.X(elecs);
        ys = D_SuMAP.Su.Y(elecs);
        %plot(xs,ys,'black.'); hold on;
        
        mstyle = '.';
        for j = 1:length(xs)
            %plot(xs(j),ys(j),'.','Color',col_elec_border,'MarkerSize',size_elec); hold on;
            plot(xs(j),ys(j),mstyle,'Color',col_elecs{ii},'MarkerSize',factor_elec_border*size_elec); hold on;
        end
    end
    
    daspect([1 1 1]);
    set(gca,'visible','off');
end


subplot(3,1,2)
x = 1:n_eg;
% bar(x,Coh_avg(1,:),'grouped','FaceColor',col_bar); hold on;
% er = errorbar(x,Coh_avg(1,:),Coh_std(1,:),Coh_std(1,:)); hold on;
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% bar(x,Coh_avg(2,:),'grouped','FaceColor',col_bar_mSu); hold on;
% er2 = errorbar(x,Coh_avg(2,:),Coh_std(2,:),Coh_std(2,:)); hold on;
% er2.Color = [0 0 0];
% er2.LineStyle = 'none';

% Raw data
% ba = bar(x,Coh_avg','grouped'); hold on;
% ba(2).FaceColor = col_bar;
% ba(1).FaceColor = col_bar_mSu;
% er = errorbar(x - ba(1).BarWidth/bar_center_adjust,Coh_avg(1,:),Coh_std(1,:),Coh_std(1,:),...
%     'CapSize',size_cap); hold on;
% er2 = errorbar(x + ba(2).BarWidth/bar_center_adjust,Coh_avg(2,:),Coh_std(2,:),Coh_std(2,:),...
%     'CapSize',size_cap); hold on;


% percentages
min1 = min(Coh_avg(1,:));
min2 = min(Coh_avg(2,:));

if (trig_divide_avg_coh)
    Coh_avg(1,:) = Coh_avg(1,:) / min1;
    Coh_avg(2,:) = Coh_avg(2,:) / min2;
    Coh_std(1,:) = Coh_std(1,:) / min1;
    Coh_std(2,:) = Coh_std(2,:) / min2;
end

ba = bar(x,Coh_avg','grouped'); hold on;
ba(2).FaceColor = col_bar;
ba(1).FaceColor = col_bar_mSu;
er = errorbar(x - ba(1).BarWidth/bar_center_adjust,Coh_avg(1,:),Coh_std(1,:),Coh_std(1,:),...
    'CapSize',size_cap); hold on;
er2 = errorbar(x + ba(2).BarWidth/bar_center_adjust,Coh_avg(2,:),Coh_std(2,:),Coh_std(2,:),...
    'CapSize',size_cap); hold on;

er.LineStyle = 'none';
er.Color = [0 0 0];
er2.LineStyle = 'none';
er2.Color = [0 0 0];
bar_xoffset = ba(1).BarWidth/bar_center_adjust;

% statistics
ax = gca;
if (trig_divide_avg_coh)
    ax.YLim(1) = 0;
end

if (trig_divide_avg_coh)
    %ax.YLim(2) = 4;
    sig_count = 0.05/min1;
    sig_yoffset = 0.02/min1;
    sigbar_ycap = 0.01/min1;
else
    sig_count = 0.05;
    sig_yoffset = 0.05;
    sigbar_ycap = 0.01;
end
fontsize_pval = 6;
fontsize_annot = 8;
fontsize_nsig = 8;
ymin = ax.YLim(1);
ymax = ax.YLim(2);
for i = 1:2
    Coh_effect = [];
    fprintf('%i:\n',i)
    cohs = Coh(i,:);
    for j = 1:(length(cohs)-1)
        for k = (j+1):length(cohs)
            x_stats = cohs{j};
            y_stats = cohs{k};
            sat_fwd = is_pos(j) && (~ is_pos(k));
            sat_bak = is_pos(k) && (~ is_pos(j));
            if (sat_fwd || sat_bak)
                [p,hy,stats] = ranksum(x_stats,y_stats);
                if (sat_fwd)
                    coh_avg_effect = Coh_avg(i,j) - Coh_avg(i,k);
                else
                    coh_avg_effect = Coh_avg(i,k) - Coh_avg(i,j);
                end
                
                % Bonferroni correct
                n_comparisons = sum(is_pos) * sum(~is_pos);
                p = p * (n_comparisons);
                hy = p < 0.05;
                
                fprintf('[%.4d] %i %i\tpass: %i\tranksum: %.6f\teffect: %.6f\n',p,j,k,hy,stats.ranksum,coh_avg_effect); % stats.zval
                if (hy && (coh_avg_effect > 0))
                    xoffset = 2*(i - 1.5)*bar_xoffset;
                    sigbar_x = [j k] + xoffset;
                    sigbar_y = [ymax ymax]+sig_count;
                    line(sigbar_x,sigbar_y,'Color','black'); hold on;
                    line(sigbar_x(1)*ones(1,2),[sigbar_y(1),sigbar_y(1)-sigbar_ycap],'Color','black'); hold on;
                    line(sigbar_x(2)*ones(1,2),[sigbar_y(1),sigbar_y(1)-sigbar_ycap],'Color','black'); hold on;
                    text(mean(sigbar_x),mean(sigbar_y),sprintf('p=%.1d',p),...
                        'HorizontalAlignment','center','VerticalAlignment','bottom',...
                        'FontSize',fontsize_pval);
                    sig_count = sig_count + sig_yoffset;
                    Coh_effect = [Coh_effect coh_avg_effect];
                end
            end
            %return
        end
    end
    
    if (i == 1)
        minc = min1;
    elseif (i == 2)
        minc = min2;
    end
    fprintf('[*] Average effect size: %.6f (%.2f %%)\n',mean(Coh_effect),100*mean(Coh_effect)/minc);
    fprintf('[*] Max effect size: %.6f (%.2f %%)\n',max(Coh_effect),100*max(Coh_effect)/minc);
    fprintf('[*] Min effect size: %.6f (%.2f %%)\n',min(Coh_effect),100*min(Coh_effect)/minc);
end

% positive/negative example line
col_posneg_line = 0.5*[1 1 1];
x_line = find(~is_pos,1) - 0.5;
plot([x_line x_line],[ymin ymax],'--','Color',col_posneg_line); hold on;
if (trig_divide_avg_coh)
    ylabel('Coherence (fraction)')
else
    ylabel('Coherence')
end
box off;
set(gca,'TickDir','out');

%xannot = cell(size(Annot));
if (trig_divide_avg_coh)
    ax.YLim(2) = ax.YLim(2) + 1.5;
    line2_yoffset = -0.1/0.15;
    line_yoffset = -0.05/0.15;
else
    line2_yoffset = -0.1;
    line_yoffset = -0.05;
end
for i = 1:length(Annot)
    %xannot{i} = [Annot{i}, [Roi1{i}{1},' & ',roi2txt]];
    if (length(Roi2{i}) <= 1)
        roi2txt = Roi2{i}{1};
    else
        roi2txt = '';
        for i2 = 1:length(Roi2{i})
            if (i2 > 1)
                roi2txt = [roi2txt, ', ', Roi2{i}{i2}];
            else
                roi2txt = [roi2txt, Roi2{i}{i2}];
            end
        end
    end
    
    line2 = [Roi1{i}{1},' & ',replace(roi2txt,'_','/')];
    text(i,line2_yoffset,line2,'HorizontalAlignment','center','FontSize',fontsize_annot);
    text(i,line_yoffset,Annot{i},'HorizontalAlignment','center','FontSize',fontsize_annot);
end
%xticklabels(Annot);
xticks([]);
%legend({'Human','Macaque'})
legend({'Macaque','Human'},'Location','Best')

subplot(3,1,3)

ba = bar(x,Fracsig','grouped'); hold on;
ba(2).FaceColor = col_bar;
ba(1).FaceColor = col_bar_mSu;
% er = errorbar(x - ba(1).BarWidth/bar_center_adjust,Coh_avg(1,:),Coh_std(1,:),Coh_std(1,:),...
%     'CapSize',size_cap); hold on;
% er2 = errorbar(x + ba(2).BarWidth/bar_center_adjust,Coh_avg(2,:),Coh_std(2,:),Coh_std(2,:),...
%     'CapSize',size_cap); hold on;
% er.LineStyle = 'none';
% er.Color = [0 0 0];
% er2.LineStyle = 'none';
% er2.Color = [0 0 0];

% positive/negative example line
ax = gca;
ymin = ax.YLim(1);
ymax = ax.YLim(2);
col_posneg_line = 0.5*[1 1 1];
x_line = find(~is_pos,1) - 0.5;
plot([x_line x_line],[ymin ymax],'--','Color',col_posneg_line); hold on;


% ba(2).BarWidth/bar_center_adjust
for i = 1:n_eg
    text(i - bar_xoffset,Fracsig(1,i),...
        sprintf('%i/%i',Npairs_sig(1,i),Npairs(1,i)),...
        'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontSize',fontsize_nsig);
    text(i + bar_xoffset,Fracsig(2,i),...
        sprintf('%i/%i',Npairs_sig(2,i),Npairs(2,i)),...
        'HorizontalAlignment','center','VerticalAlignment','bottom',...
        'FontSize',fontsize_nsig);
end

ylabel('Fraction of interactions')
box off;
set(gca,'TickDir','out');
xticks([]);
%legend({'Human','Macaque'})
legend({'Macaque','Human'},'Location','Best')

% Save
pname = 'macaque_validation';
for i = 1:length(Roi1)
    roi1 = Roi1{i};
    roi2 = Roi2{i};
    for j = 1:length(roi1)
        pname = [pname,'_',roi1{j},'-',roi2{j}];
    end
end

if (trig_divide_avg_coh)
    print(h,sprintf('figures/T7d2/%s_fraction',pname),'-depsc');
    print(h,sprintf('figures/T7d2/%s_fraction',pname),'-dpng');
else
    print(h,sprintf('figures/T7d2/%s',pname),'-depsc');
    print(h,sprintf('figures/T7d2/%s',pname),'-dpng');
end
%close(h);



% July 23, 2019
% 
% figure_t7d1
% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T7d1’: File exists
% [*] Example Zhou-Desimone, V4 & TEad, TEpd
% 	n_pairs: 66
% 	n_unique_subs: 8
% 	avg mag: 0.147301, std: 0.005760
% 	nosig avg mag: 0.328500, std: 0.224745
% 	frac sig: 0.045455 (3 of 66)
% 	distance: 37.16 mm, std: 9.07 mm
% 	n_pairs_mSu: 15
% 	avg mag mSu: 0.151323, std: 0.010874
% 	nosig avg mag: NaN, std: NaN
% 	frac sig mSu: 1.000000 (15 of 15)
% 	distance: 59.02 mm?, std: 22.96 mm?
% length(coh_all) == n_pairs
% length(coh_all_mSu) == n_pairs_mSu
% [*] Example Gregoriou-Desimone, V4 & 8l, 8m
% 	n_pairs: 56
% 	n_unique_subs: 2
% 	avg mag: NaN, std: NaN
% 	nosig avg mag: 0.432920, std: 0.152223
% 	frac sig: 0.000000 (0 of 56)
% 	distance: 81.79 mm, std: 9.62 mm
% 	n_pairs_mSu: 15
% 	avg mag mSu: NaN, std: NaN
% 	nosig avg mag: 0.136924, std: 0.003696
% 	frac sig mSu: 0.000000 (0 of 15)
% 	distance: 90.58 mm?, std: 14.61 mm?
% length(coh_all) == n_pairs
% length(coh_all_mSu) == n_pairs_mSu
% [*] Example Bastos-Fries, V1 & V4
% 	n_pairs: 248
% 	n_unique_subs: 6
% 	avg mag: 0.175235, std: 0.014740
% 	nosig avg mag: 0.300564, std: 0.186477
% 	frac sig: 0.431452 (107 of 248)
% 	distance: 30.90 mm, std: 7.46 mm
% 	n_pairs_mSu: 44
% 	avg mag mSu: 0.155349, std: 0.007004
% 	nosig avg mag: NaN, std: NaN
% 	frac sig mSu: 1.000000 (44 of 44)
% 	distance: 37.07 mm?, std: 8.74 mm?
% length(coh_all) == n_pairs
% length(coh_all_mSu) == n_pairs_mSu
% [*] Example Bastos-Fries, V1 & TEO
% 	n_pairs: 32
% 	n_unique_subs: 2
% 	avg mag: 0.189168, std: 0.009524
% 	nosig avg mag: 0.214511, std: 0.125758
% 	frac sig: 0.187500 (6 of 32)
% 	distance: 38.22 mm, std: 9.34 mm
% 	n_pairs_mSu: 54
% 	avg mag mSu: 0.151640, std: 0.006959
% 	nosig avg mag: 0.144061, std: 0.001081
% 	frac sig mSu: 0.944444 (51 of 54)
% 	distance: 49.19 mm?, std: 10.88 mm?
% length(coh_all) == n_pairs
% length(coh_all_mSu) == n_pairs_mSu
% [*] Example Markov-Kennedy, V4 & 9
% 	n_pairs: 2
% 	n_unique_subs: 1
% 	avg mag: NaN, std: NaN
% 	nosig avg mag: 0.171336, std: 0.009398
% 	frac sig: 0.000000 (0 of 2)
% 	distance: 99.02 mm, std: 3.39 mm
% 	n_pairs_mSu: 10
% 	avg mag mSu: NaN, std: NaN
% 	nosig avg mag: 0.137903, std: 0.002810
% 	frac sig mSu: 0.000000 (0 of 10)
% 	distance: 121.26 mm?, std: 16.56 mm?
% length(coh_all) == n_pairs
% length(coh_all_mSu) == n_pairs_mSu
% [*] Example Markov-Kennedy, V1 & 7B
% 	n_pairs: 17
% 	n_unique_subs: 2
% 	avg mag: NaN, std: NaN
% 	nosig avg mag: 0.163018, std: 0.020873
% 	frac sig: 0.000000 (0 of 17)
% 	distance: 69.97 mm, std: 6.82 mm
% 	n_pairs_mSu: 27
% 	avg mag mSu: 0.143386, std: 0.003852
% 	nosig avg mag: 0.147859, std: 0.007965
% 	frac sig mSu: 0.074074 (2 of 27)
% 	distance: 64.76 mm?, std: 9.08 mm?
% length(coh_all) == n_pairs
% length(coh_all_mSu) == n_pairs_mSu
% 1:
% [2.7718e-03] 1 5	pass: 1	ranksum: 260.000000	effect: 0.013420
% [2.7573e+00] 1 6	pass: 0	ranksum: 359.000000	effect: 0.003796
% [2.8805e+00] 2 5	pass: 0	ranksum: 178.000000	effect: -0.000979
% [1.3359e-04] 2 6	pass: 1	ranksum: 158.000000	effect: -0.010604
% [8.1556e-06] 3 5	pass: 1	ranksum: 1430.000000	effect: 0.017446
% [1.1208e-03] 3 6	pass: 1	ranksum: 1906.000000	effect: 0.007822
% [1.5397e-05] 4 5	pass: 1	ranksum: 2013.000000	effect: 0.013316
% [3.0105e-01] 4 6	pass: 0	ranksum: 2422.000000	effect: 0.003692
% [*] Average effect size: 0.013001 (9.50 %)
% [*] Max effect size: 0.017446 (12.74 %)
% [*] Min effect size: 0.007822 (5.71 %)
% 2:
% [5.7746e+00] 1 5	pass: 0	ranksum: 1701.000000	effect: 0.147627
% [3.5869e+00] 1 6	pass: 0	ranksum: 2197.000000	effect: 0.155946
% [1.4665e-01] 2 5	pass: 0	ranksum: 1482.000000	effect: 0.261584
% [7.4709e-09] 2 6	pass: 1	ranksum: 2260.000000	effect: 0.269902
% [3.5047e+00] 3 5	pass: 0	ranksum: 28516.000000	effect: 0.072645
% [2.9652e-04] 3 6	pass: 1	ranksum: 31425.000000	effect: 0.080964
% [7.3006e+00] 4 5	pass: 0	ranksum: 562.000000	effect: 0.038423
% [1.5141e+00] 4 6	pass: 0	ranksum: 863.000000	effect: 0.046742
% [*] Average effect size: 0.175433 (107.62 %)
% [*] Max effect size: 0.269902 (165.57 %)
% [*] Min effect size: 0.080964 (49.67 %)