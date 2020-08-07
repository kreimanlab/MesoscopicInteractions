close all
clear
rng shuffle;


metric = 'pcBroadband';
n_perm = 10000;
perm_alpha_sec = 12;

% Fast i/o definitions
dir_artLp = '/media/klab/KLAB101/h5_notch20/art_nosz';
dir_resLp = '/media/klab/KLAB101/results/coh_w10';
dir_corLp = '/media/klab/internal/data/coreg';
dir_cacheLp = './cache';
subjects_dirLp = '/mnt/cuenap_ssd/coregistration';

% Slow i/o definitions
dir_h5Lp = '/media/klab/KLAB101/h5_notch20';

metricsp = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};

% Patients
Subjects = {'mSu'};
sidp = Subjects{1};
iMp = 1;

% Triggers
system('mkdir figures');
fig_fmt = '-depsc';
trig_t8 = true;
trig_title = false;
trig_save_fig = false;
x_div_secs = 2;
PLOT_W_PIX = 1080;
PLOT_H_PIX = 1080;

fn_artLp = sprintf('%s/%s_art.h5',dir_artLp,sidp);
fn_distLp = sprintf('%s/%s_dists-%s-%i.mat',dir_resLp,sidp,metric,n_perm);
fn_graphLp = sprintf('%s/%s_graph-%s.h5',dir_resLp,sidp,metric);
fn_permLp = sprintf('/mnt/cuenap2/data/results/coh_w10/%s_perm-%s-%i.h5',sidp,metric,n_perm);
fn_h5Lp = sprintf('%s/%s.h5',dir_h5Lp,sidp);
fn_coregLp = sprintf('%s/%s/label/all_parcellation.mat',dir_corLp,sidp);
fn_cacheLp = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheLp,sidp,iMp);

% Load mSu cache
load(fn_cacheLp);

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


% Pick areas to plot
% C.AtlROIs{2}.LH.struct_names
% (PFC, area 9/46) and the anterior cingulate cortex (ACC, areas 24c and 32).
atl = 7;
delta_thresh_4_plot = 0.4; % if no time sample is givem, plot with delta > this thresh, 


% Buschman Miller 2012   
%roi_1s = {'9_46d_M132','9_46v_M132'};
%roi_2s = {'24c_M132','32_M132'};
% NO COVERAGE


% Gregoriou Desimone 2009
%V4-8m negative no coverage
% prefix = 'GD';
% roi_1s = {'V4_M132'};
% %roi_2s = {'8l_M132'}; % [!] Fraction of covered pairs significantly consistent in time: 0.700
% roi_2s = {'8m_M132','8l_M132'}; % FEF subdivisions found in Markov Kennedy
% bchan1_const = NaN; %72; % V4 - 72, 64, 65
% bchan2_const = NaN; %24; % 8l - 24, 24, 24
% r_samp_const = 6152501;

% V4 - 8l, same as V4 - 8m and 8l, no coverage 8m
%V4_M132-->8l_M132: 1186 neurons
%8l_M132-->V4_M132: 186 neurons
%[!] Fraction of covered pairs significantly consistent in time: 0.700
%	Numbers: 14 of 20
%[!] Average Broadband coherence: 0.256

% Bastos Fries 2015
prefix = 'BF';
roi_1s = {'V1_M132'};
roi_2s = {'V4_M132'};
%roi_2s = {'TEO_M132'};
%roi_2s = {'TEO_M132','DP_M132'};
bchan1_const = NaN; %78; % V1 - 107, 78, 79
bchan2_const = NaN; %87; % TEO - 102, 87, 87
r_samp_const = 1250001;
% [!] Fraction of covered pairs significantly consistent in time: 1.000

%roi_2s = {'TEO_M132','DP_M132'};
%[!] Fraction of covered pairs significantly consistent in time: 1.000
%	Numbers: 117 of 117
%[!] Average Broadband coherence: 0.276

%roi_2s = {'TEO_M132'};
%[!] Fraction of covered pairs significantly consistent in time: 1.000
%	Numbers: 117 of 117
%[!] Average Broadband coherence: 0.276

%roi_2s = {'V4_M132'};
%[!] Fraction of covered pairs significantly consistent in time: 1.000
%	Numbers: 116 of 116
%[!] Average Broadband coherence: 0.281

% MK confirmed negative
% [!] Fraction of covered pairs significantly consistent in time: 0.292
% prefix = 'MK';
% roi_1s = {'V1_M132'};
% roi_2s = {'ProM_M132'};
% bchan1_const = NaN; %93;
% bchan2_const = NaN; %58;
% r_samp_const = 1427501;
%[!] Fraction of covered pairs significantly consistent in time: 0.292
%	Numbers: 7 of 24
%[!] Average Broadband coherence: 0.258


ct_thresh_perm_count = [];
ct_thresh_perm_count_mag = [];

for roi_idx1p = 1:length(roi_1s)
    for roi_idx2p = 1:length(roi_2s)

        roi_1 = roi_1s{roi_idx1p}; %'parsorbitalis';
        roi_2 = roi_2s{roi_idx2p}; %'inferiorparietal';
        
        % Find anatomical connections
        roi_1_idx = 0;
        roi_2_idx = 0;
        for i = 1:length(D_anat.MK_labels)
            if (contains(roi_1,D_anat.MK_labels{i}))
                roi_1_idx = i;
            end
            if (contains(roi_2,D_anat.MK_labels{i}))
                roi_2_idx = i;
            end
        end
        
        

        if (trig_t8)
            system('mkdir figures/T8');
        end
        
        chan_labels = h5readatt(fn_h5Lp,'/h5eeg/eeg','labels');
        bchan_labels = cell(ecog.n_bchan,1);
        for ii0 = 1:ecog.n_bchan
            bchan_labels{ii0} = sprintf('%s-%s',chan_labels{ecog.bip(ii0,1)},chan_labels{ecog.bip(ii0,2)});
        end
        % Build bipolar roi list
    %     has_roi_1 = false(ecog.n_bchan,length(atl_labels));
    %     has_roi_2 = false(ecog.n_bchan,length(atl_labels));
    %     for ii = 1:ecog.n_bchan
    %         has_roi_1(ii,:) = strcmp(atl_labels,roi_1);
    %         has_roi_2(ii,:) = strcmp(atl_labels,roi_2);
    %     end

        % Load graphs
        art_idx = h5read(fn_artLp,'/art_idx');
        art_idxB = (art_idx == 1);
        frac_art = h5readatt(fn_artLp,'/art_idx','frac_art');
        w_art = h5readatt(fn_artLp,'/art_idx','w');
        % trim last time sample of graph
        R = h5read(fn_graphLp,'/R',[1 1],size(art_idx));
        w = double(h5read(fn_graphLp,'/w'));
        [n_comb,n_graph] = size(R);

        % Find suitable bipolar pair
        countp = 1;
        for ii1p = 1:(ecog.n_bchan-1)
            for ii2p = (ii1p+1):ecog.n_bchan
                %fprintf('chan1:%i ii1:%i chan2:%i ii2:%i\n',chan1(count)+1,ii1,chan2(count)+1,ii2)

                % check if bchan pairs map onto specified roi pair
    %             has_roi_f = strcmp(atl_labels{ii1},roi_1) & strcmp(atl_labels{ii2},roi_2);
    %             has_roi_r = strcmp(atl_labels{ii1},roi_2) & strcmp(atl_labels{ii2},roi_1);
                has_roi_f = ( strcmp(atl_labels{ecog.bip(ii1p,1)},roi_1) & strcmp(atl_labels{ecog.bip(ii2p,1)},roi_2) )...
                        | ( strcmp(atl_labels{ecog.bip(ii1p,1)},roi_1) & strcmp(atl_labels{ecog.bip(ii2p,2)},roi_2) )...
                        | ( strcmp(atl_labels{ecog.bip(ii1p,2)},roi_1) & strcmp(atl_labels{ecog.bip(ii2p,1)},roi_2) )...
                        | ( strcmp(atl_labels{ecog.bip(ii1p,2)},roi_1) & strcmp(atl_labels{ecog.bip(ii2p,2)},roi_2) );
                has_roi_r = ( strcmp(atl_labels{ecog.bip(ii1p,1)},roi_2) & strcmp(atl_labels{ecog.bip(ii2p,1)},roi_1) )...
                        | ( strcmp(atl_labels{ecog.bip(ii1p,1)},roi_2) & strcmp(atl_labels{ecog.bip(ii2p,2)},roi_1) )...
                        | ( strcmp(atl_labels{ecog.bip(ii1p,2)},roi_2) & strcmp(atl_labels{ecog.bip(ii2p,1)},roi_1) )...
                        | ( strcmp(atl_labels{ecog.bip(ii1p,2)},roi_2) & strcmp(atl_labels{ecog.bip(ii2p,2)},roi_1) );

                if ((~isnan(bchan1_const)) && (~isnan(bchan2_const)))
                    % if channel pair pick bypass
                    cond_bpair_bypass_f = ((ii1p == bchan1_const) && (ii2p == bchan2_const));
                    cond_bpair_bypass_r = ((ii1p == bchan2_const) && (ii2p == bchan1_const));
                    cond_bpair_bypass = cond_bpair_bypass_f || cond_bpair_bypass_r;
                else
                    % if either bchan1_const or bchan2_const is nan, do all
                    cond_bpair_bypass = true;
                end
                    
                if (((has_roi_f || has_roi_r) && (Dmat(ii1p,ii2p) > dist_thresh)) && cond_bpair_bypass)
                %if ((has_roi_f || has_roi_r) && (Dmat(ii1p,ii2p) > dist_thresh))
                    
                    % ct threshold based on permutations
                    p_val_ct_2 = 0.01/n_comb;
                    %ct_thresh_perm = (binoinv([p_val_ct_2 (1-p_val_ct_2)],n_graph,p_val/bonf_seconds)/n_graph);
                    % ct - 99% CI for 3600 seconds (1 hr) because we
                    % only have 1 hour for monkey control data
                    ct_thresh_perm = binoinv([0.01/n_comb (1-0.01/n_comb)],(1*3600/w),p_val/bonf_seconds)/(1*3600/w);
                    % if in exploratory mode, keep track of number of
                    % covered pairs that pass ct threshold
                    if ((isnan(bchan1_const)) && (isnan(bchan2_const)))
                        ct_thresh_perm_count = [ct_thresh_perm_count; ct(countp) > ct_thresh_perm(2)];
                        if (ct(countp) > ct_thresh_perm(2))
                            ct_thresh_perm_count_mag = [ct_thresh_perm_count_mag; mag(countp)];
                        else
                            ct_thresh_perm_count_mag = [ct_thresh_perm_count_mag; NaN];
                        end
                    end
                    
                    fprintf('%i\n',countp)
                    
                    if (has_roi_r)
                        pname = sprintf('%s__%i_%s-%s_%s__%i_%s-%s_%s__%imm_ct%i_ctt%i-%i_mag%i',...
                            sidp,ii2p,atl_labels{ecog.bip(ii2p,1)},atl_labels{ecog.bip(ii2p,2)},bchan_labels{ii2p},ii1p,atl_labels{ecog.bip(ii1p,1)},atl_labels{ecog.bip(ii1p,2)},...
                            bchan_labels{ii1p},round(Dmat(ii1p,ii2p)),round(1000*ct(countp)),...
                            ceil(1000*ct_thresh_perm(1)),ceil(1000*ct_thresh_perm(2)),round(1000*mag(countp)));
                    else
                        pname = sprintf('%s__%i_%s-%s_%s__%i_%s-%s_%s__%imm_ct%i_ctt%i-%i_mag%i',...
                            sidp,ii1p,atl_labels{ecog.bip(ii1p,1)},atl_labels{ecog.bip(ii1p,2)},bchan_labels{ii1p},ii2p,atl_labels{ecog.bip(ii2p,1)},atl_labels{ecog.bip(ii2p,2)},...
                            bchan_labels{ii2p},round(Dmat(ii1p,ii2p)),round(1000*ct(countp)),...
                            ceil(1000*ct_thresh_perm(1)),ceil(1000*ct_thresh_perm(2)),round(1000*mag(countp)));
                    end
                    % append coh_thresh
                    pname = sprintf('%s_%s_coh%it%i',prefix,pname,iMp,round(1000*coh_thresh(countp)));
                    
                    
                    % Pick time sample to plot
                    Rc = R(countp,:);
                    Risart = art_idx(countp,:);
%                     Ridx = 1:length(Rc);
%                     Ridx = Ridx(randperm(length(Ridx)));
%                     r_idx = NaN;
%                     for i3 = 1:length(Ridx)
%                         if ( (Rc(Ridx(i3)) > coh_thresh(countp)) && (~Risart(Ridx(i3))) )
%                             r_idx = Ridx(i3);
%                             break;
%                         end
%                     end
%                     if (isnan(r_idx)) %if not found
%                         fprintf(2,'W> no time sample found, defaulting to first');
%                         r_idx = 1;
%                     end
%                     r_samp = (r_idx-1)*round(ecog.fs)*w + 1;


                    if (isnan(r_samp_const))
                        % read delta
                        fn_graphLp_delta = sprintf('%s/%s_graph-%s.h5',dir_resLp,sidp,'pcDelta');
                        coh_val_delta = h5read(fn_graphLp_delta,'/R',[countp 1],[1 length(Rc)]);

                        Ridx = 1:length(Rc);
                        idx_rp = randperm(length(Ridx));
                        Ridx = Ridx(idx_rp);
                        coh_val_delta_shuf = coh_val_delta(idx_rp);
                        r_idx = NaN;
                        for i3 = 1:length(Ridx)
                            % Threshold random time sample to plot by coherence values
                            if ( ( (Rc(Ridx(i3)) > coh_thresh(countp)) && (~Risart(Ridx(i3))) )&&(coh_val_delta_shuf(i3) > delta_thresh_4_plot) )
                            % nonsignificant condition
                            %if ( ( (Rc(Ridx(i3)) < coh_thresh(countp)) && (~Risart(Ridx(i3))) )&&(coh_val_delta_shuf(i3) < 0.6) )
                                r_idx = Ridx(i3);
                                break;
                            end
                        end
                        if (isnan(r_idx)) %if not found600
                            fprintf(2,'W> no time sample found, defaulting to first');
                            r_idx = 1;
                        end
                        r_samp = (r_idx-1)*round(ecog.fs)*w + 1;

                    else
                        % bypassround_y_nearest = 100;
                        r_samp = r_samp_const;
                    end

                    r_samp_n = r_samp;

                    % convert r_idx, r_idx_n
                    r_idx = (r_samp - 1)/(round(ecog.fs)*w) + 1;
                    r_idx_n = (r_samp_n - 1)/(round(ecog.fs)*w) + 1;


%                     % Pick time shifted sample to plot
%                     r_idx_n = round(length(Rc)/2);
%                     for i3 = 1:length(Ridx)
%                         if ( (abs(Ridx(i3) - r_idx) > perm_alpha_sec*w) && (~Risart(Ridx(i3))) )
%                             r_idx_n = Ridx(i3);
%                             break;
%                         end
%                     end
%                     r_samp_n = (r_idx_n-1)*round(ecog.fs)*w + 1;

                    % Read timeseries
                    b1 = ecog.bip(ii1p,1:2);
                    b2 = ecog.bip(ii2p,1:2);
                    
%                     if (has_roi_r)
%                         b2t = b2;
%                         b2 = b1;
%                         b1 = b2t;
%                     end
                            
                    v_b1c1 = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(1) r_samp],[1 round(ecog.fs)*w]);
                    v_b1c2 = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(2) r_samp],[1 round(ecog.fs)*w]);
                    v_b2c1 = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(1) r_samp],[1 round(ecog.fs)*w]);
                    v_b2c2 = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(2) r_samp],[1 round(ecog.fs)*w]);
                    v0_b2c1 = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(1) r_samp_n],[1 round(ecog.fs)*w]);
                    v0_b2c2 = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(2) r_samp_n],[1 round(ecog.fs)*w]);
                    % Bipolar montage
                    v_b1 = v_b1c1 - v_b1c2;
                    v_b2 = v_b2c1 - v_b2c2;
                    v0_b2 = v0_b2c1 - v0_b2c2;
                    % Subtract mean
                    v_b1 = v_b1 - mean(v_b1);
                    v_b2 = v_b2 - mean(v_b2);
                    v0_b2 = v0_b2 - mean(v0_b2);

                    % If correlation is negative, flip sign
                    if (corr2(v_b1',v_b2') < 0)
                        v_b2 = (-1) * v_b2;
                        v0_b2 = (-1) * v0_b2;
                    end
                    
                    % Compute null coherence
                    [ coh_val0 ] = coherence( v_b1, v0_b2, round(ecog.fs) );
                    coh_val0 = coh_val0([6,1:5]); % bring broadband to front

                    % Append coherence values to pname
                    coh_vals = zeros(6,1);
                    for i3 = 1:length(metricsp)
                        fn_graphLp = sprintf('%s/%s_graph-%s.h5',dir_resLp,sidp,metricsp{i3});
                        coh_val = h5read(fn_graphLp,'/R',[countp r_idx],[1 1]);
                        coh_vals(i3) = coh_val;
                        mname = metricsp{i3}(3:4);
                        pname = sprintf('%s_%s%i',pname,mname,round(1000*coh_val));
                    end                   
                    pname = sprintf('%s_Samp%i',pname,r_samp);
                    for i3 = 1:length(metricsp)
                        mname = [metricsp{i3}(3:4),'N'];
                        pname = sprintf('%s_%s%i',pname,mname,round(1000*coh_val0(i3)));
                    end
                    pname = sprintf('%s_SampN%i',pname,r_samp_n);

                    fprintf('%s\n',pname);

                    % times
                    t = linspace(0,length(v_b1)/round(ecog.fs),length(v_b1));
                    
                    
                    % plot
                    ytxt = sprintf('IFP (\\muV)');
                    xtxt = 'Time (seconds)';

                    if (trig_t8)
                        h = figure;
                        set(h,'Position',[0 0 PLOT_W_PIX PLOT_H_PIX]);
                        subplot(4,2,1)
                        plot(t,v_b1,'black-','LineWidth',0.5); hold on;
                        ylabel({ytxt;replace(roi_1(1:(end-5)),'_','/')});
                        xlabel(xtxt);
                        if (trig_title)
                            title(sprintf('Coherence: %.2f (\\delta:%.2f \\theta:%.2f \\alpha:%.2f \\beta:%.2f \\gamma:%.2f )',...
                                coh_vals(1),coh_vals(2),coh_vals(3),coh_vals(4),coh_vals(5),coh_vals(6) ));
                        end
                        xticks(min(t):x_div_secs:max(t));
                        set(gca,'TickDir','out');
                        box off;
                        subplot(4,2,3)
                        plot(t,v_b2,'black-','LineWidth',0.5); hold on;
                        ylabel({ytxt;replace(roi_2(1:(end-5)),'_','/')});
                        xlabel(xtxt);
                        xticks(min(t):x_div_secs:max(t));
                        set(gca,'TickDir','out');
                        box off;
                        neurons_f = D_anat.AdjMKneurons(roi_1_idx,roi_2_idx);
                        neurons_r = D_anat.AdjMKneurons(roi_2_idx,roi_1_idx);
                        fprintf('%s-->%s: %i neurons\n',roi_1,roi_2,neurons_f)
                        fprintf('%s-->%s: %i neurons\n',roi_2,roi_1,neurons_r)
                        if (trig_save_fig)
                            print(h,sprintf('figures/T8/%s_src-%s_%i_src-%s_%i',pname,roi_1,neurons_f,roi_2,neurons_r),fig_fmt);
                        end
                        close(h);
                        
                        
                        
                        % plot only brain
                        h = figure;
                        set(h,'Position',[0 0 PLOT_W_PIX PLOT_H_PIX]);
%                         subplot(4,2,5)
%                         plot(t,v_b1,'black-','LineWidth',0.5); hold on;
%                         ylabel({ytxt;replace(roi_1(1:(end-5)),'_','/')});
%                         xlabel(xtxt);
%                         xticks(min(t):x_div_secs:max(t));
%                         set(gca,'TickDir','out');
%                         if (trig_title)
%                             title(sprintf('Coherence: %.2f (\\delta:%.2f \\theta:%.2f \\alpha:%.2f \\beta:%.2f \\gamma:%.2f )',...
%                                 coh_val0(1),coh_val0(2),coh_val0(3),coh_val0(4),coh_val0(5),coh_val0(6) ));
%                         end
%                         box off;
%                         subplot(4,2,7)
%                         %tshift = w*abs(r_idx-r_idx_n);
%                         tshift = w*(r_idx_n - r_idx);
%                         plot(t,v0_b2,'black-','LineWidth',0.5); hold on;
%                         ylabel({ytxt;replace(roi_2(1:(end-5)),'_','/')});
%                         if (abs(tshift) > 3600)
%                             xlabel({xtxt;sprintf('Shifted by %.1f hours',tshift/3600)});
%                         elseif (abs(tshift) > 60)
%                             xlabel({xtxt;sprintf('Shifted by %.1f minutes',tshift/60)});
%                         else
%                             xlabel({xtxt;sprintf('Shifted by %.1f seconds',tshift)});
%                         end
%                         xticks(min(t):x_div_secs:max(t));
%                         set(gca,'TickDir','out');
%                         box off;
                        
                        
                        % reverse labels
                        if (has_roi_r)
                            b2t = b2;
                            b2 = b1;
                            b1 = b2t;
                        end
                            
                        marker_color = 0.99*[1 1 1];
                        marker_outline_color = [0 0 0];
                        hsub = subplot(4,2,2);
                        image(Su.I);hold on
                        for i6 = 1:128
                            if ((i6 == b1(1)) || (i6 == b1(2)) ) 
                                plot(Su.X(i6)-crop_pix_h,Su.Y(i6)-crop_pix,'.','MarkerSize',30,'color', marker_outline_color);hold on
                                plot(Su.X(i6)-crop_pix_h,Su.Y(i6)-crop_pix,'.','MarkerSize',20,'color', marker_color);hold on
                            else
                                %plot(Su.X(i6),Su.Y(i6),'.','MarkerSize',10,'color', 0.9*[1 1 1]);hold on
                            end
                        end
                        daspect([1,1,1])
                        axis equal
                        set(gca,'Visible','off')
                        
                        hsub = subplot(4,2,4);
                        image(Su.I);hold on
                        for i6 = 1:128
                            if ((i6 == b2(1)) || (i6 == b2(2)) ) 
                                plot(Su.X(i6)-crop_pix_h,Su.Y(i6)-crop_pix,'.','MarkerSize',30,'color', marker_outline_color);hold on
                                plot(Su.X(i6)-crop_pix_h,Su.Y(i6)-crop_pix,'.','MarkerSize',20,'color', marker_color);hold on
                            else
                                %plot(Su.X(i6),Su.Y(i6),'.','MarkerSize',10,'color', 0.9*[1 1 1]);hold on
                            end
                        end
                        daspect([1,1,1])
                        axis equal
                        set(gca,'Visible','off')
                        
%                         hsub = subplot(4,2,6);
%                         image(Su.I);hold on
%                         for i6 = 1:128
%                             if ((i6 == b1(1)) || (i6 == b1(2)) ) 
%                                 plot(Su.X(i6),Su.Y(i6),'.','MarkerSize',30,'color', marker_outline_color);hold on
%                                 plot(Su.X(i6),Su.Y(i6),'.','MarkerSize',20,'color', marker_color);hold on
%                             else
%                                 %plot(Su.X(i6),Su.Y(i6),'.','MarkerSize',10,'color', 0.9*[1 1 1]);hold on
%                             end
%                         end
%                         axis equal
%                         set(gca,'Visible','off')
%                         
%                         hsub = subplot(4,2,8);
%                         image(Su.I);hold on
%                         for i6 = 1:128
%                             if ((i6 == b2(1)) || (i6 == b2(2)) ) 
%                                 plot(Su.X(i6),Su.Y(i6),'.','MarkerSize',30,'color', marker_outline_color);hold on
%                                 plot(Su.X(i6),Su.Y(i6),'.','MarkerSize',20,'color', marker_color);hold on
%                             else
%                                 %plot(Su.X(i6),Su.Y(i6),'.','MarkerSize',10,'color', 0.9*[1 1 1]);hold on
%                             end
%                         end
%                         axis equal
%                         set(gca,'Visible','off')
                        
%                         neurons_f = D_anat.AdjMKneurons(roi_1_idx,roi_2_idx);
%                         neurons_r = D_anat.AdjMKneurons(roi_2_idx,roi_1_idx);
%                         fprintf('%s-->%s: %i neurons\n',roi_1,roi_2,neurons_f)
%                         fprintf('%s-->%s: %i neurons\n',roi_2,roi_1,neurons_r)
                        
                        if (trig_save_fig)
                            print(h,sprintf('figures/T8/%s_src-%s_%i_src-%s_%i_brain',pname,roi_1,neurons_f,roi_2,neurons_r),fig_fmt);
                        end
                        close(h);
                        %return
                    end
                    
                    
                end
                
                countp = countp + 1;
            end
        end
        
        
        
        
        
    end
end


% if valid, eport fraction of covered pairs that passed significance test
if (~isempty(ct_thresh_perm_count))
    frac_ct_thresh = sum(ct_thresh_perm_count) / length(ct_thresh_perm_count);
    fprintf('[!] Fraction of covered pairs significantly consistent in time: %.3f\n',frac_ct_thresh);
    fprintf('\tNumbers: %i of %i\n',sum(ct_thresh_perm_count),length(ct_thresh_perm_count))
    
    mag_roi = nanmean(ct_thresh_perm_count_mag);
    mag_roi_std = nanstd(ct_thresh_perm_count_mag);
    fprintf('[!] Average %s coherence: %.3f =- %.3f\n',metric(3:end),mag_roi,mag_roi_std)
end
