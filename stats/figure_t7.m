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
fig_fmt = '-depsc';
trig_t1 = true;
trig_title = false;
trig_save_fig = false;
x_div_secs = 2;
PLOT_W_PIX = 1080;
PLOT_H_PIX = 1080;
delta_thresh_4_plot = 0; % if no time sample is givem, plot with delta > this thresh, 

% Pick areas to plot
% C.AtlROIs{2}.LH.struct_names
atlp = 7;

for ix = 1:4 %1:4
    ct_thresh_perm_count = [];
    ct_thresh_perm_count_mag = [];

    switch(ix)
        case(1)
            % Buschman Miller 2012            
            % (PFC, area 9/46) and the anterior cingulate cortex (ACC, areas 24c and 32).
            prefix = 'BM';
            % MK confirmed connection for all comparisons
            roi_1s = {'9_46d_M132','9_46v_M132'};
            roi_2s = {'24c_M132','32_M132'};
            %121130001 0.490 m43 69 1 ct10 - too fuzzy maybe
            %67855001 0.265 m43 69 17 ct12
            %36505001 0.261 m43 70 1 ct29 mag291 (max ct) (*)
            %23305001 0.274 m43 70 1 ct29 mag291 (max ct)
            %[!] Fraction of covered pairs significantly consistent in time: 0.875
            %       Numbers: 7 of 8
            %[!] Average Broadband coherence: 0.288
            %Subjectsp = {'m00043'};
            bchan1_const = NaN; % 70
            bchan2_const = NaN; % 1
            r_samp_const = 1;
        case(2)
            % Gregoriou Desimone 2009
            prefix = 'GD';
            roi_1s = {'V4_M132'};
            % MK confirmed negative, V4 - 8m
            %roi_2s = {'8l_M132'};
            roi_2s = {'8m_M132','8l_M132'}; % FEF subdivisions found in Markov Kennedy
            % samp: 
            %43932501 0.263 m68 50 32 ct11 (*)
            %2440001 0.262 m68 50 32 ct11
            %2342501 0.260 m68 50 40 ct18
            %        0.259 m38 15 45 ct39 (max ct)
            %[!] Fraction of covered pairs significantly consistent in
            %       time: 0.481 (0.190 V4-8m, 0.829 V4-8l)
            %       time: 0.829 Numbers: 29 of 35 (8l)
            %       time: 0.190 Numbers: 8 of 42 (8m)
            %       time: 0.481 Numbers: 37 of 77 (8l, 8m)

            % V4 - 8m
            %[!] Fraction of covered pairs significantly consistent in time: 0.190
            %	Numbers: 8 of 42
            %[!] Average Broadband coherence: 0.279
            
            % V4 - 8l
            %[!] Fraction of covered pairs significantly consistent in time: 0.829
            %	Numbers: 29 of 35
            %[!] Average Broadband coherence: 0.277

            % V4 - 8m and 8l
            %[!] Fraction of covered pairs significantly consistent in time: 0.481
            %	Numbers: 37 of 77
            %[!] Average Broadband coherence: 0.277
            
            
            %Subjectsp = {'m00068'};
            bchan1_const = NaN; % 50
            bchan2_const = NaN; % 32
            r_samp_const = 1;
        case(3)
            % Bastos Fries 2015
            prefix = 'BF';
            roi_1s = {'V1_M132'};
            % MK confirmed negative: V1 - 8m, V1 - 7A, V1 -> TEO
            % MK confirmed positive: TEO -> V1, V1 <-> DP
            % %roi_2s = {'V1_M132','V2_M132','8l_M132','V4_M132','TEO_M132','DP_M132','8m_M132','7A_M132'};
            %roi_2s = {'TEO_M132','DP_M132'};
            roi_2s = {'V4_M132'};
            % [!] Fraction of covered pairs significantly consistent in time: 1.000
            %[!] Fraction of covered pairs significantly consistent in time: 0.792
            %       Numbers: 156 of 197
            %[!] Average Broadband coherence: 0.291
            
            %roi_2s = {'V2_M132','8l_M132','V4_M132','TEO_M132','DP_M132','8m_M132','7A_M132'};
            
            %19475001 0.272 m48 21 3 ct71 (*)
            %120462501 0.305 m79 59 33 ct3
            %24687501 0.270 m48 27 3 ct58
            %         0.288 m83 41 66 ct106
            %         0.268 m83 38 88 ct123 (39mm)
            %         0.277 m83 35 88 ct143 (30mm)
            %         0.296 m83 34 87 ct124 (33mm)
            %         0.266 m83 31 88 ct128 (41mm)
            %Subjectsp = {'m00048'};

            %roi_2s = {'TEO_M132','DP_M132'};
            %[!] Fraction of covered pairs significantly consistent in time: 1.000
            %	Numbers: 112 of 112
            %[!] Average Broadband coherence: 0.296
            
            %roi_2s = {'TEO_M132'};
            %[!] Fraction of covered pairs significantly consistent in time: 1.000
            %        Numbers: 44 of 44
            %[!] Average Broadband coherence: 0.306
            
            %roi_2s = {'V4_M132'};
            %[!] Fraction of covered pairs significantly consistent in time: 1.000
            %        Numbers: 509 of 509
            %[!] Average Broadband coherence: 0.300
            
            bchan1_const = NaN; % 21
            bchan2_const = NaN; % 3
            r_samp_const = 1;
        case(4)
            % Markov Kennedy 2014 confirmed negatives
            prefix = 'MK';
            roi_1s = {'V1_M132'};
            roi_2s = {'ProM_M132'}; % prefrontal, inferior
            %roi_2s = {'2_M132'}; % somatosensory
            %roi_2s = {'3_M132'}; % somatosensory
            %roi_2s = {'10_M132'}; % F2 no coverage
            %[!] Fraction of covered pairs significantly consistent in
            %   time: 0.000 (V1 - ProM)
            %   time: 0.220 (V1 - 2)
            %   time: 0.578 (V1 - 3)
            %   time: 0.714 (V1 - F1)
            %   time: 0.704 (V1 - 10)
            
            % V1 - ProM
            %[!] Fraction of covered pairs significantly consistent in time: 0.000
            %    Numbers: 0 of 26
            %[!] Average Broadband coherence: NaN
            
            %Subjectsp = {'m00095'};
            bchan1_const = NaN; % 79 V1, V1
            bchan2_const = NaN; % 14 ProM, ProM
            r_samp_const = 47272501;
    end
    
for roi_idx1 = 1:length(roi_1s)
    for roi_idx2 = 1:length(roi_2s)

        roi_1 = roi_1s{roi_idx1}; %'parsorbitalis';
        roi_2 = roi_2s{roi_idx2}; %'inferiorparietal';
        %roi_1 = 'parsopercularis';
        %roi_2 = 'supramarginal';
        %roi_1 = find(strcmp(rois,'parsorbitalis'),1);
        %roi_2 = find(strcmp(rois,'inferiorparietal'),1);
        
        % Don't turn any of these on
        trig_t2 = false;
        trig_t3 = false;
        trig_t4 = false;
        trig_t6 = false;

        if (trig_t1)
            system('mkdir figures/T7');
        end
        if (trig_t2)
            system('mkdir figures/T2');
        end
        if (trig_t3)
            system('mkdir figures/T3');
        end
        if (trig_t4)
            system('mkdir figures/T4');
        end
        if (trig_t6)
            system('mkdir figures/T6');
        end

        %i_s = randi([1 length(Subjects)]);
        for i_s = 1:length(Subjectsp)
        % Main loop
        %for iM = 1:length(metrics)
        for iMp = 1
            metricp = metricsp{iMp};
            sidp = Subjectsp{i_s};
            fn_artLp = sprintf('%s/%s_art.h5',dir_artLp,sidp);
            fn_distLp = sprintf('%s/%s_dists-%s-%i.mat',dir_resLp,sidp,metricp,n_perm);
            fn_graphLp = sprintf('%s/%s_graph-%s.h5',dir_resLp,sidp,metricp);
            fn_permLp = sprintf('/mnt/cuenap2/data/results/coh_w10/%s_perm-%s-%i.h5',sidp,metricp,n_perm);
            fn_h5Lp = sprintf('%s/%s.h5',dir_h5Lp,sidp);
            fn_coregLp = sprintf('%s/%s/label/all_parcellation.mat',dir_corLp,sidp);
            fn_cacheLp = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheLp,sidp,iMp);

            % Check if files exist
            ckf = {fn_artLp,fn_distLp,fn_graphLp,fn_h5Lp,fn_coregLp,fn_cacheLp};
            f_exist = true;
            for j = 1:length(ckf)
                if (~exist(ckf{j},'file'))
                    fprintf(2,'W> File not found: %s\n',ckf{j});
                    f_exist = false;
                    %return
                end
            end

            if (f_exist)
            %ecog = H5eeg(fn_h5);
            load(fn_cacheLp);
            
            % === OVERRIDE ATLAS TO M132 ==================================
            atlp = 7;
            atl_labels = C.AtlLabels{atlp};
            
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
            for ii1 = 1:(ecog.n_bchan-1)
                for ii2 = (ii1+1):ecog.n_bchan
                    %fprintf('chan1:%i ii1:%i chan2:%i ii2:%i\n',chan1(count)+1,ii1,chan2(count)+1,ii2)

                    % check if bchan pairs map onto specified roi pair
        %             has_roi_f = strcmp(atl_labels{ii1},roi_1) & strcmp(atl_labels{ii2},roi_2);
        %             has_roi_r = strcmp(atl_labels{ii1},roi_2) & strcmp(atl_labels{ii2},roi_1);
                    has_roi_f = ( strcmp(atl_labels{ecog.bip(ii1,1)},roi_1) & strcmp(atl_labels{ecog.bip(ii2,1)},roi_2) )...
                            | ( strcmp(atl_labels{ecog.bip(ii1,1)},roi_1) & strcmp(atl_labels{ecog.bip(ii2,2)},roi_2) )...
                            | ( strcmp(atl_labels{ecog.bip(ii1,2)},roi_1) & strcmp(atl_labels{ecog.bip(ii2,1)},roi_2) )...
                            | ( strcmp(atl_labels{ecog.bip(ii1,2)},roi_1) & strcmp(atl_labels{ecog.bip(ii2,2)},roi_2) );
                    has_roi_r = ( strcmp(atl_labels{ecog.bip(ii1,1)},roi_2) & strcmp(atl_labels{ecog.bip(ii2,1)},roi_1) )...
                            | ( strcmp(atl_labels{ecog.bip(ii1,1)},roi_2) & strcmp(atl_labels{ecog.bip(ii2,2)},roi_1) )...
                            | ( strcmp(atl_labels{ecog.bip(ii1,2)},roi_2) & strcmp(atl_labels{ecog.bip(ii2,1)},roi_1) )...
                            | ( strcmp(atl_labels{ecog.bip(ii1,2)},roi_2) & strcmp(atl_labels{ecog.bip(ii2,2)},roi_1) );
                   
                    if ((~isnan(bchan1_const)) && (~isnan(bchan2_const)))
                        % if channel pair pick bypass
                        cond_bpair_bypass_f = ((ii1 == bchan1_const) && (ii2 == bchan2_const));
                        cond_bpair_bypass_r = ((ii1 == bchan2_const) && (ii2 == bchan1_const));
                        cond_bpair_bypass = cond_bpair_bypass_f || cond_bpair_bypass_r;
                    else
                        % if either bchan1_const or bchan2_const is nan, do all
                        cond_bpair_bypass = true;
                    end

                    if (((has_roi_f || has_roi_r) && (Dmat(ii1,ii2) > dist_thresh)) && cond_bpair_bypass)
                    %if ((has_roi_f || has_roi_r) && (Dmat(ii1,ii2) > dist_thresh))
                    %if ( (Dmat(ii1,ii2) > dist_thresh) && cond_bpair_bypass )
                    
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
                        
                        
                        % reverse labels
                        if (has_roi_r)
                            pname = sprintf('%s__%i_%s-%s_%s__%i_%s-%s_%s__%imm_ct%i_ctt%i-%i_mag%i',...
                                sidp,ii2,atl_labels{ecog.bip(ii2,1)},atl_labels{ecog.bip(ii2,2)},bchan_labels{ii2},ii1,atl_labels{ecog.bip(ii1,1)},atl_labels{ecog.bip(ii1,2)},...
                                bchan_labels{ii1},round(Dmat(ii2,ii1)),round(1000*ct(countp)),...
                                ceil(1000*ct_thresh_perm(1)),ceil(1000*ct_thresh_perm(2)),round(1000*mag(countp)));
                        else
                            pname = sprintf('%s__%i_%s-%s_%s__%i_%s-%s_%s__%imm_ct%i_ctt%i-%i_mag%i',...
                                sidp,ii1,atl_labels{ecog.bip(ii1,1)},atl_labels{ecog.bip(ii1,2)},bchan_labels{ii1},ii2,atl_labels{ecog.bip(ii2,1)},atl_labels{ecog.bip(ii2,2)},...
                                bchan_labels{ii2},round(Dmat(ii1,ii2)),round(1000*ct(countp)),...
                                ceil(1000*ct_thresh_perm(1)),ceil(1000*ct_thresh_perm(2)),round(1000*mag(countp)));
                        end
                        % append coh_thresh
                        pname = sprintf('%s_%s_coh%it%i',prefix,pname,iMp,round(1000*coh_thresh(countp)));
                        %return
                        %fprintf('%s\n',pname);
                        
                        % Plot T1
                        % check passes ct threshold
                        %if (trig_t1 && ( ct(count) > ct_thresh ))
                        if ((trig_t1 || trig_t2 || trig_t3) && ( ct(countp) > 0 ))
                            % Pick time sample to plot
                            Rc = R(countp,:);
                            Risart = art_idx(countp,:);
                            
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
                                    %if ( ( (Rc(Ridx(i3)) > coh_thresh(countp)) && (~Risart(Ridx(i3))) )&&(coh_val_delta_shuf(i3) > delta_thresh_4_plot) )
                                    % nonsignificant condition
                                    if ( ( (Rc(Ridx(i3)) < coh_thresh(countp)) && (~Risart(Ridx(i3))) )&&(coh_val_delta_shuf(i3) < 0.5) )
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
%                             
%                             
%                             Ridx = 1:length(Rc);
%                             Ridx = Ridx(randperm(length(Ridx)));
%                             r_idx = NaN;
%                             for i3 = 1:length(Ridx)
%                                 if ( (Rc(Ridx(i3)) > coh_thresh(countp)) && (~Risart(Ridx(i3))) )
%                                     r_idx = Ridx(i3);
%                                     break;
%                                 end
%                             end
%                             if (isnan(r_idx)) %if not found
%                                 fprintf(2,'W> no time sample found, defaulting to first');
%                                 r_idx = 1;
%                             end
%                             r_samp = (r_idx-1)*round(ecog.fs)*w + 1;

%                             % Pick time shifted sample to plot
%                             r_idx_n = round(length(Rc)/2);
%                             for i3 = 1:length(Ridx)
%                                 if ( (abs(Ridx(i3) - r_idx) > perm_alpha_sec*w) && (~Risart(Ridx(i3))) )
%                                     r_idx_n = Ridx(i3);
%                                     break;
%                                 end
%                             end
%                             r_samp_n = (r_idx_n-1)*round(ecog.fs)*w + 1;
                            r_samp_n = r_samp;
                            
                            % convert r_idx, r_idx_n
                            r_idx = (r_samp - 1)/(round(ecog.fs)*w) + 1;
                            r_idx_n = (r_samp_n - 1)/(round(ecog.fs)*w) + 1;

                            % Read timeseries
                            b1 = ecog.bip(ii1,1:2);
                            b2 = ecog.bip(ii2,1:2);
                            
                            if (has_roi_r)
                                b2t = b2;
                                b2 = b1;
                                b1 = b2t;
                            end
                            
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

                            if (trig_t2 || trig_t3)
                                Fs = round(ecog.fs);
                                PLI_S = 56;      % Power line interference frequency                          
                                PLI_E = 64;                                                                   
                                PL2_S = (Fs-180)-2;     % Power line interference frequency second band       
                                PL2_E = (Fs-180)+2;                                                           
                                PL3_S = 117;      % Power line interference frequency third band              
                                PL3_E = 123; 
                                %
                                THZ_S = 17;
                                THZ_E = 23;
                                %
                                DEL_S = 0.5;     % Delta wave                                                 
                                DEL_E = 3;                                                                    
                                THE_S = 3;       % Theta wave                                                 
                                THE_E = 8;                                                                    
                                ALP_S = 8;       % Alpha wave                                                 
                                ALP_E = 12;                                                                   
                                BET_S = 12;      % Beta wave                                                  
                                BET_E = 30;                                                                   
                                GAM_S = 30;      % Gamma wave                                                 
                                GAM_E = 100;                                                                  
                                BRO_S = 0.5;     % Broadband                                                  
                                BRO_E = 125; 
                            end

                            % plot
                            ytxt = sprintf('IFP (\\muV)');
                            xtxt = 'Time (seconds)';

                            if (trig_t1)
                                h = figure('visible','off');
                                set(h,'Position',[0 0 PLOT_W_PIX PLOT_H_PIX]);
                                subplot(4,2,1)
                                plot(t,v_b1,'black-','LineWidth',0.5); hold on;
                                %ylabel({ytxt;replace(roi_1(1:(end-5)),'_','/')});
                                ylabel(ytxt);
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
                                %ylabel({ytxt;replace(roi_2(1:(end-5)),'_','/')});
                                ylabel(ytxt);
                                xlabel(xtxt);
                                xticks(min(t):x_div_secs:max(t));
                                set(gca,'TickDir','out');
                                box off;
                                
                                if (trig_save_fig)
                                    print(h,sprintf('figures/T7/%s',pname),fig_fmt);
                                end
                                close(h);
%                                 
%                                 subplot(4,2,5)
%                                 plot(t,v_b1,'black-','LineWidth',0.5); hold on;
%                                 ylabel({ytxt;replace(roi_1(1:(end-5)),'_','/')});
%                                 xlabel(xtxt);
%                                 xticks(min(t):x_div_secs:max(t));
%                                 set(gca,'TickDir','out');
%                                 if (trig_title)
%                                     title(sprintf('Coherence: %.2f (\\delta:%.2f \\theta:%.2f \\alpha:%.2f \\beta:%.2f \\gamma:%.2f )',...
%                                         coh_val0(1),coh_val0(2),coh_val0(3),coh_val0(4),coh_val0(5),coh_val0(6) ));
%                                 end
%                                 box off;
%                                 subplot(4,2,7)
%                                 %tshift = w*abs(r_idx-r_idx_n);
%                                 tshift = w*(r_idx_n - r_idx);
%                                 plot(t,v0_b2,'black-','LineWidth',0.5); hold on;
%                                 ylabel({ytxt;replace(roi_2(1:(end-5)),'_','/')});
%                                 if (abs(tshift) > 3600)
%                                     xlabel({xtxt;sprintf('Shifted by %.1f hours',tshift/3600)});
%                                 elseif (abs(tshift) > 60)
%                                     xlabel({xtxt;sprintf('Shifted by %.1f minutes',tshift/60)});
%                                 else
%                                     xlabel({xtxt;sprintf('Shifted by %.1f seconds',tshift)});
%                                 end
%                                 xticks(min(t):x_div_secs:max(t));
%                                 set(gca,'TickDir','out');
%                                 box off;
                                

                                % Plot brain separately
                                
                                h = figure('visible','off');
                                set(h,'Position',[0 0 PLOT_W_PIX PLOT_H_PIX]);
                                hsub = subplot(4,2,2);
                                if (trig_save_fig)
                                    hb = brainplot_one(sidp,[b1(1),b1(2)],'');
                                    copyobj(hb,hsub);
                                end
                                if (trig_title)
                                    title(sidp);
                                end
                                hsub = subplot(4,2,4);
                                if (trig_save_fig)
                                    hb = brainplot_one(sidp,[b2(1),b2(2)],'');
                                    copyobj(hb,hsub);
                                end
                                if (trig_title)
                                    title(sidp);
                                end
%                                 hsub = subplot(4,2,6);
%                                 hb = brainplot_one(sid,[b1(1),b1(2)],'');
%                                 copyobj(hb,hsub);
%                                 if (trig_title)
%                                     title(sid);
%                                 end
%                                 hsub = subplot(4,2,8);
%                                 hb = brainplot_one(sid,[b2(1),b2(2)],'');
%                                 copyobj(hb,hsub);
%                                 if (trig_title)
%                                     title(sid);
%                                 end

                                %return

                                if (trig_save_fig)
                                    print(h,sprintf('figures/T7/%s_brain',pname),fig_fmt);
                                end
                                close(h);
                            end


                            if (trig_t2)
                                % backup
                                v_b1_b = v_b1;
                                v_b2_b = v_b2;
                                v0_b2_b = v0_b2;
                                forder = 4;

                                for i4 = 1:length(metricsp)
                                    % restore from backup
                                    v_b1 = v_b1_b;
                                    v_b2 = v_b2_b;
                                    v0_b2 = v0_b2_b;

                                    switch metricsp{i4}
                                        case 'pcBroadband'
                                            [b,a] = butter(forder,[BRO_S BRO_E-1]/(Fs/2),'bandpass');
                                        case 'pcDelta'
                                            [b,a] = butter(forder,[DEL_S DEL_E]/(Fs/2),'bandpass');
                                        case 'pcTheta'
                                            [b,a] = butter(forder,[THE_S THE_E]/(Fs/2),'bandpass');
                                        case 'pcAlpha'
                                            [b,a] = butter(forder,[ALP_S ALP_E]/(Fs/2),'bandpass');
                                        case 'pcBeta'
                                            [b,a] = butter(forder,[BET_S BET_E]/(Fs/2),'bandpass');
                                        case 'pcGamma'
                                            [b,a] = butter(forder,[GAM_S GAM_E]/(Fs/2),'bandpass');
                                    end
                                    v_b1 = filtfilt(b,a,double(v_b1));
                                    v_b2 = filtfilt(b,a,double(v_b2));
                                    v0_b2 = filtfilt(b,a,double(v0_b2));
        %                             v_b1 = FilterData(v_b1,Fs,'notch',f);

                                    h = figure('visible','off');
                                    set(h,'Position',[0 0 1080 1080])
                                    subplot(4,2,1)
                                    plot(t,v_b1,'black-','LineWidth',0.5); hold on;
                                    ylabel({ytxt;roi_1});
                                    xlabel(xtxt);
                                    title(sprintf('Coherence: %.2f (\\delta:%.2f \\theta:%.2f \\alpha:%.2f \\beta:%.2f \\gamma:%.2f )',...
                                        coh_vals(1),coh_vals(2),coh_vals(3),coh_vals(4),coh_vals(5),coh_vals(6) ));
                                    xticks(min(t):x_div_secs:max(t));
                                    box off;
                                    subplot(4,2,3)
                                    plot(t,v_b2,'black-','LineWidth',0.5); hold on;
                                    ylabel({ytxt;roi_2});
                                    xlabel(xtxt);
                                    xticks(min(t):x_div_secs:max(t));

                                    box off;
                                    subplot(4,2,5)
                                    plot(t,v_b1,'black-','LineWidth',0.5); hold on;
                                    ylabel({ytxt;roi_1});
                                    xlabel(xtxt);
                                    xticks(min(t):x_div_secs:max(t));
                                    title(sprintf('Coherence: %.2f (\\delta:%.2f \\theta:%.2f \\alpha:%.2f \\beta:%.2f \\gamma:%.2f )',...
                                        coh_val0(1),coh_val0(2),coh_val0(3),coh_val0(4),coh_val0(5),coh_val0(6) ));
                                    box off;
                                    subplot(4,2,7)
                                    %tshift = w*abs(r_idx-r_idx_n);
                                    tshift = w*(r_idx_n - r_idx);
                                    plot(t,v0_b2,'black-','LineWidth',0.5); hold on;
                                    ylabel({ytxt;roi_2});
                                    if (abs(tshift) > 3600)
                                        xlabel({xtxt;sprintf('Shifted by %.1f hours',tshift/3600)});
                                    elseif (abs(tshift) > 60)
                                        xlabel({xtxt;sprintf('Shifted by %.1f minutes',tshift/60)});
                                    else
                                        xlabel({xtxt;sprintf('Shifted by %.1f seconds',tshift)});
                                    end
                                    xticks(min(t):x_div_secs:max(t));
                                    box off;

                                    if (trig_save_fig)
                                        hsub = subplot(4,2,2);
                                        hb = brainplot_one(sidp,[b1(1),b1(2)],'');
                                        copyobj(hb,hsub);
                                        title(sidp);
                                        hsub = subplot(4,2,4);
                                        hb = brainplot_one(sidp,[b2(1),b2(2)],'');
                                        copyobj(hb,hsub);
                                        title(sidp);
                                        hsub = subplot(4,2,6);
                                        hb = brainplot_one(sidp,[b1(1),b1(2)],'');
                                        copyobj(hb,hsub);
                                        title(sidp);
                                        hsub = subplot(4,2,8);
                                        hb = brainplot_one(sidp,[b2(1),b2(2)],'');
                                        copyobj(hb,hsub);
                                        title(sidp);
                                    end
                                        

                                    %return

                                    if (trig_save_fig)
                                        print(h,sprintf('figures/T2/%s_FILT-%s',pname,metricsp{i4}(3:4)),fig_fmt);
                                    end
                                    close(h);
                                end

                                %return
                            end


                            if (trig_t3)


                                [cxx,f] = mscohere(v_b1', v_b2',hamming(Fs/1),[],[],Fs);
                                cxx = sqrt(cxx);
                                [cxx0,f] = mscohere(v_b1', v0_b2',hamming(Fs/1),[],[],Fs);
                                cxx0 = sqrt(cxx0);

                                mask_del = ((f > DEL_S) & (f < DEL_E));
                                mask_the = ((f >= THE_S) & (f < THE_E));
                                mask_alp = ((f >= ALP_S) & (f < ALP_E));
                                mask_bet = ((f >= BET_S) & (f < THZ_S)) | ((f > THZ_E) & (f < BET_E));
                                mask_gam = ((f >= GAM_S) & (f < PLI_S)) | ((f > PLI_E) & (f < PL2_S)) | ((f > PL2_E) & (f < GAM_E));
                                mask_bro = ((f > BRO_S) & (f < THZ_S)) | ((f > THZ_E) & (f < PLI_S)) | ((f > PLI_E) & (f < PL2_S)) | ((f > PL2_E) & (f < PL3_S)) | ((f > PL3_E) & (f < BRO_E));

                                color_divider = ones(1,3)*0.3;
                                mean_line_width = 4;
                                line_fontsize = 20;
                                h = figure('visible','off');
                                subplot(2,1,1)
                                set(h,'Position',[0 0 1080 1080])
                                cxxp = cxx;
                                cxxp(~mask_bro) = nan;
                                plot(f,cxxp,'-','LineWidth',0.5,'color',0.3*[1 1 1]); hold on;
                                plot([PLI_S PLI_S],[0 1],'--','color',color_divider); hold on;
                                plot([PLI_E PLI_E],[0 1],'--','color',color_divider); hold on;
                                plot([PL2_S PL2_S],[0 1],'--','color',color_divider); hold on;
                                plot([PL2_E PL2_E],[0 1],'--','color',color_divider); hold on;
                                plot([PL3_S PL3_S],[0 1],'--','color',color_divider); hold on;
                                plot([PL3_E PL3_E],[0 1],'--','color',color_divider); hold on;
                                plot([THZ_S THZ_S],[0 1],'--','color',color_divider); hold on;
                                plot([THZ_E THZ_E],[0 1],'--','color',color_divider); hold on;
                                coh_del = mean(cxx(mask_del));
                                plot([DEL_S DEL_E],[coh_del coh_del],'black-','LineWidth',mean_line_width); hold on;
                                text(mean([DEL_S,DEL_E]),double(coh_del),sprintf('\\delta'),'HorizontalAlignment',...
                                    'center','VerticalAlignment','bottom','FontWeight','bold','FontSize',line_fontsize)
                                coh_the = mean(cxx(mask_the));
                                plot([THE_S THE_E],[coh_the coh_the],'black-','LineWidth',mean_line_width); hold on;
                                text(mean([THE_S,THE_E]),double(coh_the),sprintf('\\theta'),'HorizontalAlignment',...
                                    'center','VerticalAlignment','bottom','FontWeight','bold','FontSize',line_fontsize)
                                coh_alp = mean(cxx(mask_alp));
                                plot([ALP_S ALP_E],[coh_alp coh_alp],'black-','LineWidth',mean_line_width); hold on;
                                text(mean([ALP_S,ALP_E]),double(coh_alp),sprintf('\\alpha'),'HorizontalAlignment',...
                                    'center','VerticalAlignment','bottom','FontWeight','bold','FontSize',line_fontsize)
                                coh_bet = mean(cxx(mask_bet));
                                plot([BET_S BET_E],[coh_bet coh_bet],'black-','LineWidth',mean_line_width); hold on;
                                text(mean([BET_S,BET_E]),double(coh_bet),sprintf('\\beta'),'HorizontalAlignment',...
                                    'center','VerticalAlignment','bottom','FontWeight','bold','FontSize',line_fontsize)
                                coh_gam = mean(cxx(mask_gam));
                                plot([GAM_S GAM_E],[coh_gam coh_gam],'black-','LineWidth',mean_line_width); hold on;
                                text(mean([GAM_S,GAM_E]),double(coh_gam),sprintf('\\gamma'),'HorizontalAlignment',...
                                    'center','VerticalAlignment','bottom','FontWeight','bold','FontSize',line_fontsize)
                                xlabel('Frequency (Hz)')
                                ylabel('Magnitude Coherence')
                                yticks(0:0.25:1)
                                box off;
                                axis tight;
                                subplot(2,1,2)
                                set(h,'Position',[0 0 1080 1080])
                                cxx0p = cxx0;
                                cxx0p(~mask_bro) = nan;
                                plot(f,cxx0p,'-','LineWidth',0.5,'color',0.3*[1 1 1]); hold on;
                                plot([PLI_S PLI_S],[0 1],'--','color',color_divider); hold on;
                                plot([PLI_E PLI_E],[0 1],'--','color',color_divider); hold on;
                                plot([PL2_S PL2_S],[0 1],'--','color',color_divider); hold on;
                                plot([PL2_E PL2_E],[0 1],'--','color',color_divider); hold on;
                                plot([PL3_S PL3_S],[0 1],'--','color',color_divider); hold on;
                                plot([PL3_E PL3_E],[0 1],'--','color',color_divider); hold on;
                                plot([THZ_S THZ_S],[0 1],'--','color',color_divider); hold on;
                                plot([THZ_E THZ_E],[0 1],'--','color',color_divider); hold on;
                                coh_del = mean(cxx0(mask_del));
                                plot([DEL_S DEL_E],[coh_del coh_del],'black-','LineWidth',mean_line_width); hold on;
                                text(mean([DEL_S,DEL_E]),double(coh_del),sprintf('\\delta'),'HorizontalAlignment',...
                                    'center','VerticalAlignment','bottom','FontWeight','bold','FontSize',line_fontsize)
                                coh_the = mean(cxx0(mask_the));
                                plot([THE_S THE_E],[coh_the coh_the],'black-','LineWidth',mean_line_width); hold on;
                                text(mean([THE_S,THE_E]),double(coh_the),sprintf('\\theta'),'HorizontalAlignment',...
                                    'center','VerticalAlignment','bottom','FontWeight','bold','FontSize',line_fontsize)
                                coh_alp = mean(cxx0(mask_alp));
                                plot([ALP_S ALP_E],[coh_alp coh_alp],'black-','LineWidth',mean_line_width); hold on;
                                text(mean([ALP_S,ALP_E]),double(coh_alp),sprintf('\\alpha'),'HorizontalAlignment',...
                                    'center','VerticalAlignment','bottom','FontWeight','bold','FontSize',line_fontsize)
                                coh_bet = mean(cxx0(mask_bet));
                                plot([BET_S BET_E],[coh_bet coh_bet],'black-','LineWidth',mean_line_width); hold on;
                                text(mean([BET_S,BET_E]),double(coh_bet),sprintf('\\beta'),'HorizontalAlignment',...
                                    'center','VerticalAlignment','bottom','FontWeight','bold','FontSize',line_fontsize)
                                coh_gam = mean(cxx0(mask_gam));
                                plot([GAM_S GAM_E],[coh_gam coh_gam],'black-','LineWidth',mean_line_width); hold on;
                                text(mean([GAM_S,GAM_E]),double(coh_gam),sprintf('\\gamma'),'HorizontalAlignment',...
                                    'center','VerticalAlignment','bottom','FontWeight','bold','FontSize',line_fontsize)
                                xlabel('Frequency (Hz)')
                                tshift = w*(r_idx_n - r_idx);
                                if (abs(tshift) > 3600)
                                    shift_txt = sprintf('(shifted by %.1f hours)',tshift/3600);
                                elseif (abs(tshift) > 60)
                                    shift_txt = sprintf('(shifted by %.1f minutes)',tshift/60);
                                else
                                    shift_txt = sprintf('(shifted by %.1f seconds)',tshift);
                                end
                                ylabel({'Magnitude Coherence';shift_txt})
                                yticks(0:0.25:1)
                                box off;
                                axis tight;

                                %return

                                if (trig_save_fig)
                                    print(h,sprintf('figures/T3/%s',pname),fig_fmt);
                                end
                                close(h);
                            end



                        end


                        if (trig_t4)
                            h = figure('visible','off');
                            set(h,'Position',[0 0 0.5*1080 0.5*1080])

                            Rc = R(countp,:);
                            Risart = art_idx(countp,:);

                            % Pick artifact to plot
                            Ridx = 1:length(Rc);
                            Ridx = Ridx(randperm(length(Ridx)));
                            r_idx = NaN;
                            for i3 = 1:length(Ridx)
                                if ( Risart(Ridx(i3)) )
                                    r_idx = Ridx(i3);
                                    break;
                                end
                            end
                            if (isnan(r_idx)) %if not found
                                fprintf(2,'W> no time sample found, defaulting to first');
                                r_idx = 1;
                            end
                            r_samp = (r_idx-1)*round(ecog.fs)*w + 1;                    
                            r_samp_1s = (r_idx-1)*w + 1;

                            % Load 1 second resolution artifacts
                            a1_b1 = h5read(fn_artLp,'/artifacts',[ii1 r_samp_1s],[1 w]);
                            a1_b2 = h5read(fn_artLp,'/artifacts',[ii2 r_samp_1s],[1 w]);

                            % Read timeseries
                            b1 = ecog.bip(ii1,1:2);
                            b2 = ecog.bip(ii2,1:2);
                            v_b1c1 = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(1) r_samp],[1 round(ecog.fs)*w]);
                            v_b1c2 = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(2) r_samp],[1 round(ecog.fs)*w]);
                            v_b2c1 = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(1) r_samp],[1 round(ecog.fs)*w]);
                            v_b2c2 = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(2) r_samp],[1 round(ecog.fs)*w]);
                            % Bipolar montage
                            v_b1 = v_b1c1 - v_b1c2;
                            v_b2 = v_b2c1 - v_b2c2;
                            % Subtract mean
                            v_b1 = v_b1 - mean(v_b1);
                            v_b2 = v_b2 - mean(v_b2);

                            % times
                            secs = length(v_b1)/round(ecog.fs);
                            t = linspace(0,secs,length(v_b1));

                            color1 = [1 0 0]; % V < 10 uV
                            color2 = [0 0 1]; % V > 2000 uV
                            color3 = [0 1 0]; % dV/dt > 100 uV/ms
                            color4 = [1 0 1]; % horiz
                            color5 = [0 1 1]; % vert
                            subplot(2,1,1)
                            plot(t,v_b1,'black'); hold on;
                            cond1 = repelem(a1_b1,round(ecog.fs)) == 1;
                            v_b1p = v_b1;
                            v_b1p(~cond1) = NaN;
                            plot(t,v_b1p,'-','color',color1); hold on;
                            cond2 = repelem(a1_b1,round(ecog.fs)) == 2;
                            v_b1p = v_b1;
                            v_b1p(~cond2) = NaN;
                            plot(t,v_b1p,'-','color',color2); hold on;
                            cond3 = repelem(a1_b1,round(ecog.fs)) == 3;
                            v_b1p = v_b1;
                            v_b1p(~cond3) = NaN;
                            plot(t,v_b1p,'-','color',color3); hold on;
                            cond4 = repelem(a1_b1,round(ecog.fs)) == 4;
                            v_b1p = v_b1;
                            v_b1p(~cond4) = NaN;
                            plot(t,v_b1p,'-','color',color4); hold on;
                            cond5 = repelem(a1_b1,round(ecog.fs)) == 5;
                            v_b1p = v_b1;
                            v_b1p(~cond5) = NaN;
                            plot(t,v_b1p,'-','color',color5); hold on;
                            xticks(0:(secs/10):secs);
                            xlabel('Time (seconds)')
                            ylabel({sprintf('IFP \\muV'),roi_1})
                            box off;

                            subplot(2,1,2)
                            plot(t,v_b2,'black'); hold on;
                            cond1 = repelem(a1_b2,round(ecog.fs)) == 1;
                            v_b2p = v_b2;
                            v_b2p(~cond1) = NaN;
                            plot(t,v_b2p,'-','color',color1); hold on;
                            cond2 = repelem(a1_b2,round(ecog.fs)) == 2;
                            v_b2p = v_b2;
                            v_b2p(~cond2) = NaN;
                            plot(t,v_b2p,'-','color',color2); hold on;
                            cond3 = repelem(a1_b2,round(ecog.fs)) == 3;
                            v_b2p = v_b2;
                            v_b2p(~cond3) = NaN;
                            plot(t,v_b2p,'-','color',color3); hold on;
                            cond4 = repelem(a1_b2,round(ecog.fs)) == 4;
                            v_b2p = v_b2;
                            v_b2p(~cond4) = NaN;
                            plot(t,v_b2p,'-','color',color4); hold on;
                            cond5 = repelem(a1_b2,round(ecog.fs)) == 5;
                            v_b2p = v_b2;
                            v_b2p(~cond5) = NaN;
                            plot(t,v_b2p,'-','color',color5); hold on;
                            xticks(0:(secs/10):secs);
                            xlabel('Time (seconds)')
                            ylabel({sprintf('IFP \\muV'),roi_2})
                            box off;
                            %return
                            if (trig_save_fig)
                                print(h,sprintf('figures/T4/%s_Samp%i',pname,r_samp),fig_fmt);
                            end
                            close(h);
                        end


                        if (trig_t6)
                            x = linspace(0,1,250);
                            h = figure('visible','off');
                            set(h,'Position',[0 0 1080 1080])

                            for i5 = 1:length(metricsp)

                                metricp = metricsp{i5};
                                subplot(length(metricsp)/2,2,i5)
                                fn_artLp = sprintf('%s/%s_art.h5',dir_artLp,sidp);
                                fn_distLp = sprintf('%s/%s_dists-%s-%i.mat',dir_resLp,sidp,metricp,n_perm);
                                fn_graphLp = sprintf('%s/%s_graph-%s.h5',dir_resLp,sidp,metricp);
                                fn_permLp = sprintf('/mnt/cuenap2/data/results/coh_w10/%s_perm-%s-%i.h5',sidp,metricp,n_perm);

                                load(fn_distLp);
                                Distr = d{countp};
                                Rc = h5read(fn_graphLp,'/R',[countp 1],[1 n_graph]);
                                Risart = h5read(fn_artLp,'/art_idx',[countp 1],[1 n_graph]) ~= 0;
                                Rc(Risart) = NaN;
                                Rperm = h5read(fn_permLp,'/R',[countp 1],[1 n_perm]);

                                [counts,~] = hist(Rc,x);
                                plot(x,counts/trapz(x,counts),'-','color',[1 1 1]*0); hold on;
                                [countsP,~] = hist(Rperm,x);
                                plot(x,countsP/trapz(x,countsP),'-','color',[1 0 0]*0.3); hold on;
                                plot(x,pdf(d{countp},x),'--','color',[1 0 0]*0.3);
                                box off;
                                xlabel(sprintf('Coherence'));
                                ylabel({'pdf';sprintf('%s',metricp(3:end))});


                            end

                            %return

                            if (trig_save_fig)
                                print(h,sprintf('figures/T6/%s',pname),fig_fmt);
                            end
                            close(h);

                        end

                    end
                    countp = countp + 1;
                end
            end
            
            end % if f_exist

        end

        end

    end
end

% if valid, report fraction of covered pairs that passed significance test

if (~isempty(ct_thresh_perm_count))
    frac_ct_thresh = sum(ct_thresh_perm_count) / length(ct_thresh_perm_count);
    fprintf('[!] Fraction of covered pairs significantly consistent in time: %.3f\n',frac_ct_thresh);
    fprintf('\tNumbers: %i of %i\n',sum(ct_thresh_perm_count),length(ct_thresh_perm_count))
    
    mag_roi = nanmean(ct_thresh_perm_count_mag);
    mag_roi_std = nanstd(ct_thresh_perm_count_mag);
    fprintf('[!] Average %s coherence: %.3f =- %.3f\n',metric(3:end),mag_roi,mag_roi_std)
    %return
end

% Clear loop indices
clear i;
clear j;

end

% Print finish message
fprintf('Done.\n')
