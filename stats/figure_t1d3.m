% Feb 19, 2019
% bchan1: 35
% bchan2: 84
% n_bchan: 91
% n_pairs: 4095
% n_pass_dist: 3193
% n_pd_pairs_bchan1: 73
% n_pd_pairs_bchan2: 58
% n_pd_pairs_pos_bchan1: 33
% n_pd_pairs_pos_bchan2: 19
% n_pd_pairs_neg_both: 14
% n_pd_pairs_pos_both: 22
% n_pd_pairs_pos_b1: 0
% n_pd_pairs_pos_b2: 8
% n_pd: 44
% Done.


close all;
clear;
rng shuffle;

metricpp = 'pcBroadband';
n_perm = 10000;
perm_alpha_sec = 12;
cp_thresh_override = 0.05;
p_val = 0.01;

% Fast i/o definitions
dir_resLpp =  '/media/jerry/untitled/results/coh_w10'; %'/home/jerry/data/results/coh_w10'; %'/media/jerry/KLAB101/results/coh_w10';
dir_corLpp = '/media/jerry/internal/data/coreg';
dir_cacheLpp = './cache';
%subjects_dirL = '/mnt/cuenap_ssd/coregistration';

% Slow i/o definitions
dir_h5L = '/media/jerry/untitled/h5_notch20';%'/media/klab/internal/data/h5_notch20'; %'/media/klab/KLAB101/h5_notch20';
dir_artLpp = sprintf('%s/art_nosz',dir_h5L);
%dir_h5L = '/mnt/cuenap/data/h5_notch20';

metricspp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'}; % ,'pcDelta'

% Patients
Subjectspp = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124',...
    'mSu'};

% Exclude monkey
Subjectspp = Subjectspp(1:(end-1));

% bypass
%[*] cost: -0.213723722	m00100	30	80	7315001
Subjectspp = {'m00005'};
%r_samp_const = 106865001; % set to nan to pick randomly
%r_samp_const = nan;%90545001;
%r_samp_const = 22755001; 
r_samp_const = 72056321;

% Second way
% bchan1_const = NaN; % 84 - 35
% bchan2_const = 84;
% %bchan2_const = [44, 55, 73]; % NaN;
% roi_1 = 'superiortemporal';
% roi_2 = 'parsopercularis';

% One way
bchan1_const = 35; % 84 - 35
bchan2_const = 84; % NaN;
roi_1 = 'superiortemporal';
roi_2 = 'parsopercularis';
roi_1_rois = {'parsorbitalis','parstriangularis','parsopercularis'};
roi_2_chans = {'PT5','PT6','PT7'};

% bchan1_const = 13; % 27 13
% bchan2_const = NaN; % 46 58
%bchan2_const = [41,58]; %bchan1_const = 30
%bchan2_const = [94]; %bchan1_const = 80
%bchan2_const = [58,53,57,61]; %m00001


% roi_1 = 'parstriangularis';
% roi_2 = 'inferiorparietal';

% === m00084 ===================
% - 30 precentral
% - 33 supramarginal
% - 49 middletemporal
% ==============================

% Subjectspp = {'m00061'};
% r_samp_const = 79817501;
% bchan1_const = 27;
% bchan2_const = NaN;
% roi_1 = 'parsopercularis';
% roi_2 = 'inferiorparietal';

% Subjectspp = {'m00047'};
% r_samp_const = 64875001;
% bchan1_const = 54;
% bchan2_const = NaN;
% roi_1 = 'parstriangularis';
% roi_2 = 'inferiorparietal';

%Strengths = [0,0.25,0.5,0.75,1];
%Strengths = [0.1,0.9,1];
Strengths = [1];
%Strengths = linspace(0,1,200);
%r_samp_const = NaN;

% Pick areas to plot
% C.AtlROIs{2}.LH.struct_names
% (PFC, area 9/46) and the anterior cingulate cortex (ACC, areas 24c and 32).

%roi_1 = 'parsorbitalis';

% roi_1 = 'parstriangularis';
% roi_2 = 'inferiorparietal';

%roi_1 = 'parsopercularis';
%roi_2 = 'supramarginal';
%roi_1 = find(strcmp(rois,'parsorbitalis'),1);
%roi_2 = find(strcmp(rois,'inferiorparietal'),1);

% Triggers
system('mkdir figures');
fig_fmt = '-dpng';
trig_title = false; % show titles
trig_both_elec = true; % bipolar
trig_eps = true;
trig_brain = true;
trig_t1 = true;
trig_t2 = false;
trig_t3 = false;
trig_t4 = false;
trig_t6 = false;
trig_t9 = false;
trig_t9d1 = false;

% time shift
trig_samp_offset_override = true;
samp_offset_override = 600; %seconds
round_y_nearest = 100; % uV

% plot 
PLOT_W_PIX = 1080;
PLOT_H_PIX = 1080;
x_div_secs = 2;
col_faint = 0.5*[1 1 1];

if (trig_t1)
    %system('mkdir figures/T1');
    system('mkdir figures/T1d3');
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
if (trig_t9)
    system('mkdir figures/T9');
end
if (trig_t9d1)
    system('mkdir figures/T9d1');
end

%i_s = randi([1 length(Subjects)]);
for i_str = 1:length(Strengths)
    strength = Strengths(i_str);
for i_s = 1:length(Subjectspp)
% Main loop
for iMpp = 1:1%length(metricspp)
%for iM = 1
    metricpp = metricspp{iMpp};
    sidpp = Subjectspp{i_s};
    fn_artLpp = sprintf('%s/%s_art.h5',dir_artLpp,sidpp);
    fn_distLpp = sprintf('%s/%s_dists-%s-%i.mat',dir_resLpp,sidpp,metricpp,n_perm);
    fn_graphLpp = sprintf('%s/%s_graph-%s.h5',dir_resLpp,sidpp,metricpp);
    fn_permLpp = sprintf('/mnt/cuenap2/data/results/coh_w10/%s_perm-%s-%i.h5',sidpp,metricpp,n_perm);
    fn_h5Lpp = sprintf('%s/%s.h5',dir_h5L,sidpp);
    fn_coregLpp = sprintf('%s/%s/label/all_parcellation.mat',dir_corLpp,sidpp);
    fn_cacheLpp = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheLpp,sidpp,iMpp);
    
    % Check if files exist
    ckf = {fn_artLpp,fn_distLpp,fn_graphLpp,fn_h5Lpp,fn_coregLpp,fn_cacheLpp};
    for j = 1:length(ckf)
        if (~exist(ckf{j},'file'))
            fprintf(2,'W> File not found: %s\n',ckf{j});
            return
        end
    end
    
    %ecog = H5eeg(fn_h5);
    load(fn_cacheLpp);
    
    % --- Apply threshold override --------------------------------------------
    cp_thresh = cp_thresh_override;
    %return
    
    chan_labels = h5readatt(fn_h5Lpp,'/h5eeg/eeg','labels');
    bchan_labels = cell(ecog.n_bchan,1);
    for ii0 = 1:ecog.n_bchan
        bchan_labels{ii0} = sprintf('%s-%s',chan_labels{ecog.bip(ii0,1)},chan_labels{ecog.bip(ii0,2)});
    end
    % Build bipolar roi list
%     has_roi_1 = false(ecog.n_bchan,length(atl_labels));
%     has_roi_2 = false(ecog.n_bchan,length(atl_labels));
%     for ii = 1:ecog.n_bchan600
%         has_roi_1(ii,:) = strcmp(atl_labels,roi_1);
%         has_roi_2(ii,:) = strcmp(atl_labels,roi_2);
%     end

    % Load graphs
    art_idx = h5read(fn_artLpp,'/art_idx');
    art_idxB = (art_idx == 1);
    frac_art = h5readatt(fn_artLpp,'/art_idx','frac_art');
    w_art = h5readatt(fn_artLpp,'/art_idx','w');
    % trim last time sample of graph
    R = h5read(fn_graphLpp,'/R',[1 1],size(art_idx));
    w = double(h5read(fn_graphLpp,'/w'));
    [~,n_graph] = size(R);
    
    % Find suitable bipolar pair
    count = 1;
    n_pass_dist_thresh = 0;
    n_corr_bchan1 = 0;
    n_corr_bchan2 = 0;
    n_corr_bchan1_pos = 0;
    n_corr_bchan2_pos = 0;
    n_not_coh_with_both = 0;
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
            if ((~isnan(bchan1_const)) && all(~isnan(bchan2_const)))
                % if channel pair pick bypass
                if (length(bchan2_const) == 1)
                    % Only 1 chan2 specified
                    %cond_bpair_bypass = (ii1 == bchan1_const) && (ii2 == bchan2_const);
                    % both forward and backward
                    found_fwd = (ii1 == bchan1_const);
                    found_bak = (ii2 == bchan1_const);
                    cond_bpair_bypass = found_fwd || found_bak;
                else
                    cond_bpair_bypass = (ii1 == bchan1_const) && any(ii2 == bchan2_const);
                end
                t1d3_cond = true;
            elseif (isnan(bchan1_const))
                cond_bpair_bypass = true;
                t1d3_cond = (ii2 == bchan2_const);
            else
                % if either bchan1_const or bchan2_const is nan, do all
                cond_bpair_bypass = true;
                % pick pairs that are not significant
                t1d3_cond = (ii1 == bchan1_const) || (ii2 == bchan1_const);% && (ct(count) < ct_thresh);
            end
            
            if ( ((Dmat(ii1,ii2) > dist_thresh) && cond_bpair_bypass ) && t1d3_cond )
                pname = sprintf('%s__%i_%s_%s__%i_%s_%s__%imm_ct%i_mag%i_coh%it%i_str%i',...
                    sidpp,ii1,atl_labels{ecog.bip(ii1,1)},bchan_labels{ii1},ii2,atl_labels{ecog.bip(ii2,1)},...
                    bchan_labels{ii2},round(Dmat(ii1,ii2)),round(1000*ct(count)),round(1000*mag(count)),...
                    iMpp,round(1000*coh_thresh(count)),round(1000*strength)); % (has_roi_f || has_roi_r) &&
                %fprintf('%s\n',pname);
                
                
                
                
                % Plot T1
                % check passes ct threshold
                %if (trig_t1 && ( ct(count) > ct_thresh ))
                if ((trig_t1 || trig_t2 || trig_t3) && ( ct(count) > 0 ))
                    
                    % Pick time sample to plot
                    Rc = R(count,:);
                    Risart = art_idx(count,:);
                    if (isnan(r_samp_const))
                        Ridx = 1:length(Rc);
                        Ridx = Ridx(randperm(length(Ridx)));
                        %r_idx = NaN;
                        
                        % Build list of just significant interactions
                        r_idx = [];
                        r_idx_cohval = [];
                        for i3 = 1:length(Ridx)
                            if ( (Rc(Ridx(i3)) > coh_thresh(count)) && (~Risart(Ridx(i3))) )
                                r_idx = [r_idx; Ridx(i3)];
                                r_idx_cohval = [r_idx_cohval; Rc(Ridx(i3))];
                                %break;
                            end
                        end
                        fprintf('%s\n',pname)
                        fprintf('\t[*] found %i significant interactions\n',length(r_idx));
                        
                        % sort by coherence threshold
                        [~,sidx] = sort(r_idx_cohval);
                        
                        % significance indexed based on strength
                        %r_idx = r_idx(sidx(end));
                        sidx2 = sidx( floor(strength*(length(sidx)-1)) +1 );
                        r_idx = r_idx( sidx2 );
                        fprintf('\t[*] most significant graph index: %i, coh: %.3f\n',r_idx,r_idx_cohval(sidx2))
                        
                        % empty case - hopefully won't happen
                        if (isempty(r_idx)) %if not found600
                            fprintf(2,'W> no time sample found, defaulting to first');
                            r_idx = 1;
                        end
                        r_samp = (r_idx-1)*round(ecog.fs)*w + 1;
                        
                        
%                         if (trig_samp_offset_override)
%                             samp_offset = round(ecog.fs)*(samp_offset_override);
%                             r_samp_n = r_samp + samp_offset;
%                             r_idx_n = (r_samp_n - 1)/(round(ecog.fs)*w) + 1;
%                         else
%                             % Pick time shifted sample to plot
%                             r_idx_n = round(length(Rc)/2);
%                             for i3 = 1:length(Ridx)
%                                 if ( (abs(Ridx(i3) - r_idx) > perm_alpha_sec*w) && (~Risart(Ridx(i3))) )
%                                     r_idx_n = Ridx(i3);
%                                     break;
%                                 end
%                             end
%                             r_samp_n = (r_idx_n-1)*round(ecog.fs)*w + 1;
%                             samp_offset = r_samp_n - r_samp;
%                         end
                    
                    else
                        % bypassround_y_nearest = 100;
                        r_samp = r_samp_const;
                    end
                    
                    r_samp_n = r_samp;
%                     sat = false;
%                     first_run = true;
%                     while(~sat)
%                         %rand_sign = 1 - randi([0,1])*2; % -1 or 1 randomly
%                         rand_sign = 1;
%                         if (trig_samp_offset_override && first_run)
%                             samp_offset = round(ecog.fs)*(samp_offset_override);
%                             first_run = false;
%                         else
%                             samp_offset = randi([perm_alpha_sec*round(ecog.fs),ecog.n_samples]) * rand_sign;
%                         end
%                         r_samp_n = r_samp + samp_offset;
%                         %r_samp_n_end = r_samp_n + w*round(ecog.fs) - 1;
%                         is_valid = (r_samp_n >= 1) && (r_samp_n <= (ecog.n_samples - round(ecog.fs)*w));
%                         is_art = true;
%                         if (is_valid)
%                             coeff = 1/(round(ecog.fs)*w);
%                             %is_art = any(Risart( ceil(r_samp_n*coeff):ceil(r_samp_n_end*coeff) ));
%                             is_art = Risart(ceil(r_samp_n*coeff));
%                         end
%                         sat = is_valid && (~ is_art);
%                     end

                    % convert r_idx, r_idx_n
                    r_idx = (r_samp - 1)/(round(ecog.fs)*w) + 1;
                    r_idx_n = (r_samp_n - 1)/(round(ecog.fs)*w) + 1;
                    
                    % Read timeseries
                    b1 = ecog.bip(ii1,1:2);
                    b2 = ecog.bip(ii2,1:2);
                    v_b1c1 = h5read(fn_h5Lpp,'/h5eeg/eeg',[b1(1) r_samp],[1 round(ecog.fs)*w]);
                    v_b1c2 = h5read(fn_h5Lpp,'/h5eeg/eeg',[b1(2) r_samp],[1 round(ecog.fs)*w]);
                    v_b2c1 = h5read(fn_h5Lpp,'/h5eeg/eeg',[b2(1) r_samp],[1 round(ecog.fs)*w]);
                    v_b2c2 = h5read(fn_h5Lpp,'/h5eeg/eeg',[b2(2) r_samp],[1 round(ecog.fs)*w]);
                    v0_b1c1 = h5read(fn_h5Lpp,'/h5eeg/eeg',[b1(1) r_samp_n],[1 round(ecog.fs)*w]);
                    v0_b1c2 = h5read(fn_h5Lpp,'/h5eeg/eeg',[b1(2) r_samp_n],[1 round(ecog.fs)*w]);
                    v0_b2c1 = h5read(fn_h5Lpp,'/h5eeg/eeg',[b2(1) r_samp_n],[1 round(ecog.fs)*w]);
                    v0_b2c2 = h5read(fn_h5Lpp,'/h5eeg/eeg',[b2(2) r_samp_n],[1 round(ecog.fs)*w]);
                    % Bipolar montage
                    v_b1 = v_b1c1 - v_b1c2;
                    v_b2 = v_b2c1 - v_b2c2;
                    v0_b1 = v0_b1c1 - v0_b1c2;
                    v0_b2 = v0_b2c1 - v0_b2c2;
                    % Subtract mean
                    v_b1 = v_b1 - mean(v_b1);
                    v_b2 = v_b2 - mean(v_b2);
                    v0_b1 = v0_b1 - mean(v_b1);
                    v0_b2 = v0_b2 - mean(v_b2);
                    
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
                    for i3 = 1:length(metricspp)
                        fn_graphLpp = sprintf('%s/%s_graph-%s.h5',dir_resLpp,sidpp,metricspp{i3});
                        coh_val = h5read(fn_graphLpp,'/R',[count r_idx],[1 1]);
                        coh_vals(i3) = coh_val;
                        mname = metricspp{i3}(3:4);
                        pname = sprintf('%s_%s%i',pname,mname,round(1000*coh_val));
                    end                   
                    pname = sprintf('%s_Samp%i',pname,r_samp);
                    for i3 = 1:length(metricspp)
                        mname = [metricspp{i3}(3:4),'N'];
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
                    
                    % Pick pair only if it is not significant
                    plot_check = false;
                    if (trig_t1 && (coh_vals(1) < coh_thresh(count)))
                        
                        % Check to see if it is interaction with bchan 2
                        if (found_fwd)
                            bchan1_n = bchan2_const;
                            bchan2_n = ii2;
                        elseif (found_bak)
                            bchan1_n = ii1;
                            bchan2_n = bchan2_const;
                        end
                        count2 = 1;
                        for i_t1 = 1:(ecog.n_bchan - 1)
                            for i_t2 = (i_t1 + 1):ecog.n_bchan
                                is_bpair = ((i_t1 == bchan1_n) && (i_t2 == bchan2_n));
                                is_bpair2 = ((i_t2 == bchan1_n) && (i_t1 == bchan2_n));
                                if (is_bpair || is_bpair2)
                                    
                                    % Finally, check if any of electrodes
                                    % are in Brocas or Wernickes
                                    %roi_1_rois = {'parsorbitalis','parstriangularis','parsopercularis'};
                                    %roi_2_chans = {'PT5','PT6','PT7'};
%                                     in_roi_b1 = any(strcmp(roi_1_rois,C.AtlLabels{atl}{i_t1}));
%                                     in_roi_b2 = any(strcmp(roi_1_rois,C.AtlLabels{atl}{i_t2}));
%                                     in_chan_b1c1 = any(strcmp(roi_2_chans,C.EleLabels{ecog.bip(i_t1,1)}));
%                                     in_chan_b1c2 = any(strcmp(roi_2_chans,C.EleLabels{ecog.bip(i_t1,2)}));
%                                     in_chan_b2c1 = any(strcmp(roi_2_chans,C.EleLabels{ecog.bip(i_t2,1)}));
%                                     in_chan_b2c2 = any(strcmp(roi_2_chans,C.EleLabels{ecog.bip(i_t2,2)}));
%                                     cond_region = (~in_roi_b1) && (~in_roi_b2) && (~in_chan_b1c1) && (~in_chan_b1c2) && (~in_chan_b2c1) && (~in_chan_b2c2);
                                    %cond_region = true;
                                    %if (cond_region || true)
                                    fn_graphLpp = sprintf('%s/%s_graph-%s.h5',dir_resLpp,sidpp,metricspp{metrici});
                                    coh_val = h5read(fn_graphLpp,'/R',[count2 r_idx],[1 1]);
                                    if (coh_val < coh_thresh(count2))
                                        fprintf('[debug] Found negative interaction: %i %i coh: %.3f, thresh: %.3f\n',...
                                           i_t1,i_t2,coh_val,coh_thresh(count2));
                                        %fprintf('r_idx: %i\n',r_idx)
                                        plot_check = true;
                                        n_not_coh_with_both = n_not_coh_with_both + 1;
                                    end
                                    %end                                    
                                end
                                count2 = count2 + 1;
                            end
                        end
                        
                        
                        if (plot_check)
                            
                            h = figure('visible','off');
                            set(h,'Position',[0 0 PLOT_W_PIX PLOT_H_PIX]);
                            subplot(4,2,1)
                            plot(t,v_b1,'black-','LineWidth',0.5); hold on;
                            %ylabel({ytxt;roi_1});
                            ylabel(ytxt);
                            xlabel(xtxt);
                            if (trig_title)
                                title(sprintf('Coherence: %.2f (\\delta:%.2f \\theta:%.2f \\alpha:%.2f \\beta:%.2f \\gamma:%.2f )',...
                                    coh_vals(1),coh_vals(2),coh_vals(3),coh_vals(4),coh_vals(5),coh_vals(6) ));
                            end
                            xticks(min(t):x_div_secs:max(t));
                            set(gca,'TickDir','out')
                            box off;
                            subplot(4,2,3)
                            plot(t,v_b2,'black-','LineWidth',0.5); hold on;
                            %ylabel({ytxt;roi_2});
                            ylabel(ytxt);
                            xlabel(xtxt);
                            xticks(min(t):x_div_secs:max(t));
                            set(gca,'TickDir','out')
                            box off;

                            % reverse labels
                            if (has_roi_r)
                                b2t = b2;
                                b2 = b1;
                                b1 = b2t;
                            end

    %                         hsub = subplot(4,2,2);
    %                         hb = brainplot_one(sidpp,[b1(1),b1(2)],'');
    %                         copyobj(hb,hsub);
    %                         if (trig_title)
    %                             title(sidpp);
    %                         end
    %                         hsub = subplot(4,2,4);
    %                         hb = brainplot_one(sidpp,[b2(1),b2(2)],'');
    %                         copyobj(hb,hsub);
    %                         if (trig_title)
    %                             title(sidpp);
    %                         end
    %                         % ------- Bottom right 2 --------------------------
    %                         hsub = subplot(4,2,6);
    % %                         hb = brainplot_one(sid,[b1(1),b1(2)],'');
    % %                         copyobj(hb,hsub);
    % %                         if (trig_title)
    % %                             title(sid);
    % %                         end
    %                         t2 = t + samp_offset/round(ecog.fs);
    %                         plot(t2,v0_b1,'black-','LineWidth',0.5,'Color',col_faint); hold on;
    %                         xticks(min(t2):x_div_secs:max(t2));
    %                         set(gca,'TickDir','out');
    %                         ylim(yl1);
    %                         yticks(yt1);
    %                         set(gca,'Ycolor','none');
    %                         box off;
    % 
    %                         hsub = subplot(4,2,8);
    % %                         hb = brainplot_one(sid,[b2(1),b2(2)],'');
    % %                         copyobj(hb,hsub);
    % %                         if (trig_title)
    % %                             title(sid);
    % %                         end
    %                         plot(t2,v0_b2,'black-','LineWidth',0.5); hold on;
    %                         xticks(min(t2):x_div_secs:max(t2));
    %                         set(gca,'TickDir','out');
    %                         ylim(yl2);
    %                         yticks(yt2);
    %                         set(gca,'Ycolor','none');
    %                         box off;



                            print(h,sprintf('figures/T1d3/%s',pname),fig_fmt);
                            if (trig_eps)
                                print(h,sprintf('figures/T1d3/%s',pname),'-depsc');
                            end
                            close(h);



                            % Plot brain separately so .eps can import easily
                            if (trig_brain)
                                h = figure('visible','off');
                                set(h,'Position',[0 0 PLOT_W_PIX PLOT_H_PIX]);

                                hsub = subplot(4,2,2);
                                if (trig_both_elec)
                                    hb = brainplot_one(sidpp,[b1(1),b1(2)],'');
                                else
                                    hb = brainplot_one(sidpp,[b1(1)],'');
                                end
                                copyobj(hb,hsub);
                                if (trig_title)
                                    title(sidpp);
                                end
                                hsub = subplot(4,2,4);
                                if (trig_both_elec)
                                    hb = brainplot_one(sidpp,[b2(1),b2(2)],'');
                                else
                                    hb = brainplot_one(sidpp,[b2(1)],'');
                                end
                                copyobj(hb,hsub);
                                if (trig_title)
                                    title(sidpp);
                                end
                                print(h,sprintf('figures/T1d3/%s_brain',pname),fig_fmt);
                                if (trig_eps)
                                    print(h,sprintf('figures/T1d3/%s_brain',pname),'-depsc');
                                end
                                close(h);
                            end
                        
                        end
                        
                        
%                         % --- Print Figure S04 ---------------------------
%                         
%                         h = figure;
%                         set(h,'Position',[0 0 PLOT_W_PIX PLOT_H_PIX]);
%                         
%                         subplot(4,2,5)
%                         plot(t,v_b1,'black-','LineWidth',0.5); hold on;
%                         %ylabel({ytxt;roi_1});
%                         ylabel(ytxt);
%                         xlabel(xtxt);
%                         xticks(min(t):x_div_secs:max(t));
%                         set(gca,'TickDir','out')
%                         if (trig_title)
%                             title(sprintf('Coherence: %.2f (\\delta:%.2f \\theta:%.2f \\alpha:%.2f \\beta:%.2f \\gamma:%.2f )',...
%                                 coh_val0(1),coh_val0(2),coh_val0(3),coh_val0(4),coh_val0(5),coh_val0(6) ));
%                         end
%                         ymax1 = ceil(max([v_b1,v0_b1])/round_y_nearest) * round_y_nearest;
%                         ymin1 = floor(min([v_b1,v0_b1])/round_y_nearest) * round_y_nearest;
%                         ylim([ymin1 ymax1]);
%                         yl1 = ylim;
%                         yt1 = yticks;
%                         box off;
%                         
%                         subplot(4,2,7)
%                         plot(t,v_b2,'black-','LineWidth',0.5,'Color',col_faint); hold on;
%                         %ylabel({ytxt;roi_2});
%                         ylabel(ytxt);
%                         xlabel(xtxt);
%                         xticks(min(t):x_div_secs:max(t));
%                         set(gca,'TickDir','out')
%                         ymax2 = ceil(max([v_b2,v0_b2])/round_y_nearest) * round_y_nearest;
%                         ymin2 = floor(min([v_b2,v0_b2])/round_y_nearest) * round_y_nearest;
%                         ylim([ymin2 ymax2]);
%                         yl2 = ylim;
%                         yt2 = yticks;
%                         box off;
%                         
%                         hsub = subplot(4,2,6);
% %                         hb = brainplot_one(sid,[b1(1),b1(2)],'');
% %                         copyobj(hb,hsub);
% %                         if (trig_title)
% %                             title(sid);
% %                         end
%                         t2 = t + samp_offset/round(ecog.fs);
%                         plot(t2,v0_b1,'black-','LineWidth',0.5,'Color',col_faint); hold on;
%                         xticks(min(t2):x_div_secs:max(t2));
%                         set(gca,'TickDir','out');
%                         ylim(yl1);
%                         yticks(yt1);
%                         set(gca,'Ycolor','none');
%                         box off;
% 
%                         hsub = subplot(4,2,8);
% %                         hb = brainplot_one(sid,[b2(1),b2(2)],'');
% %                         copyobj(hb,hsub);
% %                         if (trig_title)
% %                             title(sid);
% %                         end
%                         plot(t2,v0_b2,'black-','LineWidth',0.5); hold on;
%                         xticks(min(t2):x_div_secs:max(t2));
%                         set(gca,'TickDir','out');
%                         ylim(yl2);
%                         yticks(yt2);
%                         set(gca,'Ycolor','none');
%                         box off;
%                         
%                         print(h,sprintf('figures/T1d1/%s',pname),fig_fmt);
%                         if (trig_eps)
%                             print(h,sprintf('figures/T1d1/%s',pname),'-depsc');
%                         end
%                         close(h);
                        %return
                    end
                    
                    
                    if (trig_t2)
                        % backup
                        v_b1_b = v_b1;
                        v_b2_b = v_b2;
                        v0_b2_b = v0_b2;
                        forder = 4;
                        
                        for i4 = 1:length(metricspp)
                            % restore from backup
                            v_b1 = v_b1_b;
                            v_b2 = v_b2_b;
                            v0_b2 = v0_b2_b;
                            
                            switch metricspp{i4}
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
                            if (trig_title)
                                title(sprintf('Coherence: %.2f (\\delta:%.2f \\theta:%.2f \\alpha:%.2f \\beta:%.2f \\gamma:%.2f )',...
                                    coh_vals(1),coh_vals(2),coh_vals(3),coh_vals(4),coh_vals(5),coh_vals(6) ));
                            end
                            xticks(min(t):x_div_secs:max(t));
                            set(gca,'TickDir','out')
                            box off;
                            subplot(4,2,3)
                            plot(t,v_b2,'black-','LineWidth',0.5); hold on;
                            ylabel({ytxt;roi_2});
                            xlabel(xtxt);
                            xticks(min(t):x_div_secs:max(t));
                            set(gca,'TickDir','out')

                            box off;
%                             subplot(4,2,5)
%                             plot(t,v_b1,'black-','LineWidth',0.5); hold on;
%                             ylabel({ytxt;roi_1});
%                             xlabel(xtxt);
%                             xticks(min(t):x_div_secs:max(t));
%                             set(gca,'TickDir','out')
%                             if (trig_title)
%                                 title(sprintf('Coherence: %.2f (\\delta:%.2f \\theta:%.2f \\alpha:%.2f \\beta:%.2f \\gamma:%.2f )',...
%                                     coh_val0(1),coh_val0(2),coh_val0(3),coh_val0(4),coh_val0(5),coh_val0(6) ));
%                             end
%                             box off;
%                             subplot(4,2,7)
%                             plot(t,v0_b2,'black-','LineWidth',0.5); hold on;
%                             ylabel({ytxt;roi_2});
%                             xlabel(xtxt);
%                             xticks(min(t):x_div_secs:max(t));
%                             set(gca,'TickDir','out')
%                             box off;
%                             
                            % reverse labels
                            if (has_roi_r)
                                b2t = b2;
                                b2 = b1;
                                b1 = b2t;
                            end
% 
%                             hsub = subplot(4,2,2);
%                             hb = brainplot_one(sid,[b1(1),b1(2)],'');
%                             copyobj(hb,hsub);
%                             if (trig_title)
%                                 title(sid);
%                             end
%                             hsub = subplot(4,2,4);
%                             hb = brainplot_one(sid,[b2(1),b2(2)],'');
%                             copyobj(hb,hsub);
%                             if (trig_title)
%                                 title(sid);
%                             end

%                             hsub = subplot(4,2,6);
%                             hb = brainplot_one(sid,[b1(1),b1(2)],'');
%                             copyobj(hb,hsub);
%                             if (trig_title)
%                                 title(sid);
%                             end
%                             hsub = subplot(4,2,8);
%                             hb = brainplot_one(sid,[b2(1),b2(2)],'');
%                             copyobj(hb,hsub);
%                             if (trig_title)
%                                 title(sid);
%                             end

                            %return

                            print(h,sprintf('figures/T2/%s_FILT-%s',pname,metricspp{i4}(3:4)),fig_fmt);
                            if (trig_eps)
                                print(h,sprintf('figures/T2/%s_FILT-%s',pname,metricspp{i4}(3:4)),'-depsc');
                            end
                            close(h);
                        end
                        
                        %return
                    end
                    
                    
                    if (trig_t3)
                        

                        [cxx,f] = mscohere(v_b1', v_b2',hamming(Fs/1),[],[],Fs);fn_cacheLpp = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheLpp,sidpp,iMpp);
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
                        %axis tight;
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
                        %axis tight;
                        
                        %return
                        
                        print(h,sprintf('figures/T3/%s',pname),fig_fmt);
                        if (trig_eps)
                            print(h,sprintf('figures/T3/%s',pname),'-depsc');
                        end
                        close(h);
                    end
                    
                    
                    
                end
                
                
                if (trig_t4)
                    h = figure('visible','off');
                    set(h,'Position',[0 0 0.5*1080 0.5*1080])
                    
                    Rc = R(count,:);
                    Risart = art_idx(count,:);
                    
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
                    a1_b1 = h5read(fn_artLpp,'/artifacts',[ii1 r_samp_1s],[1 w]);
                    a1_b2 = h5read(fn_artLpp,'/artifacts',[ii2 r_samp_1s],[1 w]);
                    
                    % Read timeseries
                    b1 = ecog.bip(ii1,1:2);
                    b2 = ecog.bip(ii2,1:2);
                    v_b1c1 = h5read(fn_h5Lpp,'/h5eeg/eeg',[b1(1) r_samp],[1 round(ecog.fs)*w]);
                    v_b1c2 = h5read(fn_h5Lpp,'/h5eeg/eeg',[b1(2) r_samp],[1 round(ecog.fs)*w]);
                    v_b2c1 = h5read(fn_h5Lpp,'/h5eeg/eeg',[b2(1) r_samp],[1 round(ecog.fs)*w]);
                    v_b2c2 = h5read(fn_h5Lpp,'/h5eeg/eeg',[b2(2) r_samp],[1 round(ecog.fs)*w]);
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
                    
                    print(h,sprintf('figures/T4/%s_Samp%i',pname,r_samp),fig_fmt);
                    if (trig_eps)
                        print(h,sprintf('figures/T4/%s_Samp%i',pname,r_samp),'-depsc');
                    end
                    close(h);
                end
                
                
                if (trig_t6)
                    x = linspace(0,1,250);
                    h = figure('visible','off');
                    set(h,'Position',[0 0 1080 1080])
                    
                    for i5 = 1:length(metricspp)
                        
                        metricpp = metricspp{i5};
                        subplot(length(metricspp)/2,2,i5)
                        fn_artLpp = sprintf('%s/%s_art.h5',dir_artLpp,sidpp);
                        fn_distLpp = sprintf('%s/%s_dists-%s-%i.mat',dir_resLpp,sidpp,metricpp,n_perm);
                        fn_graphLpp = sprintf('%s/%s_graph-%s.h5',dir_resLpp,sidpp,metricpp);
                        %fn_permL = sprintf('/mnt/cuenap2/data/results/coh_w10/%s_perm-%s-%i.h5',sid,metric,n_perm);
                        fn_permLpp = sprintf('%s/%s_perm-%s-%i.h5',dir_resLpp,sidpp,metricpp,n_perm);
                        
                        load(fn_distLpp);
                        Distr = d{count};
                        Rc = h5read(fn_graphLpp,'/R',[count 1],[1 n_graph]);
                        Risart = h5read(fn_artLpp,'/art_idx',[count 1],[1 n_graph]) ~= 0;
                        Rc(Risart) = NaN;
                        Rperm = h5read(fn_permLpp,'/R',[count 1],[1 n_perm]);

                        [counts,~] = hist(Rc,x);
                        plot(x,counts/trapz(x,counts),'-','color',[1 1 1]*0); hold on;
                        [countsP,~] = hist(Rperm,x);
                        plot(x,countsP/trapz(x,countsP),'-','color',[1 0 0]*0.3); hold on;
                        plot(x,pdf(d{count},x),'--','color',[1 0 0]*0.3);
                        box off;
                        xlabel(sprintf('Coherence'));
                        ylabel({'pdf';sprintf('%s',metricpp(3:end))});
                        
                    
                    end
                    
                    %return
                    
                    print(h,sprintf('figures/T6/%s',pname),fig_fmt);
                    if (trig_eps)
                        print(h,sprintf('figures/T6/%s',pname),'-depsc');
                    end
                    close(h);
                    
                end
                
                N_9d1 = [];
                if (trig_t9d1)
                    N_9d1 = [N_9d1 2];
                end
                if (trig_t9)
                    N_9d1 = [N_9d1 1];
                end
                for i_9d1 = N_9d1
                    for iMp = 1:length(metricspp)

                        metricpp = metricspp{iMp};
                        
                        h = figure('visible','off');
                        set(h,'Position',[0 0 0.5*1080 0.5*1080])
                        
                        if (i_9d1 == 2)
                            t_plot_hrs = n_graph * w / 3600;
                            Tmin_fac = (1/3600);
                        elseif (i_9d1 == 1)
                            t_plot_hrs = 1;
                            Tmin_fac = 1;
                        end

                        % Reload coh_thresh
                        fn_cacheLp = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheLpp,sidpp,iMp);
                        load(fn_cacheLp,'coh_thresh');
                        
                        % Load artifacts
                        Risart = (art_idx(count,:) ~= 0);
                        
                        % Load graph
                        fn_graphLpp = sprintf('%s/%s_graph-%s.h5',dir_resLpp,sidpp,metricpp);
                        Rc = h5read(fn_graphLpp,'/R',[count 1],[1 n_graph]);
                        
                        % Sample null
                        fn_distLpp = sprintf('%s/%s_dists-%s-%i.mat',dir_resLpp,sidpp,metricpp,n_perm);
                        D = load(fn_distLpp);
                        Rn =random(D.d{count},[1 n_graph]);
%                         
%                         fn_cache2L = sprintf('%s/%s_R_null.mat',dir_cacheL,sid);
%                         N = load(fn_cache2L);
%                         Rn = N.R_null(count,:);

                        cond_plot = 1:(t_plot_hrs*3600 / w);
                        T_min = linspace(0,max(cond_plot)*w/(60),length(cond_plot));
                        T_min = Tmin_fac * T_min;
                        
                        % Subsample to plot with least artifacts and high coherence
                        if (exist('r_samp','var') && exist('r_samp_n','var'))
                            % If figure t1 was plotted, choose sample containing time segment used for that
                            if (i_9d1 == 1)
                                t_plot_off =  round((t_plot_hrs*3600 / w)/2);
                                cond_offset = ( floor((r_samp - 1) / (w * round(ecog.fs))) ) - t_plot_off;
                                offset_fit = (cond_offset >= 0) & ((length(cond_plot) + cond_offset) < n_graph);
                                if (offset_fit)
                                    cond_plot = cond_plot + cond_offset;
                                end
                            end
                            % assert: r_idx == ( floor((r_samp - 1) / (w * round(ecog.fs))) + 1 )
                            cond_plot0 = cond_plot; % duplicate time range to plot for null
                        else
                            art_frac = zeros(n_graph - length(cond_plot) + 1,1);
                            for i3 = 1:(n_graph - length(cond_plot) + 1) 
                                r_mean = mean(Rc(cond_plot + i3 - 1));
                                seg_isa = Risart(cond_plot + i3 - 1);
                                art_frac(i3) = sum(seg_isa)/length(seg_isa) - 1*r_mean;
                            end
                            [~,min_idx] = min(art_frac);

                            if (~isempty(min_idx))
                                cond_plot = cond_plot + min_idx - 1;
                            end
                            cond_plot0 = cond_plot; % duplicate time range to plot for null
                        end
                        
                        R_plt = Rc(cond_plot);
                        Rn_plt = Rn(cond_plot0);
                        Risa_plt = Risart(cond_plot);
                        Rnisa_plt = Risart(cond_plot0);

                        % Remove artifacts
                        R_plt(Risa_plt) = NaN;
                        Rn_plt(Rnisa_plt) = NaN;
                        f_art = sum(Risa_plt)/numel(Risa_plt);

                        color_div = 0.33*[1 1 1];
                        subplot(2,1,1)
                        plot(T_min,R_plt,'black.','MarkerSize',5);hold on;
                        plot([min(T_min) max(T_min)],[coh_thresh(count) coh_thresh(count)],'--','color',color_div)
                        if (i_9d1 == 2)
                            plot(T_min(r_idx),R_plt(r_idx),'redo','MarkerSize',5); hold on;
                        else
                            if (offset_fit)
                                plot(T_min(t_plot_off+1),R_plt(t_plot_off+1),'redo','MarkerSize',5); hold on;
                            end
                        end
                        
                        axis([min(T_min) max(T_min) 0 1]); 
                        xlabel('Time (minutes)')
                        ylabel({sprintf('%s Coherence',metricpp(3:end));sprintf('%s-%s',roi_1,roi_2)})
                        box off;

                        subplot(2,1,2)
                        plot(T_min,Rn_plt,'black.','MarkerSize',5); hold on;
                        plot([min(T_min) max(T_min)],[coh_thresh(count) coh_thresh(count)],'--','color',color_div)
%                         if (offset_fit)
%                             plot(T_min(cond_offset+1),Rn_plt(cond_offset+1),'blacko','MarkerSize',10);
%                         end
                        
                        axis([min(T_min) max(T_min) 0 1])
                        xlabel('Time (minutes)')
                        ylabel({sprintf('Time-shifted %s Coherence',metricpp(3:end));sprintf('%s-%s',roi_1,roi_2)})
                        box off;

                        %return

                        ct_plt = sum(R_plt(~Risa_plt) > coh_thresh(count)) / length(R_plt(~Risa_plt));
                        ct0_plt = sum(Rn_plt(~Rnisa_plt) > coh_thresh(count)) / length(Rn_plt(~Rnisa_plt));
                        
                        if (i_9d1 == 1)
                            dir_t9 = 'T9';
                        elseif (i_9d1 == 2)
                            dir_t9 = 'T9d1';
                        end
                        print(h,sprintf('figures/%s/%s_%s_coht%i_ct%i_ctN%i_fart%i',...
                            dir_t9,pname,metricpp(3:4),round(1000*coh_thresh(count)),round(1000*ct_plt),round(1000*ct0_plt),round(1000*f_art)),fig_fmt);
                        if (trig_eps)
                            print(h,sprintf('figures/%s/%s_%s_coht%i_ct%i_ctN%i_fart%i',...
                                dir_t9,pname,metricpp(3:4),round(1000*coh_thresh(count)),round(1000*ct_plt),round(1000*ct0_plt),round(1000*f_art)),'-depsc');
                        end
                        close(h);
                    end
                end
                
               
                
            end
            count = count + 1;
        end
    end
    
end
end
end
% 
% fprintf('Total number of pairs passing distance threshold: %i\n',n_pass_dist_thresh);
% fprintf('Total number of pairs cohered to bchan1 (%i): %i\n',bchan1_const,n_corr_bchan1)
% fprintf('Total number of pairs cohered to bchan2 (%i): %i\n',bchan2_const,n_corr_bchan2)
% fprintf('Total number of pairs not cohered to either bip %i or %i: %i\n',bchan1_const,bchan2_const,n_not_coh_with_both);
% fprintf('Fraction w.r.t. bchan1 pairs: %.6f\n',n_not_coh_with_both/n_corr_bchan1);
% fprintf('Fraction w.r.t. all pairs: %.6f\n',n_not_coh_with_both/n_pass_dist_thresh);
% fprintf('N_bchan: %i\n',ecog.n_bchan)

% Print coherence with other electrode
% bchan2_const = [44, 55, 73]; % NaN;
% bchan2_2 = 84;
% c = 1;
% for ia = 1:(ecog.n_bchan-1)
%     for ib = (ia + 1):ecog.n_bchan
%         has_a = any(ia == bchan2_const) || (ia == bchan2_2);
%         has_b = any(ib == bchan2_const) || (ib == bchan2_2);
%         if (has_a && has_b)
%             coh_val = h5read(fn_graphLpp,'/R',[c r_idx],[1 1]);
%             fprintf('%i_%i_Br%i_coht%i\n',ia,ib,round(1000*coh_val),round(1000*coh_thresh(c)))
%         end
%         c = c + 1;
%     end
% end

%return

% count numbers
fn_graphLpp = sprintf('%s/%s_graph-%s.h5',dir_resLpp,sidpp,metricspp{metrici});
n_pass_dist = 0;
n_pairs = 0;
n_pd_pairs_bchan1 = 0;
n_pd_pairs_bchan2 = 0;
n_pd_pairs_pos_bchan1 = 0;
n_pd_pairs_pos_bchan2 = 0;
n_pd_pairs_neg_bchan1 = 0;
n_pd_pairs_neg_bchan2 = 0;
for ii1 = 1:(ecog.n_bchan-1)
    for ii2 = (ii1+1):ecog.n_bchan
        count3 = n_pairs + 1;
        if (Dmat(ii1,ii2) > dist_thresh)
            n_pass_dist = n_pass_dist + 1;
            
            has_b1 = (ii1 == bchan1_const) || (ii2 == bchan1_const);
            has_b2 = (ii1 == bchan2_const) || (ii2 == bchan2_const);
            
            if (has_b1 && (~has_b2))
                n_pd_pairs_bchan1 = n_pd_pairs_bchan1 + 1;
                coh_val = h5read(fn_graphLpp,'/R',[count3 r_idx],[1 1]);
                if (coh_val < coh_thresh(count3))
                    n_pd_pairs_pos_bchan1 = n_pd_pairs_pos_bchan1 + 1;
                else
                    n_pd_pairs_neg_bchan1 = n_pd_pairs_neg_bchan1 + 1;
                end
            elseif (has_b2 && (~has_b1))
                n_pd_pairs_bchan2 = n_pd_pairs_bchan2 + 1;
                coh_val = h5read(fn_graphLpp,'/R',[count3 r_idx],[1 1]);
                if (coh_val < coh_thresh(count3))
                    n_pd_pairs_pos_bchan2 = n_pd_pairs_pos_bchan2 + 1;
                else
                    n_pd_pairs_neg_bchan2 = n_pd_pairs_neg_bchan2 + 1;
                end
            end
        end
        n_pairs = n_pairs + 1;
    end
end

% count number of electrodes not significant to either 1 or 2
n_pd_pairs_neg_both = 0;
n_pd_pairs_pos_both = 0;
n_pd_pairs_pos_b1 = 0;
n_pd_pairs_pos_b2 = 0;
n_pd = 0;
for iii = 1:ecog.n_bchan
    not_1_2 = (iii ~= bchan1_const) && (iii ~= bchan2_const);
    isdist = (Dmat(iii,bchan1_const) > dist_thresh) && (Dmat(iii,bchan2_const) > dist_thresh);
    
    if (not_1_2 && isdist)
        n_pd = n_pd + 1;
        % get coherence between iii and bchan1
        ilower = iii;
        ihigher = bchan1_const;
        if (ilower > ihigher)
            ilower_sav = ilower;
            ilower = ihigher;
            ihigher = ilower_sav;
        end
        count_i = find( ((chan1+1) == ilower) & ((chan2+1) == ihigher) );
        coh_val = h5read(fn_graphLpp,'/R',[count_i r_idx],[1 1]);
        neg_iii_bchan1 = coh_val < coh_thresh(count_i);

        % get coherence between iii and bchan2
        ilower = iii;
        ihigher = bchan2_const;
        if (ilower > ihigher)
            ilower_sav = ilower;
            ilower = ihigher;
            ihigher = ilower_sav;
        end
        count_i = find( ((chan1+1) == ilower) & ((chan2+1) == ihigher) );
        coh_val = h5read(fn_graphLpp,'/R',[count_i r_idx],[1 1]);
        neg_iii_bchan2 = coh_val < coh_thresh(count_i);

        if (neg_iii_bchan1 && neg_iii_bchan2)
            n_pd_pairs_neg_both = n_pd_pairs_neg_both + 1;
        elseif ((~neg_iii_bchan1) && (neg_iii_bchan2))
            n_pd_pairs_pos_b1 = n_pd_pairs_pos_b1 + 1;
        elseif ((~neg_iii_bchan2) && (neg_iii_bchan1))
            n_pd_pairs_pos_b2 = n_pd_pairs_pos_b2 + 1;
        elseif ((~neg_iii_bchan1) && (~neg_iii_bchan2))
            n_pd_pairs_pos_both = n_pd_pairs_pos_both + 1;
        end
    end
end

fprintf('bchan1: %i\n',bchan1_const)
fprintf('bchan2: %i\n',bchan2_const)
fprintf('n_bchan: %i\n',ecog.n_bchan)
fprintf('n_pairs: %i\n',n_pairs)
fprintf('n_pass_dist: %i\n',n_pass_dist)
fprintf('n_pd_pairs_bchan1: %i\n',n_pd_pairs_bchan1)
fprintf('n_pd_pairs_bchan2: %i\n',n_pd_pairs_bchan2)
fprintf('n_pd_pairs_pos_bchan1: %i\n',n_pd_pairs_pos_bchan1)
fprintf('n_pd_pairs_pos_bchan2: %i\n',n_pd_pairs_pos_bchan2)

fprintf('n_pd_pairs_neg_both: %i\n',n_pd_pairs_neg_both);
fprintf('n_pd_pairs_pos_both: %i\n',n_pd_pairs_pos_both);
fprintf('n_pd_pairs_pos_b1: %i\n',n_pd_pairs_pos_b1);
fprintf('n_pd_pairs_pos_b2: %i\n',n_pd_pairs_pos_b2);
fprintf('n_pd: %i\n',n_pd);


% for ii1 = 1:(ecog.n_bchan-1)
%     for ii2 = (ii1+1):ecog.n_bchan
%         if (Dmat(ii1,ii2) > dist_thresh)
%             n_pass_dist_thresh = n_pass_dist_thresh + 1;
%             if ((ii1 == bchan2_const) || (ii2 == bchan2_const))
%                 n_corr_bchan2 = n_corr_bchan2 + 1;
%                 % count number of positives
%                 fn_graphLpp = sprintf('%s/%s_graph-%s.h5',dir_resLpp,sidpp,metricspp{metrici});
%                 ctmp = 1;
%                 n_corr_bchan2_pos = 1;
%                 n_corr_bchan2_neg = 1;
%                 for it1 = 1:(ecog.n_bchan-1)
%                     for it2 = (it1 + 1):ecog.n_bchan
%                         coh_val = h5read(fn_graphLpp,'/R',[ctmp r_idx],[1 1]);
% 
%                         if ( (Dmat(it1,it2) > dist_thresh) && (coh_val > coh_thresh(ctmp)) )
%                             n_corr_bchan2_pos = n_corr_bchan2_pos + 1;
%                         elseif ( (Dmat(it1,it2) > dist_thresh) && (coh_val > coh_thresh(ctmp)) )
%                             n_corr_bchan2_neg = n_corr_bchan2_neg + 1;
%                         end
%                         ctmp = ctmp + 1;
%                     end
%                 end
%             end
%             if ((ii1 == bchan1_const) || (ii2 == bchan1_const))
%                 n_corr_bchan1 = n_corr_bchan1 + 1;
%                 % count number of positives
%                 fn_graphLpp = sprintf('%s/%s_graph-%s.h5',dir_resLpp,sidpp,metricspp{metrici});
%                 ctmp = 1;
%                 n_corr_bchan1_pos = 1;
%                 n_corr_bchan1_neg = 1;
%                 for it1 = 1:(ecog.n_bchan-1)
%                     for it2 = (it1 + 1):ecog.n_bchan
%                         coh_val = h5read(fn_graphLpp,'/R',[ctmp r_idx],[1 1]);
% 
%                         if ( (Dmat(it1,it2) > dist_thresh) && (coh_val > coh_thresh(ctmp)) )
%                             n_corr_bchan1_pos = n_corr_bchan1_pos + 1;
%                         elseif ( (Dmat(it1,it2) > dist_thresh) && (coh_val > coh_thresh(ctmp)) )
%                             n_corr_bchan1_neg = n_corr_bchan1_neg + 1;
%                         end
%                         ctmp = ctmp + 1;
%                     end
%                 end
%             end
% 
%         end
%     end
% end

% Clear loop indices
clear i;
clear j;

% Print finish message
fprintf('Done.\n')
