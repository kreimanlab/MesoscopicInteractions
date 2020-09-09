% this script is only for finding the optimal example to show.

close all;
clear;
rng shuffle;

metricp = 'pcBroadband';
n_perm = 10000;
perm_alpha_sec = 12;
cp_thresh_override = 0.05;
p_val = 0.01;

% Fast i/o definitions
dir_artLp = '/media/klab/KLAB101/h5_notch20/art_nosz';
dir_resLp = '/media/klab/KLAB101/results/coh_w10'; %dir_resLp = '/media/klab/internal/data/results/coh_w10';
dir_corLp = '/media/klab/internal/data/coreg';
dir_cacheLp = './cache';
dir_stamp = '/media/klab/internal/data/stamps';
dir_video = '/media/klab/internal/data/videos';
%subjects_dirLp = '/mnt/cuenap_ssd/coregistration';

% Slow i/o definitions
dir_h5Lp = '/media/klab/KLAB101/h5_notch20';
%dir_h5Lp = '/mnt/cuenap/data/h5_notch20';

%metricsp = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};
metricsp = {'pcBroadband'}; % debug

% Patients
Subjectsp = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48',...
    'mSu'};

% Exclude monkey
Subjectsp = Subjectsp(1:(end-1));

% bypass
%Subjectsp = {'sub41'};

%r_samp_const = 106865001; % set to nan to pick randomly
%r_samp_const = 90545001; % strong t1, t9 good but next to art
r_samp_const = NaN; % 22755001; % good t9, magnitude of t1 is only 0.339 (DEFAULT)
%r_samp_const = 22547501; % bad t1, but strong t9 1 hr
%r_samp_const = 90835001; % good t9, but not as good t1, high amplitude
%r_samp_const = 91430001; % t9 missing right half
%bchan1_const = 13;
%bchan2_const = 58;

bchan1_const = NaN; % 27
bchan2_const = NaN; % 46

%NaN
%r_samp_const = NaN;

% Pick areas to plot
% C.AtlROIs{2}.LH.struct_names
% (PFC, area 9/46) and the anterior cingulate cortex (ACC, areas 24c and 32).

%roi_1 = 'parstriangularis';
roi_1 = 'parsopercularis';
roi_2 = 'superiortemporal';
%roi_2 = 'middletemporal';
%roi_2 = 'inferiorparietal';

% --- parsopercularis - middletemporal
% sub1 50 l/r 
% sub2 66 l/r
% sub5 6  l/r vid
% sub10 5  r/l 
% sub13 2  l/r
% sub15 3  l/r vid
% sub17 1  l/r vid
% sub34 24 l/r vid
% sub37 2  l/r vid
% sub38 1  l/r vid
% sub41 30 l/r vid
% sub42 16 l/r vid

n_opti = 1;

% -- new calculations ----------------------------------------------------
% parsopercularis, inferiorparietal
%Subjectsp = {'sub34'}; % dominant hemi and video
%Subjectsp = {'sub1','sub3','sub24','sub34'};
%[*] cost: -0.676649978	sub34	27	54	79115001
%[*] cost: -0.649445107	sub34	27	54	80445001
%[*] cost: 0.093009810	sub34	13	54	74092501
%[*] cost: -0.354445984	sub34	27	54	74210001

% parsopercularis, superiortemporal
%Subjectsp = {'sub1','sub2','sub3','sub5','sub7','sub8','sub10','sub11','sub14','sub15','sub23','sub34','sub37','sub38','sub39','sub41','sub42'};
% 3n,19,25n,30,32,61,73,75,84,95
%Subjectsp = {'sub5','sub14','sub15','sub34','sub37','sub38','sub41','sub42'}; % dominant hemi and video
%Subjectsp = {'sub15'};
%[*] cost: 0.564330871	sub42	20	53	23932501
%[*] cost: 0.539396019	sub42	20	33	28727501
%[*] cost: 0.572074498	sub42	20	47	47595001
%[*] cost: 0.546873493	sub42	20	41	54105001
%[*] cost: 0.560554637	sub42	20	41	61880001
%[*] cost: 0.584691258	sub41	27	61	96155001
%[*] cost: 0.497183125	sub41	20	36	71732501
%[*] cost: 0.288932666	sub41	20	54	64352501
%[*] cost: 0.597722641	sub41	20	41	72342501
%[*] cost: 0.814054461	sub41	22	46	78242501
%[*] cost: 0.204371957	sub38	14	98	29332501
%[*] cost: 0.224854736	sub38	14	98	32395001
%[*] cost: 0.236381678	sub38	14	98	38712501
%[*] cost: 0.261469712	sub37	7	21	3327501
%[*] cost: 0.318063100	sub37	7	21	5400001
%[*] cost: 0.314823946	sub37	7	21	8427501
%[*] cost: -0.155321119	sub34	13	39	74105001
%[*] cost: 0.134045247	sub34	27	36	76422501
%[*] cost: -0.320159630	sub34	21	36	80770001
%[*] cost: 0.837716711	sub15	55	60	23692501
%[*] cost: 0.761445999	sub15	28	60	51460001
%[*] cost: 0.849530846	sub15	28	60	66437501
%[*] cost: 0.835929138	sub15	55	60	124770001
%[*] cost: 0.239254323	sub14	29	31	31370001
%[*] cost: 0.206529823	sub14	29	31	66517501
%[*] cost: 0.236791518	sub14	29	31	73600001
%[*] cost: -0.139612804	sub5	14	28	55607501
%[*] cost: 0.215532971	sub5	5	14	57387501

% parstriangularis, superiortemporal
%Subjectsp = {'sub1','sub2','sub3','sub5','sub6','sub7','sub8','sub10','sub11','sub15','sub17','sub23','sub28','sub41'};
% 3n,19,21n,25n,32,35,53,84
%Subjectsp = {'sub5','sub15','sub17','sub28','sub41'};% dominant hemi and video
%Subjectsp = {'sub5'};
%[*] cost: 0.463727641	sub41	3	61	66562501
%[*] cost: 0.507391921	sub41	20	41	77422501
%[*] cost: 0.275485745	sub41	4	42	92715001
%[*] cost: 0.587059769	sub28	2	40	104900001
%[*] cost: 0.586993876	sub28	2	40	112397501
%[*] cost: 0.576190564	sub28	2	40	120232501
%[*] cost: 0.523224774	sub17	77	88	522501
%[*] cost: 0.528096888	sub17	77	88	15890001
%[*] cost: 0.789399004	sub17	77	90	8185001
%[*] cost: 0.793726301	sub17	77	90	19687501
%[*] cost: 0.821242760	sub15	26	61	30217501
%[*] cost: 0.531650747	sub15	26	61	96310001
%[*] cost: 0.425622409	sub15	26	60	103780001
%[*] cost: 0.786905825	sub15	28	60	120367501
%[*] cost: -0.136224830	sub5	28	56	1700001
%[*] cost: -0.325221711	sub5	28	56	5937501
%[*] cost: 0.305643355	sub5	41	56	73312501
%[*] cost: -0.064330439	sub5	28	54	12265001
%[*] cost: 0.371564081	sub5	7	20	60667501

% Triggers
system('mkdir figures');
fig_fmt = '-depsc';
trig_title = false; % show titles
trig_eps = false;
trig_t1 = true;
trig_t2 = false; %
trig_t3 = false; %
trig_t4 = false; %default: false
trig_t6 = false; %
trig_t9 = true; %
trig_t9d1 = true; %
t9d1_offset_minutes = 25; % if 0, t9d1 centers the chosen time point
trig_save_fig = false; % used for debug


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
    system('mkdir figures/T1');
    system('mkdir figures/T1d1');
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

% not optimized:
% t1d3 - separate function
% t1d2 - separate function
% t2 - filtered
% t3 - coherence
% t4 - artifacts
% t6 - distributions

% optimization constraints
% t1
%   - broadband coherence must be significant
%   - no artifacts
%   - not in the very beginning or end of data
%
% t9
%   - CT > 0.1
%   - patient with video


% optimization variables
v_sub = ones(n_opti,1);
v_bchans = zeros(n_opti,2);
v_samp = zeros(n_opti,1);
c_bro_t1 = zeros(n_opti,1);
c_brot_t1 = zeros(n_opti,1);
c_delta_t1 = zeros(n_opti,1);
% c_bro_t1d2 = zeros(n_opti,1);
% c_delta_t1d2 = zeros(n_opti,1);
c_days_t9 = zeros(n_opti,1);
c_mag_t9 = zeros(n_opti,1);
c_magt_t9 = zeros(n_opti,1);
c_ct_t9 = zeros(n_opti,1);
c_art_t9 = zeros(n_opti,1);

%[*] cost: -0.043560	sub1	21	58	83470001
%[*] cost: 0.011246898	sub1	14	51	82335001
%[*] cost: -0.065562469	sub1	4	57	83825001
%[*] cost: -0.164624247	sub45	30	80	14870001
%[*] cost: -0.165981444	sub45	30	80	10840001
%[*] cost: -0.170712146	sub45	30	80	2245001
%[*] cost: -0.174177681	sub45	30	79	3407501
%[*] cost: -0.213723722	sub45	30	80	7315001
%[*] cost: 0.660226596	sub41	20	66	89857501
%[*] cost: 0.627504541	sub41	20	66	83155001
%[*] cost: 0.613476766	sub41	20	66	60247501
%[*] cost: 0.136753701	sub24	26	54	64410001
%[*] cost: 0.269184195	sub24	26	54	86635001

% untested:
%[*] cost: -26116.387557812	sub24	26	54	64875001 - line 64873354 - 02/07/2010 04:34:46
%[*] cost: -17351.387557812	sub24	26	54	84207501 - line 84205854 - 02/07/2010 15:19:11
%[*] cost: 0.232203059	sub24	26	54	86870001 - line 86868354 - 02/07/2010 16:47:56
%[*] cost: 0.241688453	sub24	26	54	86840001 - line 86838354 - 02/07/2010 16:46:56
%[*] cost: 0.265921794	sub24	26	54	86665001 - line 86663354 - 02/07/2010 16:41:06

% untested 'parahippocampal' <-> 'precentral'
%[*] cost: -0.523349863	sub34	27	62	80437501
%

%Subjectsp = {'sub1','sub2','sub15','sub36','sub41','sub46'};
%Subjectsp = {'sub45','sub41','sub24'}; % sub1
%Subjectsp = {'sub24'};

% parsopercularis
% Subjectsp = {'sub1','sub2','sub3','sub7','sub10','sub11','sub13',...
%     'sub18','sub20','sub21','sub24','sub32','sub34','sub35','sub39',...
% 	'sub41','sub45','sub46'};

% 'parahippocampal' <-> 'precentral';
%Subjectsp = {'sub1','sub5','sub34','sub42','sub46'};


for iO = 1:n_opti

    i_s = randi([1 length(Subjectsp)]);
    iMp = 1;
    
    % --------------------------------------------------------------------
    for i_s = i_s
    %for i_s = 1:length(Subjectsp)
    %     
    % % Main loop
    % for iMp = 1:1%length(metricsp)
    %for iM = 1
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
        for j = 1:length(ckf)
            if (~exist(ckf{j},'file'))
                fprintf(2,'W> File not found: %s\n',ckf{j});
                return
            end
        end

        %ecog = H5eeg(fn_h5);
        load(fn_cacheLp);

        % --- Apply threshold override --------------------------------------------
        cp_thresh = cp_thresh_override;
        %return

        chan_labels = h5readatt(fn_h5Lp,'/h5eeg/eeg','labels');
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
        art_idx = h5read(fn_artLp,'/art_idx');
        art_idxB = (art_idx == 1);
        frac_art = h5readatt(fn_artLp,'/art_idx','frac_art');
        w_art = h5readatt(fn_artLp,'/art_idx','w');
        % trim last time sample of graph
        RM = h5read(fn_graphLp,'/R',[1 1],size(art_idx));
        fn_graphLp_delta = sprintf('%s/%s_graph-%s.h5',dir_resLp,sidp,metricp);
        RM_delta = h5read(fn_graphLp_delta,'/R',[1 1],size(art_idx));
        w = double(h5read(fn_graphLp,'/w'));
        [n_comb,n_graph] = size(RM);

        % --- optimize I/O for t9 ---
        % Load graph for t9
    %     fn_graphLp = sprintf('%s/%s_graph-%s.h5',dir_resLp,sidp,metricp);
    %     Rc_t9 = h5read(fn_graphLp,'/R',[countp 1],[1 n_graph]);
        % Sample null
        fn_distLp = sprintf('%s/%s_dists-%s-%i.mat',dir_resLp,sidp,metricp,n_perm);
        D_t9 = load(fn_distLp);

        % Build list of all significant bipolar pairs
        sig_counts = [];
        sig_bchans = [];
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
                    cond_bpair_bypass = (ii1 == bchan1_const) && (ii2 == bchan2_const);
                else
                    % if either bchan1_const or bchan2_const is nan, do all
                    cond_bpair_bypass = true;
                end
                if ((((has_roi_f || has_roi_r) && (Dmat(ii1,ii2) > dist_thresh)) )  && (ct(countp) > ct_thresh) ) % cond_bpair_bypass = true;
                    sig_counts = [sig_counts; countp];
                    sig_bchans = [sig_bchans; [ii1,ii2]];
                end
                %fprintf('[%i] %i %i\n',ecog.n_bchan,ii1,ii2);
                countp = countp + 1;
            end
        end

        fprintf('[%i/%i] %s found %i significant bipolar pairs\n',iO,n_opti,sidp,length(sig_counts))

    end
    %return
    
    % randomly choose among significant
    if (isempty(sig_counts))
        sig_counts = [NaN];
        sig_bchans = [NaN];
    end
    randidx = randi([1 length(sig_counts)]);
    countp = sig_counts(randidx);
    ii1O = sig_bchans(randidx,1);
    ii2O = sig_bchans(randidx,2);
    
    

    
    % Find suitable bipolar pair
%     countp = 1;
%     for ii1 = 1:(ecog.n_bchan-1)
%         for ii2 = (ii1+1):ecog.n_bchan
    for ii1 = ii1O
        for ii2 = ii2O
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
                cond_bpair_bypass = (ii1 == bchan1_const) && (ii2 == bchan2_const);
            else
                % if either bchan1_const or bchan2_const is nan, do all
                cond_bpair_bypass = true;
            end
            if (((has_roi_f || has_roi_r) && (Dmat(ii1,ii2) > dist_thresh)) && cond_bpair_bypass)
                % ct threshold based on permutations
                p_val_ct_2 = 0.01/n_comb;
                ct_thresh_perm = (binoinv([p_val_ct_2 (1-p_val_ct_2)],n_graph,p_val/bonf_seconds)/n_graph);
                        
                pname = sprintf('%s__%i_%s_%s__%i_%s_%s__%imm_ct%i_ctt%i-%i_mag%i_coh%it%i',...
                    sidp,ii1,atl_labels{ecog.bip(ii1,1)},bchan_labels{ii1},ii2,atl_labels{ecog.bip(ii2,1)},...
                    bchan_labels{ii2},round(Dmat(ii1,ii2)),round(1000*ct(countp)),...
                    ceil(1000*ct_thresh_perm(1)),ceil(1000*ct_thresh_perm(2)),round(1000*mag(countp)),...
                    iMp,round(1000*coh_thresh(countp)));
                %fprintf('%s\n',pname);
                
                % Plot T1
                % check passes ct threshold
                %if (trig_t1 && ( ct(count) > ct_thresh ))
                if ((trig_t1 || trig_t2 || trig_t3) && ( ct(countp) > 0 ))
                    
                    % Pick time sample to plot
                    Rc = RM(countp,:);
                    Risart = art_idx(countp,:);
                    if (isnan(r_samp_const))
                        Ridx = 1:length(Rc);
                        Ridx = Ridx(randperm(length(Ridx)));
                        r_idx = NaN;
                        for i3 = 1:length(Ridx)
                            if ( (Rc(Ridx(i3)) > coh_thresh(countp)) && (~Risart(Ridx(i3))) )
                                r_idx = Ridx(i3);
                                break;
                            end
                        end
                        if (isnan(r_idx)) %if not found600
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
                    
                    
                    
                    % save optimization variables
                    v_sub(iO) = i_s;
                    v_bchans(iO,:) = sig_bchans(randidx,:);
                    v_samp(iO) = r_samp;
                    
                    
                    sat = false;
                    first_run = true;
                    while(~sat)
                        %rand_sign = 1 - randi([0,1])*2; % -1 or 1 randomly
                        rand_sign = 1;
                        if (trig_samp_offset_override && first_run)
                            samp_offset = round(ecog.fs)*(samp_offset_override);
                            first_run = false;
                        else
                            samp_offset = randi([perm_alpha_sec*round(ecog.fs),ecog.n_samples]) * rand_sign;
                        end
                        r_samp_n = r_samp + samp_offset;
                        %r_samp_n_end = r_samp_n + w*round(ecog.fs) - 1;
                        is_valid = (r_samp_n >= 1) && (r_samp_n <= (ecog.n_samples - round(ecog.fs)*w));
                        is_art = true;
                        if (is_valid)
                            coeff = 1/(round(ecog.fs)*w);
                            %is_art = any(Risart( ceil(r_samp_n*coeff):ceil(r_samp_n_end*coeff) ));
                            is_art = Risart(ceil(r_samp_n*coeff));
                        end
                        sat = is_valid && (~ is_art);
                    end

                    % convert r_idx, r_idx_n
                    r_idx = (r_samp - 1)/(round(ecog.fs)*w) + 1;
                    r_idx_n = (r_samp_n - 1)/(round(ecog.fs)*w) + 1;
                    
                    % Show sample time
%                     ecog = H5eeg(fn_h5Lp);
%                     try
%                         realTime = ecog.getTimeFromSample(dir_stamp,r_samp);
%                     catch
%                         realTime = '';
%                     end
%                     fprintf('Sample number: %i, Time: %s\n',r_samp,realTime)
                    
                    % Read timeseries
                    b1 = ecog.bip(ii1,1:2);
                    b2 = ecog.bip(ii2,1:2);
                    v_b1c1 = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(1) r_samp],[1 round(ecog.fs)*w]);
                    v_b1c2 = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(2) r_samp],[1 round(ecog.fs)*w]);
                    v_b2c1 = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(1) r_samp],[1 round(ecog.fs)*w]);
                    v_b2c2 = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(2) r_samp],[1 round(ecog.fs)*w]);
                    v0_b1c1 = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(1) r_samp_n],[1 round(ecog.fs)*w]);
                    v0_b1c2 = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(2) r_samp_n],[1 round(ecog.fs)*w]);
                    v0_b2c1 = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(1) r_samp_n],[1 round(ecog.fs)*w]);
                    v0_b2c2 = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(2) r_samp_n],[1 round(ecog.fs)*w]);
                    
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
                        
                        % save cost variables
                        c_bro_t1(iO) = coh_vals(1);
                        c_brot_t1(iO) = coh_thresh(countp);
                        c_delta_t1(iO) = coh_vals(2);
                        
                        if (trig_save_fig)
                        
                            h = figure;
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
    %                         hb = brainplot_one(sidp,[b1(1),b1(2)],'');
    %                         copyobj(hb,hsub);
    %                         if (trig_title)
    %                             title(sidp);
    %                         end
    %                         hsub = subplot(4,2,4);
    %                         hb = brainplot_one(sidp,[b2(1),b2(2)],'');
    %                         copyobj(hb,hsub);
    %                         if (trig_title)
    %                             title(sidp);
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

                            print(h,sprintf('figures/T1/%s',pname),fig_fmt);
                            if (trig_eps)
                                print(h,sprintf('figures/T1/%s',pname),'-depsc');
                            end
                            close(h);
                            
                            

                            % --- Print Figure S04 ---------------------------

                            h = figure;
                            set(h,'Position',[0 0 PLOT_W_PIX PLOT_H_PIX]);

                            subplot(4,2,5)
                            plot(t,v_b1,'black-','LineWidth',0.5); hold on;
                            %ylabel({ytxt;roi_1});
                            ylabel(ytxt);
                            xlabel(xtxt);
                            xticks(min(t):x_div_secs:max(t));
                            set(gca,'TickDir','out')
                            if (trig_title)
                                title(sprintf('Coherence: %.2f (\\delta:%.2f \\theta:%.2f \\alpha:%.2f \\beta:%.2f \\gamma:%.2f )',...
                                    coh_val0(1),coh_val0(2),coh_val0(3),coh_val0(4),coh_val0(5),coh_val0(6) ));
                            end
                            ymax1 = ceil(max([v_b1,v0_b1])/round_y_nearest) * round_y_nearest;
                            ymin1 = floor(min([v_b1,v0_b1])/round_y_nearest) * round_y_nearest;
                            ylim([ymin1 ymax1]);
                            yl1 = ylim;
                            yt1 = yticks;
                            box off;

                            subplot(4,2,7)
                            plot(t,v_b2,'black-','LineWidth',0.5,'Color',col_faint); hold on;
                            %ylabel({ytxt;roi_2});
                            ylabel(ytxt);
                            xlabel(xtxt);
                            xticks(min(t):x_div_secs:max(t));
                            set(gca,'TickDir','out')
                            ymax2 = ceil(max([v_b2,v0_b2])/round_y_nearest) * round_y_nearest;
                            ymin2 = floor(min([v_b2,v0_b2])/round_y_nearest) * round_y_nearest;
                            ylim([ymin2 ymax2]);
                            yl2 = ylim;
                            yt2 = yticks;
                            box off;

                            hsub = subplot(4,2,6);
    %                         hb = brainplot_one(sid,[b1(1),b1(2)],'');
    %                         copyobj(hb,hsub);
    %                         if (trig_title)
    %                             title(sid);
    %                         end
                            t2 = t + samp_offset/round(ecog.fs);
                            plot(t2,v0_b1,'black-','LineWidth',0.5,'Color',col_faint); hold on;
                            xticks(min(t2):x_div_secs:max(t2));
                            set(gca,'TickDir','out');
                            ylim(yl1);
                            yticks(yt1);
                            set(gca,'Ycolor','none');
                            box off;

                            hsub = subplot(4,2,8);
    %                         hb = brainplot_one(sid,[b2(1),b2(2)],'');
    %                         copyobj(hb,hsub);
    %                         if (trig_title)
    %                             title(sid);
    %                         end
                            plot(t2,v0_b2,'black-','LineWidth',0.5); hold on;
                            xticks(min(t2):x_div_secs:max(t2));
                            set(gca,'TickDir','out');
                            ylim(yl2);
                            yticks(yt2);
                            set(gca,'Ycolor','none');
                            box off;

                            print(h,sprintf('figures/T1d1/%s',pname),fig_fmt);
                            if (trig_eps)
                                print(h,sprintf('figures/T1d1/%s',pname),'-depsc');
                            end
                            close(h);
                            %return
                        end
                    end
                    
                    
                    if (trig_t2)
                        % backup
                        v_b1_b = v_b1;
                        v_b2_b = v_b2;
                        v0_b2_b = v0_b2;
                        forder = 4;
                        
                        
                        
                        if (trig_save_fig)
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
                                v_b1_filt = filtfilt(b,a,double(v_b1));
                                v_b2_filt = filtfilt(b,a,double(v_b2));
                                v0_b2_filt = filtfilt(b,a,double(v0_b2));
    %                             v_b1 = FilterData(v_b1,Fs,'notch',f);

                                h = figure;
                                set(h,'Position',[0 0 1080 1080])
                                subplot(4,2,1)
                                plot(t,v_b1_filt,'black-','LineWidth',0.5); hold on;
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
                                plot(t,v_b2_filt,'black-','LineWidth',0.5); hold on;
                                %ylabel({ytxt;roi_2});
                                ylabel(ytxt);
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
    %                             if (has_roi_r)
    %                                 b2t = b2;
    %                                 b2 = b1;
    %                                 b1 = b2t;
    %                             end
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

                                print(h,sprintf('figures/T2/%s_FILT-%s',pname,metricsp{i4}(3:4)),fig_fmt);
                                if (trig_eps)
                                    print(h,sprintf('figures/T2/%s_FILT-%s',pname,metricsp{i4}(3:4)),'-depsc');
                                end
                                close(h);
                            end
                        end
                        
                        %return
                    end
                    
                    
                    if (trig_t3 && trig_save_fig)
                        
                        [cxx,f] = mscohere(v_b1', v_b2',hamming(Fs/1),[],[],Fs);fn_cacheLp = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheLp,sidp,iMp);
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
                        h = figure;
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
                        
                        print(h,sprintf('figures/T3/%s',pname),fig_fmt);
                        if (trig_eps)
                            print(h,sprintf('figures/T3/%s',pname),'-depsc');
                        end
                        close(h);
                    end
                    
                    
                    
                end
                
                
                if (trig_t4)
                    
                    
                    Rc = RM(countp,:);
                    Risart = art_idx(countp,:);
                    
                    % Pick artifact to plot
                    Ridx = 1:length(Rc);
                    Ridx = Ridx(randperm(length(Ridx)));
                    r_idx_a = NaN;
                    for i3 = 1:length(Ridx)
                        if ( Risart(Ridx(i3)) )
                            r_idx_a = Ridx(i3);
                            break;
                        end
                    end
                    if (isnan(r_idx_a)) %if not found
                        fprintf(2,'W> no time sample found, defaulting to first');
                        r_idx_a = 1;
                    end
                    r_samp_a = (r_idx_a-1)*round(ecog.fs)*w + 1;                    
                    r_samp_1s = (r_idx_a-1)*w + 1;
                    
                    % Load 1 second resolution artifacts
                    a1_b1 = h5read(fn_artLp,'/artifacts',[ii1 r_samp_1s],[1 w]);
                    a1_b2 = h5read(fn_artLp,'/artifacts',[ii2 r_samp_1s],[1 w]);
                    
                    % Read timeseries
                    b1 = ecog.bip(ii1,1:2);
                    b2 = ecog.bip(ii2,1:2);
                    v_b1c1 = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(1) r_samp_a],[1 round(ecog.fs)*w]);
                    v_b1c2 = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(2) r_samp_a],[1 round(ecog.fs)*w]);
                    v_b2c1 = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(1) r_samp_a],[1 round(ecog.fs)*w]);
                    v_b2c2 = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(2) r_samp_a],[1 round(ecog.fs)*w]);
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
                    
                    if (trig_save_fig)

                        h = figure;
                        set(h,'Position',[0 0 0.5*1080 0.5*1080])
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
                end
                
                
                if (trig_t6 && trig_save_fig)
                    x = linspace(0,1,250);
                    h = figure;
                    set(h,'Position',[0 0 1080 1080])
                    
                    for i5 = 1:length(metricsp)
                        
                        metricp = metricsp{i5};
                        subplot(length(metricsp)/2,2,i5)
                        fn_artLp = sprintf('%s/%s_art.h5',dir_artLp,sidp);
                        fn_distLp = sprintf('%s/%s_dists-%s-%i.mat',dir_resLp,sidp,metricp,n_perm);
                        fn_graphLp = sprintf('%s/%s_graph-%s.h5',dir_resLp,sidp,metricp);
                        %fn_permL = sprintf('/mnt/cuenap2/data/results/coh_w10/%s_perm-%s-%i.h5',sid,metric,n_perm);
                        fn_permLp = sprintf('%s/%s_perm-%s-%i.h5',dir_resLp,sidp,metricp,n_perm);
                        
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
                        set(gca,'TickDir','out')
                        
                    
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
                    for iMp = 1:length(metricsp)

                        metricp = metricsp{iMp};
                        
                        if (i_9d1 == 1)
                            t_plot_hrs = n_graph * w / 3600;
                            Tmin_fac = (1/60);
                            
                            % save parameters for opti
                            c_days_t9(iO) = t_plot_hrs / 24;
                            c_mag_t9(iO) = mag(countp);
                            c_magt_t9(iO) = coh_thresh(countp);
                            c_ct_t9(iO) = ct(countp);
                            c_art_t9(iO) = frac_art;
                            
                        elseif (i_9d1 == 2)
                            t_plot_hrs = 1;
                            Tmin_fac = 1;
                        end

                        % Reload coh_thresh
                        fn_cacheLp = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheLp,sidp,iMp);
                        load(fn_cacheLp,'coh_thresh');
                        
                        Rc_t9 = RM(countp,:); %h5read(fn_graphLp,'/R',[countp 1],[1 n_graph]);
                        
                        % Load artifacts
                        Risart = (art_idx(countp,:) ~= 0);
                        
%                         % Load graph
%                         fn_graphLp = sprintf('%s/%s_graph-%s.h5',dir_resLp,sidp,metricp);
%                         Rc_t9 = h5read(fn_graphLp,'/R',[countp 1],[1 n_graph]);
                        
                        
%                         
%                         fn_cache2L = sprintf('%s/%s_R_null.mat',dir_cacheL,sid);
%                         N = load(fn_cache2L);
%                         Rn = N.R_null(count,:);

                        cond_plot = 1:(t_plot_hrs*3600 / w);
                        T_min = linspace(0,max(cond_plot)*w/(60),length(cond_plot));
                        T_min = Tmin_fac * T_min;
                        
                        %return
                        % Subsample to plot with least artifacts and high coherence
                        if (exist('r_samp','var') && exist('r_samp_n','var'))
                            %fprintf('[*] figure t1 plot trigger detected in figure t9 plot\n')
                            % If figure t1 was plotted, choose sample containing time segment used for that
                            if (i_9d1 == 2)
                                % plot 1-hr
                                t_plot_off =  round((t_plot_hrs*3600 / w)/2);
                                cond_offset = ( floor((r_samp - 1) / (w * round(ecog.fs))) ) - t_plot_off;
                                offset_fits = (cond_offset >= 0) & ((length(cond_plot) + cond_offset) < n_graph);
                                if (offset_fits)
                                    cond_plot = cond_plot + cond_offset + round(t9d1_offset_minutes * 60 / w);
                                end
                            end
                            % assert: r_idx == ( floor((r_samp - 1) / (w * round(ecog.fs))) + 1 )
                            cond_plot0 = cond_plot; % duplicate time range to plot for null
                        else
                            art_frac = zeros(n_graph - length(cond_plot) + 1,1);
                            for i3 = 1:(n_graph - length(cond_plot) + 1) 
                                r_mean = mean(Rc_t9(cond_plot + i3 - 1));
                                seg_isa = Risart(cond_plot + i3 - 1);
                                art_frac(i3) = sum(seg_isa)/length(seg_isa) - 1*r_mean;
                            end
                            [~,min_idx] = min(art_frac);

                            if (~isempty(min_idx))
                                cond_plot = cond_plot + min_idx - 1;
                            end
                            cond_plot0 = cond_plot; % duplicate time range to plot for null
                        end
                        
                        Rn_t9 =random(D_t9.d{countp},[1 n_graph]);
                        
                        if (max(cond_plot) > length(Rc_t9))
                            cond_plot = cond_plot(cond_plot <= length(Rc_t9));
                        end
                        if (max(cond_plot0) > length(Rn_t9))
                            cond_plot0 = cond_plot0(cond_plot0 <= length(Rn_t9));
                        end
                        
                        R_plt = Rc_t9(cond_plot);
                        Rn_plt = Rn_t9(cond_plot0);
                        Risa_plt = Risart(cond_plot);
                        Rnisa_plt = Risart(cond_plot0);

                        % Remove artifacts
                        R_plt(Risa_plt) = NaN;
                        Rn_plt(Rnisa_plt) = NaN;
                        f_art = sum(Risa_plt)/numel(Risa_plt);

                        color_div = 0.33*[1 1 1];
                        % Initialize figure
                        if (trig_save_fig)
                        h = figure;
                            if (i_9d1 == 1) % all
                                msize_main = 3;
                                msize_mark = 3;
                                set(h,'Position',[0 0 1*1080 0.5*1080])
                            elseif (i_9d1 == 2) % 1 hour
                                msize_main = 5;
                                msize_mark = 5;
                                set(h,'Position',[0 0 0.5*1080 0.5*1080])
                            end
                            
                            subplot(2,1,1)
                            plot(T_min,R_plt,'black.','MarkerSize',msize_main);hold on;
                            plot([min(T_min) max(T_min)],[coh_thresh(countp) coh_thresh(countp)],'--','color',color_div)
                        
                        
                        
                            if (i_9d1 == 1)
                                % When plotting everything, simply index by
                                % r_idx
                                plot(T_min(r_idx),R_plt(r_idx),'redo','MarkerSize',msize_mark); hold on;
                                xlabel('Time (hours)')
                            elseif(i_9d1 == 2)
                                % When plotting just 1 hr, use index offset 

                                if (offset_fits)
                                    t_plot_off = t_plot_off - round(t9d1_offset_minutes * 60 / w);
                                    plot(T_min(t_plot_off+1),R_plt(t_plot_off+1),'redo','MarkerSize',msize_mark); hold on;
                                    xlabel('Time (minutes)')
                                end
                            end
                            axis([min(T_min) max(T_min) 0 1]); 
                            set(gca,'TickDir','out')


                            %ylabel({sprintf('%s Coherence',metricp(3:end));sprintf('%s-%s',roi_1,roi_2)})
                            ylabel({sprintf('%s Coherence',metricp(3:end))})
                            box off;


                            subplot(2,1,2)
                            plot(T_min,Rn_plt,'black.','MarkerSize',msize_main); hold on;
                            plot([min(T_min) max(T_min)],[coh_thresh(countp) coh_thresh(countp)],'--','color',color_div)
    %                         if (offset_fit)
    %                             plot(T_min(cond_offset+1),Rn_plt(cond_offset+1),'blacko','MarkerSize',10);
    %                         end
                            axis([min(T_min) max(T_min) 0 1])
                            if (i_9d1 == 1) % all
                                xlabel('Time (hours)')
                            elseif (i_9d1 == 2) % 1 hour
                                xlabel('Time (minutes)')
                            end
                            %ylabel({sprintf('Time-shifted %s Coherence',metricp(3:end));sprintf('%s-%s',roi_1,roi_2)})
                            ylabel({sprintf('%s Coherence',metricp(3:end))})
                            set(gca,'TickDir','out')
                            box off;
                        
                        end
                        %return

                        ct_plt = sum(R_plt(~Risa_plt) > coh_thresh(countp)) / length(R_plt(~Risa_plt));
                        ct0_plt = sum(Rn_plt(~Rnisa_plt) > coh_thresh(countp)) / length(Rn_plt(~Rnisa_plt));
                        
                        if (i_9d1 == 1)
                            dir_t9 = 'T9';
                        elseif (i_9d1 == 2)
                            dir_t9 = 'T9d1';
                        end
                        
                        if (trig_save_fig)
                            print(h,sprintf('figures/%s/%s_%s_coh%it%i_ct%i_ctN%i_fart%i',...
                                dir_t9,pname,metricp(3:4),iMp,round(1000*coh_thresh(countp)),round(1000*ct_plt),round(1000*ct0_plt),round(1000*f_art)),fig_fmt);
                            if (trig_eps)
                                print(h,sprintf('figures/%s/%s_%s_coh%it%i_ct%i_ctN%i_fart%i',...
                                    dir_t9,pname,metricp(3:4),iMp,round(1000*coh_thresh(countp)),round(1000*ct_plt),round(1000*ct0_plt),round(1000*f_art)),'-depsc');
                            end
                            close(h);
                        end
                        
                    end
                end
                
               
                
            end
            countp = countp + 1;
        end
    end
end
% end
% end

% Clear loop indices
clear i;
clear j;

% build cost function
% % optimization variables
% v_sub = zeros(n_opti,1);
% v_bchans = zeros(n_opti,2);
% v_samp = zeros(n_opti,1);
% c_bro_t1 = zeros(n_opti,1);
% c_brot_t1 = zeros(n_opti,1);
% c_delta_t1 = zeros(n_opti,1);
% % c_bro_t1d2 = zeros(n_opti,1);
% % c_delta_t1d2 = zeros(n_opti,1);
% c_days_t9 = zeros(n_opti,1);
% c_mag_t9 = zeros(n_opti,1);
% c_magt_t9 = zeros(n_opti,1);
% c_ct_t9 = zeros(n_opti,1);
% c_art_t9 = zeros(n_opti,1);

cost = (c_brot_t1 - c_bro_t1) + 1.5*(-c_delta_t1) + (c_days_t9 / 7) + (-c_mag_t9) + (-c_ct_t9) + 2*(c_art_t9);
[cost_s,sIdx] = sort(cost);
min_cost = cost_s(1);
min_sub = v_sub(sIdx(1));
min_bchans = v_bchans(sIdx(1),:);
min_samp = v_samp(sIdx(1));

fprintf('[*] cost: %.9f',min_cost)
fprintf('\t%s\t%i\t%i\t%i\n',Subjectsp{min_sub},min_bchans(1),min_bchans(2),min_samp)

% Print finish message
fprintf('Done.\n')
