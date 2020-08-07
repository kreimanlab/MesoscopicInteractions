% Feb 19, 2019
% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T1’: File exists
% mkdir: cannot create directory ‘figures/T1d1’: File exists
% mkdir: cannot create directory ‘figures/T2’: File exists
% mkdir: cannot create directory ‘figures/T3’: File exists
% mkdir: cannot create directory ‘figures/T6’: File exists
% mkdir: cannot create directory ‘figures/T9’: File exists
% mkdir: cannot create directory ‘figures/T9d1’: File exists
% /media/klab/internal/data/stamps/m00005_00000000_stamps.txt
% --- H5eeg.getTimeFromSample BEGIN ---
% [!] Corrected stamp 74645095 -> 73347934, error: -1850
% --- H5eeg.getTimeFromSample EXIT_SUCCESS ---
% [E] video file not found for: m00005, 00000000
% Sample number: 72056321, Time: 05/09/2007 01:14:53
% , Video: 
% Br364	172
% BrN58	172
% m00005__35_superiortemporal_PT7-PT8__84_parsopercularis_RP55-RP56__29mm_ct405_ctt0-1_mag221-49_coh1t172_Br364_Samp72056321_BrN58_SampN72102401
% Time-shifted coherence: 0.053
% Time-shifted gamma coherence: NaN
% [*] figure t1 plot trigger detected in figure t9 plot
% Starting parallel pool (parpool) using the 'local' profile ...
% connected to 6 workers.
% /media/klab/internal/data/stamps/m00005_00000000_stamps.txt
% --- H5eeg.getTimeFromSample BEGIN ---
% --- H5eeg.getTimeFromSample EXIT_SUCCESS ---
% Total number of days: 5.524
% Number of positive coherences: 19339
% [*] figure t1 plot trigger detected in figure t9 plot
% /media/klab/internal/data/stamps/m00005_00000000_stamps.txt
% --- H5eeg.getTimeFromSample BEGIN ---
% --- H5eeg.getTimeFromSample EXIT_SUCCESS ---
% Done.

close all;
clear;
rng shuffle;

COLOR_MIDNIGHT_LINE = 0.5*[1 1 1];

metricp = 'pcBroadband';
perm_alpha_sec = 12;
n_perm = 10000;
cp_thresh_override = 0.05;
p_val = 0.01;

% Path definitions
dir_h5Lp = '/media/klab/internal/data/h5_notch20'; %'/media/klab/KLAB101/h5_notch20';
dir_artLp = sprintf('%s/art_nosz',dir_h5Lp); %dir_artLp = '/media/klab/KLAB101/h5_notch20/art_nosz';
dir_resLp = '/home/jerry/data/results/coh_w10'; %dir_resLp = '/media/klab/KLAB101/results/coh_w10'; 
dir_corLp = '/media/klab/internal/data/coreg';
dir_cacheLp = './cache';
dir_stamp = '/media/klab/internal/data/stamps';
dir_video = '/media/klab/internal/data/videos';
%subjects_dirLp = '/mnt/cuenap_ssd/coregistration';

%dir_h5Lp = '/mnt/cuenap/data/h5_notch20';

% OLD DO NOT USE
%metricsp = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};

metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
%metricsp = {'pcBroadband','pcGamma'};
%metricsp = {'pcBroadband'};

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

% m00047
% 64410001 - stamp - 128821664
% 87692501 - stamp - 175386656
%[*] cost: -26116.387557812	m00047	26	54	64875001 - stamp 129751664
%[*] cost: -17351.387557812	m00047	26	54	84207501 - stamp 168416656
%[*] cost: 0.232203059	m00047	26	54	86870001 - stamp 173741664 - line
    %173740017 - 02/09/2010 17:03:39 - _2895.avi 28 seconds 
%[*] cost: 0.241688453	m00047	26	54	86840001 - stamp 173681664
%[*] cost: 0.265921794	m00047	26	54	86665001 - stamp 173331664
%[*] cost: 0.261862211	m00047	26	54	86545001

% m00061
%[*] cost: -0.423532230	m00061	27	54	79817501 - stamp 44116036
%[*] cost: -0.488127392	m00061	27	54	79137501
%[*] cost: -0.499706310	m00061	27	54	79770001
%[*] cost: 0.158448989	m00061	20	54	8037501
%[*] cost: 0.144253090	m00061	21	54	76060001
%[*] cost: -0.128401051	m00061	13	54	78460001

% bypass
%Subjectsp = {'m00100'}; '79ba6df8', rhanded, rh
%Subjectsp = {'m00084'}; % rhanded, lh
%Subjectsp = {'m00001'}; % rhanded, lh
%Subjectsp = {'m00061'}; % rhanded, lh 'files':  'dcdce1b0' -novid, 'dd69e165' -vid, 'c5036be4' -vid, '70706eef' -vid
%Subjectsp = {'m00047'}; % a289baf2 only, rhanded, rh

%r_samp_const = 106865001; % set to nan to pick randomly
%r_samp_const = 90545001; % strong t1, t9 good but next to art
%r_samp_const = 22547501; % bad t1, but strong t9 1 hr
%r_samp_const = 90835001; % good t9, but not as good t1, high amplitude
%r_samp_const = 91430001; % t9 missing right half

% working examples
%75120001 - 27 parsopercularis, 46 superiortemporal, m00084 - small number of artifacts just before
% old example
%22755001 - 13 parstriangularis, 58 inferiorparietal, m00001

%----------------------------------- OLD ---------------------------------
%r_samp_const = 22755001;

% Candidates:
% 75120001 - 27, 46, m00084
% 75132501 - low frequency, not as strong but doable
% 30102501 - high frequency, not ideal
% Failed:
%31512501 37867501 39132501 53297501 83542501 98377501
%77837501; % 22755001; % good t9, magnitude of t1 is only 0.339 (DEFAULT)
% bchan1_const = 13; % 13; % m00001 13
% bchan2_const = 58; % 58; % m00001 58


%roi_1 = 'parstriangularis'; %
%roi_1 = 'superiortemporal';
%roi_1 = 'parsopercularis';
%roi_1 = 'inferiortemporal';
%roi_2 = 'middletemporal';
%roi_2 = 'superiortemporal';
%roi_2 = 'inferiorparietal'; %

%Subjectsp = {'m00001'};
%Subjectsp = {'m00084'};
%-----------------------------------

% new new - Jan 10, 2019
Subjectsp = {'m00005'};
roi_1 = 'superiortemporal';
roi_2 = 'parsopercularis';
bchan1_const = 35;
bchan2_const = 84;
r_samp_const = 72056321;


% m00084 - 27 - 46 
% 78515001 - 01/29/2012 07:14:30 - m00084_8f2a6492_0686.avi,43 - SLEEPING
% 68527501 presses button 2 seconds before
% 18762501 - sleep? higher frequency coherence
% 38300001 -  01/26/2012 09:01:01 - m00084_501329d1_0027.avi,61 - CONVERSATION
% 37070001 - talks, but looks like discharges during speech
% 4877501 - 01/24/2012 13:31:02 - m00084_ac5bb0aa_0162.avi,42 - speaks during the last second, raw looks spiky and low voltage

% m00084 - 21 - 46
% 59087501 - wake from sleep
% 59720001 - presses button, makes a vocal noise
% 77837501 - 01/29/2012 06:29:20 - m00084_8f2a6492_0663.avi,93 - reposition stuffed animal in front of nurse

% -- new calculations -----------------------------------------------------
% parsopercularis, inferiorparietal
%Subjectsp = {'m00061'}; % dominant hemi and video
%Subjectsp = {'m00001','m00005','m00047','m00061'};
%[*] cost: -0.676649978	m00061	27	54	79115001 - t9 end cluster
%[*] cost: -0.649445107	m00061	27	54	80445001
%[*] cost: 0.093009810	m00061	13	54	74092501 - 80mm middle of art
%[*] cost: -0.354445984	m00061	27	54	74210001 
%[*] cost: -0.182848577	m00061	21	54	76445001 - 65mm
%[*] cost: -0.539553795	m00061	20	54	79240001 - 75mm

% parsopercularis, superiortemporal
%Subjectsp = {'m00001','m00003','m00005','m00019','m00022','m00023','m00025','m00026','m00030','m00032','m00045','m00061','m00073','m00075','m00079','m00084','m00095'};
% 3n,19,25n,30,32,61,73,75,84,95
%Subjectsp = {'m00019','m00030','m00032','m00061','m00073','m00075','m00084','m00095'}; % dominant hemi and video
%Subjectsp = {'m00095'};
%[*] cost: 0.564330871	m00095	20	53	23932501 - weak interaction
%[*] cost: 0.539396019	m00095	20	33	28727501 - too many arts
%[*] cost: 0.572074498	m00095	20	47	47595001 - looks ok, could be stronger
%[*] cost: 0.546873493	m00095	20	41	54105001 - looks ok, 26mm
%[*] cost: 0.560554637	m00095	20	41	61880001 - 
%[*] cost: 0.584691258	m00084	27	61	96155001 - arts, big change at 50%
%[*] cost: 0.497183125	m00084	20	36	71732501 - big change at 50%
%[*] cost: 0.288932666	m00084	20	54	64352501 - arts, big change
%[*] cost: 0.597722641	m00084	20	41	72342501 - change more gradual, t9 looks good
%[*] cost: 0.814054461	m00084	22	46	78242501
%[*] cost: 0.204371957	m00075	14	98	29332501 - 18 mm
%[*] cost: 0.224854736	m00075	14	98	32395001
%[*] cost: 0.236381678	m00075	14	98	38712501
%[*] cost: 0.261469712	m00073	7	21	3327501 - 21 mm
%[*] cost: 0.318063100	m00073	7	21	5400001
%[*] cost: 0.314823946	m00073	7	21	8427501
%[*] cost: -0.155321119	m00061	13	39	74105001 - 60 mm - ct59 - 10/29/2010
    %00:31:44 - 70706eef - _0543.avi - 1:38 - putting on blindfold going to bed
%[*] cost: -0.112810579	m00061	13	39	72645001
%[*] cost: 0.134045247	m00061	27	36	76422501 - t9 not good
%[*] cost: -0.320159630	m00061	21	36	80770001 - t9 not good
%[*] cost: -0.426196223	m00061	21	48	80590001 - 40 mm not best placement
%[*] cost: 0.079583935	m00061	27	39	22232501 - 52mm ct60 - 10/26/2010 14:52:13 - dd69e165 - _0196.avi - -1 sec - reading during stim
%[*] cost: -0.358221876	m00061	21	39	80067501
%[*] cost: 0.837716711	m00032	55	60	23692501 - 22mm
%[*] cost: 0.761445999	m00032	28	60	51460001 - 30mm
%[*] cost: 0.849530846	m00032	28	60	66437501
%[*] cost: 0.835929138	m00032	55	60	124770001
%[*] cost: 0.239254323	m00030	29	31	31370001 - 20mm
%[*] cost: 0.206529823	m00030	29	31	66517501
%[*] cost: 0.236791518	m00030	29	31	73600001
% m00030	29	31  72342501 - 10/28/2008 22:04:19 0a9d761f - _0443.avi 23 sec - talking - elec too close together
%[*] cost: -0.139612804	m00019	14	28	55607501 - 21mm
%[*] cost: 0.215532971	m00019	5	14	57387501 - 23mm
%[*] cost: 0.362740277	m00019	14	41	53862501 - 35mm
%[*] cost: 0.155424700	m00019	14	27	28620001 - 25mm

% parstriangularis, superiortemporal
%Subjectsp = {'m00001','m00003','m00005','m00019','m00021','m00022','m00023','m00025','m00026','m00032','m00035','m00045','m00053','m00084'};
% 3n,19,21n,25n,32,35,53,84
%Subjectsp = {'m00019','m00032','m00035','m00053','m00084'};% dominant hemi and video
%Subjectsp = {'m00035'};
%[*] cost: 0.463727641	m00084	3	61	66562501 - arts
%[*] cost: 0.507391921	m00084	20	41	77422501
%[*] cost: 0.275485745	m00084	4	42	92715001
%[*] cost: 0.587059769	m00053	2	40	104900001 - arts
%[*] cost: 0.586993876	m00053	2	40	112397501
%[*] cost: 0.576190564	m00053	2	40	120232501
%[*] cost: 0.523224774	m00035	77	88	522501
%[*] cost: 0.528096888	m00035	77	88	15890001
%[*] cost: 0.789399004	m00035	77	90	8185001
%[*] cost: 0.793726301	m00035	77	90	19687501
%[*] cost: 0.821242760	m00032	26	61	30217501
%[*] cost: 0.531650747	m00032	26	61	96310001
%[*] cost: 0.425622409	m00032	26	60	103780001
%[*] cost: 0.786905825	m00032	28	60	120367501
%[*] cost: -0.136224830	m00019	28	56	1700001
%[*] cost: -0.325221711	m00019	28	56	5937501
%[*] cost: 0.305643355	m00019	41	56	73312501
%[*] cost: -0.064330439	m00019	28	54	12265001
%[*] cost: 0.371564081	m00019	7	20	60667501
%--------------------------------------------------------------------------

%bchan1_const = 30;
%bchan2_const = 80;
%r_samp_const = NaN;
%[*] cost: 0.011246898	m00001	14	51	82335001
%[*] cost: -0.065562469	m00001	4	57	83825001
%[*] cost: -0.064870114	m00001	4	51	83425001
%[*] cost: -0.025408129	m00001	21	50	116170001
%[*] cost: -0.165981444	m00100	30	80	10840001
%[*] cost: -0.213723722	m00100	30	80	7315001 % Sample 7315001, stamp 14630936 
%   time: 05/10/2013 23:35:28 - estimated time: 10-May-2013 19:31:38 -
%   00:01:27.957 into video 0243.avi
%[*] cost: 0.660226596	m00084	20	66	89857501
%[*] cost: 0.627504541	m00084	20	66	83155001
%[*] cost: 0.613476766	m00084	20	66	60247501
%[*] cost: 0.136753701	m00047	26	54	64410001
%[*] cost: 0.269184195	m00047	26	54	86635001
%[*] cost: 0.261862211	m00047	26	54	86545001
%[*] cost: 0.091736786	m00047	26	54	87692501 % Sample 87692501, stamp 175386656 estimated time: 07-Feb-2010 17:15:24


% Pick areas to plot
% C.AtlROIs{2}.LH.struct_names
% (PFC, area 9/46) and the anterior cingulate cortex (ACC, areas 24c and 32).

%roi_1 = 'parsorbitalis';
%roi_1 = 'parsopercularis';


% roi_1 = 'parahippocampal';
% roi_2 = 'lateralorbitofrontal';

%roi_1 = 'middletemporal';
%roi_2 = 'precentral';

%roi_2 = 'supramarginal';
%roi_1 = find(strcmp(rois,'parsorbitalis'),1);
%roi_2 = find(strcmp(rois,'inferiorparietal'),1);

% Triggers
system('mkdir figures');
fig_fmt = '-depsc';
trig_title = false; % show titles
trig_both_elec = false; % plot both electrodes in bipolar pair
trig_eps = false;
trig_t1 = true;
trig_t2 = true; %
trig_t3 = true; %
trig_t4 = false; %default: false
trig_t6 = true; %
trig_t9 = true; %
trig_t9d1 = true; %
t9d1_offset_minutes = 25; % if 0, t9d1 centers the chosen time point
trig_plot_filebreak_x = false; % whether to indicate file breaks with black x

% time shift
trig_samp_offset_override = true;
samp_offset_override = 180; %seconds
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

%i_s = randi([1 length(Subjects)]);
for i_s = 1:length(Subjectsp)
% Main loop
for iMp = 1%1:1 %length(metricsp)
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
    
%     % print time
%     datetime_v = h5readatt(fn_h5Lp,'/h5eeg','datetime');
%     files = h5readatt(fn_h5Lp,'/h5eeg','files');
%     ds_factors = h5readatt(fn_h5Lp,'/h5eeg','ds_factors');
%     n_samples = h5readatt(fn_h5Lp,'/h5eeg','n_samples');
%     fileIdx = find(cumsum(n_samples) > r_samp_const);
%     etime = h5read(fn_h5Lp,'/h5eeg/aux',[1 r_samp_const],[1 1]);
%     start_etime = h5read(fn_h5Lp,'/h5eeg/aux',[1 1],[1 1]);
%     try
%         est_time = datestr(datenum(replace(datetime_v{fileIdx},'T',' ')) + seconds(r_samp_const*(1/(ecog.fs*ds_factors(fileIdx)))));
%         fprintf('Sample %i, stamp %i estimated time: %s\n',r_samp_const,etime,est_time);
%     catch
%     end
    
    
    
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
                cond_bpair_bypass = (ii1 == bchan1_const) && (ii2 == bchan2_const);
            else
                % if either bchan1_const or bchan2_const is nan, do all
                cond_bpair_bypass = true;
            end
            if (((has_roi_f || has_roi_r) && (Dmat(ii1,ii2) > dist_thresh)) && cond_bpair_bypass)
                % ct threshold based on permutations
                p_val_ct_2 = 0.01/n_comb;
                ct_thresh_perm = (binoinv([p_val_ct_2 (1-p_val_ct_2)],n_graph,p_val/bonf_seconds)/n_graph);
                        
                pname = sprintf('%s__%i_%s_%s__%i_%s_%s__%imm_ct%i_ctt%i-%i_mag%i-%i_coh%it%i',...
                    sidp,ii1,atl_labels{ecog.bip(ii1,1)},bchan_labels{ii1},ii2,atl_labels{ecog.bip(ii2,1)},...
                    bchan_labels{ii2},round(Dmat(ii1,ii2)),round(1000*ct(countp)),...
                    ceil(1000*ct_thresh_perm(1)),ceil(1000*ct_thresh_perm(2)),round(1000*mag(countp)),...
                    round(1000*mag_std(countp)),iMp,round(1000*coh_thresh(countp)));
                %fprintf('%s\n',pname);
                
                % Plot T1
                % check passes ct threshold
                %if (trig_t1 && ( ct(count) > ct_thresh ))
                if ((trig_t1 || trig_t2 || trig_t3) && ( ct(countp) > 0 ))
                    
                    % Pick time sample to plot
                    Rc = R(countp,:);
                    Risart = art_idx(countp,:);
                    if (isnan(r_samp_const))
                        Ridx = 1:length(Rc);
                        Ridx = Ridx(randperm(length(Ridx)));
                        r_idx = NaN;
                        for i3 = 1:length(Ridx)
                            if ( (Rc(Ridx(i3)) > coh_thresh(countp)) && (~Risart(Ridx(i3))) )
                            % -----------------------------------------------------------
                            %if ( (Rc(Ridx(i3)) > 0.35) && (Rc(Ridx(i3)) < 0.55) && (~Risart(Ridx(i3))) )
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
                    ecog = H5eeg(fn_h5Lp);
                    try
                        realTime = ecog.getTimeFromSample(dir_stamp,r_samp);
                        vidStr = ecog.getVideo(dir_video,r_samp,realTime);
                    catch
                        realTime = '';
                        vidStr = '';
                    end
                    fprintf('Sample number: %i, Time: %s, Video: %s\n',r_samp,realTime,vidStr);
                    
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
                    v0_b1 = v0_b1 - mean(v0_b1);
                    v0_b2 = v0_b2 - mean(v0_b2);
                    
                    % If correlation is negative, flip sign
                    if (corr2(v_b1',v_b2') < 0)
                        v_b2 = (-1) * v_b2;
                        v0_b2 = (-1) * v0_b2;
                    end
                    
                    if (has_roi_r)
                        b2t = b2;
                        b2 = b1;
                        b1 = b2t;
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
                        %if ((i3 == 1) || (i3 == 5))
                        if (i3 == 1)
                            pname = sprintf('%s_%s%i',pname,mname,round(1000*coh_val));
                        end
                        CaT = load(sprintf('%s/xsub_out_%s_%i',dir_cacheLp,sid,i3),'coh_thresh');
                        fprintf('%s%i\t%i\n',mname,round(1000*coh_val),round(1000*CaT.coh_thresh(countp)));
                    end                   
                    pname = sprintf('%s_Samp%i',pname,r_samp);
                    for i3 = 1:length(metricsp)
                        mname = [metricsp{i3}(3:4),'N'];
                        %if ((i3 == 1) || (i3 == 5))
                        if (i3 == 1)
                            pname = sprintf('%s_%s%i',pname,mname,round(1000*coh_val0(i3)));
                        end
                        CaT = load(sprintf('%s/xsub_out_%s_%i',dir_cacheLp,sid,i3),'coh_thresh');
                        fprintf('%s%i\t%i\n',mname,round(1000*coh_val0(i3)),round(1000*CaT.coh_thresh(countp)));
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
%                         if (has_roi_r)
%                             b2t = b2;
%                             b2 = b1;
%                             b1 = b2t;
%                         end
                        

                        h2 = figure('visible','off');
                        hsub = subplot(4,2,2);
                        if (trig_both_elec)
                            hb = brainplot_one(sidp,[b1(1),b1(2)],'');
                        else
                            hb = brainplot_one(sidp,[b1(1)],'');
                        end
                        copyobj(hb,hsub);
                        if (trig_title)
                            title(sidp);
                        end
                        hsub = subplot(4,2,4);
                        if (trig_both_elec)
                            hb = brainplot_one(sidp,[b2(1),b2(2)],'');
                        else
                            hb = brainplot_one(sidp,[b2(1)],'');
                        end
                        copyobj(hb,hsub);
                        if (trig_title)
                            title(sidp);
                        end
                        print(h2,sprintf('figures/T1/%s_brains',pname),'-dpng','-r300');
                        close(h2);

                        

                        print(h,sprintf('figures/T1/%s',pname),fig_fmt);
                        if (trig_eps)
                            print(h,sprintf('figures/T1/%s',pname),'-depsc');
                        end
                        close(h);
                        
                        % --- Print Figure S04 ---------------------------
                        sec_before = 2; % space before
                        sec_after = 2; % space after
                        t_plot_offset_sec = 1; % space between time offset
                        soffset = round(sec_before*(ecog.fs));
                        loffset = round((sec_before + sec_after)*(ecog.fs));
                        v_b1c1p = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(1) r_samp-soffset],[1 round(ecog.fs)*w+loffset]);
                        v_b1c2p = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(2) r_samp-soffset],[1 round(ecog.fs)*w+loffset]);
                        v_b2c1p = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(1) r_samp-soffset],[1 round(ecog.fs)*w+loffset]);
                        v_b2c2p = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(2) r_samp-soffset],[1 round(ecog.fs)*w+loffset]);
                        v0_b1c1p = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(1) r_samp_n-soffset],[1 round(ecog.fs)*w+loffset]);
                        v0_b1c2p = h5read(fn_h5Lp,'/h5eeg/eeg',[b1(2) r_samp_n-soffset],[1 round(ecog.fs)*w+loffset]);
                        v0_b2c1p = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(1) r_samp_n-soffset],[1 round(ecog.fs)*w+loffset]);
                        v0_b2c2p = h5read(fn_h5Lp,'/h5eeg/eeg',[b2(2) r_samp_n-soffset],[1 round(ecog.fs)*w+loffset]);
                        % Bipolar montage
                        v_b1p = v_b1c1p - v_b1c2p;
                        v_b2p = v_b2c1p - v_b2c2p;
                        v0_b1p = v0_b1c1p - v0_b1c2p;
                        v0_b2p = v0_b2c1p - v0_b2c2p;
                        % Subtract mean
                        v_b1p = v_b1p - mean(v_b1p);
                        v_b2p = v_b2p - mean(v_b2p);
                        v0_b1p = v0_b1p - mean(v0_b1p);
                        v0_b2p = v0_b2p - mean(v0_b2p);
                        % If correlation is negative, flip sign
                        if (corr2(v_b1p',v_b2p') < 0)
                            v_b2p = (-1) * v_b2p;
                            v0_b2p = (-1) * v0_b2p;
                        end
                        
                        
                        y_div_uv = 500;
                        t = linspace(0,length(v_b1p)/round(ecog.fs),length(v_b1p));
                        t2 = t + samp_offset/round(ecog.fs);
                        vmax = max([v_b1p, v0_b1p, v_b2p, v0_b2p]);
                        vmin = min([v_b1p, v0_b1p, v_b2p, v0_b2p]);
                        % extend to nearest division
                        ynearest = 100;
                        vmin = floor(vmin/ynearest)*ynearest;
                        vmax = ceil(vmax/ynearest)*ynearest;
                        
                        % Filter
                        % Gaussian filter
                        filter_sec = 5;
                        gain_adjust = 2;
                        wf = gausswin(filter_sec*round(ecog.fs));
                        wf = gain_adjust * wf/sum(wf);
                        col_filt = 0.7*[1 1 1];
                        v_b1p_f = filter(wf,1,v_b1p);
                        v0_b1p_f = filter(wf,1,v0_b1p);
                        v_b2p_f = filter(wf,1,v_b2p);
                        v0_b2p_f = filter(wf,1,v0_b2p);
                        
                        vmax1 = max([v_b1p, v0_b1p, v_b1p_f, v0_b1p_f                                ]);
                        vmin1 = min([v_b1p, v0_b1p, v_b1p_f, v0_b1p_f                                ]);
                        vmax2 = max([                                v_b2p, v0_b2p, v_b2p_f, v0_b2p_f]);
                        vmin2 = min([                                v_b2p, v0_b2p, v_b2p_f, v0_b2p_f]);
                        vmin1 = floor(vmin1/ynearest)*ynearest;
                        vmax1 = ceil(vmax1/ynearest)*ynearest;
                        vmin2 = floor(vmin2/ynearest)*ynearest;
                        vmax2 = ceil(vmax2/ynearest)*ynearest;
                        
                        %h = figure('visible','off');
                        h = figure('visible','on');
                        
                        
                        
                        % tick locations
                        sTicks = [sec_before sec_before+w 2*sec_before+w+sec_after+t_plot_offset_sec 2*(sec_before+w)+sec_after+t_plot_offset_sec];
                        sTicksLab = {'0',num2str(w),num2str(samp_offset/round(ecog.fs)),num2str(samp_offset/round(ecog.fs) + w)};
                        xtxtp = sec_before+w+sec_after+t_plot_offset_sec/2;
                        set(h,'Position',[0 0 PLOT_W_PIX round(PLOT_H_PIX/2)]);
                        
                        subplot(2,1,1)
                        vmin = vmin1;
                        vmax = vmax1;
                        hold all;
                        
                        % filtered wave
                        
                        %plot(t,(v_b1p_f/(max(v_b1p_f)-min(v_b1p_f)))*(max(v_b1p)-min(v_b1p)),'-','LineWidth',0.5,'Color',col_filt); hold on;
                        plot(t,v_b1p_f,'-','LineWidth',0.5,'Color',col_filt); hold on;
                        
                        % IFP
                        plot(t,v_b1p,'black-','LineWidth',0.5); hold on;
                        
                        t2_p = linspace(max(t)+t_plot_offset_sec,max(t) + (max(t)-min(t))+t_plot_offset_sec,length(t));
                        % filtered wave
                        
                        %plot(t2_p,(v0_b1p_f/(max(v0_b1p_f)-min(v0_b1p_f)))*(max(v0_b1p)-min(v0_b1p)),'-','LineWidth',0.5,'Color',col_filt); hold on; 
                        plot(t2_p,v0_b1p_f,'-','LineWidth',0.5,'Color',col_filt); hold on; 
                        
                        % IFP
                        plot(t2_p,v0_b1p,'black-','LineWidth',0.5); hold on; %,'Color',col_faint
                        
                        plot([sTicks(1) sTicks(1)],[vmin vmax],'black:','Color',col_faint); hold on;
                        plot([sTicks(2) sTicks(2)],[vmin vmax],'black:','Color',col_faint); hold on;
                        %plot([sTicks(3) sTicks(3)],[vmin vmax],'black:','Color',col_faint); hold on;
                        %plot([sTicks(4) sTicks(4)],[vmin vmax],'black:','Color',col_faint); hold on;
                        text(xtxtp,double(vmin),'//','HorizontalAlignment','Center','VerticalAlignment','middle','FontSize',12);
                        t_ticks = [min(t):x_div_secs:max(t), min(t2):x_div_secs:max(t2)];
                        t_ticks_p = [min(t):x_div_secs:max(t), min(t2_p):x_div_secs:max(t2_p)];
                        ylabel(ytxt);
                        xlabel(xtxt);
                        xticks(t_ticks_p(1:2));
                        xticklabels(t_ticks(1:2));
                        %yticks(round(vmin:y_div_uv:vmax));
                        set(gca,'TickDir','out');
                        box off;
                        axis([min(t) max(t2_p) vmin vmax]);
                        % Arrowhead params
                        arrow_x_w = 0.04;
                        arrow_y_w = 0.05;
                        % Arrows
                        axp = get(gca,'Position');
                        xs=axp(1);
                        xe=axp(1)+axp(3)+arrow_x_w;
                        ys=axp(2);
                        ye=axp(2)+axp(4)+arrow_y_w;
                        annotation('arrow', [xs xe],[ys ys]);
                        xticks(sTicks);
                        xticklabels(sTicksLab);
                        %set(gca,'XTick',[]);
                        
                        subplot(2,1,2)
                        vmin = vmin2;
                        vmax = vmax2;
                        %plot(t,v_b2,'black-','LineWidth',0.5,'Color',col_faint); hold on;
                        %plot(t2,v0_b2,'black-','LineWidth',0.5); hold on;
                        hold all;
                        
                        %plot(t,v_b2p,'black-','LineWidth',0.5); hold on;
                        
                        %plot(t,(v_b2p_f/(max(v_b2p_f)-min(v_b2p_f)))*(max(v_b2p)-min(v_b2p)),'-','LineWidth',0.5,'Color',col_filt); hold on;
                        plot(t,v_b2p_f,'-','LineWidth',0.5,'Color',col_filt); hold on;
                        
                        % IFP
                        plot(t,v_b2p,'black-','LineWidth',0.5); hold on;
                        
                        t2_p = linspace(max(t)+t_plot_offset_sec,max(t) + (max(t)-min(t))+t_plot_offset_sec,length(t));
                        %plot(t2_p,v0_b2p,'black-','LineWidth',0.5); hold on; %,'Color',col_faint
                        % filtered wave
                        
                        %plot(t2_p,(v0_b2p_f/(max(v0_b2p_f)-min(v0_b2p_f)))*(max(v0_b2p)-min(v0_b2p)),'-','LineWidth',0.5,'Color',col_filt); hold on; 
                        plot(t2_p,v0_b2p_f,'-','LineWidth',0.5,'Color',col_filt); hold on; 
                        
                        % IFP
                        plot(t2_p,v0_b2p,'black-','LineWidth',0.5); hold on; %,'Color',col_faint
                        
                        
                        plot([sTicks(1) sTicks(1)],[vmin vmax],'black:','Color',col_faint); hold on;
                        plot([sTicks(2) sTicks(2)],[vmin vmax],'black:','Color',col_faint); hold on;
                        plot([sTicks(3) sTicks(3)],[vmin vmax],'black:','Color',col_faint); hold on;
                        plot([sTicks(4) sTicks(4)],[vmin vmax],'black:','Color',col_faint); hold on;
                        text(xtxtp,double(vmin),'//','HorizontalAlignment','Center','VerticalAlignment','middle','FontSize',12);
                        t_ticks = [min(t):x_div_secs:max(t), min(t2):x_div_secs:max(t2)];
                        t_ticks_p = [min(t):x_div_secs:max(t), min(t2_p):x_div_secs:max(t2_p)];
                        ylabel(ytxt);
                        xlabel(xtxt);
                        xticks(t_ticks_p);
                        xticklabels(t_ticks);
                        %yticks(round(vmin:y_div_uv:vmax));
                        set(gca,'TickDir','out');
                        box off;
                        axis([min(t) max(t2_p) vmin vmax]);
                        % Arrows
                        axp = get(gca,'Position');
                        xs=axp(1);
                        xe=axp(1)+axp(3)+arrow_x_w;
                        ys=axp(2);
                        ye=axp(2)+axp(4)+arrow_y_w;
                        annotation('arrow', [xs xe],[ys ys]);
                        xticks(sTicks);
                        xticklabels(sTicksLab);
                        %set(gca,'XTick',[]);
                        
                        %return
                        
                        % --- OLD timeshift plot -------------------------
%                         set(h,'Position',[0 0 PLOT_W_PIX PLOT_H_PIX]);
%                         subplot(4,2,5)
%                         % top-left
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
                        
                        % --- OLD timeshift plot -------------------------
                        
                        % time-shifted coherence
                        coh_ts = coherence( v0_b2, v_b1, round(ecog.fs) );
                        coh_ts_intime = coherence( v_b2, v_b1, round(ecog.fs) );
                        fprintf('In-time coherence: %.3f\n',coh_ts_intime(1));
                        fprintf('Time-shifted coherence: %.3f\n',coh_ts(1));
                        %fprintf('Time-shifted gamma coherence: %.3f\n',coh_ts(5));
                        
                        % Apply moving mean filter
                        k_filter = filter_sec*round(ecog.fs);
                        coh_ts_f = coherence( movmean(v0_b2,k_filter), movmean(v_b1,k_filter), round(ecog.fs) );
                        coh_ts_f_intime = coherence( movmean(v_b2,k_filter), movmean(v_b1,k_filter), round(ecog.fs) );
                        fprintf('Movmean (%i s) In-time coherence: %.3f\n',filter_sec,coh_ts_f_intime(1));
                        fprintf('Movmean (%i s) Time-shifted coherence: %.3f\n',filter_sec,coh_ts_f(1));
                        
                        % Apply gaussian filter
                        coh_ts_f = coherence( filter(wf,1,v0_b2), filter(wf,1,v_b1), round(ecog.fs) );
                        coh_ts_f_intime = coherence( filter(wf,1,v_b2), filter(wf,1,v_b1), round(ecog.fs) );
                        fprintf('Gaussian (%i s) In-time coherence: %.3f\n',filter_sec,coh_ts_f_intime(1));
                        fprintf('Gaussian (%i s) Time-shifted coherence: %.3f\n',filter_sec,coh_ts_f(1));
                        
%                         In-time coherence: 0.405
%                         Time-shifted coherence: 0.053
%                         Movmean (5s) In-time coherence: 0.643
%                         Movmean (5s) Time-shifted coherence: 0.554
%                         Gaussian (5s) In-time coherence: 0.736
%                         Gaussian (5s) Time-shifted coherence: 0.613
                        
                        %figure; subplot(2,1,1); plot(v_b2); axis tight; subplot(2,1,2); plot(v_b1); axis tight;
                        
                        %return
                        
                        print(h,sprintf('figures/T1d1/%s_cohTS%i',pname,round(1000*coh_ts(1))),fig_fmt);
                        if (trig_eps)
                            print(h,sprintf('figures/T1d1/%s_coh_TS%i',pname,round(1000*coh_ts(1))),'-depsc');
                        end
                        close(h);
                        %return
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
                            t = linspace(0,length(v_b1)/round(ecog.fs),length(v_b1));
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
                        
                        %return
                    end
                    
                    
                    if (trig_t3)
                        

                        [cxx,f] = mscohere(v_b1', v_b2',hamming(round(Fs/5)),[],[],Fs);fn_cacheLp = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheLp,sidp,iMp);
                        cxx = sqrt(cxx);
                        [cxx0,f] = mscohere(v_b1', v0_b2',hamming(round(Fs/5)),[],[],Fs);
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
                        set(gca,'TickDir','out');
                        
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
                        set(gca,'TickDir','out');
                        
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
                    
                    Rc = R(countp,:);
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
                    
                    for i5 = 1:length(metricsp)
                        
                        metricp = metricsp{i5};
                        subplot(ceil(length(metricsp)/2),2,i5)
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
                    N_9d1 = [1 N_9d1];
                end
                for i_9d1 = N_9d1
                    for iMp = 1:length(metricsp)

                        metricp = metricsp{iMp};
                        
                        if (i_9d1 == 1)
                            t_plot_hrs = n_graph * w / 3600;
                            %t_plot_hrs = 24;%n_graph * w / 3600; %doesnt work
                            Tmin_fac = (1/60);
                        elseif (i_9d1 == 2)
                            t_plot_hrs = 24; % default - 1 hr
                            Tmin_fac = (1/60); % default - 1
                        end

                        % Reload coh_thresh
                        fn_cacheLp = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheLp,sidp,iMp);
                        load(fn_cacheLp,'coh_thresh');
                        
                        % Load artifacts
                        Risart = (art_idx(countp,:) ~= 0);
                        
                        % Load graph
                        fn_graphLp = sprintf('%s/%s_graph-%s.h5',dir_resLp,sidp,metricp);
                        Rc = h5read(fn_graphLp,'/R',[countp 1],[1 n_graph]);
                        
                        % Sample null
                        fn_distLp = sprintf('%s/%s_dists-%s-%i.mat',dir_resLp,sidp,metricp,n_perm);
                        D = load(fn_distLp);
                        Rn =random(D.d{countp},[1 n_graph]);
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
                            fprintf('[*] figure t1 plot trigger detected in figure t9 plot\n')
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
                        
                        
                        % Calculate midnight --------------------------
                        n_samples_f = h5readatt(ecog.filename,'/h5eeg','n_samples');
                        starts_f = cumsum([1;n_samples_f(1:(end))]);
                        realTime_f = cell(length(starts_f),1);
                        starts_g = ceil(starts_f/(round(ecog.fs)*w));
                        mnight_g = [];
                        parfor i_f = 1:(length(starts_f)-1)
                            dstr = ecog.getTimeFromSample(dir_stamp,starts_f(i_f));
                            sat = true;
                            nday = 1;
                            while(sat)
                                mnight_offset_g = seconds((datetime(dstr(1:10),'format','MM/dd/uuuu') + days(nday)) - datetime(dstr,'format','MM/dd/uuuu HH:mm:ss'))/w;
                                mng = (starts_g(i_f)+mnight_offset_g);
                                if (mng < starts_g(i_f+1))
                                    mnight_g = [mnight_g; round(mng)];
                                    nday = nday + 1;
                                else
                                    break;
                                end
                            end
                            realTime_f{i_f} = datetime(dstr,'format','MM/dd/uuuu HH:mm:ss');
                        end
                        
                        

                        
                        % Initialize figure
                        h = figure('visible','off');
                        if (i_9d1 == 1) % all
                            msize_main = 3;
                            msize_mark = 3;
                            set(h,'Position',[0 0 1*1080 0.5*1080])
                        elseif (i_9d1 == 2) % 1 hour
                            msize_main = 5;
                            msize_mark = 5;
                            %set(h,'Position',[0 0 0.5*1080 0.5*1080])
                            set(h,'Position',[0 0 1*1080 0.5*1080])
                        end
                        
                        color_div = 0.33*[1 1 1];
                        subplot(2,1,1)
                        plot(T_min,R_plt,'black.','MarkerSize',msize_main);hold on;
                        plot([min(T_min) max(T_min)],[coh_thresh(countp) coh_thresh(countp)],'--','color',color_div)
                        R_plt_max = nanmax(R_plt);
                        R_plt_min = nanmin(R_plt);
                        if (i_9d1 == 1)
                            % When plotting everything, simply index by
                            % r_idx
                            try
                                plot(T_min(r_idx),R_plt(r_idx),'redo','MarkerSize',msize_mark); hold all;
                                
                                % plot time markers
                                starts_g = starts_g(1:(end-1));
                                if (trig_plot_filebreak_x)
                                    plot(T_min(starts_g),R_plt_max*ones(size(starts_g)),'blackx');
                                end
                                % midnight
                                %return
                                %plot(T_min(mnight_g),R_plt_max*ones(size(mnight_g)),'redx');
                                for itmp = 1:length(mnight_g)
                                    ttmp = T_min(mnight_g(itmp));
                                    plot([ttmp ttmp],[0 R_plt_max],':','Color',COLOR_MIDNIGHT_LINE);
                                end
                                
                                
                                xlabel('Time (hours)')
                                cax_max = max(R_plt);
                            catch
                                fprintf('[*] figure t9 did not plot red circle.\n');
                                cax_max = 1;
                            end
                            cax_max_sav = cax_max;
                        elseif(i_9d1 == 2)
                            % When plotting just 1 hr, use index offset 
                            
                            if (offset_fits)
                                t_plot_off = t_plot_off - round(t9d1_offset_minutes * 60 / w);
                                plot(T_min(t_plot_off+1),R_plt(t_plot_off+1),'redo','MarkerSize',msize_mark); hold all;
                                
%                                 % plot time markers
%                                 plot(starts_g,ones(size(starts_g)),'black.');
%                                 plot(mnight_g,2*ones(size(mnight_g)),'red.');
                                
                                if (Tmin_fac == 1)
                                    xlabel('Time (minutes)');
                                elseif (Tmin_fac == (1/60))
                                    xlabel('Time (hours)');
                                    if (mod(t_plot_hrs,12) == 0)
                                        xticks(linspace(0,t_plot_hrs,9));
                                    end
                                else
                                    xlabel('Time');
                                end
                                
                                
                                %cax_max = double(max(R_plt(t_plot_off+1)));
                                if (exist('cax_max_sav','var'))
                                    cax_max = cax_max_sav;
                                else
                                    cax_max = max(R_plt);
                                end
                                %axis([min(T_min) max(T_min) min(Rn_plt) max(R_plt(t_plot_off+1))]); 
                                
                                % t9d1 print CT of segment in name
                                %return
                                Rtmp = R_plt(~isnan(R_plt));
                                ct_t9d1 = sum(Rtmp > coh_thresh(countp))/length(Rtmp);
                                n_t9d1 = length(Rtmp);
                                pname = sprintf('%s_9d1-ct%i_n%i',pname,round(1000*ct_t9d1),n_t9d1);
                            end
                        end
                        
                        % round up cax_max to nearest 0.1
                        cax_nearest = 0.1;
                        cax_max = ceil(cax_max/cax_nearest)*cax_nearest;
                        
                        %cax_max = cax_max + 0.1;
                        if (cax_max > 1)
                            cax_max = 1;
                        elseif (cax_max <= 0)
                            cax_max = 0;
                        end
                        if (min(Rn_plt) > 0)
                            if (min(Rn_plt) < cax_max)
                                axis([min(T_min) max(T_min) 0 cax_max]);
                            else
                                axis([min(T_min) max(T_min) 0 cax_max]);
                            end
                        else
                            axis([min(T_min) max(T_min) 0 cax_max]); 
                        end
                        
                        set(gca,'TickDir','out')
                        
                        
                        %ylabel({sprintf('%s Coherence',metricp(3:end));sprintf('%s-%s',roi_1,roi_2)})
                        ylabel({sprintf('%s Coherence',metricp(3:end))})
                        box off; 
                        if (i_9d1 == 1)
                            % When plotting everything, tick every 24 hours
                            xticks(0:24:max(T_min));
                        end

                        
                        subplot(2,1,2)
                        plot(T_min,Rn_plt,'black.','MarkerSize',msize_main); hold on;
                        plot([min(T_min) max(T_min)],[coh_thresh(countp) coh_thresh(countp)],'--','color',color_div);
                        if ((i_9d1 == 2) && (offset_fits))
                        %if ((i_9d1 == 1) || ((i_9d1 == 2) && (offset_fits)))
                            plot(T_min(t_plot_off+1),Rn_plt(t_plot_off+1),'redo','MarkerSize',msize_mark); hold all;
                        elseif (i_9d1 == 1)
                            plot(T_min(r_idx),Rn_plt(r_idx),'redo','MarkerSize',msize_mark); hold all;
                        end
%                         if (offset_fit)
%                             plot(T_min(cond_offset+1),Rn_plt(cond_offset+1),'blacko','MarkerSize',10);
%                         end
                        if (cax_max > min(Rn_plt))
                            axis([min(T_min) max(T_min) 0 cax_max]);
                        else
                            axis([min(T_min) max(T_min) 0 cax_max]);
                        end
                        if (i_9d1 == 1) % all
                            xlabel('Time (hours)')
                        elseif (i_9d1 == 2) % 1 hour
                            xlabel('Time (minutes)')
                        end
                        %ylabel({sprintf('Time-shifted %s Coherence',metricp(3:end));sprintf('%s-%s',roi_1,roi_2)})
                        ylabel({sprintf('%s Coherence',metricp(3:end))})
                        set(gca,'TickDir','out')
                        box off;
                        if (i_9d1 == 1)
                            % When plotting everything, tick every 24 hours
                            xticks(0:24:max(T_min));
                        end
                        
%                         if (i_9d1 == 1)
%                             return
%                         end

                        ct_plt = sum(R_plt(~Risa_plt) > coh_thresh(countp)) / length(R_plt(~Risa_plt));
                        ct0_plt = sum(Rn_plt(~Rnisa_plt) > coh_thresh(countp)) / length(Rn_plt(~Rnisa_plt));
                        
                        if (i_9d1 == 1)
                            dir_t9 = 'T9';
                            
                            n_days = (n_graph * w) / (3600*24);
                            fprintf('Total number of days: %.3f\n',n_days)
                            fprintf('Number of positive coherences: %i\n',round(ct_plt*n_graph))
                            %return
                        elseif (i_9d1 == 2)
                            dir_t9 = 'T9d1';
                        end
                        fn_pname = sprintf('figures/%s/%s_%s_T%i_ctN%i_fart%i',...
                            dir_t9,pname,metricp(3:4),round(1000*coh_thresh(countp)),round(1000*ct0_plt),round(1000*f_art));
                        print(h,fn_pname,fig_fmt);
                        if (trig_eps)
                            print(h,fn_pname,'-depsc');
                        end
                        close(h);
                    end
                end
                
               
                
            end
            countp = countp + 1;
        end
    end
    
end
end

% Clear loop indices
clear i;
clear j;

% Print finish message
fprintf('Done.\n')
