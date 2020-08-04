close all;
clear;
rng shuffle;

metricp = 'pcBroadband';
n_perm = 10000;
perm_alpha_sec = 12;
cp_thresh_override = 0.05;
p_val = 0.01;

% Fast i/o definitions
%dir_artLp = '/media/klab/KLAB101/h5_notch20/art_nosz';
%dir_resLp = '/media/klab/KLAB101/results/coh_w10'; %dir_resLp = '/media/klab/internal/data/results/coh_w10';
%dir_corLp = '/media/klab/internal/data/coreg';
dir_cacheLp = './cache';
%dir_stamp = '/media/klab/internal/data/stamps';
%dir_video = '/media/klab/internal/data/videos';
%subjects_dirLp = '/mnt/cuenap_ssd/coregistration';

% Slow i/o definitions
%dir_h5Lp = '/media/jerry/KLAB101/h5_notch20';
%dir_h5Lp = '/media/jerry/internal/data/h5_notch20';

[~,hname] = system('hostname');

% Paths
if contains(hname,'kraken')
    dir_artLp = '/media/klab/internal/data/h5_notch20/art_nosz';
    dir_h5Lp = '/mnt/cuenap/data/h5_notch20';
elseif contains(hname,'hopperu')
    dir_artLp = '/mnt/cuenap/data/h5_notch20/art_nosz';
    dir_h5Lp = '/mnt/cuenap/data/h5_notch20';
end

%metricsp = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};
metricsp = {'pcBroadband','pcGamma'};

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

% Artifact definitions from v2/
x1_t1 = 10; % uV
x1_t2 = 2000; % uV
x2_t = 400; % uV


% Show art stats
fn_ca = sprintf('%s/%s',dir_cacheLp,'figure_t4_artstat.mat');
if (exist(fn_ca,'file'))
    fprintf('[*] Found cached art stats.\n')
    Ca = load(fn_ca);
    Frac = Ca.Frac;
    FracN = Ca.FracN;
    art_uids = Ca.art_uids;
else
    art_uids = 0:6; %unique(art(:));
    Frac = zeros(length(art_uids),length(Subjectsp));
    FracN = Frac;
    for i = 1:length(Subjectsp)
        sid = Subjectsp{i};
        fn_art = sprintf('%s/%s_art.h5',dir_artLp,sid);
        fn_h5 = sprintf('%s/%s.h5',dir_h5Lp,sid);
        ecog = H5eeg(fn_h5);
        art = h5read(fn_art,'/artifacts');
        [n_artbchan,n_art] = size(art);

        % Report art fractions

        for j2 = 1:length(art_uids)
            frac = sum(sum(art == art_uids(j2))) / (n_art * n_artbchan);
            fprintf('%s\t%i\t%.3f%%\n',sid,art_uids(j2),100*frac);
            Frac(j2,i) = frac;
            FracN(j2,i) = n_art*n_artbchan;
        end
    end
    save(fn_ca,'FracN','Frac','art_uids');
end

FracW = FracN / sum(FracN(1,:));
Frac_art = sum(FracW.*Frac,2);
Frac_art_N = sum(FracN(1,:));

fprintf('[!] Total artifacts:\n')
ArtKey = {'No Art','LoSNR','HiMag','HiDiff','Elec','Time','SzStim'};

sIdx = [1 3 4 2 5 6 7];
ArtKey = ArtKey(sIdx);
Frac_art = Frac_art(sIdx);
Frac = Frac(sIdx,:);
FracW = FracW(sIdx,:);

for i = 1:length(Frac_art)
    fprintf('%s\t\t%.3f%%\t\t%i\n',ArtKey{i},100*Frac_art(i),round(Frac_art(i)*Frac_art_N))
end
fprintf('Total number of segments: %i (%i patient-electrode-days)\n',Frac_art_N, round(Frac_art_N / (24*3600)))
fprintf('Total artifacts only: %.2f%%\n',100*sum(Frac_art(2:6)));

fprintf('\n');

% Variation across subjects
Frac_mean = mean(Frac,2);
Frac_std = zeros(length(Frac_art),1);
for j = 1:length(Frac_art)
    Frac_std(j) = std(Frac(j,:));
end
Frac_median = median(Frac,2);
Frac_min = min(Frac')';
Frac_max = max(Frac')';
%tot = 0;
for i = 1:length(Frac_art)
    fprintf('%s\tmean: %.6f, std: %.6f, max: %.6f\n',ArtKey{i},Frac_mean(i),Frac_std(i),Frac_max(i));
    
    %tot = tot + Frac_mean(i);
end
tot = sum(Frac(2:6,:));
fprintf('Total artifacts mean: %.2f%%\n',100*mean(sum(Frac(2:6,:))));
fprintf('Total artifacts std: %.2f%%\n',100*std(tot));
fprintf('Total artifacts max: %.2f%%\n',100*max(tot));


%%
% Plot t4

% voltage threshold for showing saturation
vrange = 800; % 2.5 * 1e3;

% Find best artifact
n_samp = 2;
Artids = 1:3;

sid_const = 'm00083';
Bchans_i = NaN; 
Samps_j = NaN; 

Bchans_i = [65       69       96        148       37       56       ]; 
Samps_j =  [1298945  2175233  86966017  10280961  2123777  64310017 ]; 

% --- June 18, 2019 ---
% Artifact 1 - small magnitude
% 70       69       65        
% 1292545  2175233  1298945  
% Artifact 2 - large magnitude
% 115       28       13       78       13       153       64
% 89255937  3256321  2389249  2178049  1071105  25284097  908801
% Artifact 3 - large difference
% 37       108      120      71       12      
% 2123777  1241089  1240321  1234945  839937  

%Bchans_i = [105,12,15];
%Samps_j = [95708673,53505,1240321];

% top 2
%m00083__x-7_bchan-105_samp-95708673_x-12648_bchan-12_samp-53505_x-1391_bchan-45_samp-1240065
% last
%m00083__x-6_bchan-105_samp-90326273_x-11930_bchan-62_samp-883201_x-1311_bchan-15_samp-1240321

chosen_already = all(~isnan([Bchans_i, Samps_j]));
if (chosen_already)
    n_samp = 1;
end

for i = 1:1 %length(Subjectsp)
    %sid = Subjectsp{i};
    sid = sid_const; %'m00083';
    fn_art = sprintf('%s/%s_art.h5',dir_artLp,sid);
    fn_h5 = sprintf('%s/%s.h5',dir_h5Lp,sid);
    Ca = load(sprintf('cache/xsub_out_%s_1.mat',sid));
    ecog = Ca.ecog;
    %ecog = H5eeg(fn_h5);
    art = h5read(fn_art,'/artifacts');
    [n_artbchan,n_art] = size(art);
    
%     % Report art fractions
%     art_uids = unique(art(:));
%     for j2 = 1:length(art_uids)
%         frac = sum(sum(art == art_uids(j2))) / (n_art * n_artbchan);
%         fprintf('\t%i\t%.3f%%\n',art_uids(j2),100*frac);
%     end
%     
%     return
    
    h = figure('visible','off');
    set(h,'PaperUnits','inches');
    set(h,'PaperPosition',[0 0 8.5 11]);
    ofstr = [sid,'_'];
    for j = 1:(2*length(Artids))
        art_id = Artids(ceil(j/2));
        Xs = zeros(n_samp,1);
        Bchans = zeros(n_samp,1);
        Samps = zeros(n_samp,1);
        Vs = cell(n_samp,1);
        Vfs = cell(n_samp,1);
        As = cell(n_samp,1);
        
        % extract all artifacts for art_id
        [n_i,n_j] = size(art);
        Is = [];
        Js = [];
        parfor i1 = 1:n_i
            for j1 = 1:n_j
                if (art_id == art(i1,j1))
                    Is = [Is; i1];
                    Js = [Js; j1];
                end
            end
        end
        pIdx = randperm(length(Is));
        Is = Is(pIdx);
        Js = Js(pIdx);
        
        for k = 1:n_samp
        %for k = 1
            % row number
            if (chosen_already)
                idx_i = Bchans_i(j);
            else
                idx_i = Is(k); %randi([1 ecog.n_bchan]);
            end
            
            % artifact index
            if (chosen_already)
                idx_j = Samps_j(j);
            else
                idx_j = Js(k); %randi([1 n_art]);
            end
            
            if (~chosen_already)
                % convert artifact index to real time index
                idx_ja = idx_j;
                idx_j = ((idx_j - 1)*round(ecog.fs)) + 1;
            else
                idx_ja = ((idx_j - 1)/round(ecog.fs)) + 1;
            end
            
            % read v
            b1c1 = ecog.bip(idx_i,1);
            b1c2 = ecog.bip(idx_i,2);
            v1 = h5read(fn_h5,'/h5eeg/eeg',[b1c1 idx_j],[1 round(ecog.fs)]);
            v2 = h5read(fn_h5,'/h5eeg/eeg',[b1c2 idx_j],[1 round(ecog.fs)]);
            v = v1 - v2;
            v = v - mean(v);
            
            % full v
            sec_offset = 5;
            sec_length = 10;
            vf1 = h5read(fn_h5,'/h5eeg/eeg',[b1c1 idx_j-sec_offset*round(ecog.fs)],[1 sec_length*round(ecog.fs)]);
            vf2 = h5read(fn_h5,'/h5eeg/eeg',[b1c2 idx_j-sec_offset*round(ecog.fs)],[1 sec_length*round(ecog.fs)]);
            vf = vf1 - vf2;
            vf = vf - mean(vf);
            
            % calculate features
            if (art_id < 3)
                x = max(v) - min(v);
            else
                x = max(abs(diff(v)));
            end
            
            a = h5read(fn_art,'/artifacts',[idx_i,idx_ja-sec_offset],[1 sec_length]);
            
            % store
            Xs(k) = x;
            Bchans(k) = idx_i;
            Samps(k) = idx_j;
            Vs{k} = v;
            Vfs{k} = vf;
            As{k} = a;
            
        end
        
        
        if (art_id > 1)
            [~,mIdx] = max(Xs);
        else
            [~,mIdx] = min(Xs);
        end
        
        
        subplot(2*length(Artids),1,j);
        T = linspace(0,sec_length,length(Vfs{mIdx}));
        V = Vfs{mIdx};
        plot(T,V,'black-','LineWidth',1); hold all;
        %tIdx = false(1,length(T));
        T2 = T;
        V2 = V;
        for ii = 1:length(T2)
            if (( As{mIdx}( floor((ii-1)/(round(ecog.fs))) + 1 ) ) > 0)
                V2(ii) = V(ii);
            else
                V2(ii) = NaN;
            end
        end
        %tIdx = (T >= sec_offset) & (T < (sec_offset + 1));
        %plot(T(tIdx),Vfs{mIdx}(tIdx),'black-','LineWidth',2);
        
        % Artifact line spec
        plot(T2,V2,'-','LineWidth',1,'Color',[0 0 1]);
        
        % plot saturated points
        idx_satHi = (V>=vrange);
        idx_satLo = (V<=(-1*vrange));
        V_SATH = vrange*ones(size(T));
        V_SATL = (-1)*vrange*ones(size(T));
        plot(T(idx_satHi),V_SATH(idx_satHi),'red.');
        plot(T(idx_satLo),V_SATL(idx_satLo),'red.');
        
        set(gca,'box','off');
        set(gca,'TickDir','out');
        ylabel(sprintf('IFP (\\muV)'));
        if (j == (2*length(Artids)))
            xlabel('Time (s)');
        end
        axis([0 sec_length (-1)*vrange vrange]);
        
        %ofstr = [ofstr,sprintf('_x-%i_bchan-%i_samp-%i',round(Xs(k)),Bchans(k),Samps(k))];
        ofstr = [ofstr,sprintf('_x-%i_bchan-%i_samp-%i',round(Xs(mIdx)),Bchans(mIdx),Samps(mIdx))];% mIdx
        
    end
    
    print(h,sprintf('figures/T4/%s',ofstr),'-dpng');
    print(h,sprintf('figures/T4/%s',ofstr),'-depsc');
end


