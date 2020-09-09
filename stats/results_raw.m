close all;
clear;

% === Figure 2 analysis === frequency band: 1
% [!] Found region: 110 of 150 for bchan2
% 	post-cluster region: 30
% [!] Found region: 114 of 150 for bchan1
% 	post-cluster region: 54
% [*] Mean coherence between regions: 0.2231
% [*] Number of pairs: 32.0000
% [*] Number of usubs: 37.0000
% [*] Consistency across subjects (pairs): 0.0625 (2.00)
% [*] Consistency across subjects (usubs): 0.0270 (1.00)
% [*] Average distance: 19.3258 mm
% [!] Checking the number of subjects..

Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

dir_cache = './cache';
subjects_dir = '/media/klab/internal/data/coreg';

for iM = [1]

    Fsig = [];
    Nbip = nan(length(Subjects),1);
    Ncov = nan(length(Subjects),1);
    Nsig = nan(length(Subjects),1);
    Hours = nan(length(Subjects),1);
    for iSub = 1:length(Subjects)
        sid = Subjects{iSub};
        Ca = load(sprintf('%s/xsub_out_%s_%i_atl2.mat',dir_cache,sid,iM));
        Adj = Ca.AdjMag;
        n_Adj = length(Adj);
        A = [];
        for i = 1:(n_Adj-1)
            for j = (i+1):n_Adj
                A = [A;Adj(i,j)];
            end
        end
        
        % Counts
        Nbip(iSub) = length(A);
        Ncov(iSub) = sum(~isnan(A));
        Nsig(iSub) = sum(A>0);
        Hours(iSub) = (Ca.ecog.n_samples / (Ca.ecog.fs)) * (1/3600);
        
        f_sig = sum((~isnan(A)) & ((A~=0))) / sum(~isnan(A));
        Fsig(iSub) = f_sig;
    end
end