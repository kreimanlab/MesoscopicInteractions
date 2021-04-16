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

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

dir_cache = './cache';
subjects_dir = '/media/klab/internal/data/coreg';

for iM = [1]

    Fsig = [];
    Nbip = nan(length(Subjects),1);
    Ncov = nan(length(Subjects),1);
    Nsig = nan(length(Subjects),1);
    Ngraph = nan(length(Subjects),1);
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
        Ngraph(iSub) = Ca.n_graph;
        Hours(iSub) = (Ca.ecog.n_samples / (Ca.ecog.fs)) * (1/3600);
        
        f_sig = sum((~isnan(A)) & ((A~=0))) / sum(~isnan(A));
        Fsig(iSub) = f_sig;
    end
end