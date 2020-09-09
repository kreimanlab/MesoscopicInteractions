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
% 	total: 37
% ------------------------------------------------
% (*) xsub method rois: 54, 30
% (*) post-clustering rois: 36, 45
% (*) coherence: 0.000000
% 	n_pairs: 13
% 	n_usub: 4
% [*] === Left-Right lateralization analysis === Frequency band: 1
% 	Subject: sub1


sid_int_const = 5; 
sid_const = 'sub3';
bchan1_const = 35;
bchan2_const = 84;

Ca = load(sprintf('cache/xsub_out_%s'));