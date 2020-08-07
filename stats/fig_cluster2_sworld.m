close all;
clear;
% The Watts-strogatz definition for small world network is:
%   L >= L_random but C >> C_random
%           ws    - Watts-Strogatz 1998
%           sigma - Humphries-Gurney,2008 (sigma > 1 for small world network)
%           omega - Telesford-Laurienti,2011 (0 to 1: 1 is most small world)
%           I     - Neal,2017 (0 to 1: 1 is most small world)

n_MC = 12;

% Load human adjacency matrix
iM = 1;
fn_ca3 = sprintf('./cache/fig_cluster2_reduce_%i_new.mat',iM);
%fn_ca3 = './cache/fig_cluster2_fill_1.mat';
Ca3 = load(fn_ca3);
Ahs = Ca3.A; %Ca3.A2;
n_Ahs = length(Ahs);


fprintf('[*] Computing small world index for humans in 150 parcellation..\n')
[ws,sigma,omega,I] = swi(Ahs,n_MC);
fprintf('\tn = %i\n',n_Ahs);
fprintf('\t%s\n',ws)
fprintf('\tsigma = %.6f\n',sigma);
fprintf('\tomega = %.6f\n',omega);
fprintf('\tSWI = %.6f\n',I);

% Load individual adjacency matrices
Ca = load('./cache/xsub_out_all_1.mat');
iM = 1;
for i = 1:length(Ca.Subjects)
    sid = Ca.Subjects{i};
    Cas = load(sprintf('./cache/xsub_out_%s_%i.mat',sid,iM));
end

% fprintf('[*] Calculating Watts-Strogatz metrics for H sapiens..\n')
% [Lhs,Chs] = watts_strogatz(Ahs);
% fprintf('\tn = %i\n',length(Ahs));
% fprintf('\tL = %.6f\n\tC = %.6f\n',Lhs,Chs);
% 
% % Compute random graph metrics
% fprintf('[*] Computing random graph metrics..\n')
% n_MC = 12;
% [Lhs_rand,Chs_rand] = watts_strogatz_random(Ahs, n_MC);
% 
% p_Lhs = 1 - sum((Lhs > sort(Lhs_rand)) | (Lhs < sort(Lhs_rand)))/length(Lhs_rand);
% p_Chs = 1 - sum((Chs > sort(Chs_rand)) | (Chs < sort(Chs_rand)))/length(Chs_rand);
% fprintf('\tL_rand = %.6f +- %.6f\n\tC_rand = %.6f +- %.6f\n',...
%     mean(Lhs_rand),std(Lhs_rand),mean(Chs_rand),std(Chs_rand));
% fprintf('[*] Permutation test with n: %i, 2-sided p_L = %.3d, p_C = %.3d\n',...
%     n_MC,p_Lhs,p_Chs);
% 
% % small world metric
% C = Chs;
% Cr = mean(Chs_rand);
% L = Lhs;
% Lr = mean(Lhs_rand);
% sigma = (C/Cr)/(L/Lr); % if sigma > 1 then small world
% fprintf('\tsigma: %.6f\n',sigma);






% Load c elegans adjacency matrix
% L.R. Varshney, B.L. Chen, E. Paniagua, D.H. Hall and D.B. Chklovskii (2011)
fn_ce = './connectome/NeuronConnect.csv';
f_ce = fopen(fn_ce,'r');
Dce = textscan(f_ce,'%s','Delimiter','\n');
Dce = Dce{1};
fclose(f_ce);
Ece = {};
j = 1;
for i = 2:length(Dce)
    line = strsplit(Dce{i},',');
    if (~strcmp(line{3},'NMJ'))
        Ece{j,1} = line{1};
        Ece{j,2} = line{2};
        Ece{j,3} = line{3};
        Ece{j,4} = str2double(line{4});
        j = j + 1;
    end
end
Ace_labels = unique({Ece{:,1}; Ece{:,2}});
n_Ace = length(Ace_labels);
Ace = zeros(n_Ace,n_Ace);
% build adjacency from edges
[n_Ece,~] = size(Ece);
for i = 1:n_Ece
    neur1 = Ece{i,1};
    neur2 = Ece{i,2};
    idx_neur1 = find(strcmp(Ace_labels,neur1),1);
    idx_neur2 = find(strcmp(Ace_labels,neur2),1);
    Ace(idx_neur1,idx_neur2) = Ece{i,4};
    Ace(idx_neur2,idx_neur1) = Ece{i,4};
end

fprintf('[*] Computing small world index for C elegans..\n')
[ws,sigma,omega,I] = swi(Ace,n_MC);
fprintf('\tn = %i\n',n_Ace);
fprintf('\t%s\n',ws)
fprintf('\tsigma = %.6f\n',sigma);
fprintf('\tomega = %.6f\n',omega);
fprintf('\tSWI = %.6f\n',I);


% fprintf('\n[*] Calculating Watts-Strogatz metrics for C elegans..\n')
% [Lce,Cce] = watts_strogatz(Ace);
% fprintf('\tn = %i\n',length(Ace));
% fprintf('\tL = %.6f\n\tC = %.6f\n',Lce,Cce);
% 
% % Compute random graph metrics
% fprintf('[*] Computing random graph metrics..\n')
% %n_MC = 12;
% [Lce_rand,Cce_rand] = watts_strogatz_random(Ahs, n_MC);
% 
% p_Lce = 1 - sum((Lce > sort(Lce_rand)) | (Lce < sort(Lce_rand)))/length(Lce_rand);
% p_Cce = 1 - sum((Cce > sort(Cce_rand)) | (Cce < sort(Cce_rand)))/length(Cce_rand);
% fprintf('\tL_rand = %.6f +- %.6f\n\tC_rand = %.6f +- %.6f\n',...
%     mean(Lce_rand),std(Lce_rand),mean(Cce_rand),std(Cce_rand));
% fprintf('[*] Permutation test with n: %i, 2-sided p_L = %.3d, p_C = %.3d\n',...
%     n_MC,p_Lce,p_Cce);
% 
% % small world metric
% C = Cce;
% Cr = mean(Cce_rand);
% L = Lce;
% Lr = mean(Lce_rand);
% sigma = (C/Cr)/(L/Lr); % if sigma > 1 then small world
% fprintf('\tsigma: %.6f\n',sigma);





% Macaque
D_anat = load('../coreg/AdjMKFV.mat');
D = load('../coreg/AdjSu.mat');

% convert to symmetric
AdjMKneurons = D_anat.AdjMKneurons;
AdjMKneurons = AdjMKneurons + AdjMKneurons';
AdjMKneuronst = AdjMKneurons;
AdjMKneuronst(AdjMKneuronst == 0) = NaN;
AdjMK = log10(AdjMKneuronst);
AdjMK_bin = ~isnan(AdjMK);
AdjMK_bin = double(AdjMK_bin);

Amaca = AdjMK;
Amacf = D.AdjMag;
isn = all(isnan(Amacf));
Amacf(isn,:) = [];
Amacf(:,isn) = [];
Amacf(isnan(Amacf)) = 0;

fprintf('[*] Computing small world index for Markov-Kennedy..\n')
[ws,sigma,omega,I] = swi(Amaca,n_MC);
fprintf('\tn = %i\n',length(Amaca));
fprintf('\t%s\n',ws);
fprintf('\tsigma = %.6f\n',sigma);
fprintf('\tomega = %.6f\n',omega);
fprintf('\tSWI = %.6f\n',I);

fprintf('[*] Computing small world index for macaque functional interactions..\n')
[ws,sigma,omega,I] = swi(Amacf,n_MC);
fprintf('\tn = %i\n',length(Amacf));
fprintf('\t%s\n',ws);
fprintf('\tsigma = %.6f\n',sigma);
fprintf('\tomega = %.6f\n',omega);
fprintf('\tSWI = %.6f\n',I);

% 
% %Compute
% fprintf('\n[*] Calculating Watts-Strogatz metrics for Markov-Kennedy..\n')
% [Lmaca,Cmaca] = watts_strogatz(Amaca);
% fprintf('\tn = %i\n',length(Amaca));
% fprintf('\tL = %.6f\n\tC = %.6f\n',Lmaca,Cmaca);
% 
% % Compute random graph metrics
% fprintf('[*] Computing random graph metrics..\n')
% %n_MC = 12;
% [Lmaca_rand,Cmaca_rand] = watts_strogatz_random(Amaca, n_MC);
% 
% p_Lmaca = 1 - sum((Lmaca > sort(Lmaca_rand)) | (Lmaca < sort(Lmaca_rand)))/length(Lmaca_rand);
% p_Cmaca = 1 - sum((Cmaca > sort(Cmaca_rand)) | (Cmaca < sort(Cmaca_rand)))/length(Cmaca_rand);
% fprintf('\tL_rand = %.6f +- %.6f\n\tC_rand = %.6f +- %.6f\n',...
%     mean(Lmaca_rand),std(Lmaca_rand),mean(Cmaca_rand),std(Cmaca_rand));
% fprintf('[*] Permutation test with n: %i, 2-sided p_L = %.3d, p_C = %.3d\n',...
%     n_MC,p_Lmaca,p_Cmaca);
% 
% % small world metric
% C = Cmaca;
% Cr = mean(Cmaca_rand);
% L = Lmaca;
% Lr = mean(Lmaca_rand);
% sigma = (C/Cr)/(L/Lr); % if sigma > 1 then small world
% fprintf('\tsigma: %.6f\n',sigma);
% 
% 
% 
% %Compute
% fprintf('\n[*] Calculating Watts-Strogatz metrics for mSu..\n')
% [Lmacf,Cmacf] = watts_strogatz(Amacf);
% fprintf('\tn = %i\n',length(Amacf));
% fprintf('\tL = %.6f\n\tC = %.6f\n',Lmacf,Cmacf);
% 
% % Compute random graph metrics
% fprintf('[*] Computing random graph metrics..\n')
% %n_MC = 12;
% [Lmacf_rand,Cmacf_rand] = watts_strogatz_random(Amacf, n_MC);
% 
% p_Lmacf = 1 - sum((Lmacf > sort(Lmacf_rand)) | (Lmacf < sort(Lmacf_rand)))/length(Lmacf_rand);
% p_Cmacf = 1 - sum((Cmacf > sort(Cmacf_rand)) | (Cmacf < sort(Cmacf_rand)))/length(Cmacf_rand);
% fprintf('\tL_rand = %.6f +- %.6f\n\tC_rand = %.6f +- %.6f\n',...
%     mean(Lmacf_rand),std(Lmacf_rand),mean(Cmacf_rand),std(Cmacf_rand));
% fprintf('[*] Permutation test with n: %i, 2-sided p_L = %.3d, p_C = %.3d\n',...
%     n_MC,p_Lmacf,p_Cmacf);
% 
% % small world metric
% C = Cmacf;
% Cr = mean(Cmacf_rand);
% L = Lmacf;
% Lr = mean(Lmacf_rand);
% sigma = (C/Cr)/(L/Lr); % if sigma > 1 then small world
% fprintf('\tsigma: %.6f\n',sigma);
% 
