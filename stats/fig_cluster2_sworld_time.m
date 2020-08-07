close all;
clear;
% The Watts-strogatz definition for small world network is:
%   L >= L_random but C >> C_random
%           ws    - Watts-Strogatz 1998
%           sigma - Humphries-Gurney,2008 (sigma > 1 for small world network)
%           omega - Telesford-Laurienti,2011 (0 to 1: 1 is most small world)
%           I     - Neal,2017 (0 to 1: 1 is most small world)

n_MC = 2;

% Load human adjacency matrix
fn_ca3 = './cache/fig_cluster2_fill_1.mat';
dir_cache = './cache';
dir_res = '/home/jerry/data/results/coh_w10';
dir_art = '/media/jerry/internal/data/h5_notch20/art_nosz';
%Ca3 = load(fn_ca3);
%Ahs = Ca3.A2;

% Load 150-adjacency matrix
fn_caT14 = './cache/figure_t14_1_150.mat';
CaT14 = load(fn_caT14);
Ahs = CaT14.Adj_plt2;
    
n_Ahs = length(Ahs);



fprintf('[*] Computing small world index for humans..\n')
[ws,sigma,omega,I] = swi(Ahs,n_MC);
fprintf('\tn = %i\n',n_Ahs);
fprintf('\t%s\n',ws)
fprintf('\tsigma = %.6f\n',sigma);
fprintf('\tomega = %.6f\n',omega);
fprintf('\tSWI = %.6f\n',I);

% Load individual adjacency matrices
Ca = load('./cache/xsub_out_all_1.mat');
iM = 1;
Sigmas = [];
Omegas = [];
Is = [];
Sidints = [];
for i = 1:length(Ca.Subjects)
    sid = Ca.Subjects{i};
    sidint = str2double(sid(2:end));
    Cas = load(sprintf('./cache/xsub_out_%s_%i.mat',sid,iM));
    A = Cas.AdjMag;
    A(isnan(A)) = nanmean(A(:));
    
    %fprintf('[*] Computing small world index for %s..\n',sid)
    [ws,sigma,omega,I] = swi(A,n_MC);
    fprintf('%s\t%s\t%.6f\t%.6f\t%.6d\n',sid,ws,sigma,omega,I)
    Sigmas = [Sigmas; sigma];
    Omegas = [Omegas; omega];
    Is = [Is; I];
    Sidints = [Sidints; sidint];
end

n_range = 0;
ri = (1+n_range):(length(Sigmas)-n_range);
sigma = nanmean(Sigmas);
sigma_std = nanstd(Sigmas);
sigma_range = [min(Sigmas(ri)) max(Sigmas(ri))];
omega = nanmean(Omegas);
omega_std = nanstd(Omegas);
omega_range = [min(Omegas(ri)) max(Omegas(ri))];
Is(isinf(Is)) = NaN;
Is(Is>1) = NaN;
Is(Is<0) = NaN;
I = nanmean(Is);
I_std = nanstd(Is);
I_range = [min(Is(ri)) max(Is(ri))];
fprintf('\n')
fprintf('sigma:\t%.6f +- %.6f\t[%.6f ~ %.6f: %.2f%%]\n',sigma,sigma_std,sigma_range(1),sigma_range(2),100*(1-2*n_range/length(Sigmas)))
fprintf('omega:\t%.6f +- %.6f\t[%.6f ~ %.6f: %.2f%%]\n',omega,omega_std,omega_range(1),omega_range(2),100*(1-2*n_range/length(Omegas)))
fprintf('SWI:\t%.6f +- %.6f\t[%.6f ~ %.6f: %.2f%%]\n',I,I_std,I_range(1),I_range(2),100*(1-2*n_range/length(Is)));


%%
Subjects = Ca.Subjects;
SMI_all = {};
for i = 1:length(Ca.Subjects)
    sid = Subjects{i};
    Ca = load(sprintf('%s/xsub_out_%s_%i',dir_cache,sid,iM));
    R = h5read(sprintf('%s/%s_graph-%s.h5',dir_res,sid,Ca.metrics{iM}),'/R');
    art_idx = h5read(sprintf('%s/%s_art.h5',dir_art,sid),'/art_idx');
    R(art_idx>0) = NaN;
    [n_comb,n_graph] = size(R);
    for j = 1:n_comb
        R(j,:) = fillmissing(R(j,:),'linear');
    end

    SMIs = zeros(3,n_graph);
    for j = 1:n_graph
        n_bchan = Ca.ecog.n_bchan;
        A = zeros(n_bchan,n_bchan);
        kk = 1;
        for ii = 1:(n_bchan-1)
            for jj = (ii+1):n_bchan
                r = R(kk,j);
                if (r > Ca.coh_thresh(kk))
                    rr = r;
                else
                    rr = 0;
                end
                A(ii,jj) = rr;
                A(jj,ii) = rr;
                kk = kk + 1;
            end
        end
        
        % fill missing for A
        A(Ca.Dmat < 0) = nanmean(A(:));
        
        [ws,sigma,omega,I] = swi(A,n_MC);
        if (I > 1)
            I = 1;
        end
        if (I < 0)
            I = 0;
        end
        SMIs(1,j) = sigma;
        SMIs(2,j) = omega;
        SMIs(3,j) = I;
        
        fprintf('%s: %i of %i\n',sid,j,n_graph)
    end
    
    SMIs(any(art_idx)) = NaN;
    
    SMI_all{i} = SMIs;
    %SMI_all = [SMI_all; SMIs];
    %return
end

save('./cache/SMI_all','SMI_all');
fprintf('Done.\n')