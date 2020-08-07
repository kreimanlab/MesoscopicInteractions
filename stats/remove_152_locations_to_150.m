close all;
clear;

load('cache/fig_cluster2_reduce.mat');
n_loc_before = length(A);
rem = [5; 29];
A(rem,:) = [];
A(:,rem) = [];
n_loc_after = length(A);
Ad(rem,:) = [];
Ad(:,rem) = [];
E(rem,:) = [];
E_sphere(rem,:) = [];
Es(rem,:) = [];
Es_sphere(rem,:) = [];
% cc = 1;
% ss_samesid_old = ss_samesid;
% ss_samesid = [];
% for i = 1:(n_loc_before-1)
%     for j = (i+1):n_loc_before
%         if (~(any(i == rem) || any(j == rem)))
%             ss_samesid = [ss_samesid; ss_samesid_old(cc,:)];
%         end
%         cc = cc + 1;
%     end
% end
% clear cc;
% clear ss_samesid_old;
%return
save('cache/fig_cluster2_reduce_clean.mat');