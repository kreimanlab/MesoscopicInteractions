close all;
clear;

metric = 'pcBroadband';
n_perm = 10000;

% Fast i/o definitions
dir_art = '/media/klab/internal/data/h5_notch20/art';
dir_res = '/media/klab/internal/data/results/coh_w10';
dir_cor = '/media/klab/internal/data/coreg';

% Slow i/o definitions
dir_h5 = '/media/klab/KLAB101/h5_notch20';

% Patients
Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124',...
    'mSu'};

% Exclude monkey
%Subjects = Subjects(1:(end-1));

% Main loop
for i = 1:length(Subjects)
    sid = Subjects{i};
    fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
    fn_dist = sprintf('%s/%s_dists-%s-%i.mat',dir_res,sid,metric,n_perm);
    fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
    fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
    fn_coreg = sprintf('%s/%s/label/all_parcellation.mat',dir_cor,sid);
%     
%     % Check if files exist
%     ckf = {fn_art,fn_dist,fn_graph,fn_h5,fn_coreg};
%     for j = 1:length(ckf)
%         if (~exist(ckf{j},'file'))
%             fprintf(2,'E> File not found: %s\n',ckf{j});
%             return
%         end
%     end
    
    frac_art = h5readatt(fn_art,'/art_idx','frac_art');
    n_arts = h5readatt(fn_art,'/artifacts','n_arts');
    n_bchan = h5readatt(fn_art,'/artifacts','n_bchan');
    h5i = h5info(fn_art);
    h5i = {h5i.Datasets.Dataspace};
    h5i = h5i{1}.Size;
    fprintf('%i\t%s\tn_bchan: %i\tn_comb: %i\tn_times: %i\tartifacts: %.2f %%\n',i,sid,n_bchan,h5i(1),h5i(2),100*frac_art)
    %return
end
