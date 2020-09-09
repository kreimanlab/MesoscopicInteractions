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

metrics = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};

% Patients
Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48',...
    'mSu'};

% Exclude monkey
Subjects = Subjects(1:(end-1));

% Main loop
for i = 1:length(Subjects)
    sid = Subjects{i};
    fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
    fn_dist = sprintf('%s/%s_dists-%s-%i.mat',dir_res,sid,metric,n_perm);
    fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
    fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
    fn_coreg = sprintf('%s/%s/label/all_parcellation.mat',dir_cor,sid);
    
    % Check if files exist
    ckf = {fn_art,fn_dist,fn_graph,fn_h5,fn_coreg};
    for j = 1:length(ckf)
        if (~exist(ckf{j},'file'))
            fprintf(2,'E> File not found: %s\n',ckf{j});
            return
        end
    end
end

% Clear loop indices
clear i;
clear j;

% Print finish message
fprintf('Done.\n')