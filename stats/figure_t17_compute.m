close all;
clear;
rng('shuffle');

%/media/klab/internal/data/coreg/fsaverage_sym/label/all_surf_ielvis_m00037.label
dir_artLp = '/media/klab/internal/data/h5_notch20/art';
dir_corLp = '/media/klab/internal/data/coreg';
%dir_resLp = '/home/jerry/data/results/coh_w10';
setenv('SUBJECTS_DIR',dir_corLp);
dir_cacheLp = './cache';
dir_h5Lp = '/media/klab/internal/data/h5_notch20';
mkdir('figures/T17');

metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

% Colors
BRAIN_MESH_ALPHA = 1;
COLOR_MEDIAL_WALL = 0.3*[1 1 1];
COLOR_NOSIG = 0.5*[1 1 1];
SURFACE_TYPE = 'pial';
paper_width = 8.5;
paper_height = 6.5;

% Options
trig_show_roi_label = false;
trig_blot_nosig = true;

% --- register t14 adjacency matrix to 150 parcellation -------------------
cluster_i = load('cache/fig_cluster3_cluster_i.mat');
cluster_i = cluster_i.cluster_i;
CaT14 = load(sprintf('./cache/figure_t14_%i_150',1));
cluster_i2 = zeros(size(cluster_i));
for i = 1:length(cluster_i)
    idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(i))));
    cluster_i2(i) = idx;
end
% cluster_i = cluster_i2;
% cluster_i2: converts Es indices to final adjacency matrix indices
% -------------------------------------------------------------------------


% load average surface
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_corLp,'fsaverage_sym','r',SURFACE_TYPE));

% load atlas definition
CaC1 = load('cache/fig_cluster2_reduce_new_annot_vertex_color2.mat');
CaC2 = load('cache/fig_cluster2_reduce_annot');


n_compute = 32;
% 'K','Ssumd','Res','A2km';
compute_K = cell(n_compute,1);
compute_Ssumd = cell(n_compute,1);
compute_Res = cell(n_compute,1);
compute_A2km = cell(n_compute,1);
for compute = 1:n_compute
    tic;
    
    for i = 1 %1:length(metricsp)
        metric = metricsp{i};

        % Load adjacency matrix
        CaT14 = load(sprintf('%s/figure_t14_%i_150.mat',dir_cacheLp,i));
        CaR = load(sprintf('%s/fig_cluster2_reduce_%i_new.mat',dir_cacheLp,i));

        % Define clustering variables
        cl_adj = CaT14.Adj_plt2;

        % Fill missing data
        %cl_adj = fillmissing(cl_adj,'nearest');
        cl_adj(isnan(cl_adj)) = 0;

        % nosig
        cl_adj_nosig = all(isnan(cl_adj) | (cl_adj == 0));


        %======================================================================
        % k-means clustering
        fn_kmeans = sprintf('cache/figure_t17_%i_kmeans_compute-%i',i,compute);

        if (false) %(exist([fn_kmeans,'.mat'],'file'))
            load([fn_kmeans,'.mat']);
        else
            K = 2:75;
            Ssumd = zeros(length(K),1);
            Res = cell(length(K),4);
            A2km = cl_adj;
            %A2km(isnan(A2km)) = 0; % set diagonals to zeros;
            for ki = 1:length(K)
                k = K(ki);
                stream = RandStream('mlfg6331_64','seed','shuffle');  % Random number stream
                %stream = RandStream('mt19937ar','Seed','shuffle');
                options = statset('UseParallel',1,'UseSubstreams',1,'Streams',stream);
                %[idx,C,sumd,D] = kmeans(A2km,k,'Options',options,'MaxIter',10000,'Display','off','Replicates',6*200);
                [idx,C,sumd,D] = kmedoids(A2km,k,'distance',@nanDist4clustering,'Options',options,'OnlinePhase','on','Replicates',6);
                %return
                
                Res{ki,1} = idx;
                Res{ki,2} = C;
                Res{ki,3} = sumd;
                Res{ki,4} = D;

                % cross-val metric
                if (length(K) > 1)
                    cv = sum(sumd); 
                    Ssumd(ki) = cv; %sum(sumd);
                end
            end
            save(fn_kmeans,'K','Ssumd','Res','A2km');
            fprintf('[*] %s\n',fn_kmeans);
        end
        %======================================================================

        compute_K{compute} = K;
        compute_Ssumd{compute} = Ssumd;
        compute_Res{compute} = Res;
        compute_A2km{compute} = A2km;

    end
    t_single = toc;
    
    fprintf('[*] Computed stage %i of %i (ETA: %.2f hrs)\n',compute,n_compute,t_single*(n_compute-compute)/(60*60))
end

fn_kmeans_compute = sprintf('cache/figure_t17_%i_kmedoids_compute-all_k-75',i);
save(fn_kmeans_compute);