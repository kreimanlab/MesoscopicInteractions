close all;
clear;

Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};


[~,host] = system('hostname');
if (contains(host,'kraken'))
    %dir_res = '/home/jerry/data/results/coh_w10';
    dir_cor = '/media/klab/internal/data/coreg';
else
    dir_cor = '/n/groups/kreiman/jerry/data/coreg';
end


C1 = load(sprintf('%s/%s/label/all_parcellation.mat',dir_cor,Subjects{1}));

% --- Params ---
atl = 2;
trig_reset = true; % whether to start permutation from beginning
trig_plt = false;
N_PERM = 1;
N_SIGMA = 0.05; % normal distribution sigma
E_MU = 0.04; % exponential mu, proportional to temperature

iM = 1;
metric = metrics{iM};
ca14_fn = (sprintf('./cache/figure_t14_%i_atl%i_%s',iM,atl,C1.AtlNames{atl}));
Ca14 = load(ca14_fn);
Adj_plt2 = Ca14.Adj_plt2;
cluster_i = Ca14.cluster_i;
CaA = load(sprintf('./cache/xsub_out_all_%i_atl%i.mat',iM,atl));


% proportion of work per subject
trig_subject_combs = false;
if (trig_subject_combs)
    Scombs = zeros(1,length(Subjects));
    for iSub = 1:length(Subjects)
        sid = Subjects{iSub};
        Ca = load(sprintf('./cache/xsub_out_%s_%i_atl%i.mat',sid,iM,atl));
        Scombs(iSub) = Ca.ecog.n_comb;
    end
    cores = ceil(Scombs / min(Scombs));
    disp(cores);
    return;
end


for iSub = 1:length(Subjects)
    sid = Subjects{iSub};
    Ca = load(sprintf('./cache/xsub_out_%s_%i_atl%i.mat',sid,iM,atl));
    CaBE = load(sprintf('brainexport/%s_6all.mat',sid));
    CaBE_faces = CaBE.fa_pia + 1;
    CaBE_verts = CaBE.v_pia;
    CaBE_vf = CaBE.vf_pia;
    CaBE_paths = CaBE.Paths;

    % plot surface
    if (trig_plt)
        p = patch('faces',(CaBE_faces),'vertices',CaBE_verts); %'EdgeColor',0.1*[1 1 1] 
        p.LineStyle = 'none';
        p.FaceColor = [1 1 1]*0.9;
        p.FaceAlpha = 0.3;
        daspect([1 1 1]);
        brainlight;
        view(90,0);
    end

    n_bpairs = length(CaBE_paths);

    
    fn_tmp = sprintf('./cache_wdist/paths_sub-%i',iSub);
    if (trig_reset)
        %mkdir("cache_wdist");
        Path_iters = cell(1,n_bpairs);
    else
        CaTmp = load(fn_tmp);
        Path_iters = CaTmp.Path_iters;
    end
    
    dist_history_bip = [];
    dist_history_ibp = [];
    parfor ibp = 1:n_bpairs
        path = CaBE_paths{ibp};

        % get path size
        [n_pts,~] = size(path);
        n_segs = n_pts - 1;

        % initialize experimental path
        if (trig_reset)
            path_iter = path;
        else
            path_iter = Path_iters{ibp};
        end

        dist_path = pathsum(path_iter);
        fprintf('[s:%i %i/%i O] %.12f mm\n',iSub,ibp,n_bpairs,pathsum(path));
        fprintf('[s:%i %i/%i I] %.12f mm\n',iSub,ibp,n_bpairs,pathsum(path_iter));

        if (trig_plt)
            n_perm = 0;
        else
            n_perm = N_PERM;
        end
        
        dist_history = zeros(n_perm,1);
        for ip = 1:n_perm
            % randomly shift path
            path_new = path_iter + N_SIGMA*normrnd(0,1,size(path));
            path_new(1,:) = path_iter(1,:);
            path_new(end,:) = path_iter(end,:);

            % Check path points to make sure they are inside
            path_isin = inpolyhedron(CaBE_faces,CaBE_verts,path_new);

            % pull points outside inside
            for ipt = 1:n_pts
                if (~ path_isin(ipt))
                    ID = nearestNeighbor(CaBE_vf,path_new(ipt,:));
                    c1 = CaBE_vf.Points(ID,:); % new point
                    c2 = path_new(ipt,:); % old point
                    if (trig_plt)
                        hold on;
                        plot3([c1(1) c2(1)],[c1(2) c2(2)],[c1(3) c2(3)],'-red','LineWidth',2);
                    end
                    % set new point
                    path_new(ipt,:) = c1;
                end
            end

            dist_path_new = pathsum(path_new);
            cond_keep = false;
            expmu = E_MU; % 0.015 exponential distribution mean
            if (dist_path_new < dist_path)
                cond_keep = true;
            else
                if (rand() <= (1 - expcdf(abs(dist_path_new-dist_path),expmu)))
                    cond_keep = true;
                end
            end

            if (cond_keep)
                %fprintf('[*] %.12f mm\n',dist_path_new);
                fprintf('[s:%i %i/%i *] %.12f mm\n',iSub,ibp,n_bpairs,dist_path_new);
                path_iter = path_new;
                dist_path = dist_path_new;
            end
            
            dist_history(ip) = dist_path;
        end


        % plot path
        if (trig_plt)
            cc = viridis(n_segs);
            for iseg = 1:n_segs
                coord1 = path(iseg,:);
                coord2 = path(iseg+1,:);
                hold on;
                pl = plot3([coord1(1) coord2(1)],[coord1(2) coord2(2)],[coord1(3) coord2(3)],'-');
                pl.Color = cc(iseg,:);
                pl.LineWidth = 1;
            end
        end

        % plot new path
        if (trig_plt)
            cc = magma(n_segs);
            for iseg = 1:n_segs
                coord1 = path_iter(iseg,:);
                coord2 = path_iter(iseg+1,:);
                hold on;
                pl = plot3([coord1(1) coord2(1)],[coord1(2) coord2(2)],[coord1(3) coord2(3)],'-');
                pl.Color = cc(iseg,:);
                pl.LineWidth = 1;
            end
        end

        % plot electrodes
        if (trig_plt)
            plot3(path(1,1),path(2,2),path(3,3),'black.','MarkerSize',30);
            plot3(path(end,1),path(end,2),path(end,3),'black.','MarkerSize',30);
        end
        
        % Save
        Path_iters{ibp} = path_iter;
        dist_history_bip = [dist_history_bip, dist_history];
        dist_history_ibp = [dist_history_ibp, ibp];
        %fprintf('progress: %i of %i\n',length(dist_history_ibp),n_bpairs);
    end
    
    if (~trig_plt)
        if (trig_reset)
            n_perms = [];
            n_sigmas = [];
            e_mus = [];
            dist_histories = {};
            n_perms = [n_perms; N_PERM];
            n_sigmas = [n_sigmas; N_SIGMA];
            e_mus = [e_mus; E_MU];
            dist_histories = [dist_histories, {dist_history_bip; dist_history_ibp}];
            %return
        else
            n_perms = [CaTmp.n_perms; N_PERM];
            n_sigmas = [CaTmp.n_sigmas; N_SIGMA];
            e_mus = [CaTmp.e_mus; E_MU];
            dist_histories = [CaTmp.dist_histories, {dist_history_bip; dist_history_ibp}];
        end
        
        save_func(fn_tmp,Path_iters,n_perms,n_sigmas,e_mus,dist_histories);
        %save(fn_tmp,'Path_iters','n_perms','n_sigmas','e_mus','dist_histories');
        %path_iter = Path_iters{ibp};
    end

end

function [] = save_func(fn_tmp,Path_iters,n_perms,n_sigmas,e_mus,dist_histories)
    save(fn_tmp,'Path_iters','n_perms','n_sigmas','e_mus','dist_histories');
end

function [d] = pathsum(path)
    d = sum(sqrt(sum(diff(path).^2,2)));
end