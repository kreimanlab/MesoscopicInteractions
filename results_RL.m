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

    fprintf('=== Figure 2 analysis === frequency band: %i\n',iM);
    % Load 150 parcellation cache
    Cai = load('cache/fig_cluster3_cluster_i.mat');
    cluster_i = Cai.cluster_i;
    %Ca5 = load('cache/fig_cluster2_reduce_annot.mat');
    Ca5 = load(sprintf('/home/jerry/data/stats/cache/fig_cluster2_reduce_%i_new.mat',iM));
    fn_ca1 = sprintf('/home/jerry/data/stats/cache/fig_cluster2_A_%i.mat',iM);
    Ca1 = load(fn_ca1);
    
    % Figure 2 interaction
    sid_int_const = 5; %Subjectsp = {'sub3'};
    bchan1_const = 35;
    bchan2_const = 84;
    
    % Search through 150 areas for two bip electrodes
    b1_r150 = NaN;
    b2_r150 = NaN;
    for i = 1:length(Ca5.Es)
        E = Ca5.Es{i};
        
        % Note: format of E from fig_cluster2.m
        %   E = [E; [sidint, j, coord1]];
        
        [n_E,~] = size(E);
        
        for j = 1:n_E
            e_sid_int = E(j,1);
            e_bchan = E(j,2);
            if (e_sid_int == sid_int_const)
                if (e_bchan == bchan1_const)
                    b1_r150 = i;
                    fprintf('[!] Found region: %i of 150 for bchan1\n',i);
                    fprintf('\tpost-cluster region: %i\n',cluster_i(i));
                end
                if (e_bchan == bchan2_const)
                    b2_r150 = i;
                    fprintf('[!] Found region: %i of 150 for bchan2\n',i);
                    fprintf('\tpost-cluster region: %i\n',cluster_i(i));
                end
            end
        end
        %return
    end
    fprintf('[*] Mean coherence between regions: %.4f\n',Ca5.A_nocov(b1_r150,b2_r150));
    fprintf('[*] Number of pairs: %.4f\n',Ca5.A_npairs(b1_r150,b2_r150));
    fprintf('[*] Number of usubs: %.4f\n',Ca5.A_nusubs(b1_r150,b2_r150));
    fprintf('[*] Consistency across subjects (pairs): %.4f (%.2f)\n',Ca5.Acp(b1_r150,b2_r150),Ca5.Acp(b1_r150,b2_r150)*Ca5.A_npairs(b1_r150,b2_r150));
    fprintf('[*] Consistency across subjects (usubs): %.4f (%.2f)\n',Ca5.Acs(b1_r150,b2_r150),Ca5.Acs(b1_r150,b2_r150)*Ca5.A_nusubs(b1_r150,b2_r150));
    fprintf('[*] Average distance: %.4f mm\n',Ca5.Ad(b1_r150,b2_r150));
    
    % Check the number of pairs
    fprintf('[!] Checking the number of subjects..\n');
    u_sid_ints = [];
    for i = [b1_r150,b2_r150]
        E = Ca5.Es{i};
        [n_E,~] = size(E);
        for j = 1:n_E
            e_sid_int = E(j,1);
            e_bchan = E(j,2);
        end
        u_sid_ints = [u_sid_ints; E(:,1)];
    end
    n_usubs_chk = unique(u_sid_ints);
    disp(n_usubs_chk');
    %fprintf('\tsubject: %i\n',n_usubs_chk);
    fprintf('\ttotal: %i\n',length(n_usubs_chk));
    return
    
    
    fprintf('[*] === Left-Right lateralization analysis === Frequency band: %i\n',iM);
    for i = 1:length(Subjects)
        sid = Subjects{i};
        fprintf('\tSubject: %s\n',sid);
        
        % Load cache
        Ca = load(sprintf('%s/xsub_out_%s_%i.mat',dir_cache,sid,iM));
        
        % Load hemisphere
        fn_hemi = sprintf('%s/%s/elec_recon/%s.electrodeNames',subjects_dir,sid,sid);
        He = textscan(fopen(fn_hemi,'r'),'%s','HeaderLines',2);
        He = reshape(He{1},3,[]);
        hemis = He(3,:)';
        
        % load DK cache
        fn_ca4 = sprintf('./cache/figure_t14_%i',iM);
        Ca4 = load(fn_ca4);
        
        return
    end
end