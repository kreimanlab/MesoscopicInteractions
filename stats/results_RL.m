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

    fprintf('=== Figure 2 analysis === frequency band: %i\n',iM);
    % Load 150 parcellation cache
    Cai = load('cache/fig_cluster3_cluster_i.mat');
    cluster_i = Cai.cluster_i;
    %Ca5 = load('cache/fig_cluster2_reduce_annot.mat');
    Ca5 = load(sprintf('/home/jerry/data/stats/cache/fig_cluster2_reduce_%i_new.mat',iM));
    fn_ca1 = sprintf('/home/jerry/data/stats/cache/fig_cluster2_A_%i.mat',iM);
    Ca1 = load(fn_ca1);
    
    % Figure 2 interaction
    sid_int_const = 5; %Subjectsp = {'m00005'};
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