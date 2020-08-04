close all;
clear;


% BROADBAND
% >> fig_cluster2_sworld_circle_150
% [*] Computing small world index for humans..
% L: 1.81
% C: 0.54
% Ll: 2.35
% Cl: 0.73
% Lr: 1.73
% Cr: 0.27
% 	n = 150
% 	(0.544958 >> 0.274678 and 1.812897 =~ 1.725280)
% 	sigma = 1.888102
% 	omega = 0.206519
% 	SWI = 0.509122


% GAMMA
% [*] Computing small world index for humans..
% L: 1.86
% C: 0.49
% Ll: 2.63
% Cl: 0.73
% Lr: 1.76
% Cr: 0.24
% 	n = 150
% 	(0.492845 >> 0.235812 and 1.855614 =~ 1.763527)
% 	sigma = 1.986278
% 	omega = 0.273614
% 	SWI = 0.466807


%     
% The Watts-strogatz definition for small world network is:
%   L >= L_random but C >> C_random
%           ws    - Watts-Strogatz 1998 (C >> Cr,L >= Lr)
%           sigma - Humphries-Gurney,2008 (sigma > 1 for small world network)
%           omega - Telesford-Laurienti,2011 (0 to 1: 1 is most small world)
%           I     - Neal,2017 (0 to 1: 1 is most small world)



iM = 1;
for iM = 1;%1:5 %[1 5]

    % Number of random initializations for small world index (default: 12)
    n_MC = 12;

    % Atlas index
    atl = 2; % 2 - Desikan-Killiany

    % Load human adjacency matrix
    fn_ca3 = sprintf('./cache/fig_cluster2_reduce_%i_new.mat',iM);
    Ca3 = load(fn_ca3);
    
    % Load 150-adjacency matrix
    fn_caT14 = './cache/figure_t14_1_150.mat';
    CaT14 = load(fn_caT14);
    Ahs = CaT14.Adj_plt2;
    
    %Ahs = Ca3.A; %Ca3.A2;
    n_Ahs = length(Ahs);
    

%     % Load DK
%     fn_ca = sprintf('./cache/xsub_out_all_%i.mat',iM);
%     Ca = load(fn_ca);
%     Ahs = nan(Ca.n_rois,Ca.n_rois);
%     An = Ahs;
%     Ausub = Ahs;
%     Ad = Ahs;
%     for i = 1:Ca.n_rois
%         for j = 1:Ca.n_rois
%             if (~isempty(Ca.AdjAtl{i,j}))
%                 lt = Ca.AdjAtl{i,j};
%                 li = isnan(lt) | (lt == 0);
%                 Ahs(i,j) = mean(lt(~li));
%                 An(i,j) = length(Ca.AdjAtlN{i,j});
%                 Ausub(i,j) = length(unique(Ca.AdjAtl_sid{i,j}));
%                 Ad(i,j) = mean(Ca.adjct_dist{i,j});
%             end
%         end
%     end
% 
%     % Remove unknown node
%     rois = Ca.rois;
%     unk_i = strcmpi(rois,'UNKNOWN');
%     rois = rois(~unk_i);
%     Ahs = Ahs(~unk_i,~unk_i);
%     An = An(~unk_i,~unk_i);
%     Ausub = Ausub(~unk_i,~unk_i);
%     Ad = Ad(~unk_i,~unk_i);
% 
% 
%     % load cache
%     fn_ca4 = sprintf('./cache/figure_t14_%i',iM);
%     Ca4 = load(fn_ca4);
%     rois2 = Ca4.rois_plt(Ca4.cluster_i);
%     Ahs = Ca4.Adj_plt2(Ca4.cluster_i,Ca4.cluster_i);
% 
%     % 
%     % 
%     % % Label nodes
%     % rois2 = cell(size(rois));
%     % for i = 1:length(rois2)
%     %     rois2{i} = convertRoiDK(rois{i});
%     % end
% 

    Ca4 = load('cache/fig_cluster3_cluster_i.mat');
    cluster_i = Ca4.cluster_i;
    rois2 = cell(length(Ahs),1);
    for i = 1:length(rois2)
        %rois2{i} = sprintf('%i',cluster_i(i));
        rois2{i} = sprintf('%i',i);
    end
    
% 
%     % Build node colors
%     NodeColors = Ca.C.AtlROIs{atl}.LH.table(:,1:3) / (2^8 - 1);
%     NodeRois = Ca.C.AtlROIs{atl}.LH.struct_names;
%     NodeRoisDK = cell(size(NodeRois));
%     for i = 1:length(NodeRois)
%         NodeRoisDK{i} = convertRoiDK(NodeRois{i});
%     end
%     rois2_col = zeros(length(rois2),3);
%     for i = 1:length(rois2)
%         % Find roi
%         sIdx = strcmpi(NodeRoisDK,rois2{i});
%         rois2_col(i,:) = NodeColors(sIdx,:);
%     end

    % make graph
    Agr = Ahs;
    
%     % ----------------------------------------------------------------
%     % Remove NaNs
%     A = Agr;
%     cond_pass = true;
%     comp_idx = 1:length(A);
%     while (cond_pass)
%         [~,A_b] = max(sum(isnan(A)));
%         A(A_b,:) = [];
%         A(:,A_b) = [];
%         comp_idx(A_b) = [];
%         cond_pass = (sum(isnan(A(:))) ~= 0);
%     end
%     Agr = A;
%     rois2 = rois2(comp_idx);
%     %Ahs_dist = Ahs_dist(comp_idx,comp_idx);
%     % ----------------------------------------------------------------
    
    NodeTable = table(rois2,'VariableNames',{'Name'});
    Agr(isnan(Agr)) = 0;
    G = graph(Agr,NodeTable,'upper','omitselfloops');

%     % Get rois2 original names
%     rois2_raw = cell(size(rois2));
%     for i = 1:length(rois2)
%         sIdx = strcmpi(NodeRoisDK,rois2{i});
%         rois2_raw{i} = NodeRois{sIdx};
%     end
% 
%     % 3-letter roi name
%     def_rois_short;
%     rois2_ab = cell(size(rois2));
%     for i = 1:length(rois2)
%         rois2_ab{i} = rois_short{strcmp(NodeRois,rois2_raw{i})};
%     end
    
    rois2_ab = rois2;

    % % network plot
    % h = figure;
    % Lcolor = [0 0 0];
    % hg = plot(G,'Layout','force','EdgeColor',Lcolor,'NodeColor',rois2_col,'NodeLabel',[]); %'Layout','force','EdgeLabel',G.Edges.Weight
    % %layout(hg,'force','UseGravity',true);
    % %layout(hg,'force','WeightEffect','direct');
    % axis off;


    % --- Random --------------------------------------------------------------
    n_nodes = length(Agr);
    n_edges = sum(sum(triu(Agr)>0));
    n_eposs = nchoosek(n_nodes,2);
    Agr2 = triu(Agr);
    Agr2 = Agr2(Agr2 > 0);
    Esamp = zeros(n_eposs,1);
    Esamp(1:n_edges) = Agr2;
    
    
    % --- random network sampling -----------------------------------------

    % Load previously initialized network
    Eca = load('fig_cluster2_sworld_circle_150_E');
    E = Eca.E;

    % Sample a new network
%     E = Esamp(randperm(n_eposs));
%     save('fig_cluster2_sworld_circle_150_E','E');

    % Cluster the random network by coherence
    fprintf('[*] Clustering random network..\n');
    Y = E';
    roi_dist_random = squareform(Y);
    Z = linkage(Y,'average');
    cluster_i_rand = optimalleaforder(Z,Y); % ,'transformation','inverse'
    roi_dist_random = roi_dist_random(cluster_i_rand,cluster_i_rand);
    %roi_dist(isinf(roi_dist)) = ;
    clash = nansum(nansum(triu(roi_dist_random,1) - triu(roi_dist_random,2)));
    fprintf('cluster clash: %.12f mm\n',clash)

    % ---------------------------------------------------------------------

    A_rand = zeros(n_nodes,n_nodes);
    kk = 1;
    for ii = 1:(n_nodes-1)
        for jj = (ii+1):n_nodes
            e = E(kk);
            A_rand(ii,jj) = e;
            A_rand(jj,ii) = e;
            kk = kk + 1;
        end
    end
    
    % Apply random clustering
    A_rand = A_rand(cluster_i_rand,cluster_i_rand);
    

    % --- Lattice -------------------------------------------------------------
    n_nodes = length(Agr);
    n_edges = sum(sum(triu(Agr)>0));
    A_lattice = zeros(n_nodes,n_nodes);
    k = 1; % edge count
    n = 1; % node position count
    nn = 1; % nearest neighbor count
    while (true)
        % break when number of edges run out
        if (k > n_edges)
            break;
        end

        % set lattice edges
        n1 = n;
        n2 = n + nn;
        % wrap around
        if (n2 > n_nodes)
            n2 = n2 - n_nodes;
        end
        A_lattice(n1,n2) = Agr2(k); %Agr(n1,n2); %1;
        A_lattice(n2,n1) = Agr2(k); %Agr(n2,n1); %1;

        % increment nearest neighbor count on full circle
        if (n == n_nodes)
            n = 0;
            nn = nn + 1;
        end

        k = k + 1;
        n = n + 1;
    end


    % circular plot
    col_circ = repmat(0.1*[1 1 1],length(rois2),1);
    h = figure('Position',round(1*[0, 0, 8*(6/3), 4]*100));
    set(h,'PaperUnits','inches');
    set(h,'PaperPosition',2*[0, 0, 8*(6/3), 4]);
    addpath('circularGraph');
    [ha, pos] = tight_subplot(1,4,[.01 .01],[.01 .01],[.01 .01]);
    
    % Dropout ROI labels
    for i = 1:length(rois2_ab)
        % TODO
    end
    
    axes(ha(1))
    circularGraphJW150(Agr,'Label',rois2_ab,'Colormap',col_circ);
    %circularGraphJW150(Agr(cluster_i,cluster_i),'Label',rois2_ab,'Colormap',col_circ(cluster_i,:)); %rois2_col
    axes(ha(2))
    circularGraphJW150(A_rand,'Label',rois2_ab,'Colormap',col_circ); %rois2_col
    axes(ha(3))
    circularGraphJW150(A_lattice,'Label',rois2_ab,'Colormap',col_circ); %rois2_col


    % Colorbar
    axes(ha(4))
    cc = colormap(inferno(2^8));
    col_var_min = min(Agr2);
    col_var_max = max(Agr2);
    cb = colorbar;
    colormap(cc);
    caxis([col_var_min col_var_max])
    set(cb,'YTick',[col_var_min,col_var_max]);
    set(cb,'YTickLabel',{sprintf('%.2f',col_var_min),sprintf('%.2f',col_var_max)});
    set(cb,'TickLength',0);
    set(cb,'Position',[0.8 0.1 0.01 0.4]);
    set(cb,'FontSize',10);

    axis off;
    % circularGraphJW(Agr,'Label',rois2_ab,'Colormap',col_circ); %rois2_col

    % Expand space
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];

    % save
    print(h,sprintf('./figures/fig_cluster2_sworld_circle_150_metric-%i',iM),'-depsc','-r400');
    print(h,sprintf('./figures/fig_cluster2_sworld_circle_150_metric-%i',iM),'-dpng','-r400');
    close(h);

    % Small world index calculations
    trig_swi = true;

    if (trig_swi)

        %Ahs = Ca3.A2;
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
        %iM = 1;
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

    end

end

% == Freq band 1 - Broadband ==
% fig_cluster2_sworld_circle_150
% [*] Clustering random network..
% cluster clash: 0.000000000000 mm
% [*] Computing small world index for humans..
% L: 2.42
% C: 0.44
% Ll: 5.24
% Cl: 0.70
% Lr: 2.06
% Cr: 0.10
% 	n = 150
% 	(0.435876 >> 0.103649 and 2.420218 =~ 2.062722)
% 	sigma = 3.584140
% 	omega = 0.229197
% 	SWI = 0.494839
% [*] Computing small world index for C elegans..
% L: 2.44
% C: 0.34
% Ll: 9.11
% Cl: 0.70
% Lr: 2.31
% Cr: 0.06
% 	n = 281
% 	(0.335113 >> 0.057862 and 2.436001 =~ 2.307355)
% 	sigma = 5.485725
% 	omega = 0.469172
% 	SWI = 0.422907
% [*] Computing small world index for Markov-Kennedy..
% L: 1.66
% C: 0.74
% Ll: 1.99
% Cl: 0.72
% Lr: 1.66
% Cr: 0.34
% 	n = 91
% 	(0.741518 >> 0.338230 and 1.662271 =~ 1.662047)
% 	sigma = 2.192056
% 	omega = -0.023230
% 	SWI = 1.042583
% [*] Computing small world index for macaque functional interactions..
% L: 1.74
% C: 0.59
% Ll: 2.17
% Cl: 0.72
% Lr: 1.69
% Cr: 0.30
% 	n = 86
% 	(0.586906 >> 0.304734 and 1.735494 =~ 1.694802)
% 	sigma = 1.880801
% 	omega = 0.161480
% 	SWI = 0.620964
% [*] Clustering random network..
% cluster clash: 0.000000000000 mm


% == Freq band 2 ==
% [*] Computing small world index for humans..
% L: 2.42
% C: 0.44
% Ll: 5.24
% Cl: 0.70
% Lr: 2.06
% Cr: 0.10
% 	n = 150
% 	(0.435876 >> 0.103568 and 2.420218 =~ 2.062215)
% 	sigma = 3.586041
% 	omega = 0.228987
% 	SWI = 0.494813
% [*] Computing small world index for C elegans..
% L: 2.44
% C: 0.34
% Ll: 9.11
% Cl: 0.70
% Lr: 2.31
% Cr: 0.06
% 	n = 281
% 	(0.335113 >> 0.058269 and 2.436001 =~ 2.307422)
% 	sigma = 5.447593
% 	omega = 0.469199
% 	SWI = 0.422558
% [*] Computing small world index for Markov-Kennedy..
% L: 1.66
% C: 0.74
% Ll: 1.99
% Cl: 0.72
% Lr: 1.66
% Cr: 0.34
% 	n = 91
% 	(0.741518 >> 0.336685 and 1.662271 =~ 1.662027)
% 	sigma = 2.202089
% 	omega = -0.023242
% 	SWI = 1.042345
% [*] Computing small world index for macaque functional interactions..
% L: 1.74
% C: 0.59
% Ll: 2.17
% Cl: 0.72
% Lr: 1.69
% Cr: 0.31
% 	n = 86
% 	(0.586906 >> 0.306863 and 1.735494 =~ 1.694847)
% 	sigma = 1.867801
% 	omega = 0.161507
% 	SWI = 0.619514
% [*] Clustering random network..
% cluster clash: 0.000000000000 mm


% == Freq band 3 ==
% [*] Computing small world index for humans..
% L: 2.42
% C: 0.44
% Ll: 5.24
% Cl: 0.70
% Lr: 2.06
% Cr: 0.11
% 	n = 150
% 	(0.435876 >> 0.106616 and 2.420218 =~ 2.061350)
% 	sigma = 3.482076
% 	omega = 0.228630
% 	SWI = 0.492661
% [*] Computing small world index for C elegans..
% L: 2.44
% C: 0.34
% Ll: 9.11
% Cl: 0.70
% Lr: 2.31
% Cr: 0.06
% 	n = 281
% 	(0.335113 >> 0.057862 and 2.436001 =~ 2.307846)
% 	sigma = 5.486882
% 	omega = 0.469373
% 	SWI = 0.422937
% [*] Computing small world index for Markov-Kennedy..
% L: 1.66
% C: 0.74
% Ll: 1.99
% Cl: 0.72
% Lr: 1.66
% Cr: 0.34
% 	n = 91
% 	(0.741518 >> 0.337732 and 1.662271 =~ 1.662047)
% 	sigma = 2.195287
% 	omega = -0.023230
% 	SWI = 1.042527
% [*] Computing small world index for macaque functional interactions..
% L: 1.74
% C: 0.59
% Ll: 2.17
% Cl: 0.72
% Lr: 1.69
% Cr: 0.31
% 	n = 86
% 	(0.586906 >> 0.305532 and 1.735494 =~ 1.694756)
% 	sigma = 1.875839
% 	omega = 0.161454
% 	SWI = 0.620340
% [*] Clustering random network..
% cluster clash: 0.000000000000 mm

% == Freq band 4 ==
% [*] Computing small world index for humans..
% L: 2.42
% C: 0.44
% Ll: 5.24
% Cl: 0.70
% Lr: 2.06
% Cr: 0.10
% 	n = 150
% 	(0.435876 >> 0.103971 and 2.420218 =~ 2.059083)
% 	sigma = 3.566741
% 	omega = 0.227693
% 	SWI = 0.494061
% [*] Computing small world index for C elegans..
% L: 2.44
% C: 0.34
% Ll: 9.11
% Cl: 0.70
% Lr: 2.31
% Cr: 0.06
% 	n = 281
% 	(0.335113 >> 0.058702 and 2.436001 =~ 2.307340)
% 	sigma = 5.407238
% 	omega = 0.469165
% 	SWI = 0.422176
% [*] Computing small world index for Markov-Kennedy..
% L: 1.66
% C: 0.74
% Ll: 1.99
% Cl: 0.72
% Lr: 1.66
% Cr: 0.34
% 	n = 91
% 	(0.741518 >> 0.338941 and 1.662271 =~ 1.662047)
% 	sigma = 2.187454
% 	omega = -0.023230
% 	SWI = 1.042662
% [*] Computing small world index for macaque functional interactions..
% L: 1.74
% C: 0.59
% Ll: 2.17
% Cl: 0.72
% Lr: 1.69
% Cr: 0.30
% 	n = 86
% 	(0.586906 >> 0.304985 and 1.735494 =~ 1.694779)
% 	sigma = 1.879232
% 	omega = 0.161467
% 	SWI = 0.620758
% [*] Clustering random network..
% cluster clash: 0.000000000000 mm


% == Freq band 5 - Gamma ==
% [*] Computing small world index for humans..
% L: 2.42
% C: 0.44
% Ll: 5.24
% Cl: 0.70
% Lr: 2.07
% Cr: 0.11
% 	n = 150
% 	(0.435876 >> 0.106632 and 2.420218 =~ 2.065011)
% 	sigma = 3.487725
% 	omega = 0.230143
% 	SWI = 0.493218
% [*] Computing small world index for C elegans..
% L: 2.44
% C: 0.34
% Ll: 9.11
% Cl: 0.70
% Lr: 2.31
% Cr: 0.06
% 	n = 281
% 	(0.335113 >> 0.058054 and 2.436001 =~ 2.307048)
% 	sigma = 5.466898
% 	omega = 0.469045
% 	SWI = 0.422721
% [*] Computing small world index for Markov-Kennedy..
% L: 1.66
% C: 0.74
% Ll: 1.99
% Cl: 0.72
% Lr: 1.66
% Cr: 0.34
% 	n = 91
% 	(0.741518 >> 0.337507 and 1.662271 =~ 1.662027)
% 	sigma = 2.196721
% 	omega = -0.023242
% 	SWI = 1.042436
% [*] Computing small world index for macaque functional interactions..
% L: 1.74
% C: 0.59
% Ll: 2.17
% Cl: 0.72
% Lr: 1.69
% Cr: 0.30
% 	n = 86
% 	(0.586906 >> 0.304910 and 1.735494 =~ 1.694847)
% 	sigma = 1.879766
% 	omega = 0.161507
% 	SWI = 0.620900