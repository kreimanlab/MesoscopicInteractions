close all;
clear;

addpath(genpath('SWP'));



% The Watts-strogatz definition for small world network is:
%   L >= L_random but C >> C_random
%           ws    - Watts-Strogatz 1998 (C >> Cr,L >= Lr)
%           sigma - Humphries-Gurney,2008 (sigma > 1 for small world network)
%           omega - Telesford-Laurienti,2011 (0 to 1: 1 is most small world)
%           I     - Neal,2017 (0 to 1: 1 is most small world)



%iM = 1;
for iM = 1 %[1:5] %[1 5]
    fprintf('=== Frequency band: %i ===\n',iM);

    % Number of random initializations for small world index (default: 12)
    n_MC = 12;

    % Atlas index
    atl = 2; % 2 - Desikan-Killiany

    % Load human adjacency matrix
    %fn_ca3 = './cache/fig_cluster2_fill_1.mat';
    %Ca3 = load(fn_ca3);

    fn_ca = sprintf('./cache/xsub_out_all_%i_atl2_pCorr.mat',iM);
    Ca = load(fn_ca);
    Ahs = nan(Ca.n_rois,Ca.n_rois);
    An = Ahs;
    Ausub = Ahs;
    Ad = Ahs;
    for i = 1:Ca.n_rois
        for j = 1:Ca.n_rois
            if (~isempty(Ca.AdjAtl{i,j}))
                lt = Ca.AdjAtl{i,j};
                li = isnan(lt) | (lt == 0);
                Ahs(i,j) = mean(lt(~li));
                An(i,j) = length(Ca.AdjAtlN{i,j});
                Ausub(i,j) = length(unique(Ca.AdjAtl_sid{i,j}));
                Ad(i,j) = mean(Ca.adjct_dist{i,j});
            end
        end
    end

    % Remove unknown node
    rois = Ca.rois;
    unk_i = strcmpi(rois,'UNKNOWN');
    rois = rois(~unk_i);
    Ahs = Ahs(~unk_i,~unk_i);
    An = An(~unk_i,~unk_i);
    Ausub = Ausub(~unk_i,~unk_i);
    Ad = Ad(~unk_i,~unk_i);


    % load cache
    % figure_t14_%i_atl%i_%s
    fn_ca4 = sprintf('./cache/figure_t14_%i_atl2_Desikan-Killiany_pCorr',iM);
    Ca4 = load(fn_ca4);
    rois2 = Ca4.rois_plt(Ca4.cluster_i);
    Ahs = Ca4.Adj_plt2(Ca4.cluster_i,Ca4.cluster_i);
    %Ahs_dist = Ca4.Adj_dist(Ca4.cluster_i,Ca4.cluster_i);
    Ahs_dist = ones(size(Ca4.Adj_plt2));
    
    
    
    % ----------------------------------------------------------------
    % Remove NaNs
    %return
    Ahs_nodiag = Ahs;
    for inod = 1:length(Ahs)
        Ahs_nodiag(inod,inod) = 0;
    end
    A = Ahs_nodiag;
    
    A = Ahs;
    cond_pass = true;
    comp_idx = 1:length(A);
    while (cond_pass)
        [~,A_b] = max(sum(isnan(A)));
        A(A_b,:) = [];
        A(:,A_b) = [];
        comp_idx(A_b) = [];
        cond_pass = (sum(isnan(A(:))) ~= 0);
    end
    Ahs = A;
    rois2 = rois2(comp_idx);
    Ahs_dist = Ahs_dist(comp_idx,comp_idx);
    % ----------------------------------------------------------------
    
    
    
%     % Filter out nocoverage areas
%     cov_idx = ~ all(isnan(Ahs));
%     Ahs = Ahs(cov_idx,cov_idx);
%     rois2 = rois2(cov_idx);

    % 
    % 
    % % Label nodes
    % rois2 = cell(size(rois));
    % for i = 1:length(rois2)
    %     rois2{i} = convertRoiDK(rois{i});
    % end

    NodeTable = table(rois2,'VariableNames',{'Name'});

    % Build node colors
    NodeColors = Ca.C.AtlROIs{atl}.LH.table(:,1:3) / (2^8 - 1);
    NodeRois = Ca.C.AtlROIs{atl}.LH.struct_names;
    NodeRoisDK = cell(size(NodeRois));
    for i = 1:length(NodeRois)
        NodeRoisDK{i} = convertRoiDK(NodeRois{i});
    end
    rois2_col = zeros(length(rois2),3);
    for i = 1:length(rois2)
        % Find roi
        sIdx = strcmpi(NodeRoisDK,rois2{i});
        rois2_col(i,:) = NodeColors(sIdx,:);
    end

    
    
    %======================================================================
    % make graph
    Agr = Ahs;
    Agr(isnan(Agr)) = 0;
    G = graph(Agr,NodeTable,'upper','omitselfloops');

    % Get rois2 original names
    rois2_raw = cell(size(rois2));
    for i = 1:length(rois2)
        sIdx = strcmpi(NodeRoisDK,rois2{i});
        rois2_raw{i} = NodeRois{sIdx};
    end

    % 3-letter roi name
    def_rois_short;
    rois2_ab = cell(size(rois2));
    for i = 1:length(rois2)
        rois2_ab{i} = rois_short{strcmp(NodeRois,rois2_raw{i})};
    end

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

    %return
    trig_perform_random_init = true;
    
    if (trig_perform_random_init)
        % Perform random initialization
        E = Esamp(randperm(n_eposs));
        save('fig_cluster2_sworld_circle_E_pCorr','E');
    else
        % Load existing randomly initialized network
        Eca = load('fig_cluster2_sworld_circle_E_pCorr');
        E = Eca.E;
    end

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
    
    % Random matrix override
    A_rand = randomize_matrix(Agr);
    

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
    
    % Lattice network override
    %calculate the number of nodes
    n_n = length(Agr);  
    %compute the average degree of the unweighted network, to give
    %the approximate radius
    numb_connections = length(find(Agr>0));
    avg_deg_unw = numb_connections/n_n;
    avg_rad_unw = avg_deg_unw/2;
    avg_rad_eff = ceil(avg_rad_unw);
    A_lattice = regular_matrix_generator(Agr, avg_rad_eff);


    % circular plot
    col_circ = repmat(0.1*[1 1 1],length(rois2),1);
    h = figure('Position',round([0, 0, 8*(6/3), 4]*100));
    set(h,'PaperUnits','inches');
    set(h,'PaperPosition',[0, 0, 8*(6/3), 4]);
    addpath('circularGraph');
    [ha, pos] = tight_subplot(1,4,[.01 .01],[.01 .01],[.01 .01]);
    axes(ha(1))
    circularGraphJW(Agr,'Label',rois2_ab,'Colormap',col_circ); %rois2_col
    axes(ha(2))
    circularGraphJW(A_rand,'Label',rois2_ab,'Colormap',col_circ); %rois2_col
    axes(ha(3))
    circularGraphJW(A_lattice,'Label',rois2_ab,'Colormap',col_circ); %rois2_col


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
    print(h,sprintf('./figures/fig_cluster2_sworld_circle_metric-%i_pCorr',iM),'-depsc');
    print(h,sprintf('./figures/fig_cluster2_sworld_circle_metric-%i_pCorr',iM),'-dpng','-r900');
    close(h);
    
    
    %======================================================================
    
    
    

    % Small world index calculations
    trig_swi = false;
    D_anat = load('../coreg/AdjMKFV.mat');
    
    if (trig_swi)

        %Ahs = Ca3.A2;
        n_Ahs = length(Ahs);

        

        fprintf('[*] Computing small world index for humans..\n')
        [ws,sigmaw,omega,I] = swi(Ahs,n_MC);
        fprintf('\tn = %i\n',n_Ahs);
        fprintf('\t%s\n',ws)
        fprintf('\tsigma = %.6f\n',sigmaw);
        fprintf('\tomega = %.6f\n',omega);
        fprintf('\tSWI = %.6f\n',I);
        
        % Alternative nan interpolation
        trig_interp = true;
        if (trig_interp)
            Ahs_previous = tril(fillmissing(Ahs,'previous'));
            Ahs_previous = (Ahs_previous+Ahs_previous') - eye(size(Ahs_previous,1)).*diag(Ahs_previous);
            Ahs_nearest = tril(fillmissing(Ahs,'nearest'));
            Ahs_nearest = (Ahs_nearest+Ahs_nearest') - eye(size(Ahs_nearest,1)).*diag(Ahs_nearest);
            Ahs_linear = tril(fillmissing(Ahs,'linear'));
            Ahs_linear = (Ahs_linear+Ahs_linear') - eye(size(Ahs_linear,1)).*diag(Ahs_linear);
            
            % Distance-based imputing
            Ahs_custom = Ahs;
            n_d = length(Ahs);
            for ii = 1:(n_d-1)
                for jj = (ii+1):n_d
                    a = Ahs_custom(ii,jj);
                    if (isnan(a))
                        %return
                        a_dist = Ahs_dist(ii,:);
                        a_dist(ii) = Inf;
                        a_dist(jj) = Inf;
                        [mDist,mIdx] = min(a_dist);
                        % fprintf('[*] impute distance: %.6f mm\n',mDist);
                        % largest distance: 34.955626 mm
                        Ahs_custom(ii,jj) = Ahs_custom(ii,mIdx);
                    end
                end
            end
            Ahs_custom = tril(fillmissing(Ahs,'nearest'));
            Ahs_custom = (Ahs_custom+Ahs_custom') - eye(size(Ahs_custom,1)).*diag(Ahs_custom);
            
            %return
            
            % == custom
            fprintf('[*] Fill missing values with method: distance-based\n')
            [ws,sigmaw,omega,I] = swi(Ahs_custom,n_MC);
            fprintf('\tn = %i\n',n_Ahs);
            fprintf('\t%s\n',ws)
            fprintf('\tsigma = %.6f\n',sigmaw);
            fprintf('\tomega = %.6f\n',omega);
            fprintf('\tSWI = %.6f\n',I);
            
            % Random
            Ahs_rand = Ahs;
            n_perm = 100;
            S = cell(n_perm,4);
            Sraw = zeros(n_perm,6);
            for iperm = 1:n_perm
                Ap = Ahs_rand;
                for ip1 = 1:(length(Ap))
                    for ip2 = (ip1):length(Ap)
                        if (isnan(Ap(ip1,ip2)))
                            pap = round(rand());
                            Ap(ip1,ip2) = pap;
                            Ap(ip2,ip1) = pap;
                        end
                    end
                end
                [ws,sigmaw,omega,I] = swi(Ap,n_MC);
                [L,C,Ll,Cl,Lr,Cr] = swi_raw(Ap,n_MC);
                S{iperm,1} = ws;
                S{iperm,2} = sigmaw;
                S{iperm,3} = omega;
                S{iperm,4} = I;
                Sraw(iperm,:) = [L,C,Ll,Cl,Lr,Cr];
            end
            fprintf('[*] Fill missing values with method: random\n')
            fprintf('\tn = %i\n',n_Ahs);
            fprintf('\t%s\n',ws)
            sig = [S{:,2}];
            fprintf('\tsigma = %.6f +- %.6f (%.6f - %.6f), n=%i\n',mean(sig),std(sig),min(sig),max(sig),n_perm);
            omg = [S{:,3}];
            fprintf('\tomega = %.6f +- %.6f (%.6f - %.6f), n=%i\n',mean(omg),std(omg),min(omg),max(omg),n_perm);
            %fprintf('\tomega = %.6f\n',omega);
            isw = [S{:,4}];
            fprintf('\tSWI = %.6f +- %.6f (%.6f - %.6f), n=%i\n',mean(isw),std(isw),min(isw),max(isw),n_perm);
            %fprintf('\tSWI = %.6f\n',I);
            
            
            % == previous
%             fprintf('[*] Fill missing values with method: previous\n')
%             [ws,sigma,omega,I] = swi(Ahs_previous,n_MC);
%             fprintf('\tn = %i\n',n_Ahs);
%             fprintf('\t%s\n',ws)
%             fprintf('\tsigma = %.6f\n',sigma);
%             fprintf('\tomega = %.6f\n',omega);
%             fprintf('\tSWI = %.6f\n',I);
            
            % == nearest
%             fprintf('[*] Fill missing values with method: nearest\n')
%             [ws,sigma,omega,I] = swi(Ahs_nearest,n_MC);
%             fprintf('\tn = %i\n',n_Ahs);
%             fprintf('\t%s\n',ws)
%             fprintf('\tsigma = %.6f\n',sigma);
%             fprintf('\tomega = %.6f\n',omega);
%             fprintf('\tSWI = %.6f\n',I);
            
            % == linear
%             fprintf('[*] Fill missing values with method: linear\n')
%             [ws,sigma,omega,I] = swi(Ahs_linear,n_MC);
%             fprintf('\tn = %i\n',n_Ahs);
%             fprintf('\t%s\n',ws)
%             fprintf('\tsigma = %.6f\n',sigma);
%             fprintf('\tomega = %.6f\n',omega);
%             fprintf('\tSWI = %.6f\n',I);
        end

        % Load individual adjacency matrices
        Ca = load('./cache/xsub_out_all_1.mat');
        %iM = 1;
        for i = 1:length(Ca.Subjects)
            sid = Ca.Subjects{i};
            Cas = load(sprintf('./cache/xsub_out_%s_%i_pCorr.mat',sid,iM));
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
        
        Aswp = Ace;
        nsamps = 12;
        SWP = zeros(1,nsamps);
        DC = zeros(1,nsamps);
        DL = zeros(1,nsamps);
        for isamp = 1:nsamps
            [swp,dc,dl] = small_world_propensity(Aswp);
            SWP(isamp) = swp;
            DC(isamp) = dc;
            DL(isamp) = dl;
        end
        Devi = (4*(atan2(DL,DC)))/pi - 1;
        fprintf('[*] C elegans SWP (n=%i), weighted\n',length(Aswp));
        fprintf('\tn_samp: %i\n',nsamps);
        fprintf('\tSWP: %.6f (+- %.6f)\n',mean(SWP),std(SWP));
        fprintf('\tdeltaC: %.6f (+- %.6f)\n',mean(DC),std(DC));
        fprintf('\tdeltaL: %.6f (+- %.6f)\n',mean(DL),std(DL));
        fprintf('\tdeviance: %.6f (+- %.6f)\n',mean(Devi),std(Devi));

%         fprintf('[*] Computing small world index for C elegans..\n')
%         [ws,sigmaw,omega,I] = swi(Ace,n_MC);
%         fprintf('\tn = %i\n',n_Ace);
%         fprintf('\t%s\n',ws)
%         fprintf('\tsigma = %.6f\n',sigmaw);
%         fprintf('\tomega = %.6f\n',omega);
%         fprintf('\tSWI = %.6f\n',I);


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

%         fprintf('[*] Computing small world index for Markov-Kennedy..\n')
%         [ws,sigmaw,omega,I] = swi(Amaca,n_MC);
%         fprintf('\tn = %i\n',length(Amaca));
%         fprintf('\t%s\n',ws);
%         fprintf('\tsigma = %.6f\n',sigmaw);
%         fprintf('\tomega = %.6f\n',omega);
%         fprintf('\tSWI = %.6f\n',I);
% 
%         fprintf('[*] Computing small world index for macaque functional interactions..\n')
%         [ws,sigmaw,omega,I] = swi(Amacf,n_MC);
%         fprintf('\tn = %i\n',length(Amacf));
%         fprintf('\t%s\n',ws);
%         fprintf('\tsigma = %.6f\n',sigmaw);
%         fprintf('\tomega = %.6f\n',omega);
%         fprintf('\tSWI = %.6f\n',I);

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

    save(sprintf('cache/fig_cluster2_sworld_circle-%i_pCorr.mat',iM));
    
    %% --------------------------------------------------------------
    % Weighted SWI
    fprintf('[!] === Weighted small world index calculation ===\n');
    addpath(genpath('BCT'));
    addpath(genpath('SWP'));
    
    % Initiate graphs
    
    % === Monkey ===
    A = D_anat.AdjMKflne; %log10(AdjMKneurons);
    
    % Introduce nans
    idx_nan = all(A==0);
    rois2_mk = D_anat.MK_labels(~idx_nan);
    %A(idx_nan,idx_nan) = NaN;
    
    % Symmetricize
    A = (A + A') / 2;
    %A(idx_nan,:) = [];
    %A(:,idx_nan) = [];
    
    %A = log10(A);
    %A(isinf(A)) = 0;
    %A = threshold_proportional(A,2/3);
    %A = weight_conversion(A,'normalize');
    A = weight_conversion(A,'autofix');
    
%     % Monkey brain parameters
%     Cw = mean(clustering_coef_wu(A));
%     C = mean(clustering_coef_bu(A>0));
%     %[D,~] = distance_wei(A);
%     [Dw,~] = distance_wei(weight_conversion(A,'lengths'));
%     Lw = charpath(Dw);
%     L = charpath(distance_bin(A>0));
%     
%     n_samp = 4*60;
%     n_resamp = 60;
%     Cw_r = zeros(n_samp,1);
%     Lw_r = zeros(n_samp,1);
%     Cw_l = zeros(n_samp,1);
%     Lw_l = zeros(n_samp,1);
%     C_r = zeros(n_samp,1);
%     L_r = zeros(n_samp,1);
%     C_l = zeros(n_samp,1);
%     L_l = zeros(n_samp,1);
%     parfor isamp = 1:n_samp
%         fprintf('[%5.i/%5.i]\n',isamp,n_samp);
%         Ar = randmio_und(A,n_resamp);
%         %Ar = randomgraph(A);
%         Al = latmio_und(A,n_resamp);
%         Cw_r(isamp) = mean(clustering_coef_wu(Ar));
%         Cw_l(isamp) = mean(clustering_coef_wu(Al));
%         C_r(isamp) = mean(clustering_coef_bu(Ar>0));
%         C_l(isamp) = mean(clustering_coef_bu(Al>0));
% %         [D_r,~] = distance_wei(Ar);
% %         [D_l,~] = distance_wei(Al);
%         [Dw_r,~] = distance_wei(weight_conversion(Ar,'lengths'));
%         [Dw_l,~] = distance_wei(weight_conversion(Al,'lengths'));
%         Lw_r(isamp) = charpath(Dw_r);
%         Lw_l(isamp) = charpath(Dw_l);
%         L_r(isamp) = charpath(distance_bin(Ar));
%         L_l(isamp) = charpath(distance_bin(Al));
%     end
%     
%     gamma = C ./ C_r;
%     lambda = L ./ L_r;
%     sigma = gamma./lambda;
%     delta_C = (C_l - C)./(C_l - C_r);
%     delta_L = (L - L_r)./(L_l - L_r);
%     phi = 1 - sqrt((delta_C.^2 + delta_L.^2)./2);
%     fprintf('[*] MK FLNe:\n')
%     fprintf('Gamma: %.6f +- %.6f\n',mean(gamma),std(gamma));
%     fprintf('Lambda: %.6f +- %.6f\n',mean(lambda),std(lambda));
%     fprintf('Sigma: %.6f +- %.6f\n',mean(sigma),std(sigma));
%     fprintf('phi: %.6f +- %.6f\n',mean(phi),std(phi));
%     
%     gammaw = Cw ./ Cw_r;
%     lambdaw = Lw ./ Lw_r;
%     sigmaw = gammaw./lambdaw;
%     deltaw_C = (Cw_l - Cw)./(Cw_l - Cw_r);
%     deltaw_L = (Lw - Lw_r)./(Lw_l - Lw_r);
%     phiw = 1 - sqrt((deltaw_C.^2 + deltaw_L.^2)./2);
%     fprintf('[*] MK FLNe:\n')
%     fprintf('Gamma, weighted: %.6f +- %.6f\n',mean(gammaw),std(gammaw));
%     fprintf('Lambda, weighted: %.6f +- %.6f\n',mean(lambdaw),std(lambdaw));
%     fprintf('Sigma, weighted: %.6f +- %.6f\n',mean(sigmaw),std(sigmaw));
%     fprintf('phi, weighted: %.6f +- %.6f\n',mean(phiw),std(phiw));
%     
%     [SWP,delta_C,delta_L] = small_world_propensity(A);





    % Weighted Monkey anatomical MK, overall (raw)
    % SWP raw
    % --------------------------------------------------------------------
    A_input = A;
    nsamps = 10000;
    fprintf('[*] Calculating monkey MK-FLNe SWP with %i permutations..\n',nsamps);
    SWPraw = zeros(1,nsamps);
    DCraw = zeros(1,nsamps);
    DLraw = zeros(1,nsamps);
    Gammaraw = zeros(1,nsamps);
    Lambdaraw = zeros(1,nsamps);
    Sigmaraw = zeros(1,nsamps);
    Ln = zeros(1,nsamps);
    Lr = zeros(1,nsamps);
    Ll = zeros(1,nsamps);
    Cn = zeros(1,nsamps);
    Cr = zeros(1,nsamps);
    Cl = zeros(1,nsamps);
    Omegaraw = zeros(1,nsamps);
    SWIraw = zeros(1,nsamps);
    Qraw = zeros(1,nsamps);
    parfor isamp = 1:nsamps
        %[swp,dc,dl] = small_world_propensity(Ahs);
        [swp,dc,dl,net_clus,rand_clus,net_path,rand_path,reg_clus,reg_path,gamma,lambda,sigma] = small_world_propensity_raw(A_input);
        SWPraw(isamp) = swp;
        DCraw(isamp) = dc;
        DLraw(isamp) = dl;
        Gammaraw(isamp) = gamma;
        Lambdaraw(isamp) = lambda;
        Sigmaraw(isamp) = sigma;
        Ln(isamp) = net_path;
        Lr(isamp) = rand_path;
        Ll(isamp) = reg_path;
        Cn(isamp) = net_clus;
        Cr(isamp) = rand_clus;
        Cl(isamp) = reg_clus;
        Omegaraw(isamp) = (rand_path/net_path) - (net_clus/reg_clus);
        SWIraw(isamp) = ((net_path - reg_path)/(rand_path - reg_path)) * ((net_clus - rand_clus)/(reg_clus - rand_clus));
    end
    Deviraw = (4*(atan2(DLraw,DCraw)))/pi - 1;
    fprintf('[*] Adj, nodes=%i\n',length(A_input));
    alpha = 0.05;
    sDCraw = sort(DCraw);
    sDLraw = sort(DLraw);
    sSWPraw = sort(SWPraw);
    sDeviraw = sort(Deviraw);
    sGammaraw = sort(Gammaraw);
    sLambdaraw = sort(Lambdaraw);
    sSigmaraw = sort(Sigmaraw);
    sOmegaraw = sort(Omegaraw);
    sLn = sort(Ln);
    sLr = sort(Lr);
    sCn = sort(Cn);
    sCr = sort(Cr);
    sLl = sort(Ll);
    sCl = sort(Cl);
    sSWIraw = sort(SWIraw);
    fprintf('\tL\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sLn),std(sLn),nsamps,sLn(round((alpha/2)*nsamps)),sLn(round((1 - (alpha/2))*nsamps)));
    fprintf('\tLr\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sLr),std(sLr),nsamps,sLr(round((alpha/2)*nsamps)),sLr(round((1 - (alpha/2))*nsamps)));
    fprintf('\tLl\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sLl),std(sLl),nsamps,sLl(round((alpha/2)*nsamps)),sLl(round((1 - (alpha/2))*nsamps)));
    fprintf('\tC\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sCn),std(sCn),nsamps,sCn(round((alpha/2)*nsamps)),sCn(round((1 - (alpha/2))*nsamps)));
    fprintf('\tCr\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sCr),std(sCr),nsamps,sCr(round((alpha/2)*nsamps)),sCr(round((1 - (alpha/2))*nsamps)));
    fprintf('\tCl\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sCl),std(sCl),nsamps,sCl(round((alpha/2)*nsamps)),sCl(round((1 - (alpha/2))*nsamps)));
    
    fprintf('\tDC\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(DCraw),std(DCraw),nsamps,sDCraw(round((alpha/2)*nsamps)),sDCraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tDL\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(DLraw),std(DLraw),nsamps,sDLraw(round((alpha/2)*nsamps)),sDLraw(round((1 - (alpha/2))*nsamps)));
    
    fprintf('\tSWP\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(SWPraw),std(SWPraw),nsamps,sSWPraw(round((alpha/2)*nsamps)),sSWPraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tSWI\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(SWIraw),std(SWIraw),nsamps,sSWIraw(round((alpha/2)*nsamps)),sSWIraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tdelta\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Deviraw),std(Deviraw),nsamps,sDeviraw(round((alpha/2)*nsamps)),sDeviraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tgamma\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Gammaraw),std(Gammaraw),nsamps,sGammaraw(round((alpha/2)*nsamps)),sGammaraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tlambda\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Lambdaraw),std(Lambdaraw),nsamps,sLambdaraw(round((alpha/2)*nsamps)),sLambdaraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tsigma\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Sigmaraw),std(Sigmaraw),nsamps,sSigmaraw(round((alpha/2)*nsamps)),sSigmaraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tomega\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Omegaraw),std(Omegaraw),nsamps,sOmegaraw(round((alpha/2)*nsamps)),sOmegaraw(round((1 - (alpha/2))*nsamps)));
    % --------------------------------------------------------------------
    


    % Weighted Monkey
    ci_alpha = 0.05;
    Aswp = A;
    nsamps = 100;
    
%     [!] === Weighted small world index calculation ===
%     [*] Markov-Kennedy SWP (n=29), weighted
%         n_samp: 10000
%         SWP: 0.426996 (+- 0.060875)(95% CI: 0.342813 - 0.504378, n=10000)
%         deltaC: 0.008692 (+- 0.007812)
%         deltaL: 0.810244 (+- 0.086288)
%         deviance: 0.980934 (+- 0.101726)(95% CI: 0.959619 - 1.000000, n=10000)
    
    SWP = zeros(1,nsamps);
    DC = zeros(1,nsamps);
    DL = zeros(1,nsamps);
    for isamp = 1:nsamps
        [swp,dc,dl] = small_world_propensity(Aswp);
        SWP(isamp) = swp;
        DC(isamp) = dc;
        DL(isamp) = dl;
    end
    Devi = (4*(atan2(DL,DC)))/pi - 1;
    sSWP = sort(SWP);
    sDevi = sort(Devi);
    SWP_mk_mean = mean(SWP);
    SWP_mk_std = std(SWP);
    SWP_mk_ci_lo = sSWP(round(ci_alpha * length(sSWP)));
    SWP_mk_ci_hi = sSWP(round((1-ci_alpha) * length(sSWP)));
    Devi_mk_mean = mean(Devi);
    Devi_mk_std = std(Devi);
    Devi_mk_ci_lo = sDevi(round(ci_alpha * length(sDevi)));
    Devi_mk_ci_hi = sDevi(round((1-ci_alpha) * length(sDevi)));
    fprintf('[*] Markov-Kennedy SWP (n=%i), weighted\n',length(Aswp));
    fprintf('\tn_samp: %i\n',nsamps);
    fprintf('\tSWP: %.6f (+- %.6f)(%.0f%% CI: %.6f - %.6f, n=%i)\n',...
        mean(SWP),std(SWP),100*(1-ci_alpha),SWP_mk_ci_lo,SWP_mk_ci_hi,length(SWP));
    fprintf('\tdeltaC: %.6f (+- %.6f)\n',mean(DC),std(DC));
    fprintf('\tdeltaL: %.6f (+- %.6f)\n',mean(DL),std(DL));
    fprintf('\tdeviance: %.6f (+- %.6f)(%.0f%% CI: %.6f - %.6f, n=%i)\n',...
        mean(Devi),std(Devi),100*(1-ci_alpha),Devi_mk_ci_lo,Devi_mk_ci_hi,length(Devi));
    
    
    % Binary Monkey
%     SWPs = [];
%     n_pts = 10;
%     max_density = sum(sum(triu(A,1)>0))/nchoosek(length(A),2);
%     Thrs = linspace(0,max_density,n_pts+2);
%     Thrs = Thrs(2:(end-1));
%     for thr = Thrs
%         %Athr = threshold_absolute(A,thr);
%         Athr = threshold_proportional(A,thr);
%         Aswp = (Athr>0);
%         nsamps = 800;
%         SWP = zeros(1,nsamps);
%         DC = zeros(1,nsamps);
%         DL = zeros(1,nsamps);
%         for isamp = 1:nsamps
%             [swp,dc,dl] = small_world_propensity(Aswp);
%             SWP(isamp) = swp;
%             DC(isamp) = dc;
%             DL(isamp) = dl;
%         end
%         Devi = (4*(atan2(DL,DC)))/pi - 1;
%         fprintf('[*] Markov-Kennedy SWP (n=%i), binary\n',length(Aswp));
%         fprintf('\tdensity: %.3f\n',thr);
%         fprintf('\tn_samp: %i\n',nsamps);
%         fprintf('\tSWP: %.6f (+- %.6f)\n',mean(SWP),std(SWP));
%         fprintf('\tdeltaC: %.6f (+- %.6f)\n',mean(DC),std(DC));
%         fprintf('\tdeltaL: %.6f (+- %.6f)\n',mean(DL),std(DL));
%         fprintf('\tdeviance: %.6f (+- %.6f)\n',mean(Devi),std(Devi));
%         SWPs = [SWPs, (SWP')];
%     end
%     
%     boxplot(SWPs);
%     box off;
%     set(gca,'TickDir','out');
%     xlabel('Density (%%)');
%     ylabel('Small-World Propensity');
%     ylim([0 1]);
%     xlabtxt = {};
%     for ixt = 1:length(Thrs)
%         xlabtxt{ixt} = sprintf('%.1f',100*Thrs(ixt));
%     end
%     xticklabels(xlabtxt);
%     
    
    %A = weight_conversion(A,'normalize');
    
    
    
    
    % === Monkey ECoG ===
    fprintf('[*] --- Monkey MK anatomy, ECoG ---\n')
    ME = load(sprintf('figure_t8d1_adj_%i',iM));
    Adj_plt_MKraw = ME.Adj_plt_MKraw;
    Adj_plt_MKraw(isnan(ME.Adj_plt_mSu_red)) = NaN;
    Adj_plt_mSu_red = ME.Adj_plt_mSu_red;
    Adj_plt_mSu_red_rois = ME.Adj_plt_mSu_red_rois;
    for idiag = 1:length(Adj_plt_MKraw)
        Adj_plt_mSu_red(idiag,idiag) = 0;
        Adj_plt_MKraw(idiag,idiag) = 0;
    end
    
    % Remove NaNs
    
    % clear diagonals
    Ahs_nodiag = Adj_plt_MKraw;
    for inod = 1:length(Ahs)
        Ahs_nodiag(inod,inod) = 0;
    end
    Adj_plt_MKraw = Ahs_nodiag;
    Ahs_nodiag = Adj_plt_mSu_red;
    for inod = 1:length(Ahs)
        Ahs_nodiag(inod,inod) = 0;
    end
    Adj_plt_mSu_red = Ahs_nodiag;
    
    n_start = length(Adj_plt_MKraw);
    cond_pass = true;
    fprintf('[*] Remove NaNs, starting nodes: %i\n',n_start)
    while(cond_pass)
        count_nan = sum(isnan(Adj_plt_MKraw));
        [~,hidx] = max(count_nan);
        nidx = true(1,length(Adj_plt_MKraw));
        nidx(hidx) = false;
        Adj_plt_MKraw = Adj_plt_MKraw(nidx,nidx);
        Adj_plt_mSu_red = Adj_plt_mSu_red(nidx,nidx);
        Adj_plt_mSu_red_rois = Adj_plt_mSu_red_rois(nidx);
        cond_pass = (sum(isnan(Adj_plt_MKraw(:))) ~= 0);
    end
    n_stop = length(Adj_plt_MKraw);
    fprintf('\tfinished nodes: %i\n',n_stop)
    
    
    
    % Weighted Monkey anatomical MK, overall (raw)
    % SWP raw
    if (iM == 1)
        save('./cache/Adj_plt_mSu_red_pCorr','Adj_plt_mSu_red','Adj_plt_mSu_red_rois');
    end
    %return
    % --------------------------------------------------------------------
    A_input = Adj_plt_mSu_red;
    nsamps = 10000;
    fprintf('[*] Calculating monkey functional interactions SWP with %i permutations..\n',nsamps);
    SWPraw = zeros(1,nsamps);
    DCraw = zeros(1,nsamps);
    DLraw = zeros(1,nsamps);
    Gammaraw = zeros(1,nsamps);
    Lambdaraw = zeros(1,nsamps);
    Sigmaraw = zeros(1,nsamps);
    Ln = zeros(1,nsamps);
    Lr = zeros(1,nsamps);
    Cn = zeros(1,nsamps);
    Cr = zeros(1,nsamps);
    Omegaraw = zeros(1,nsamps);
    SWIraw = zeros(1,nsamps);
    Qraw = zeros(1,nsamps);
    parfor isamp = 1:nsamps
        %[swp,dc,dl] = small_world_propensity(Ahs);
        [swp,dc,dl,net_clus,rand_clus,net_path,rand_path,reg_clus,reg_path,gamma,lambda,sigma] = small_world_propensity_raw(A_input);
        SWPraw(isamp) = swp;
        DCraw(isamp) = dc;
        DLraw(isamp) = dl;
        Gammaraw(isamp) = gamma;
        Lambdaraw(isamp) = lambda;
        Sigmaraw(isamp) = sigma;
        Ln(isamp) = net_path;
        Lr(isamp) = rand_path;
        Ll(isamp) = reg_path;
        Cn(isamp) = net_clus;
        Cr(isamp) = rand_clus;
        Cl(isamp) = reg_clus;
        Omegaraw(isamp) = (rand_path/net_path) - (net_clus/reg_clus);
        SWIraw(isamp) = ((net_path - reg_path)/(rand_path - reg_path)) * ((net_clus - rand_clus)/(reg_clus - rand_clus));
    end
    Deviraw = (4*(atan2(DLraw,DCraw)))/pi - 1;
    fprintf('[*] Adj, nodes=%i\n',length(A_input));
    alpha = 0.05;
    sDCraw = sort(DCraw);
    sDLraw = sort(DLraw);
    sSWPraw = sort(SWPraw);
    sDeviraw = sort(Deviraw);
    sGammaraw = sort(Gammaraw);
    sLambdaraw = sort(Lambdaraw);
    sSigmaraw = sort(Sigmaraw);
    sOmegaraw = sort(Omegaraw);
    sLn = sort(Ln);
    sLr = sort(Lr);
    sCn = sort(Cn);
    sCr = sort(Cr);
    sLl = sort(Ll);
    sCl = sort(Cl);
    sSWIraw = sort(SWIraw);
    fprintf('\tL\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sLn),std(sLn),nsamps,sLn(round((alpha/2)*nsamps)),sLn(round((1 - (alpha/2))*nsamps)));
    fprintf('\tLr\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sLr),std(sLr),nsamps,sLr(round((alpha/2)*nsamps)),sLr(round((1 - (alpha/2))*nsamps)));
    fprintf('\tLl\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sLl),std(sLl),nsamps,sLl(round((alpha/2)*nsamps)),sLl(round((1 - (alpha/2))*nsamps)));
    fprintf('\tC\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sCn),std(sCn),nsamps,sCn(round((alpha/2)*nsamps)),sCn(round((1 - (alpha/2))*nsamps)));
    fprintf('\tCr\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sCr),std(sCr),nsamps,sCr(round((alpha/2)*nsamps)),sCr(round((1 - (alpha/2))*nsamps)));
    fprintf('\tCl\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sCl),std(sCl),nsamps,sCl(round((alpha/2)*nsamps)),sCl(round((1 - (alpha/2))*nsamps)));
    
    fprintf('\tDC\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(DCraw),std(DCraw),nsamps,sDCraw(round((alpha/2)*nsamps)),sDCraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tDL\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(DLraw),std(DLraw),nsamps,sDLraw(round((alpha/2)*nsamps)),sDLraw(round((1 - (alpha/2))*nsamps)));
    
    fprintf('\tSWP\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(SWPraw),std(SWPraw),nsamps,sSWPraw(round((alpha/2)*nsamps)),sSWPraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tSWI\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(SWIraw),std(SWIraw),nsamps,sSWIraw(round((alpha/2)*nsamps)),sSWIraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tdelta\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Deviraw),std(Deviraw),nsamps,sDeviraw(round((alpha/2)*nsamps)),sDeviraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tgamma\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Gammaraw),std(Gammaraw),nsamps,sGammaraw(round((alpha/2)*nsamps)),sGammaraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tlambda\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Lambdaraw),std(Lambdaraw),nsamps,sLambdaraw(round((alpha/2)*nsamps)),sLambdaraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tsigma\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Sigmaraw),std(Sigmaraw),nsamps,sSigmaraw(round((alpha/2)*nsamps)),sSigmaraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tomega\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Omegaraw),std(Omegaraw),nsamps,sOmegaraw(round((alpha/2)*nsamps)),sOmegaraw(round((1 - (alpha/2))*nsamps)));
    % --------------------------------------------------------------------
    
    
    
    
    % === Calculate SWPs ===
    %A_queue = {Adj_plt_MKraw,Adj_plt_mSu_red};
    A_queue = {Adj_plt_mSu_red};
    SC_q = cell(1,length(A_queue));
    for iq = 1:length(A_queue)
        A = weight_conversion(A_queue{iq},'autofix');
        A = weight_conversion(A,'normalize');
        Aswp = A;
        SWPs = [];
        Devis = [];
        n_pts = 12;
        max_density = sum(sum(triu(Aswp,1)>0))/nchoosek(length(Aswp),2);
        %Thrs = linspace(0,max_density,n_pts+2);
        Thrs = linspace(0.15,max_density,n_pts+2);
        Thrs = Thrs(2:(end-1));
        fn_cache_dk = sprintf('fig_cluster2_sworld_circle_swp_dk_%i_macaque-%i_pCorr',iM,iq);
        if (~exist([fn_cache_dk,'.mat'],'file'))
            n_thr = 1;
            for thr = Thrs
                nsamps = 10000;
                SWP = zeros(1,nsamps);
                DC = zeros(1,nsamps);
                DL = zeros(1,nsamps);
                parfor isamp = 1:nsamps
                    %[swp,dc,dl] = small_world_propensity(Aswp);
                    [swp,dc,dl] = small_world_propensity(threshold_proportional(Aswp,thr));
                    SWP(isamp) = swp;
                    DC(isamp) = dc;
                    DL(isamp) = dl;
                end
                Devi = (4*(atan2(DL,DC)))/pi - 1;
                fprintf('[%4i/%4i] Macaque-%i SWP (n=%i), weighted\n',n_thr,length(Thrs),iq,length(Aswp));
                fprintf('\tthresh: %.3f\n',thr);
                fprintf('\tn_samp: %i\n',nsamps);
                fprintf('\tSWP: %.6f (+- %.6f)\n',mean(SWP),std(SWP));
                fprintf('\tdeltaC: %.6f (+- %.6f)\n',mean(DC),std(DC));
                fprintf('\tdeltaL: %.6f (+- %.6f)\n',mean(DL),std(DL));
                fprintf('\tdeviance: %.6f (+- %.6f)\n',mean(Devi),std(Devi));
                SWPs = [SWPs, (SWP')];
                Devis = [Devis, (Devi')];
                n_thr = n_thr + 1;
            end
            save(fn_cache_dk,'-v7.3','-nocompression','SWPs','Devis','Thrs');
            CaSwp = load(fn_cache_dk);
            SC_q{iq} = CaSwp;
        else
            fprintf('[*] Found cache: %s, loading..\n',fn_cache_dk);
            CaSwp = load(fn_cache_dk);
            SC_q{iq} = CaSwp;
        end
    end
    
    % select
    Cmk = SC_q{1};
    
    
    % === Human ===
    % DK atlas
    % Weighted Human
    A = weight_conversion(Ahs,'autofix');
    
    
    
    
    % SWP raw
    % --------------------------------------------------------------------
    A_input = Ahs;
    nsamps = 10000;
    fprintf('[*] Calculating human SWP with %i permutations..\n',nsamps);
    SWPraw = zeros(1,nsamps);
    DCraw = zeros(1,nsamps);
    DLraw = zeros(1,nsamps);
    Gammaraw = zeros(1,nsamps);
    Lambdaraw = zeros(1,nsamps);
    Sigmaraw = zeros(1,nsamps);
    Ln = zeros(1,nsamps);
    Lr = zeros(1,nsamps);
    Cn = zeros(1,nsamps);
    Cr = zeros(1,nsamps);
    Omegaraw = zeros(1,nsamps);
    SWIraw = zeros(1,nsamps);
    Qraw = zeros(1,nsamps);
    parfor isamp = 1:nsamps
        %[swp,dc,dl] = small_world_propensity(Ahs);
        [swp,dc,dl,net_clus,rand_clus,net_path,rand_path,reg_clus,reg_path,gamma,lambda,sigma] = small_world_propensity_raw(A_input);
        SWPraw(isamp) = swp;
        DCraw(isamp) = dc;
        DLraw(isamp) = dl;
        Gammaraw(isamp) = gamma;
        Lambdaraw(isamp) = lambda;
        Sigmaraw(isamp) = sigma;
        Ln(isamp) = net_path;
        Lr(isamp) = rand_path;
        Ll(isamp) = reg_path;
        Cn(isamp) = net_clus;
        Cr(isamp) = rand_clus;
        Cl(isamp) = reg_clus;
        Omegaraw(isamp) = (rand_path/net_path) - (net_clus/reg_clus);
        SWIraw(isamp) = ((net_path - reg_path)/(rand_path - reg_path)) * ((net_clus - rand_clus)/(reg_clus - rand_clus));
    end
    Deviraw = (4*(atan2(DLraw,DCraw)))/pi - 1;
    fprintf('[*] Adj, nodes=%i\n',length(A_input));
    alpha = 0.05;
    sDCraw = sort(DCraw);
    sDLraw = sort(DLraw);
    sSWPraw = sort(SWPraw);
    sDeviraw = sort(Deviraw);
    sGammaraw = sort(Gammaraw);
    sLambdaraw = sort(Lambdaraw);
    sSigmaraw = sort(Sigmaraw);
    sOmegaraw = sort(Omegaraw);
    sLn = sort(Ln);
    sLr = sort(Lr);
    sCn = sort(Cn);
    sCr = sort(Cr);
    sLl = sort(Ll);
    sCl = sort(Cl);
    sSWIraw = sort(SWIraw);
    fprintf('\tL\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sLn),std(sLn),nsamps,sLn(round((alpha/2)*nsamps)),sLn(round((1 - (alpha/2))*nsamps)));
    fprintf('\tLr\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sLr),std(sLr),nsamps,sLr(round((alpha/2)*nsamps)),sLr(round((1 - (alpha/2))*nsamps)));
    fprintf('\tLl\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sLl),std(sLl),nsamps,sLl(round((alpha/2)*nsamps)),sLl(round((1 - (alpha/2))*nsamps)));
    fprintf('\tC\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sCn),std(sCn),nsamps,sCn(round((alpha/2)*nsamps)),sCn(round((1 - (alpha/2))*nsamps)));
    fprintf('\tCr\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sCr),std(sCr),nsamps,sCr(round((alpha/2)*nsamps)),sCr(round((1 - (alpha/2))*nsamps)));
    fprintf('\tCl\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(sCl),std(sCl),nsamps,sCl(round((alpha/2)*nsamps)),sCl(round((1 - (alpha/2))*nsamps)));
    
    fprintf('\tDC\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(DCraw),std(DCraw),nsamps,sDCraw(round((alpha/2)*nsamps)),sDCraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tDL\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(DLraw),std(DLraw),nsamps,sDLraw(round((alpha/2)*nsamps)),sDLraw(round((1 - (alpha/2))*nsamps)));
    
    fprintf('\tSWP\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(SWPraw),std(SWPraw),nsamps,sSWPraw(round((alpha/2)*nsamps)),sSWPraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tSWI\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(SWIraw),std(SWIraw),nsamps,sSWIraw(round((alpha/2)*nsamps)),sSWIraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tdelta\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Deviraw),std(Deviraw),nsamps,sDeviraw(round((alpha/2)*nsamps)),sDeviraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tgamma\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Gammaraw),std(Gammaraw),nsamps,sGammaraw(round((alpha/2)*nsamps)),sGammaraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tlambda\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Lambdaraw),std(Lambdaraw),nsamps,sLambdaraw(round((alpha/2)*nsamps)),sLambdaraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tsigma\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Sigmaraw),std(Sigmaraw),nsamps,sSigmaraw(round((alpha/2)*nsamps)),sSigmaraw(round((1 - (alpha/2))*nsamps)));
    fprintf('\tomega\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
        mean(Omegaraw),std(Omegaraw),nsamps,sOmegaraw(round((alpha/2)*nsamps)),sOmegaraw(round((1 - (alpha/2))*nsamps)));
    % --------------------------------------------------------------------
    
    
    
    
    
    Aswp = A;
    SWPs = [];
    Devis = [];
    n_pts = 12;
    max_density = sum(sum(triu(Aswp,1)>0))/nchoosek(length(Aswp),2);
    %Thrs = linspace(0,max_density,n_pts+2);
    Thrs = linspace(0.15,max_density,n_pts+2);
    Thrs = Thrs(2:(end-1));
    fn_cache_dk = sprintf('fig_cluster2_sworld_circle_swp_dk_%i_pCorr',iM);
    if (~exist([fn_cache_dk,'.mat'],'file'))
        for thr = Thrs
            nsamps = 10000;
            SWP = zeros(1,nsamps);
            DC = zeros(1,nsamps);
            DL = zeros(1,nsamps);
            parfor isamp = 1:nsamps
                %[swp,dc,dl] = small_world_propensity(Aswp);
                [swp,dc,dl] = small_world_propensity(threshold_proportional(Ahs,thr));
                SWP(isamp) = swp;
                DC(isamp) = dc;
                DL(isamp) = dl;
            end
            Devi = (4*(atan2(DL,DC)))/pi - 1;
            fprintf('[*] Human SWP (n=%i), weighted\n',length(Aswp));
            fprintf('\tthresh: %.3f\n',thr);
            fprintf('\tn_samp: %i\n',nsamps);
            fprintf('\tSWP: %.6f (+- %.6f)\n',mean(SWP),std(SWP));
            fprintf('\tdeltaC: %.6f (+- %.6f)\n',mean(DC),std(DC));
            fprintf('\tdeltaL: %.6f (+- %.6f)\n',mean(DL),std(DL));
            fprintf('\tdeviance: %.6f (+- %.6f)\n',mean(Devi),std(Devi));
            SWPs = [SWPs, (SWP')];
            Devis = [Devis, (Devi')];
        end
        save(fn_cache_dk,'-v7.3','-nocompression','SWPs','Devis','Thrs');
    else
        fprintf('[*] Found cache: %s, loading..\n',fn_cache_dk);
        load(fn_cache_dk);
    end
    
    close all;
    
    col_ci = [1 1 1]*0.8;
    col_sw = [0.8 0.1 0.1];
    thresh_sw = 0.6;
    
    
    h = figure('Position',[0,0,1000,500]);
    set(h,'PaperUnits','inches');
    set(h,'PaperPosition',[0, 0, 8, 4]);
    subplot(2,1,1);
    % -- whisker plot --
%     boxplot(SWPs,'OutlierSize',1);
%     xlabtxt = {};
%     for ixt = 1:length(Thrs)
%         xlabtxt{ixt} = sprintf('%.0f',100*Thrs(ixt));
%     end
%     xticklabels(xlabtxt);
%     xtickangle(90);

    % -- confidence plot --
    SWPs_ci_lo = zeros(1,length(Thrs));
    SWPs_ci_hi = zeros(1,length(Thrs));
    for ip = 1:length(Thrs)
        Svec = sort(SWPs(:,ip));
        SWPs_ci_lo(ip) = Svec(round(0.5*ci_alpha*length(Svec)));
        SWPs_ci_hi(ip) = Svec(round((1 - 0.5*ci_alpha)*length(Svec)));
    end
    
    % repeat for monkey
    SWPsm_ci_lo = zeros(1,length(Thrs));
    SWPsm_ci_hi = zeros(1,length(Thrs));
    for ip = 1:length(Thrs)
        Svec = sort(Cmk.SWPs(:,ip));
        SWPsm_ci_lo(ip) = Svec(round(0.5*ci_alpha*length(Svec)));
        SWPsm_ci_hi(ip) = Svec(round((1 - 0.5*ci_alpha)*length(Svec)));
    end
    
    % Threshold SW
%     x_sws = 100*Thrs((SWPs_ci_lo > thresh_sw));
%     q = fill([x_sws(1) x_sws(end) x_sws(end) x_sws(1)],[0 0 1 1],'r',...
%         'FaceColor',col_sw,'LineStyle','none','FaceAlpha',0.1);
    
    x = 100*Thrs;
    y = median(SWPs);
    hold all;
    
    % Plot Monkey SWP
    col_cim = plasma(10);
    col_cim = col_cim(3,:);
    fill([x';flipud(x')],[SWPsm_ci_lo';flipud(SWPsm_ci_hi')],col_ci,'linestyle','none','FaceAlpha',1);
    plot(x,median(Cmk.SWPs),'-','Color',[1 1 1]*0);
    
    % Plot Human SWP
    fill([x';flipud(x')],[SWPs_ci_lo';flipud(SWPs_ci_hi')],col_cim,'linestyle','none','FaceAlpha',0.4);
    plot(x,y,'-','Color',0.4*col_cim);
    
%     plot([min(x) max(x)],[1 1]*SWP_mk_mean,'-','Color',[0 0 1]);
%     plot([min(x) max(x)],[1 1]*SWP_mk_ci_lo,'--','Color',[0 0 1]);
%     plot([min(x) max(x)],[1 1]*SWP_mk_ci_hi,'--','Color',[0 0 1]);
    %plot([min(x) max(x)],[1 1]*(SWP_mk_mean + SWP_mk_std),'--','Color',[0 0 1]);
    %plot([min(x) max(x)],[1 1]*(SWP_mk_mean - SWP_mk_std),'--','Color',[0 0 1]);
    %plot(x,mean(SWPs),'-','Color',[1 1 1]*0.6);
    
    
    %plot([min(x) max(x)],[thresh_sw thresh_sw],'--','Color',[0.9 0 0]);
    %plot(x,SWPs_ci_lo,'-','Color',col_ci_bound,'LineWidth',1);
    %plot(x,SWPs_ci_hi,'-','Color',col_ci_bound,'LineWidth',1);
    
    plot([x(1) x(end)],[1 1]*thresh_sw,'--','Color',[1 1 1]*0.7);
    
    box off;
    set(gca,'TickDir','out');
    xlabel('Density (%)');
    ylabel('Small-World Propensity');
    ylim([0 1]);
    xlim([min(100*Thrs),max(100*Thrs)]);
    set(gca,'Layer','top');
    %legend({'un','deux','trois'},'Location','BestOutside');
    fprintf('Overall SWP: %.6f (+- %.6f)\n',nanmean(SWPs(:)),nanstd(SWPs(:)));
    
    
    
    %figure;
    subplot(2,1,2);
%     boxplot(Devis,'OutlierSize',1);

    % -- confidence plot --
    Devis_ci_lo = zeros(1,length(Thrs));
    Devis_ci_hi = zeros(1,length(Thrs));
    for ip = 1:length(Thrs)
        Svec = sort(Devis(:,ip));
        Devis_ci_lo(ip) = Svec(round(0.5*ci_alpha*length(Svec)));
        Devis_ci_hi(ip) = Svec(round((1 - 0.5*ci_alpha)*length(Svec)));
    end
    
    %monkey
    Devism_ci_lo = zeros(1,length(Thrs));
    Devism_ci_hi = zeros(1,length(Thrs));
    for ip = 1:length(Thrs)
        Svec = sort(Cmk.Devis(:,ip));
        Devism_ci_lo(ip) = Svec(round(0.5*ci_alpha*length(Svec)));
        Devism_ci_hi(ip) = Svec(round((1 - 0.5*ci_alpha)*length(Svec)));
    end
    
    x = 100*Thrs;
    y = median(Devis);
    
%     % Threshold SW
%     x_sws = 100*Thrs((SWPs_ci_lo > thresh_sw));
%     q = fill([x_sws(1) x_sws(end) x_sws(end) x_sws(1)],[-1 -1 1 1],'r',...
%         'FaceColor',col_sw,'LineStyle','none','FaceAlpha',0.1);
    
    hold all;
    
    % Plot Monkey SWP
    col_cim = plasma(10);
    col_cim = col_cim(3,:);
    fill([x';flipud(x')],[Devism_ci_lo';flipud(Devism_ci_hi')],col_ci,'linestyle','none','FaceAlpha',1);
    plot(x,median(Cmk.Devis),'-','Color',[1 1 1]*0);
    
    % Plot Human
    fill([x';flipud(x')],[Devis_ci_lo';flipud(Devis_ci_hi')],col_cim,'linestyle','none','FaceAlpha',0.4);
    plot(x,y,'-','Color',0.4*col_cim); 
    %plot(x,mean(Devis),'-','Color',[1 1 1]*0.6);
    plot([min(x) max(x)],[0 0],'--','Color',[1 1 1]*0.7);
    

    box off;
    set(gca,'TickDir','out');
    xlabel('Density (%)');
    ylabel('Contribution to Deviation');
%     xlabtxt = {};
%     for ixt = 1:length(Thrs)
%         xlabtxt{ixt} = sprintf('%.0f',100*Thrs(ixt));
%     end
%     xticklabels(xlabtxt);
%     xtickangle(90);
    ylim([-1 1]);
    xlim([min(100*Thrs),max(100*Thrs)]);
    set(gca,'Layer','top');
    fprintf('Overall SWP: %.6f (+- %.6f)\n',nanmean(Devis(:)),nanstd(Devis(:)));
    
    print(h,sprintf('./figures/fig_cluster2_sworld_circle_im-%i_swp_pCorr',iM),'-depsc','-r900');
    print(h,sprintf('./figures/fig_cluster2_sworld_circle_im-%i_swp_pCorr',iM),'-dpng','-r900');
    close(h);
    %return
%     
%     % Brain metrics
%     C = clustering_coef_wu(A);
%     [D,~] = distance_wei(A);
%     Lw = charpath(D);
%     Cw = mean(C);
%     
%     % Lattice metrics
%     n_resample = 20;
%     Lw_l_s = zeros(1,n_resample);
%     Cw_l_s = zeros(1,n_resample);
%     for isamp = 1:n_resample
%         Al = weight_conversion(latticegraph(A),'autofix');
%         C = clustering_coef_wu(Al);
%         %[D,~] = distance_wei(weight_conversion(Al,'lengths'));
%         [D,~] = distance_wei(Al);
%         Lw_l_s(isamp) = charpath(D);
%         Cw_l_s(isamp) = mean(C);
%     end
%     Lw_l = mean(Lw_l_s);
%     Cw_l = mean(Cw_l_s);
%     
%     % Random metrics
%     n_resample = 20;
%     Lw_r_s = zeros(1,n_resample);
%     Cw_r_s = zeros(1,n_resample);
%     for isamp = 1:n_resample
%         Ar = weight_conversion(randomgraph(A),'autofix');
%         C = clustering_coef_wu(Ar);
%         %[D,~] = distance_wei(weight_conversion(Ar,'lengths'));
%         [D,~] = distance_wei(Ar);
%         Lw_r_s(isamp) = charpath(D);
%         Cw_r_s(isamp) = mean(C);
%     end
%     Lw_r = mean(Lw_r_s);
%     Cw_r = mean(Cw_r_s);
%     
%     gammaw = Cw / Cw_r;
%     lambdaw = Lw / Lw_r;
%     sigmaw = gammaw/lambdaw;
%     deltaw_C = (Cw_l - Cw)/(Cw_l - Cw_r);
%     deltaw_L = (Lw - Lw_r)/(Lw_l - Lw_r);
%     phiw = 1 - sqrt((deltaw_C^2 + deltaw_L^2)/2);
%     fprintf('Gamma: %.6f\n',gammaw);
%     fprintf('Lambda: %.6f\n',lambdaw);
%     fprintf('Sigma: %.6f\n',sigmaw);
%     fprintf('phi: %.6f\n',phiw);
%     
%     
%     return
end

% 23 oct. 2020

% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 0.441606367176 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating monkey MK-FLNe SWP with 10000 permutations..
% [*] Adj, nodes=91
% 	L	269.980605, std 5.20e-11, n=10000, 95%CI: 2.700e+02 - 2.700e+02
% 	Lr	71.949370, std 3.15e+00, n=10000, 95%CI: 6.635e+01 - 7.857e+01
% 	Ll	280.193482, std 4.26e+00, n=10000, 95%CI: 2.716e+02 - 2.883e+02
% 	C	0.002804, std 4.58e-16, n=10000, 95%CI: 2.804e-03 - 2.804e-03
% 	Cr	0.000755, std 2.81e-05, n=10000, 95%CI: 7.030e-04 - 8.127e-04
% 	Cl	0.002825, std 5.07e-06, n=10000, 95%CI: 2.816e-03 - 2.836e-03
% 	DC	0.010164, std 2.43e-03, n=10000, 95%CI: 5.810e-03 - 1.531e-02
% 	DL	0.951275, std 1.93e-02, n=10000, 95%CI: 9.151e-01 - 9.918e-01
% 	SWP	0.327306, std 1.37e-02, n=10000, 95%CI: 2.986e-01 - 3.528e-01
% 	SWI	0.048155, std 1.93e-02, n=10000, 95%CI: 8.104e-03 - 8.387e-02
% 	delta	0.986387, std 3.28e-03, n=10000, 95%CI: 9.794e-01 - 9.923e-01
% 	gamma	3.718142, std 1.38e-01, n=10000, 95%CI: 3.450e+00 - 3.988e+00
% 	lambda	3.759403, std 1.61e-01, n=10000, 95%CI: 3.436e+00 - 4.069e+00
% 	sigma	0.990737, std 5.42e-02, n=10000, 95%CI: 8.900e-01 - 1.101e+00
% 	omega	-0.726055, std 1.18e-02, n=10000, 95%CI: -7.472e-01 - -7.010e-01
% [*] Markov-Kennedy SWP (n=91), weighted
% 	n_samp: 100
% 	SWP: 0.326849 (+- 0.012862)(95% CI: 0.303759 - 0.347161, n=100)
% 	deltaC: 0.009805 (+- 0.002113)
% 	deltaL: 0.951927 (+- 0.018190)
% 	deviance: 0.986882 (+- 0.002821)(95% CI: 0.981056 - 0.990796, n=100)
% [*] --- Monkey MK anatomy, ECoG ---
% [*] Remove NaNs, starting nodes: 39
% 	finished nodes: 18
% [*] Calculating monkey functional interactions SWP with 10000 permutations..
% [*] Adj, nodes=18
% 	L	8.815830, std 1.28e-12, n=10000, 95%CI: 8.816e+00 - 8.816e+00
% 	Lr	8.497218, std 2.95e-02, n=10000, 95%CI: 8.442e+00 - 8.556e+00
% 	Ll	8.848081, std 7.09e-03, n=10000, 95%CI: 8.834e+00 - 8.862e+00
% 	C	0.567225, std 4.83e-14, n=10000, 95%CI: 5.672e-01 - 5.672e-01
% 	Cr	0.477155, std 1.26e-02, n=10000, 95%CI: 4.550e-01 - 5.039e-01
% 	Cl	0.492266, std 5.66e-05, n=10000, 95%CI: 4.922e-01 - 4.924e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.907772, std 2.03e-02, n=10000, 95%CI: 8.670e-01 - 9.472e-01
% 	SWP	0.358108, std 1.43e-02, n=10000, 95%CI: 3.302e-01 - 3.869e-01
% 	SWI	0.420172, std 1.27e+01, n=10000, 95%CI: -3.601e+00 - 4.242e+00
% 	delta	1.000000, std 00, n=10000, 95%CI: 001 - 001
% 	gamma	1.189580, std 3.10e-02, n=10000, 95%CI: 1.126e+00 - 1.247e+00
% 	lambda	1.037508, std 3.60e-03, n=10000, 95%CI: 1.030e+00 - 1.044e+00
% 	sigma	1.146555, std 2.90e-02, n=10000, 95%CI: 1.087e+00 - 1.200e+00
% 	omega	-0.188413, std 3.35e-03, n=10000, 95%CI: -1.947e-01 - -1.817e-01
% [*] Found cache: fig_cluster2_sworld_circle_swp_dk_1_macaque-1_pCorr, loading..
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=20
% 	L	3.271470, std 5.26e-13, n=10000, 95%CI: 3.271e+00 - 3.271e+00
% 	Lr	3.068406, std 1.28e-02, n=10000, 95%CI: 3.048e+00 - 3.097e+00
% 	Ll	2.987090, std 1.40e-03, n=10000, 95%CI: 2.984e+00 - 2.990e+00
% 	C	0.711050, std 1.54e-14, n=10000, 95%CI: 7.111e-01 - 7.111e-01
% 	Cr	0.377589, std 2.10e-02, n=10000, 95%CI: 3.393e-01 - 4.234e-01
% 	Cl	0.533164, std 4.68e-05, n=10000, 95%CI: 5.331e-01 - 5.333e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	-2.578043, std 5.28e-01, n=10000, 95%CI: -3.686e+00 - -1.576e+00
% 	SWP	-0.822952, std 3.74e-01, n=10000, 95%CI: -1.606e+00 - -1.147e-01
% 	SWI	7.714067, std 1.05e+00, n=10000, 95%CI: 5.817e+00 - 9.961e+00
% 	delta	-3.000000, std 00, n=10000, 95%CI: -003 - -003
% 	gamma	1.888932, std 1.04e-01, n=10000, 95%CI: 1.679e+00 - 2.096e+00
% 	lambda	1.066197, std 4.43e-03, n=10000, 95%CI: 1.056e+00 - 1.073e+00
% 	sigma	1.771508, std 9.51e-02, n=10000, 95%CI: 1.580e+00 - 1.960e+00
% 	omega	-0.395713, std 3.92e-03, n=10000, 95%CI: -4.019e-01 - -3.869e-01
% [*] Found cache: fig_cluster2_sworld_circle_swp_dk_1_pCorr, loading..
% Overall SWP: 0.621409 (+- 0.290804)
% Overall SWP: -0.126046 (+- 0.371051)
