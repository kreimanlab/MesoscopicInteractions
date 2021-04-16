close all;
clear;

addpath(genpath('SWP'));



% The Watts-strogatz definition for small world network is:
%   L >= L_random but C >> C_random
%           ws    - Watts-Strogatz 1998 (C >> Cr,L >= Lr)
%           sigma - Humphries-Gurney,2008 (sigma > 1 for small world network)
%           omega - Telesford-Laurienti,2011 (0 to 1: 1 is most small world)
%           I     - Neal,2017 (0 to 1: 1 is most small world)



subjects_dirL = '/mnt/cuenap_ssd/coregistration';

SubjectsL = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124',...
    'mSu'};

% Exclude monkey
SubjectsL = SubjectsL(1:(end-1));

%SubjectsL = {'m00001'};

% Slow i/o definitions
%dir_h5L = '/media/klab/KLAB101/h5_notch20';
dir_cacheL = './cache';
cp_thresh_override = 0.05;
n_pairs_thresh = 0; % 20 at least this many electrode pairs to be considered
n_subs_thresh = 0; 
font_size = 10;
trig_mag_no_cp = true;
trig_saveplot = false;

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'}; % 'pcDelta'
metrics_suffix = {'0.5-125 Hz','3-8 Hz','8-12 Hz','12-30 Hz','30-100 Hz'};

T = nan(length(SubjectsL),15,5);

for iM = 1 %1:length(metrics)
    for iSub = 1:length(SubjectsL)
        fprintf('[*] Subject %i\n',iSub)
        for atlM = 2 %1:20
        
            sid = SubjectsL{iSub};
            sid_int = iSub;

            metric = metrics{iM};

            % Load human cache
            sid_sav = sid;
            iSub_sav = iSub;
            %load(sprintf('%s/xsub_out_all_%i.mat',dir_cacheL,iM));
            load(sprintf('%s/xsub_out_%s_%i_atl%i_pCorr.mat',dir_cacheL,SubjectsL{iSub},iM,atlM));
            n_rois_sav = n_rois;
            n_rois = ecog.n_bchan;
            sid = sid_sav;
            iSub = iSub_sav;

            % Calculate final functional interaction matrix
            Adj = nan(n_rois,n_rois);
            AdjMag = nan(n_rois,n_rois);
            AdjVar = nan(n_rois,n_rois);
            AdjNpairs = nan(n_rois,n_rois);
            AdjNpairs_sig = nan(n_rois,n_rois);
            AdjMag4cl = nan(n_rois,n_rois);
            AdjCP = nan(n_rois,n_rois);
            N_bchan = nan(n_rois,n_rois);
            cp_thresh = cp_thresh_override;
            Dmat = Inf(n_rois,n_rois); % set to inf to avoid removing from every instance

            count = 1;
            for i1 = 1:(n_rois-1)
                for i2 = (i1+1):n_rois
                    %AA = AdjAtl{i1,i2};
                    if (mag(count) > coh_thresh(count))
                        AA = mag(count);
                    else
                        AA = 0;
                    end
                    %AA_dist = adjct_dist{i1,i2};
                    %AA_sub = AdjAtl_sid{i1,i2};
                    %n_pairs = length(AA_sub);
                    %n_pairs = length(AA(AA_sub == iSub));
                    %n_subs = length(unique(AA_sub));

                    % index by subject
                    %AA = AA(AdjAtl_sid{i1,i2}==iSub);

                    % ROI pair coverage condition (grey)
                    %if ( (~ isempty(AA))  && (n_pairs >= n_pairs_thresh) && (n_subs >= n_subs_thresh))
                    if (~ isempty(AA))
                        frac_cp = sum(AA ~= 0)/length(AA);
                        N_bchan(i1,i2) = length(AA);
                        AdjCP(i1,i2) = frac_cp;
                        
                        N_bchan(i2,i1) = length(AA);
                        AdjCP(i2,i1) = frac_cp;
                        % if ( (frac_cp > cp_thresh) && (n_subs_ct >= n_subs_ct_thresh) )
                        if (frac_cp > cp_thresh)
                            % ROI pair significance condition (white)
                            Adj(i1,i2) = 1;
                            AdjMag(i1,i2) = mean(AA(AA ~= 0));
                            AdjVar(i1,i2) = var(AA(AA ~= 0));
                            AdjNpairs(i1,i2) = 1; %n_pairs;
                            AdjNpairs_sig(i1,i2) = length(AA(AA ~= 0));
                            
                            Adj(i2,i1) = 1;
                            AdjMag(i2,i1) = mean(AA(AA ~= 0));
                            AdjVar(i2,i1) = var(AA(AA ~= 0));
                            AdjNpairs(i2,i1) = 1; %n_pairs;
                            AdjNpairs_sig(i2,i1) = length(AA(AA ~= 0));
                        else
                            Adj(i1,i2) = 0;
                            AdjMag(i1,i2) = 0;
                            AdjVar(i1,i2) = 0;
                            
                            Adj(i2,i1) = 0;
                            AdjMag(i2,i1) = 0;
                            AdjVar(i2,i1) = 0;
                        end
                    end

                    AdjMag4cl(i1,i2) = mean(AA(AA~=0));
                    AdjMag4cl(i2,i1) = mean(AA(AA~=0));
                    
                    count = count + 1;
                end
            end

            frac_xsub_msu = sum(Adj(~isnan(Adj))) /numel(Adj(~isnan(Adj)));
            fprintf('%s - human fraction of ROIs significant: %.4f\n',metric,frac_xsub_msu)

            
            
            % --- PLOT  --------------------------------------------------
            colordef_adj;
            color_nocov = COLOR_ADJ_NOCOV;
            color_distt = color_nocov;
            color_isnan = color_nocov;
            color_not_sig0 = COLOR_ADJ_NOSIG;

            %color_distt = 0.6*[1 1 1];
            %color_isnan = 0.6*[1 1 1];
            fontsz = font_size; %10;
            %fsize = fontsz;

            % Calculate rows of nans
            Adj2 = Adj;
            ind_isnan = false(length(Adj2),1);
            for i = 1:length(ind_isnan)
                ind_isnan(i) = (sum(isnan(Adj2(i,:))) == length(ind_isnan));
            end

            for iii = 3 %1:3

                % Make figure
                if (trig_saveplot)
                    h = figure('visible','off');
                    %set(h,'Position',round(1*[0 0 0.95*1080 0.8*1080]));
                    set(h,'PaperUnits','Inches');
                    set(h,'PaperPosition',[0 0 8.5 6.8]);
                end
                %subplot(2,1,1);

                if (iii == 1)
                    if (trig_mag_no_cp)
                        Adj_plt = AdjMag;
                    else
                        Adj_plt = AdjCP;
                    end
                elseif (iii == 2)
                    Adj_plt = Adj;
                elseif (iii == 3)
                    if (trig_mag_no_cp)
                        Adj_plt = AdjMag;
                    else
                        Adj_plt = AdjCP;
                    end
                    %Adj_plt = AdjCP;
                    %Adj_plt(Adj == 0) = 0;
                    color_not_sig = color_not_sig0; %0.999*[1 1 1];
                end
                Dmat_plt = Dmat;
                rois_plt = C.EleLabels(ecog.bip(:,1));  %rois;

                % for hierarchichal clustering
                Adj_plt2_cl = AdjMag4cl;
                Adj_plt_var = AdjVar;
                Adj_plt_npairs = AdjNpairs;
                Adj_plt_npairs_sig = AdjNpairs_sig;
                
                
                rois = C.EleLabels(ecog.bip(:,1));
                % remove unknown areas
                cond_empty = false(length(rois),1);
                for ce = 1:length(cond_empty)
                    elem = rois{ce};
                    cond_empty(ce) = (all(isspace(elem)) | isempty(elem));
                end
                cond_contains = (contains(lower(rois),'unknown')) | (contains(lower(rois),'???')) | (contains(lower(rois),'wall')) | (cond_empty);
                known_idx = (~ cond_contains);
                Adj_plt = Adj_plt(known_idx,known_idx);
                Adj_plt2_cl = Adj_plt2_cl(known_idx,known_idx);
                Adj_plt_var = Adj_plt_var(known_idx,known_idx);
                Adj_plt_npairs = Adj_plt_npairs(known_idx,known_idx);
                Adj_plt_npairs_sig = Adj_plt_npairs_sig(known_idx,known_idx);
                rois_plt = rois_plt(known_idx);
                Dmat_plt = Dmat_plt(known_idx,known_idx);


                % Filter out nans nodes
                cov_idx = ~ all(isnan(Adj_plt));
                Adj_plt = Adj_plt(cov_idx,cov_idx);
                Adj_plt2_cl = Adj_plt2_cl(cov_idx,cov_idx);
                Adj_plt_var = Adj_plt_var(cov_idx,cov_idx);
                Adj_plt_npairs = Adj_plt_npairs(cov_idx,cov_idx);
                Adj_plt_npairs_sig = Adj_plt_npairs_sig(cov_idx,cov_idx);
                rois_plt = rois_plt(cov_idx);
                Dmat_plt = Dmat_plt(cov_idx,cov_idx);

                if ((length(rois_plt) <= 1) || (nanmin(Adj_plt(:)) == nanmax(Adj_plt(:))))
                    fprintf(2,'[!] No matrix to plot, skipped:\n');
                    fprintf(2,'\t%s',sprintf('./brainexport/controls_data_sub-%i_freq-%i_atl-%i\n',iSub,iM,atl))
                else




                fprintf('=== Frequency band: %i ===\n',iM);

                % Number of random initializations for small world index (default: 12)
                n_MC = 12;

                % Atlas index
                atl = 2; % 2 - Desikan-Killiany

                % Load human adjacency matrix
                %fn_ca3 = './cache/fig_cluster2_fill_1.mat';
                %Ca3 = load(fn_ca3);

                %fn_ca = sprintf('./cache/xsub_out_all_%i_atl2.mat',iM);
                fn_ca = sprintf('./cache/xsub_out_m00005_%i_atl2_pCorr.mat',iM);
                Ca = load(fn_ca);
                Ahs = nan(Ca.ecog.n_bchan,Ca.ecog.n_bchan);
                %An = Ahs;
                %Ausub = Ahs;
                %Ad = Ahs;
                for i = 1:Ca.ecog.n_bchan
                    for j = 1:Ca.ecog.n_bchan
                        if (~isempty(Ca.AdjMag(i,j)))
                            lt = Ca.AdjMag(i,j);
                            li = isnan(lt) | (lt == 0);
                            Ahs(i,j) = mean(lt(~li));
                            %An(i,j) = length(Ca.AdjAtlN{i,j});
                            %Ausub(i,j) = length(unique(Ca.AdjAtl_sid{i,j}));
                            %Ad(i,j) = mean(Ca.adjct_dist{i,j});
                        end
                    end
                end

                % Remove unknown node
                rois = Ca.C.EleLabels(Ca.ecog.bip(:,1));
                %rois = Ca.rois;
                %unk_i = strcmpi(rois,'UNKNOWN');
                %rois = rois(~unk_i);
                %Ahs = Ahs(~unk_i,~unk_i);
                %An = An(~unk_i,~unk_i);
                %Ausub = Ausub(~unk_i,~unk_i);
                %Ad = Ad(~unk_i,~unk_i);


                % load cache
                %fn_ca4 = sprintf('./cache/figure_t14_%i',iM);
                %Ca4 = load(fn_ca4);
                rois2 = rois_plt; %Ca4.rois_plt(Ca4.cluster_i);
                Ahs = Adj_plt; %Ca4.Adj_plt2(Ca4.cluster_i,Ca4.cluster_i);
                Ahs(isnan(Ahs)) = 0;
                Ahs_dist = Dmat_plt; %Ca.Dmat; %Ca4.Adj_dist(Ca4.cluster_i,Ca4.cluster_i);

%                 for ir2 = 1:length(rois2)
%                     rois2{ir2} = convertRoiDK(rois2{ir2});
%                 end


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
%                 NodeRois = Ca.C.AtlROIs{atl}.LH.struct_names;
%                 NodeRoisDK = cell(size(NodeRois));
%                 for i = 1:length(NodeRois)
%                     NodeRoisDK{i} = convertRoiDK(NodeRois{i});
%                 end
                NodeRois = rois2;
                NodeRoisDK = NodeRois;
%                 rois2_col = zeros(length(rois2),3);
%                 for i = 1:length(rois2)
%                     % Find roi
%                     sIdx = strcmpi(NodeRoisDK,rois2{i});
%                     rois2_col(i,:) = NodeColors(sIdx,:);
%                 end



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
%                 def_rois_short;
%                 rois2_ab = cell(size(rois2));
%                 for i = 1:length(rois2)
%                     rois2_ab{i} = rois_short{strcmp(NodeRois,rois2_raw{i})};
%                 end
                rois2_ab = rois2_raw;

                % % network plot
                % h = figure;
                % Lcolor = [0 0 0];
                % hg = plot(G,'Layout','force','EdgeColor',Lcolor,'NodeColor',rois2_col,'NodeLabel',[]); %'Layout','force','EdgeLabel',G.Edges.Weight
                % %layout(hg,'force','UseGravity',true);
                % %layout(hg,'force','WeightEffect','direct');
                % axis off;



                % --- Random --------------------------------------------------------------
                n_nodes = length(Agr);
                if (n_nodes >= 2)
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
                    save('fig_cluster2_sworld_circle_E_bip_pCorr','E');
                else
                    % Load existing randomly initialized network
                    Eca = load('fig_cluster2_sworld_circle_E_bip_pCorr');
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
                if (trig_saveplot)
                    col_circ = repmat(0.1*[1 1 1],length(rois2),1);
                    h = figure('Position',round([0, 0, 8*(6/3), 4]*100));
                    set(h,'PaperUnits','inches');
                    set(h,'PaperPosition',[0, 0, 8*(6/3), 4]);
                    addpath('circularGraph');
                    [ha, pos] = tight_subplot(1,4,[.01 .01],[.01 .01],[.01 .01]);
                    try
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
                        print(h,sprintf('./figures/sworld_circle/fig_cluster2_sworld_circle_metric-%i_sub%i_pCorr',iM,iSub),'-depsc');
                        print(h,sprintf('./figures/sworld_circle/fig_cluster2_sworld_circle_metric-%i_sub%i_pCorr',iM,iSub),'-dpng','-r900');
                        close(h);
                    catch e
                        rethrow(e);
                        fprintf('[!] Error: circle plot skipped.\n');
                    end
                end


                


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
%                     Ca = load('./cache/xsub_out_all_1.mat');
%                     %iM = 1;
%                     for i = 1:length(Ca.Subjects)
%                         sid = Ca.Subjects{i};
%                         % not used?
%                         Cas = load(sprintf('./cache/xsub_out_%s_%i.mat',sid,iM));
%                     end

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

                save(sprintf('cache/fig_cluster2_sworld_circle_subjects_bip-%i_pCorr.mat',iM));

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
%                 A_input = A;
%                 nsamps = 10000;
%                 fprintf('[*] Calculating monkey MK-FLNe SWP with %i permutations..\n',nsamps);
%                 SWPraw = zeros(1,nsamps);
%                 DCraw = zeros(1,nsamps);
%                 DLraw = zeros(1,nsamps);
%                 Gammaraw = zeros(1,nsamps);
%                 Lambdaraw = zeros(1,nsamps);
%                 Sigmaraw = zeros(1,nsamps);
%                 Ln = zeros(1,nsamps);
%                 Lr = zeros(1,nsamps);
%                 Ll = zeros(1,nsamps);
%                 Cn = zeros(1,nsamps);
%                 Cr = zeros(1,nsamps);
%                 Cl = zeros(1,nsamps);
%                 Omegaraw = zeros(1,nsamps);
%                 SWIraw = zeros(1,nsamps);
%                 Qraw = zeros(1,nsamps);
%                 parfor isamp = 1:nsamps
%                     %[swp,dc,dl] = small_world_propensity(Ahs);
%                     [swp,dc,dl,net_clus,rand_clus,net_path,rand_path,reg_clus,reg_path,gamma,lambda,sigma] = small_world_propensity_raw(A_input);
%                     SWPraw(isamp) = swp;
%                     DCraw(isamp) = dc;
%                     DLraw(isamp) = dl;
%                     Gammaraw(isamp) = gamma;
%                     Lambdaraw(isamp) = lambda;
%                     Sigmaraw(isamp) = sigma;
%                     Ln(isamp) = net_path;
%                     Lr(isamp) = rand_path;
%                     Ll(isamp) = reg_path;
%                     Cn(isamp) = net_clus;
%                     Cr(isamp) = rand_clus;
%                     Cl(isamp) = reg_clus;
%                     Omegaraw(isamp) = (rand_path/net_path) - (net_clus/reg_clus);
%                     SWIraw(isamp) = ((net_path - reg_path)/(rand_path - reg_path)) * ((net_clus - rand_clus)/(reg_clus - rand_clus));
%                 end
%                 Deviraw = (4*(atan2(DLraw,DCraw)))/pi - 1;
%                 fprintf('[*] Adj, nodes=%i\n',length(A_input));
%                 alpha = 0.05;
%                 sDCraw = sort(DCraw);
%                 sDLraw = sort(DLraw);
%                 sSWPraw = sort(SWPraw);
%                 sDeviraw = sort(Deviraw);
%                 sGammaraw = sort(Gammaraw);
%                 sLambdaraw = sort(Lambdaraw);
%                 sSigmaraw = sort(Sigmaraw);
%                 sOmegaraw = sort(Omegaraw);
%                 sLn = sort(Ln);
%                 sLr = sort(Lr);
%                 sCn = sort(Cn);
%                 sCr = sort(Cr);
%                 sLl = sort(Ll);
%                 sCl = sort(Cl);
%                 sSWIraw = sort(SWIraw);
%                 fprintf('\tL\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(sLn),std(sLn),nsamps,sLn(round((alpha/2)*nsamps)),sLn(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tLr\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(sLr),std(sLr),nsamps,sLr(round((alpha/2)*nsamps)),sLr(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tLl\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(sLl),std(sLl),nsamps,sLl(round((alpha/2)*nsamps)),sLl(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tC\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(sCn),std(sCn),nsamps,sCn(round((alpha/2)*nsamps)),sCn(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tCr\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(sCr),std(sCr),nsamps,sCr(round((alpha/2)*nsamps)),sCr(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tCl\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(sCl),std(sCl),nsamps,sCl(round((alpha/2)*nsamps)),sCl(round((1 - (alpha/2))*nsamps)));
% 
%                 fprintf('\tDC\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(DCraw),std(DCraw),nsamps,sDCraw(round((alpha/2)*nsamps)),sDCraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tDL\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(DLraw),std(DLraw),nsamps,sDLraw(round((alpha/2)*nsamps)),sDLraw(round((1 - (alpha/2))*nsamps)));
% 
%                 fprintf('\tSWP\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(SWPraw),std(SWPraw),nsamps,sSWPraw(round((alpha/2)*nsamps)),sSWPraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tSWI\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(SWIraw),std(SWIraw),nsamps,sSWIraw(round((alpha/2)*nsamps)),sSWIraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tdelta\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(Deviraw),std(Deviraw),nsamps,sDeviraw(round((alpha/2)*nsamps)),sDeviraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tgamma\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(Gammaraw),std(Gammaraw),nsamps,sGammaraw(round((alpha/2)*nsamps)),sGammaraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tlambda\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(Lambdaraw),std(Lambdaraw),nsamps,sLambdaraw(round((alpha/2)*nsamps)),sLambdaraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tsigma\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(Sigmaraw),std(Sigmaraw),nsamps,sSigmaraw(round((alpha/2)*nsamps)),sSigmaraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tomega\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(Omegaraw),std(Omegaraw),nsamps,sOmegaraw(round((alpha/2)*nsamps)),sOmegaraw(round((1 - (alpha/2))*nsamps)));
%                 % --------------------------------------------------------------------
% 
% 
% 
%                 % Weighted Monkey
%                 ci_alpha = 0.05;
%                 Aswp = A;
%                 nsamps = 100;
% 
%             %     [!] === Weighted small world index calculation ===
%             %     [*] Markov-Kennedy SWP (n=29), weighted
%             %         n_samp: 10000
%             %         SWP: 0.426996 (+- 0.060875)(95% CI: 0.342813 - 0.504378, n=10000)
%             %         deltaC: 0.008692 (+- 0.007812)
%             %         deltaL: 0.810244 (+- 0.086288)
%             %         deviance: 0.980934 (+- 0.101726)(95% CI: 0.959619 - 1.000000, n=10000)
% 
%                 SWP = zeros(1,nsamps);
%                 DC = zeros(1,nsamps);
%                 DL = zeros(1,nsamps);
%                 for isamp = 1:nsamps
%                     [swp,dc,dl] = small_world_propensity(Aswp);
%                     SWP(isamp) = swp;
%                     DC(isamp) = dc;
%                     DL(isamp) = dl;
%                 end
%                 Devi = (4*(atan2(DL,DC)))/pi - 1;
%                 sSWP = sort(SWP);
%                 sDevi = sort(Devi);
%                 SWP_mk_mean = mean(SWP);
%                 SWP_mk_std = std(SWP);
%                 SWP_mk_ci_lo = sSWP(round(ci_alpha * length(sSWP)));
%                 SWP_mk_ci_hi = sSWP(round((1-ci_alpha) * length(sSWP)));
%                 Devi_mk_mean = mean(Devi);
%                 Devi_mk_std = std(Devi);
%                 Devi_mk_ci_lo = sDevi(round(ci_alpha * length(sDevi)));
%                 Devi_mk_ci_hi = sDevi(round((1-ci_alpha) * length(sDevi)));
%                 fprintf('[*] Markov-Kennedy SWP (n=%i), weighted\n',length(Aswp));
%                 fprintf('\tn_samp: %i\n',nsamps);
%                 fprintf('\tSWP: %.6f (+- %.6f)(%.0f%% CI: %.6f - %.6f, n=%i)\n',...
%                     mean(SWP),std(SWP),100*(1-ci_alpha),SWP_mk_ci_lo,SWP_mk_ci_hi,length(SWP));
%                 fprintf('\tdeltaC: %.6f (+- %.6f)\n',mean(DC),std(DC));
%                 fprintf('\tdeltaL: %.6f (+- %.6f)\n',mean(DL),std(DL));
%                 fprintf('\tdeviance: %.6f (+- %.6f)(%.0f%% CI: %.6f - %.6f, n=%i)\n',...
%                     mean(Devi),std(Devi),100*(1-ci_alpha),Devi_mk_ci_lo,Devi_mk_ci_hi,length(Devi));
% 
% 
%                 % Binary Monkey
%             %     SWPs = [];
%             %     n_pts = 10;
%             %     max_density = sum(sum(triu(A,1)>0))/nchoosek(length(A),2);
%             %     Thrs = linspace(0,max_density,n_pts+2);
%             %     Thrs = Thrs(2:(end-1));
%             %     for thr = Thrs
%             %         %Athr = threshold_absolute(A,thr);
%             %         Athr = threshold_proportional(A,thr);
%             %         Aswp = (Athr>0);
%             %         nsamps = 800;
%             %         SWP = zeros(1,nsamps);
%             %         DC = zeros(1,nsamps);
%             %         DL = zeros(1,nsamps);
%             %         for isamp = 1:nsamps
%             %             [swp,dc,dl] = small_world_propensity(Aswp);
%             %             SWP(isamp) = swp;
%             %             DC(isamp) = dc;
%             %             DL(isamp) = dl;
%             %         end
%             %         Devi = (4*(atan2(DL,DC)))/pi - 1;
%             %         fprintf('[*] Markov-Kennedy SWP (n=%i), binary\n',length(Aswp));
%             %         fprintf('\tdensity: %.3f\n',thr);
%             %         fprintf('\tn_samp: %i\n',nsamps);
%             %         fprintf('\tSWP: %.6f (+- %.6f)\n',mean(SWP),std(SWP));
%             %         fprintf('\tdeltaC: %.6f (+- %.6f)\n',mean(DC),std(DC));
%             %         fprintf('\tdeltaL: %.6f (+- %.6f)\n',mean(DL),std(DL));
%             %         fprintf('\tdeviance: %.6f (+- %.6f)\n',mean(Devi),std(Devi));
%             %         SWPs = [SWPs, (SWP')];
%             %     end
%             %     
%             %     boxplot(SWPs);
%             %     box off;
%             %     set(gca,'TickDir','out');
%             %     xlabel('Density (%%)');
%             %     ylabel('Small-World Propensity');
%             %     ylim([0 1]);
%             %     xlabtxt = {};
%             %     for ixt = 1:length(Thrs)
%             %         xlabtxt{ixt} = sprintf('%.1f',100*Thrs(ixt));
%             %     end
%             %     xticklabels(xlabtxt);
%             %     
% 
%                 %A = weight_conversion(A,'normalize');
% 
% 
% 
% 
%                 % === Monkey ECoG ===
%                 fprintf('[*] --- Monkey MK anatomy, ECoG ---\n')
%                 ME = load(sprintf('figure_t8d1_adj_%i',iM));
%                 Adj_plt_MKraw = ME.Adj_plt_MKraw;
%                 Adj_plt_MKraw(isnan(ME.Adj_plt_mSu_red)) = NaN;
%                 Adj_plt_mSu_red = ME.Adj_plt_mSu_red;
%                 Adj_plt_mSu_red_rois = ME.Adj_plt_mSu_red_rois;
%                 for idiag = 1:length(Adj_plt_MKraw)
%                     Adj_plt_mSu_red(idiag,idiag) = 0;
%                     Adj_plt_MKraw(idiag,idiag) = 0;
%                 end
% 
%                 % Remove NaNs
% 
%                 % clear diagonals
%                 Ahs_nodiag = Adj_plt_MKraw;
%                 for inod = 1:length(Ahs)
%                     Ahs_nodiag(inod,inod) = 0;
%                 end
%                 Adj_plt_MKraw = Ahs_nodiag;
%                 Ahs_nodiag = Adj_plt_mSu_red;
%                 for inod = 1:length(Ahs)
%                     Ahs_nodiag(inod,inod) = 0;
%                 end
%                 Adj_plt_mSu_red = Ahs_nodiag;
% 
%                 n_start = length(Adj_plt_MKraw);
%                 cond_pass = true;
%                 fprintf('[*] Remove NaNs, starting nodes: %i\n',n_start)
%                 while(cond_pass)
%                     count_nan = sum(isnan(Adj_plt_MKraw));
%                     [~,hidx] = max(count_nan);
%                     nidx = true(1,length(Adj_plt_MKraw));
%                     nidx(hidx) = false;
%                     Adj_plt_MKraw = Adj_plt_MKraw(nidx,nidx);
%                     Adj_plt_mSu_red = Adj_plt_mSu_red(nidx,nidx);
%                     Adj_plt_mSu_red_rois = Adj_plt_mSu_red_rois(nidx);
%                     cond_pass = (sum(isnan(Adj_plt_MKraw(:))) ~= 0);
%                 end
%                 n_stop = length(Adj_plt_MKraw);
%                 fprintf('\tfinished nodes: %i\n',n_stop)
% 
% 
% 
%                 % Weighted Monkey anatomical MK, overall (raw)
%                 % SWP raw
%                 if (iM == 1)
%                     save('./cache/Adj_plt_mSu_red','Adj_plt_mSu_red','Adj_plt_mSu_red_rois');
%                 end
%                 %return
%                 % --------------------------------------------------------------------
%                 A_input = Adj_plt_mSu_red;
%                 nsamps = 10000;
%                 fprintf('[*] Calculating monkey functional interactions SWP with %i permutations..\n',nsamps);
%                 SWPraw = zeros(1,nsamps);
%                 DCraw = zeros(1,nsamps);
%                 DLraw = zeros(1,nsamps);
%                 Gammaraw = zeros(1,nsamps);
%                 Lambdaraw = zeros(1,nsamps);
%                 Sigmaraw = zeros(1,nsamps);
%                 Ln = zeros(1,nsamps);
%                 Lr = zeros(1,nsamps);
%                 Cn = zeros(1,nsamps);
%                 Cr = zeros(1,nsamps);
%                 Omegaraw = zeros(1,nsamps);
%                 SWIraw = zeros(1,nsamps);
%                 Qraw = zeros(1,nsamps);
%                 parfor isamp = 1:nsamps
%                     %[swp,dc,dl] = small_world_propensity(Ahs);
%                     [swp,dc,dl,net_clus,rand_clus,net_path,rand_path,reg_clus,reg_path,gamma,lambda,sigma] = small_world_propensity_raw(A_input);
%                     SWPraw(isamp) = swp;
%                     DCraw(isamp) = dc;
%                     DLraw(isamp) = dl;
%                     Gammaraw(isamp) = gamma;
%                     Lambdaraw(isamp) = lambda;
%                     Sigmaraw(isamp) = sigma;
%                     Ln(isamp) = net_path;
%                     Lr(isamp) = rand_path;
%                     Ll(isamp) = reg_path;
%                     Cn(isamp) = net_clus;
%                     Cr(isamp) = rand_clus;
%                     Cl(isamp) = reg_clus;
%                     Omegaraw(isamp) = (rand_path/net_path) - (net_clus/reg_clus);
%                     SWIraw(isamp) = ((net_path - reg_path)/(rand_path - reg_path)) * ((net_clus - rand_clus)/(reg_clus - rand_clus));
%                 end
%                 Deviraw = (4*(atan2(DLraw,DCraw)))/pi - 1;
%                 fprintf('[*] Adj, nodes=%i\n',length(A_input));
%                 alpha = 0.05;
%                 sDCraw = sort(DCraw);
%                 sDLraw = sort(DLraw);
%                 sSWPraw = sort(SWPraw);
%                 sDeviraw = sort(Deviraw);
%                 sGammaraw = sort(Gammaraw);
%                 sLambdaraw = sort(Lambdaraw);
%                 sSigmaraw = sort(Sigmaraw);
%                 sOmegaraw = sort(Omegaraw);
%                 sLn = sort(Ln);
%                 sLr = sort(Lr);
%                 sCn = sort(Cn);
%                 sCr = sort(Cr);
%                 sLl = sort(Ll);
%                 sCl = sort(Cl);
%                 sSWIraw = sort(SWIraw);
%                 fprintf('\tL\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(sLn),std(sLn),nsamps,sLn(round((alpha/2)*nsamps)),sLn(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tLr\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(sLr),std(sLr),nsamps,sLr(round((alpha/2)*nsamps)),sLr(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tLl\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(sLl),std(sLl),nsamps,sLl(round((alpha/2)*nsamps)),sLl(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tC\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(sCn),std(sCn),nsamps,sCn(round((alpha/2)*nsamps)),sCn(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tCr\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(sCr),std(sCr),nsamps,sCr(round((alpha/2)*nsamps)),sCr(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tCl\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(sCl),std(sCl),nsamps,sCl(round((alpha/2)*nsamps)),sCl(round((1 - (alpha/2))*nsamps)));
% 
%                 fprintf('\tDC\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(DCraw),std(DCraw),nsamps,sDCraw(round((alpha/2)*nsamps)),sDCraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tDL\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(DLraw),std(DLraw),nsamps,sDLraw(round((alpha/2)*nsamps)),sDLraw(round((1 - (alpha/2))*nsamps)));
% 
%                 fprintf('\tSWP\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(SWPraw),std(SWPraw),nsamps,sSWPraw(round((alpha/2)*nsamps)),sSWPraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tSWI\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(SWIraw),std(SWIraw),nsamps,sSWIraw(round((alpha/2)*nsamps)),sSWIraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tdelta\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(Deviraw),std(Deviraw),nsamps,sDeviraw(round((alpha/2)*nsamps)),sDeviraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tgamma\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(Gammaraw),std(Gammaraw),nsamps,sGammaraw(round((alpha/2)*nsamps)),sGammaraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tlambda\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(Lambdaraw),std(Lambdaraw),nsamps,sLambdaraw(round((alpha/2)*nsamps)),sLambdaraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tsigma\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(Sigmaraw),std(Sigmaraw),nsamps,sSigmaraw(round((alpha/2)*nsamps)),sSigmaraw(round((1 - (alpha/2))*nsamps)));
%                 fprintf('\tomega\t%.6f, std %.2d, n=%i, 95%%CI: %.3d - %.3d\n',...
%                     mean(Omegaraw),std(Omegaraw),nsamps,sOmegaraw(round((alpha/2)*nsamps)),sOmegaraw(round((1 - (alpha/2))*nsamps)));
%                 % --------------------------------------------------------------------
% 
% 
% 
% 
%                 % === Calculate SWPs ===
%                 %A_queue = {Adj_plt_MKraw,Adj_plt_mSu_red};
%                 A_queue = {Adj_plt_mSu_red};
%                 SC_q = cell(1,length(A_queue));
%                 for iq = 1:length(A_queue)
%                     A = weight_conversion(A_queue{iq},'autofix');
%                     A = weight_conversion(A,'normalize');
%                     Aswp = A;
%                     SWPs = [];
%                     Devis = [];
%                     n_pts = 120;
%                     max_density = sum(sum(triu(Aswp,1)>0))/nchoosek(length(Aswp),2);
%                     %Thrs = linspace(0,max_density,n_pts+2);
%                     Thrs = linspace(0.15,max_density,n_pts+2);
%                     Thrs = Thrs(2:(end-1));
%                     fn_cache_dk = sprintf('fig_cluster2_sworld_circle_swp_dk_%i_macaque-%i',iM,iq);
%                     if (~exist([fn_cache_dk,'.mat'],'file'))
%                         n_thr = 1;
%                         for thr = Thrs
%                             nsamps = 10000;
%                             SWP = zeros(1,nsamps);
%                             DC = zeros(1,nsamps);
%                             DL = zeros(1,nsamps);
%                             parfor isamp = 1:nsamps
%                                 %[swp,dc,dl] = small_world_propensity(Aswp);
%                                 [swp,dc,dl] = small_world_propensity(threshold_proportional(Aswp,thr));
%                                 SWP(isamp) = swp;
%                                 DC(isamp) = dc;
%                                 DL(isamp) = dl;
%                             end
%                             Devi = (4*(atan2(DL,DC)))/pi - 1;
%                             fprintf('[%4i/%4i] Macaque-%i SWP (n=%i), weighted\n',n_thr,length(Thrs),iq,length(Aswp));
%                             fprintf('\tthresh: %.3f\n',thr);
%                             fprintf('\tn_samp: %i\n',nsamps);
%                             fprintf('\tSWP: %.6f (+- %.6f)\n',mean(SWP),std(SWP));
%                             fprintf('\tdeltaC: %.6f (+- %.6f)\n',mean(DC),std(DC));
%                             fprintf('\tdeltaL: %.6f (+- %.6f)\n',mean(DL),std(DL));
%                             fprintf('\tdeviance: %.6f (+- %.6f)\n',mean(Devi),std(Devi));
%                             SWPs = [SWPs, (SWP')];
%                             Devis = [Devis, (Devi')];
%                             n_thr = n_thr + 1;
%                         end
%                         save(fn_cache_dk,'-v7.3','-nocompression','SWPs','Devis','Thrs');
%                         CaSwp = load(fn_cache_dk);
%                         SC_q{iq} = CaSwp;
%                     else
%                         fprintf('[*] Found cache: %s, loading..\n',fn_cache_dk);
%                         CaSwp = load(fn_cache_dk);
%                         SC_q{iq} = CaSwp;
%                     end
%                 end
% 
%                 % select
%                 Cmk = SC_q{1};


                % === Human ===
                % DK atlas
                % Weighted Human
                A = weight_conversion(Ahs,'autofix');




                % SWP raw
                % --------------------------------------------------------------------
                try
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
                
                
                catch e
                    rethrow(e);
                end
                % --------------------------------------------------------------------

                Xvars = cell(1,15);
                Xvars{1} = sLn;
                Xvars{2} = sLr;
                Xvars{3} = sLl;
                Xvars{4} = sCn;
                Xvars{5} = sCr;
                Xvars{6} = sCl;
                Xvars{7} = DCraw;
                Xvars{8} = DLraw;
                Xvars{9} = SWPraw;
                Xvars{10} = SWIraw;
                Xvars{11} = Deviraw;
                Xvars{12} = Gammaraw;
                Xvars{13} = Lambdaraw;
                Xvars{14} = Sigmaraw;
                Xvars{15} = Omegaraw;
                
                for ixvar = 1:length(Xvars)
                    xvar = Xvars{ixvar};
                    T(iSub,ixvar,1) = mean(xvar);
                    T(iSub,ixvar,2) = std(xvar);
                    T(iSub,ixvar,3) = nsamps;
                    T(iSub,ixvar,4) = xvar(round((alpha/2)*nsamps));
                    T(iSub,ixvar,5) = xvar(round((1 - (alpha/2))*nsamps));
                end
               % return

% 
% 
%                 Aswp = A;
%                 SWPs = [];
%                 Devis = [];
%                 n_pts = 120;
%                 max_density = sum(sum(triu(Aswp,1)>0))/nchoosek(length(Aswp),2);
%                 %Thrs = linspace(0,max_density,n_pts+2);
%                 Thrs = linspace(0.15,max_density,n_pts+2);
%                 Thrs = Thrs(2:(end-1));
%                 fn_cache_dk = sprintf('fig_cluster2_sworld_circle_swp_dk_%i',iM);
%                 if (~exist([fn_cache_dk,'.mat'],'file'))
%                     for thr = Thrs
%                         nsamps = 10000;
%                         SWP = zeros(1,nsamps);
%                         DC = zeros(1,nsamps);
%                         DL = zeros(1,nsamps);
%                         parfor isamp = 1:nsamps
%                             %[swp,dc,dl] = small_world_propensity(Aswp);
%                             [swp,dc,dl] = small_world_propensity(threshold_proportional(Ahs,thr));
%                             SWP(isamp) = swp;
%                             DC(isamp) = dc;
%                             DL(isamp) = dl;
%                         end
%                         Devi = (4*(atan2(DL,DC)))/pi - 1;
%                         fprintf('[*] Human SWP (n=%i), weighted\n',length(Aswp));
%                         fprintf('\tthresh: %.3f\n',thr);
%                         fprintf('\tn_samp: %i\n',nsamps);
%                         fprintf('\tSWP: %.6f (+- %.6f)\n',mean(SWP),std(SWP));
%                         fprintf('\tdeltaC: %.6f (+- %.6f)\n',mean(DC),std(DC));
%                         fprintf('\tdeltaL: %.6f (+- %.6f)\n',mean(DL),std(DL));
%                         fprintf('\tdeviance: %.6f (+- %.6f)\n',mean(Devi),std(Devi));
%                         SWPs = [SWPs, (SWP')];
%                         Devis = [Devis, (Devi')];
%                     end
%                     save(fn_cache_dk,'-v7.3','-nocompression','SWPs','Devis','Thrs');
%                 else
%                     fprintf('[*] Found cache: %s, loading..\n',fn_cache_dk);
%                     load(fn_cache_dk);
%                 end
% 
%                 close all;
% 
%                 col_ci = [1 1 1]*0.8;
%                 col_sw = [0.8 0.1 0.1];
%                 thresh_sw = 0.6;
% 
% 
%                 h = figure('Position',[0,0,1000,500]);
%                 set(h,'PaperUnits','inches');
%                 set(h,'PaperPosition',[0, 0, 8, 4]);
%                 subplot(2,1,1);
%                 % -- whisker plot --
%             %     boxplot(SWPs,'OutlierSize',1);
%             %     xlabtxt = {};
%             %     for ixt = 1:length(Thrs)
%             %         xlabtxt{ixt} = sprintf('%.0f',100*Thrs(ixt));
%             %     end
%             %     xticklabels(xlabtxt);
%             %     xtickangle(90);
% 
%                 % -- confidence plot --
%                 SWPs_ci_lo = zeros(1,length(Thrs));
%                 SWPs_ci_hi = zeros(1,length(Thrs));
%                 for ip = 1:length(Thrs)
%                     Svec = sort(SWPs(:,ip));
%                     SWPs_ci_lo(ip) = Svec(round(0.5*ci_alpha*length(Svec)));
%                     SWPs_ci_hi(ip) = Svec(round((1 - 0.5*ci_alpha)*length(Svec)));
%                 end
% 
%                 % repeat for monkey
%                 SWPsm_ci_lo = zeros(1,length(Thrs));
%                 SWPsm_ci_hi = zeros(1,length(Thrs));
%                 for ip = 1:length(Thrs)
%                     Svec = sort(Cmk.SWPs(:,ip));
%                     SWPsm_ci_lo(ip) = Svec(round(0.5*ci_alpha*length(Svec)));
%                     SWPsm_ci_hi(ip) = Svec(round((1 - 0.5*ci_alpha)*length(Svec)));
%                 end
% 
%                 % Threshold SW
%             %     x_sws = 100*Thrs((SWPs_ci_lo > thresh_sw));
%             %     q = fill([x_sws(1) x_sws(end) x_sws(end) x_sws(1)],[0 0 1 1],'r',...
%             %         'FaceColor',col_sw,'LineStyle','none','FaceAlpha',0.1);
% 
%                 x = 100*Thrs;
%                 y = median(SWPs);
%                 hold all;
% 
%                 % Plot Monkey SWP
%                 col_cim = plasma(10);
%                 col_cim = col_cim(3,:);
%                 fill([x';flipud(x')],[SWPsm_ci_lo';flipud(SWPsm_ci_hi')],col_ci,'linestyle','none','FaceAlpha',1);
%                 plot(x,median(Cmk.SWPs),'-','Color',[1 1 1]*0);
% 
%                 % Plot Human SWP
%                 fill([x';flipud(x')],[SWPs_ci_lo';flipud(SWPs_ci_hi')],col_cim,'linestyle','none','FaceAlpha',0.4);
%                 plot(x,y,'-','Color',0.4*col_cim);
% 
%             %     plot([min(x) max(x)],[1 1]*SWP_mk_mean,'-','Color',[0 0 1]);
%             %     plot([min(x) max(x)],[1 1]*SWP_mk_ci_lo,'--','Color',[0 0 1]);
%             %     plot([min(x) max(x)],[1 1]*SWP_mk_ci_hi,'--','Color',[0 0 1]);
%                 %plot([min(x) max(x)],[1 1]*(SWP_mk_mean + SWP_mk_std),'--','Color',[0 0 1]);
%                 %plot([min(x) max(x)],[1 1]*(SWP_mk_mean - SWP_mk_std),'--','Color',[0 0 1]);
%                 %plot(x,mean(SWPs),'-','Color',[1 1 1]*0.6);
% 
% 
%                 %plot([min(x) max(x)],[thresh_sw thresh_sw],'--','Color',[0.9 0 0]);
%                 %plot(x,SWPs_ci_lo,'-','Color',col_ci_bound,'LineWidth',1);
%                 %plot(x,SWPs_ci_hi,'-','Color',col_ci_bound,'LineWidth',1);
% 
%                 plot([x(1) x(end)],[1 1]*thresh_sw,'--','Color',[1 1 1]*0.7);
% 
%                 box off;
%                 set(gca,'TickDir','out');
%                 xlabel('Density (%)');
%                 ylabel('Small-World Propensity');
%                 ylim([0 1]);
%                 xlim([min(100*Thrs),max(100*Thrs)]);
%                 set(gca,'Layer','top');
%                 %legend({'un','deux','trois'},'Location','BestOutside');
%                 fprintf('Overall SWP: %.6f (+- %.6f)\n',nanmean(SWPs(:)),nanstd(SWPs(:)));
% 
% 
% 
%                 %figure;
%                 subplot(2,1,2);
%             %     boxplot(Devis,'OutlierSize',1);
% 
%                 % -- confidence plot --
%                 Devis_ci_lo = zeros(1,length(Thrs));
%                 Devis_ci_hi = zeros(1,length(Thrs));
%                 for ip = 1:length(Thrs)
%                     Svec = sort(Devis(:,ip));
%                     Devis_ci_lo(ip) = Svec(round(0.5*ci_alpha*length(Svec)));
%                     Devis_ci_hi(ip) = Svec(round((1 - 0.5*ci_alpha)*length(Svec)));
%                 end
% 
%                 %monkey
%                 Devism_ci_lo = zeros(1,length(Thrs));
%                 Devism_ci_hi = zeros(1,length(Thrs));
%                 for ip = 1:length(Thrs)
%                     Svec = sort(Cmk.Devis(:,ip));
%                     Devism_ci_lo(ip) = Svec(round(0.5*ci_alpha*length(Svec)));
%                     Devism_ci_hi(ip) = Svec(round((1 - 0.5*ci_alpha)*length(Svec)));
%                 end
% 
%                 x = 100*Thrs;
%                 y = median(Devis);
% 
%             %     % Threshold SW
%             %     x_sws = 100*Thrs((SWPs_ci_lo > thresh_sw));
%             %     q = fill([x_sws(1) x_sws(end) x_sws(end) x_sws(1)],[-1 -1 1 1],'r',...
%             %         'FaceColor',col_sw,'LineStyle','none','FaceAlpha',0.1);
% 
%                 hold all;
% 
%                 % Plot Monkey SWP
%                 col_cim = plasma(10);
%                 col_cim = col_cim(3,:);
%                 fill([x';flipud(x')],[Devism_ci_lo';flipud(Devism_ci_hi')],col_ci,'linestyle','none','FaceAlpha',1);
%                 plot(x,median(Cmk.Devis),'-','Color',[1 1 1]*0);
% 
%                 % Plot Human
%                 fill([x';flipud(x')],[Devis_ci_lo';flipud(Devis_ci_hi')],col_cim,'linestyle','none','FaceAlpha',0.4);
%                 plot(x,y,'-','Color',0.4*col_cim); 
%                 %plot(x,mean(Devis),'-','Color',[1 1 1]*0.6);
%                 plot([min(x) max(x)],[0 0],'--','Color',[1 1 1]*0.7);
% 
% 
%                 box off;
%                 set(gca,'TickDir','out');
%                 xlabel('Density (%)');
%                 ylabel('Contribution to Deviation');
%             %     xlabtxt = {};
%             %     for ixt = 1:length(Thrs)
%             %         xlabtxt{ixt} = sprintf('%.0f',100*Thrs(ixt));
%             %     end
%             %     xticklabels(xlabtxt);
%             %     xtickangle(90);
%                 ylim([-1 1]);
%                 xlim([min(100*Thrs),max(100*Thrs)]);
%                 set(gca,'Layer','top');
%                 fprintf('Overall SWP: %.6f (+- %.6f)\n',nanmean(Devis(:)),nanstd(Devis(:)));
% 
%                 print(h,sprintf('./figures/sworld_circle/fig_cluster2_sworld_circle_im-%i_swp',iM),'-depsc','-r900');
%                 print(h,sprintf('./figures/sworld_circle/fig_cluster2_sworld_circle_im-%i_swp',iM),'-dpng','-r900');
%                 close(h);
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
                end
            end
        end
    end
end

for ish = 1:5
    TW = squeeze(T(:,:,ish));
    writematrix(TW,'./figures/sworld_circle/output_subjects_bip.xlsx','Sheet',ish);
end



% oct. 25 2020

% fig_cluster2_sworld_circle_pCorr
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 0.000000000000 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating monkey MK-FLNe SWP with 10000 permutations..
% Starting parallel pool (parpool) using the 'local' profile ...
% Connected to the parallel pool (number of workers: 6).
% [*] Adj, nodes=91
% 	L	269.980605, std 5.20e-11, n=10000, 95%CI: 2.700e+02 - 2.700e+02
% 	Lr	71.980662, std 3.19e+00, n=10000, 95%CI: 6.631e+01 - 7.880e+01
% 	Ll	280.095229, std 4.24e+00, n=10000, 95%CI: 2.714e+02 - 2.881e+02
% 	C	0.002804, std 4.58e-16, n=10000, 95%CI: 2.804e-03 - 2.804e-03
% 	Cr	0.000755, std 2.76e-05, n=10000, 95%CI: 7.030e-04 - 8.103e-04
% 	Cl	0.002825, std 5.00e-06, n=10000, 95%CI: 2.816e-03 - 2.836e-03
% 	DC	0.010171, std 2.39e-03, n=10000, 95%CI: 5.789e-03 - 1.520e-02
% 	DL	0.951699, std 1.93e-02, n=10000, 95%CI: 9.163e-01 - 9.927e-01
% 	SWP	0.327007, std 1.36e-02, n=10000, 95%CI: 2.980e-01 - 3.520e-01
% 	SWI	0.047721, std 1.93e-02, n=10000, 95%CI: 7.197e-03 - 8.282e-02
% 	delta	0.986385, std 3.23e-03, n=10000, 95%CI: 9.797e-01 - 9.922e-01
% 	gamma	3.719839, std 1.36e-01, n=10000, 95%CI: 3.460e+00 - 3.988e+00
% 	lambda	3.757955, std 1.63e-01, n=10000, 95%CI: 3.426e+00 - 4.071e+00
% 	sigma	0.991576, std 5.36e-02, n=10000, 95%CI: 8.921e-01 - 1.101e+00
% 	omega	-0.725933, std 1.19e-02, n=10000, 95%CI: -7.470e-01 - -7.003e-01
% [*] Markov-Kennedy SWP (n=91), weighted
% 	n_samp: 100
% 	SWP: 0.328516 (+- 0.013535)(95% CI: 0.303522 - 0.347981, n=100)
% 	deltaC: 0.010032 (+- 0.002375)
% 	deltaL: 0.949566 (+- 0.019148)
% 	deviance: 0.986528 (+- 0.003255)(95% CI: 0.981133 - 0.990769, n=100)
% [*] --- Monkey MK anatomy, ECoG ---
% [*] Remove NaNs, starting nodes: 39
% 	finished nodes: 18
% [*] Calculating monkey functional interactions SWP with 10000 permutations..
% [*] Adj, nodes=18
% 	L	8.815830, std 1.28e-12, n=10000, 95%CI: 8.816e+00 - 8.816e+00
% 	Lr	8.497183, std 2.95e-02, n=10000, 95%CI: 8.440e+00 - 8.557e+00
% 	Ll	8.848119, std 7.00e-03, n=10000, 95%CI: 8.834e+00 - 8.862e+00
% 	C	0.567225, std 4.83e-14, n=10000, 95%CI: 5.672e-01 - 5.672e-01
% 	Cr	0.476979, std 1.24e-02, n=10000, 95%CI: 4.555e-01 - 5.035e-01
% 	Cl	0.492267, std 5.66e-05, n=10000, 95%CI: 4.922e-01 - 4.924e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.907697, std 1.99e-02, n=10000, 95%CI: 8.683e-01 - 9.462e-01
% 	SWP	0.358161, std 1.41e-02, n=10000, 95%CI: 3.309e-01 - 3.860e-01
% 	SWI	1.525424, std 1.24e+02, n=10000, 95%CI: -3.308e+00 - 4.031e+00
% 	delta	1.000000, std 00, n=10000, 95%CI: 001 - 001
% 	gamma	1.189995, std 3.06e-02, n=10000, 95%CI: 1.126e+00 - 1.245e+00
% 	lambda	1.037513, std 3.61e-03, n=10000, 95%CI: 1.030e+00 - 1.044e+00
% 	sigma	1.146951, std 2.85e-02, n=10000, 95%CI: 1.087e+00 - 1.199e+00
% 	omega	-0.188415, std 3.35e-03, n=10000, 95%CI: -1.949e-01 - -1.816e-01
% [*] Found cache: fig_cluster2_sworld_circle_swp_dk_1_macaque-1_pCorr, loading..
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=20
% 	L	3.271470, std 5.26e-13, n=10000, 95%CI: 3.271e+00 - 3.271e+00
% 	Lr	3.068528, std 1.27e-02, n=10000, 95%CI: 3.048e+00 - 3.098e+00
% 	Ll	2.987131, std 1.38e-03, n=10000, 95%CI: 2.984e+00 - 2.990e+00
% 	C	0.711050, std 1.54e-14, n=10000, 95%CI: 7.111e-01 - 7.111e-01
% 	Cr	0.377655, std 2.07e-02, n=10000, 95%CI: 3.397e-01 - 4.214e-01
% 	Cl	0.533165, std 4.69e-05, n=10000, 95%CI: 5.331e-01 - 5.333e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	-2.574250, std 5.32e-01, n=10000, 95%CI: -3.652e+00 - -1.563e+00
% 	SWP	-0.820270, std 3.76e-01, n=10000, 95%CI: -1.582e+00 - -1.051e-01
% 	SWI	7.705374, std 1.06e+00, n=10000, 95%CI: 5.774e+00 - 9.964e+00
% 	delta	-3.000000, std 00, n=10000, 95%CI: -003 - -003
% 	gamma	1.888414, std 1.03e-01, n=10000, 95%CI: 1.687e+00 - 2.093e+00
% 	lambda	1.066155, std 4.39e-03, n=10000, 95%CI: 1.056e+00 - 1.073e+00
% 	sigma	1.771099, std 9.35e-02, n=10000, 95%CI: 1.591e+00 - 1.958e+00
% 	omega	-0.395674, std 3.88e-03, n=10000, 95%CI: -4.018e-01 - -3.868e-01
% [*] Found cache: fig_cluster2_sworld_circle_swp_dk_1_pCorr, loading..
% Overall SWP: 0.621409 (+- 0.290804)
% Overall SWP: -0.126046 (+- 0.371051)
% figure_t20_lognormal_pCorr
% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T20_lognormal’: File exists
% [*] Log(C) kstest p=3.7032e-02, n=301
% 	min:-0.45376665, max:-0.08600742
% [*] Linear kstest p=2.4947e-03, n=301
% 	min:0.35174939, max:0.82033752
% [fusiform] n_bip_pairs: 1916
% [inferiorparietal] n_bip_pairs: 1005
% [inferiortemporal] n_bip_pairs: 3460
% [lateraloccipital] n_bip_pairs: 1041
% [lateralorbitofrontal] n_bip_pairs: 685
% [lingual] n_bip_pairs: 805
% [middletemporal] n_bip_pairs: 3775
% [parsopercularis] n_bip_pairs: 395
% [parsorbitalis] n_bip_pairs: 461
% [parstriangularis] n_bip_pairs: 350
% [postcentral] n_bip_pairs: 623
% [precentral] n_bip_pairs: 613
% [precuneus] n_bip_pairs: 287
% [rostralmiddlefrontal] n_bip_pairs: 1346
% [superiorfrontal] n_bip_pairs: 696
% [superiorparietal] n_bip_pairs: 306
% [superiortemporal] n_bip_pairs: 2161
% [supramarginal] n_bip_pairs: 1209
% [*] Log(C) kstest p=8.4181e-23, n=2387
% 	min:-0.87155125, max:-0.10260642
% [*] Linear kstest p=9.6007e-47, n=2387
% 	min:0.13441531, max:0.78957534
% figure_t20_lognormal_pCorr
% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T20_lognormal’: File exists
% [*] Log(C) kstest p=3.7032e-02, n=301
% 	min:-0.45376665, max:-0.08600742
% [*] Linear kstest p=2.4947e-03, n=301
% 	min:0.35174939, max:0.82033752
% [fusiform] n_bip_pairs: 1916
% [inferiorparietal] n_bip_pairs: 1005
% [inferiortemporal] n_bip_pairs: 3460
% [lateraloccipital] n_bip_pairs: 1041
% [lateralorbitofrontal] n_bip_pairs: 685
% [lingual] n_bip_pairs: 805
% [middletemporal] n_bip_pairs: 3775
% [parsopercularis] n_bip_pairs: 395
% [parsorbitalis] n_bip_pairs: 461
% [parstriangularis] n_bip_pairs: 350
% [postcentral] n_bip_pairs: 623
% [precentral] n_bip_pairs: 613
% [precuneus] n_bip_pairs: 287
% [rostralmiddlefrontal] n_bip_pairs: 1346
% [superiorfrontal] n_bip_pairs: 696
% [superiorparietal] n_bip_pairs: 306
% [superiortemporal] n_bip_pairs: 2161
% [supramarginal] n_bip_pairs: 1209
% [*] Log(C) kstest p=8.4181e-23, n=2387
% 	min:-0.87155125, max:-0.10260642
% [*] Linear kstest p=9.6007e-47, n=2387
% 	min:0.13441531, max:0.78957534
% IdleTimeout has been reached.
% Parallel pool using the 'local' profile is shutting down.
% fig_cluster2_sworld_circle_subjects_bip_pCorr
% [*] Subject 1
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 28.888322800398 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% Starting parallel pool (parpool) using the 'local' profile ...
% Connected to the parallel pool (number of workers: 6).
% [*] Adj, nodes=62
% 	L	1.811275, std 1.71e-14, n=10000, 95%CI: 1.811e+00 - 1.811e+00
% 	Lr	1.811275, std 3.23e-06, n=10000, 95%CI: 1.811e+00 - 1.811e+00
% 	Ll	1.830264, std 1.45e-04, n=10000, 95%CI: 1.830e+00 - 1.831e+00
% 	C	0.574483, std 7.99e-15, n=10000, 95%CI: 5.745e-01 - 5.745e-01
% 	Cr	0.573658, std 7.04e-06, n=10000, 95%CI: 5.736e-01 - 5.737e-01
% 	Cl	0.565745, std 1.24e-07, n=10000, 95%CI: 5.657e-01 - 5.657e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000014, std 1.70e-04, n=10000, 95%CI: 000 - 000
% 	SWP	0.999990, std 1.20e-04, n=10000, 95%CI: 001 - 001
% 	SWI	-0.104226, std 9.82e-04, n=10000, 95%CI: -1.060e-01 - -1.022e-01
% 	delta	-0.952200, std 3.05e-01, n=10000, 95%CI: -001 - -001
% 	gamma	1.001438, std 1.23e-05, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000000, std 1.78e-06, n=10000, 95%CI: 1.000e+00 - 001
% 	sigma	1.001438, std 1.24e-05, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.015446, std 1.79e-06, n=10000, 95%CI: -1.545e-02 - -1.545e-02
% [*] Subject 2
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 27.200641959906 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=62
% 	L	1.825695, std 1.26e-13, n=10000, 95%CI: 1.826e+00 - 1.826e+00
% 	Lr	1.825671, std 3.77e-05, n=10000, 95%CI: 1.826e+00 - 1.826e+00
% 	Ll	1.842674, std 8.03e-05, n=10000, 95%CI: 1.843e+00 - 1.843e+00
% 	C	0.592076, std 1.67e-15, n=10000, 95%CI: 5.921e-01 - 5.921e-01
% 	Cr	0.590446, std 1.20e-05, n=10000, 95%CI: 5.904e-01 - 5.905e-01
% 	Cl	0.583084, std 9.01e-08, n=10000, 95%CI: 5.831e-01 - 5.831e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.001394, std 2.20e-03, n=10000, 95%CI: 000 - 7.607e-03
% 	SWP	0.999015, std 1.56e-03, n=10000, 95%CI: 9.946e-01 - 001
% 	SWI	-0.221065, std 2.05e-03, n=10000, 95%CI: -2.248e-01 - -2.167e-01
% 	delta	0.265400, std 9.64e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.002760, std 2.04e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	lambda	1.000013, std 2.06e-05, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.002747, std 2.92e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	omega	-0.015435, std 2.06e-05, n=10000, 95%CI: -1.549e-02 - -1.542e-02
% [*] Subject 3
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 32.884566187859 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=90
% 	L	2.353063, std 2.98e-13, n=10000, 95%CI: 2.353e+00 - 2.353e+00
% 	Lr	2.353063, std 4.78e-06, n=10000, 95%CI: 2.353e+00 - 2.353e+00
% 	Ll	2.368706, std 9.92e-05, n=10000, 95%CI: 2.369e+00 - 2.369e+00
% 	C	0.554090, std 9.99e-14, n=10000, 95%CI: 5.541e-01 - 5.541e-01
% 	Cr	0.553384, std 3.11e-06, n=10000, 95%CI: 5.534e-01 - 5.534e-01
% 	Cl	0.548191, std 5.32e-08, n=10000, 95%CI: 5.482e-01 - 5.482e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000048, std 3.04e-04, n=10000, 95%CI: 000 - 6.519e-04
% 	SWP	0.999966, std 2.15e-04, n=10000, 95%CI: 9.995e-01 - 001
% 	SWI	-0.135872, std 6.81e-04, n=10000, 95%CI: -1.371e-01 - -1.344e-01
% 	delta	0.787800, std 6.16e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.001275, std 5.62e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000000, std 2.03e-06, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001275, std 5.98e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.010761, std 2.03e-06, n=10000, 95%CI: -1.077e-02 - -1.076e-02
% [*] Subject 4
% p - human fraction of ROIs significant: 0.9995
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 31.077582657337 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=62
% 	L	1.521434, std 2.50e-13, n=10000, 95%CI: 1.521e+00 - 1.521e+00
% 	Lr	1.521604, std 5.14e-05, n=10000, 95%CI: 1.521e+00 - 1.522e+00
% 	Ll	1.535422, std 1.05e-05, n=10000, 95%CI: 1.535e+00 - 1.535e+00
% 	C	0.693615, std 1.16e-13, n=10000, 95%CI: 6.936e-01 - 6.936e-01
% 	Cr	0.691828, std 1.19e-05, n=10000, 95%CI: 6.918e-01 - 6.919e-01
% 	Cl	0.684111, std 5.43e-08, n=10000, 95%CI: 6.841e-01 - 6.841e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000001, std 2.98e-05, n=10000, 95%CI: 000 - 000
% 	SWP	0.999999, std 2.10e-05, n=10000, 95%CI: 001 - 001
% 	SWI	-0.234363, std 2.05e-03, n=10000, 95%CI: -2.380e-01 - -2.301e-01
% 	delta	-0.997800, std 6.63e-02, n=10000, 95%CI: -001 - -001
% 	gamma	1.002583, std 1.72e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	lambda	0.999888, std 3.38e-05, n=10000, 95%CI: 9.998e-01 - 1.000e+00
% 	sigma	1.002694, std 3.69e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	omega	-0.013781, std 3.38e-05, n=10000, 95%CI: -1.385e-02 - -1.372e-02
% [*] Subject 5
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 22.406616508961 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=62
% 	L	2.216751, std 4.91e-13, n=10000, 95%CI: 2.217e+00 - 2.217e+00
% 	Lr	2.216495, std 1.52e-04, n=10000, 95%CI: 2.216e+00 - 2.217e+00
% 	Ll	2.235861, std 1.15e-04, n=10000, 95%CI: 2.236e+00 - 2.236e+00
% 	C	0.527666, std 5.38e-14, n=10000, 95%CI: 5.277e-01 - 5.277e-01
% 	Cr	0.526189, std 1.29e-05, n=10000, 95%CI: 5.262e-01 - 5.262e-01
% 	Cl	0.519575, std 1.32e-07, n=10000, 95%CI: 5.196e-01 - 5.196e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.013150, std 7.70e-03, n=10000, 95%CI: 1.252e-03 - 3.060e-02
% 	SWP	0.990701, std 5.44e-03, n=10000, 95%CI: 9.783e-01 - 9.991e-01
% 	SWI	-0.220463, std 2.92e-03, n=10000, 95%CI: -2.257e-01 - -2.144e-01
% 	delta	0.986000, std 1.67e-01, n=10000, 95%CI: 001 - 001
% 	gamma	1.002808, std 2.45e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	lambda	1.000115, std 6.85e-05, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.002692, std 7.32e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	omega	-0.015688, std 6.85e-05, n=10000, 95%CI: -1.584e-02 - -1.558e-02
% [*] Subject 6
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 25.002464473248 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=76
% 	L	2.531030, std 3.02e-13, n=10000, 95%CI: 2.531e+00 - 2.531e+00
% 	Lr	2.531027, std 1.56e-05, n=10000, 95%CI: 2.531e+00 - 2.531e+00
% 	Ll	2.551832, std 1.67e-04, n=10000, 95%CI: 2.552e+00 - 2.552e+00
% 	C	0.489627, std 8.37e-14, n=10000, 95%CI: 4.896e-01 - 4.896e-01
% 	Cr	0.489088, std 3.40e-06, n=10000, 95%CI: 4.891e-01 - 4.891e-01
% 	Cl	0.483820, std 9.14e-08, n=10000, 95%CI: 4.838e-01 - 4.838e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000170, std 7.46e-04, n=10000, 95%CI: 000 - 2.449e-03
% 	SWP	0.999880, std 5.28e-04, n=10000, 95%CI: 9.983e-01 - 001
% 	SWI	-0.102277, std 7.17e-04, n=10000, 95%CI: -1.035e-01 - -1.008e-01
% 	delta	0.927600, std 3.74e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.001102, std 6.96e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000001, std 6.17e-06, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001100, std 9.39e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.012002, std 6.17e-06, n=10000, 95%CI: -1.202e-02 - -1.200e-02
% [*] Subject 7
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 28.555323004723 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=76
% 	L	2.293459, std 2.10e-13, n=10000, 95%CI: 2.293e+00 - 2.293e+00
% 	Lr	2.293459, std 3.97e-07, n=10000, 95%CI: 2.293e+00 - 2.293e+00
% 	Ll	2.312543, std 1.11e-04, n=10000, 95%CI: 2.312e+00 - 2.313e+00
% 	C	0.605426, std 7.88e-14, n=10000, 95%CI: 6.054e-01 - 6.054e-01
% 	Cr	0.604979, std 3.69e-06, n=10000, 95%CI: 6.050e-01 - 6.050e-01
% 	Cl	0.598185, std 7.28e-08, n=10000, 95%CI: 5.982e-01 - 5.982e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000001, std 2.08e-05, n=10000, 95%CI: 000 - 5.339e-13
% 	SWP	1.000000, std 1.47e-05, n=10000, 95%CI: 1.000e+00 - 001
% 	SWI	-0.065909, std 5.79e-04, n=10000, 95%CI: -6.697e-02 - -6.467e-02
% 	delta	0.347800, std 9.38e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.000740, std 6.10e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000000, std 1.73e-07, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.000740, std 6.11e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.012105, std 2.13e-07, n=10000, 95%CI: -1.211e-02 - -1.210e-02
% [*] Subject 8
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 23.407025486231 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=62
% 	L	2.152632, std 2.50e-13, n=10000, 95%CI: 2.153e+00 - 2.153e+00
% 	Lr	2.152581, std 5.96e-05, n=10000, 95%CI: 2.152e+00 - 2.153e+00
% 	Ll	2.171893, std 1.14e-04, n=10000, 95%CI: 2.172e+00 - 2.172e+00
% 	C	0.579309, std 8.42e-14, n=10000, 95%CI: 5.793e-01 - 5.793e-01
% 	Cr	0.578172, std 1.17e-05, n=10000, 95%CI: 5.782e-01 - 5.782e-01
% 	Cl	0.570758, std 1.43e-07, n=10000, 95%CI: 5.708e-01 - 5.708e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.002624, std 3.06e-03, n=10000, 95%CI: 1.154e-13 - 1.042e-02
% 	SWP	0.998145, std 2.17e-03, n=10000, 95%CI: 9.926e-01 - 1.000e+00
% 	SWI	-0.152937, std 1.88e-03, n=10000, 95%CI: -1.563e-01 - -1.489e-01
% 	delta	0.990200, std 1.40e-01, n=10000, 95%CI: 001 - 001
% 	gamma	1.001966, std 2.03e-05, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	lambda	1.000024, std 2.77e-05, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001943, std 3.45e-05, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	omega	-0.015006, std 2.77e-05, n=10000, 95%CI: -1.508e-02 - -1.498e-02
% [*] Subject 9
% p - human fraction of ROIs significant: 0.9408
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 6.021188318729 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=41
% 	L	1.872644, std 4.37e-13, n=10000, 95%CI: 1.873e+00 - 1.873e+00
% 	Lr	1.871911, std 2.07e-03, n=10000, 95%CI: 1.868e+00 - 1.876e+00
% 	Ll	1.939008, std 1.62e-04, n=10000, 95%CI: 1.939e+00 - 1.939e+00
% 	C	0.542745, std 4.91e-14, n=10000, 95%CI: 5.427e-01 - 5.427e-01
% 	Cr	0.536910, std 3.53e-04, n=10000, 95%CI: 5.363e-01 - 5.377e-01
% 	Cl	0.536206, std 3.61e-07, n=10000, 95%CI: 5.362e-01 - 5.362e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.017914, std 2.05e-02, n=10000, 95%CI: 000 - 6.803e-02
% 	SWP	0.987333, std 1.45e-02, n=10000, 95%CI: 9.519e-01 - 001
% 	SWI	4.070234, std 1.05e+03, n=10000, 95%CI: -4.804e+01 - -3.098e+00
% 	delta	0.276800, std 9.61e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.010868, std 6.64e-04, n=10000, 95%CI: 1.009e+00 - 1.012e+00
% 	lambda	1.000393, std 1.10e-03, n=10000, 95%CI: 9.982e-01 - 1.003e+00
% 	sigma	1.010472, std 1.21e-03, n=10000, 95%CI: 1.008e+00 - 1.013e+00
% 	omega	-0.012585, std 1.10e-03, n=10000, 95%CI: -1.478e-02 - -1.042e-02
% [*] Subject 10
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 25.662537068129 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=60
% 	L	2.012250, std 2.38e-13, n=10000, 95%CI: 2.012e+00 - 2.012e+00
% 	Lr	2.012250, std 2.34e-13, n=10000, 95%CI: 2.012e+00 - 2.012e+00
% 	Ll	2.036944, std 9.48e-05, n=10000, 95%CI: 2.037e+00 - 2.037e+00
% 	C	0.654189, std 1.27e-13, n=10000, 95%CI: 6.542e-01 - 6.542e-01
% 	Cr	0.653719, std 4.76e-06, n=10000, 95%CI: 6.537e-01 - 6.537e-01
% 	Cl	0.644087, std 5.63e-08, n=10000, 95%CI: 6.441e-01 - 6.441e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000000, std 9.99e-14, n=10000, 95%CI: 000 - 3.598e-13
% 	SWP	1.000000, std 1.07e-13, n=10000, 95%CI: 1.000e+00 - 001
% 	SWI	-0.048813, std 5.18e-04, n=10000, 95%CI: -4.974e-02 - -4.770e-02
% 	delta	0.816800, std 5.77e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.000719, std 7.29e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000000, std 2.09e-15, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.000719, std 7.29e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.015684, std 8.88e-08, n=10000, 95%CI: -1.568e-02 - -1.568e-02
% [*] Subject 11
% p - human fraction of ROIs significant: 0.9995
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 29.111630439758 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=62
% 	L	1.641125, std 9.99e-14, n=10000, 95%CI: 1.641e+00 - 1.641e+00
% 	Lr	1.641338, std 7.66e-05, n=10000, 95%CI: 1.641e+00 - 1.641e+00
% 	Ll	1.656749, std 5.62e-05, n=10000, 95%CI: 1.657e+00 - 1.657e+00
% 	C	0.625690, std 6.63e-14, n=10000, 95%CI: 6.257e-01 - 6.257e-01
% 	Cr	0.623573, std 1.57e-05, n=10000, 95%CI: 6.235e-01 - 6.236e-01
% 	Cl	0.616183, std 6.70e-08, n=10000, 95%CI: 6.162e-01 - 6.162e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000006, std 1.08e-04, n=10000, 95%CI: 000 - 000
% 	SWP	0.999996, std 7.66e-05, n=10000, 95%CI: 001 - 001
% 	SWI	-0.290508, std 3.02e-03, n=10000, 95%CI: -2.961e-01 - -2.842e-01
% 	delta	-0.990800, std 1.35e-01, n=10000, 95%CI: -001 - -001
% 	gamma	1.003396, std 2.52e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	lambda	0.999870, std 4.66e-05, n=10000, 95%CI: 9.998e-01 - 1.000e+00
% 	sigma	1.003526, std 5.15e-05, n=10000, 95%CI: 1.003e+00 - 1.004e+00
% 	omega	-0.015300, std 4.66e-05, n=10000, 95%CI: -1.540e-02 - -1.521e-02
% [*] Subject 12
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 18.112748146057 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=55
% 	L	2.511388, std 5.96e-13, n=10000, 95%CI: 2.511e+00 - 2.511e+00
% 	Lr	2.511388, std 3.54e-06, n=10000, 95%CI: 2.511e+00 - 2.511e+00
% 	Ll	2.511388, std 5.97e-13, n=10000, 95%CI: 2.511e+00 - 2.511e+00
% 	C	0.595149, std 1.30e-13, n=10000, 95%CI: 5.951e-01 - 5.951e-01
% 	Cr	0.594056, std 1.09e-05, n=10000, 95%CI: 5.940e-01 - 5.941e-01
% 	Cl	0.594023, std 1.19e-07, n=10000, 95%CI: 5.940e-01 - 5.940e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	NaN, std NaN, n=10000, 95%CI: -2.200e+00 - NaN
% 	SWP	NaN, std NaN, n=10000, 95%CI: -5.556e-01 - NaN
% 	SWI	NaN, std NaN, n=10000, 95%CI: -3.418e+03 - Inf
% 	delta	NaN, std NaN, n=10000, 95%CI: -003 - NaN
% 	gamma	1.001840, std 1.83e-05, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	lambda	1.000000, std 1.41e-06, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001840, std 1.84e-05, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	omega	-0.001896, std 1.43e-06, n=10000, 95%CI: -1.896e-03 - -1.895e-03
% [*] Subject 13
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 20.963710784912 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=62
% 	L	2.503733, std 1.71e-13, n=10000, 95%CI: 2.504e+00 - 2.504e+00
% 	Lr	2.503713, std 4.45e-05, n=10000, 95%CI: 2.504e+00 - 2.504e+00
% 	Ll	2.526736, std 2.57e-04, n=10000, 95%CI: 2.526e+00 - 2.527e+00
% 	C	0.497300, std 4.65e-14, n=10000, 95%CI: 4.973e-01 - 4.973e-01
% 	Cr	0.496684, std 6.01e-06, n=10000, 95%CI: 4.967e-01 - 4.967e-01
% 	Cl	0.489872, std 2.21e-07, n=10000, 95%CI: 4.899e-01 - 4.899e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000857, std 1.92e-03, n=10000, 95%CI: 000 - 6.862e-03
% 	SWP	0.999394, std 1.36e-03, n=10000, 95%CI: 9.951e-01 - 001
% 	SWI	-0.090301, std 9.79e-04, n=10000, 95%CI: -9.203e-02 - -8.819e-02
% 	delta	0.861400, std 5.08e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.001240, std 1.21e-05, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000008, std 1.78e-05, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001232, std 2.16e-05, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.015170, std 1.78e-05, n=10000, 95%CI: -1.523e-02 - -1.516e-02
% [*] Subject 14
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 11.145464211702 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=34
% 	L	2.491853, std 9.33e-14, n=10000, 95%CI: 2.492e+00 - 2.492e+00
% 	Lr	2.491310, std 5.25e-04, n=10000, 95%CI: 2.490e+00 - 2.492e+00
% 	Ll	2.530719, std 9.37e-04, n=10000, 95%CI: 2.529e+00 - 2.533e+00
% 	C	0.468480, std 1.04e-13, n=10000, 95%CI: 4.685e-01 - 4.685e-01
% 	Cr	0.467473, std 2.43e-05, n=10000, 95%CI: 4.674e-01 - 4.675e-01
% 	Cl	0.455540, std 2.04e-06, n=10000, 95%CI: 4.555e-01 - 4.555e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.013610, std 1.30e-02, n=10000, 95%CI: 4.700e-14 - 4.499e-02
% 	SWP	0.990376, std 9.16e-03, n=10000, 95%CI: 9.681e-01 - 1.000e+00
% 	SWI	-0.083246, std 2.51e-03, n=10000, 95%CI: -8.769e-02 - -7.797e-02
% 	delta	0.988400, std 1.52e-01, n=10000, 95%CI: 001 - 001
% 	gamma	1.002154, std 5.21e-05, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	lambda	1.000218, std 2.11e-04, n=10000, 95%CI: 1.000e+00 - 1.001e+00
% 	sigma	1.001936, std 2.21e-04, n=10000, 95%CI: 1.001e+00 - 1.002e+00
% 	omega	-0.028624, std 2.11e-04, n=10000, 95%CI: -2.914e-02 - -2.840e-02
% [*] Subject 15
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 22.535927355289 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=62
% 	L	2.342026, std 2.56e-13, n=10000, 95%CI: 2.342e+00 - 2.342e+00
% 	Lr	2.342024, std 1.12e-05, n=10000, 95%CI: 2.342e+00 - 2.342e+00
% 	Ll	2.364824, std 1.88e-04, n=10000, 95%CI: 2.364e+00 - 2.365e+00
% 	C	0.543128, std 9.28e-14, n=10000, 95%CI: 5.431e-01 - 5.431e-01
% 	Cr	0.542402, std 6.01e-06, n=10000, 95%CI: 5.424e-01 - 5.424e-01
% 	Cl	0.534894, std 1.76e-07, n=10000, 95%CI: 5.349e-01 - 5.349e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000088, std 4.89e-04, n=10000, 95%CI: 000 - 1.346e-03
% 	SWP	0.999937, std 3.46e-04, n=10000, 95%CI: 9.990e-01 - 001
% 	SWI	-0.096708, std 8.79e-04, n=10000, 95%CI: -9.831e-02 - -9.485e-02
% 	delta	-0.131000, std 9.91e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.001339, std 1.11e-05, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000001, std 4.78e-06, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001338, std 1.21e-05, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.015393, std 4.79e-06, n=10000, 95%CI: -1.541e-02 - -1.539e-02
% [*] Subject 16
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 35.054084271193 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=111
% 	L	2.234816, std 2.27e-13, n=10000, 95%CI: 2.235e+00 - 2.235e+00
% 	Lr	2.223219, std 5.38e-04, n=10000, 95%CI: 2.222e+00 - 2.224e+00
% 	Ll	2.234816, std 2.30e-13, n=10000, 95%CI: 2.235e+00 - 2.235e+00
% 	C	0.479925, std 1.59e-14, n=10000, 95%CI: 4.799e-01 - 4.799e-01
% 	Cr	0.478118, std 7.42e-06, n=10000, 95%CI: 4.781e-01 - 4.781e-01
% 	Cl	0.478139, std 2.86e-08, n=10000, 95%CI: 4.781e-01 - 4.781e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	1.000000, std 1.84e-13, n=10000, 95%CI: 1.000e+00 - 001
% 	SWP	0.292893, std 1.29e-13, n=10000, 95%CI: 2.929e-01 - 2.929e-01
% 	SWI	-0.000000, std 1.20e-09, n=10000, 95%CI: -1.747e-10 - 7.548e-11
% 	delta	1.000000, std 00, n=10000, 95%CI: 001 - 001
% 	gamma	1.003779, std 1.56e-05, n=10000, 95%CI: 1.004e+00 - 1.004e+00
% 	lambda	1.005216, std 2.43e-04, n=10000, 95%CI: 1.005e+00 - 1.006e+00
% 	sigma	0.998571, std 2.41e-04, n=10000, 95%CI: 9.981e-01 - 9.990e-01
% 	omega	-0.008925, std 2.41e-04, n=10000, 95%CI: -9.401e-03 - -8.468e-03
% [*] Subject 17
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 35.648530185223 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=90
% 	L	2.130695, std 4.38e-13, n=10000, 95%CI: 2.131e+00 - 2.131e+00
% 	Lr	2.130695, std 2.26e-06, n=10000, 95%CI: 2.131e+00 - 2.131e+00
% 	Ll	2.145124, std 6.71e-05, n=10000, 95%CI: 2.145e+00 - 2.145e+00
% 	C	0.591704, std 8.29e-14, n=10000, 95%CI: 5.917e-01 - 5.917e-01
% 	Cr	0.590644, std 4.27e-06, n=10000, 95%CI: 5.906e-01 - 5.907e-01
% 	Cl	0.585248, std 3.72e-08, n=10000, 95%CI: 5.852e-01 - 5.852e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000018, std 1.56e-04, n=10000, 95%CI: 000 - 9.286e-13
% 	SWP	0.999988, std 1.11e-04, n=10000, 95%CI: 1.000e+00 - 001
% 	SWI	-0.196442, std 9.46e-04, n=10000, 95%CI: -1.982e-01 - -1.945e-01
% 	delta	0.375200, std 9.27e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.001795, std 7.24e-06, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	lambda	1.000000, std 1.06e-06, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001795, std 7.32e-06, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	omega	-0.011032, std 1.06e-06, n=10000, 95%CI: -1.103e-02 - -1.103e-02
% [*] Subject 18
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 30.455434381962 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=90
% 	L	2.110978, std 2.51e-13, n=10000, 95%CI: 2.111e+00 - 2.111e+00
% 	Lr	2.097635, std 6.52e-04, n=10000, 95%CI: 2.096e+00 - 2.099e+00
% 	Ll	2.120500, std 4.46e-05, n=10000, 95%CI: 2.120e+00 - 2.121e+00
% 	C	0.516595, std 1.44e-14, n=10000, 95%CI: 5.166e-01 - 5.166e-01
% 	Cr	0.513040, std 1.27e-05, n=10000, 95%CI: 5.130e-01 - 5.131e-01
% 	Cl	0.509683, std 3.81e-08, n=10000, 95%CI: 5.097e-01 - 5.097e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.583190, std 1.20e-02, n=10000, 95%CI: 5.584e-01 - 6.059e-01
% 	SWP	0.587622, std 8.46e-03, n=10000, 95%CI: 5.716e-01 - 6.051e-01
% 	SWI	-0.441312, std 1.28e-02, n=10000, 95%CI: -4.676e-01 - -4.173e-01
% 	delta	1.000000, std 00, n=10000, 95%CI: 001 - 001
% 	gamma	1.006929, std 2.49e-05, n=10000, 95%CI: 1.007e+00 - 1.007e+00
% 	lambda	1.006361, std 3.13e-04, n=10000, 95%CI: 1.006e+00 - 1.007e+00
% 	sigma	1.000565, std 3.10e-04, n=10000, 95%CI: 1.000e+00 - 1.001e+00
% 	omega	-0.019882, std 3.09e-04, n=10000, 95%CI: -2.050e-02 - -1.928e-02
% [*] Subject 19
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 41.381620347500 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=104
% 	L	2.100788, std 3.71e-13, n=10000, 95%CI: 2.101e+00 - 2.101e+00
% 	Lr	2.100788, std 1.82e-06, n=10000, 95%CI: 2.101e+00 - 2.101e+00
% 	Ll	2.113209, std 5.58e-05, n=10000, 95%CI: 2.113e+00 - 2.113e+00
% 	C	0.570248, std 8.44e-15, n=10000, 95%CI: 5.702e-01 - 5.702e-01
% 	Cr	0.569196, std 3.27e-06, n=10000, 95%CI: 5.692e-01 - 5.692e-01
% 	Cl	0.564744, std 2.41e-08, n=10000, 95%CI: 5.647e-01 - 5.647e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000019, std 1.46e-04, n=10000, 95%CI: 000 - 7.774e-05
% 	SWP	0.999987, std 1.04e-04, n=10000, 95%CI: 9.999e-01 - 001
% 	SWI	-0.236144, std 9.09e-04, n=10000, 95%CI: -2.378e-01 - -2.342e-01
% 	delta	0.508600, std 8.61e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.001847, std 5.76e-06, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	lambda	1.000000, std 8.66e-07, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001847, std 5.84e-06, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	omega	-0.009747, std 8.67e-07, n=10000, 95%CI: -9.747e-03 - -9.747e-03
% [*] Subject 20
% p - human fraction of ROIs significant: 0.9938
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 25.622929930687 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=69
% 	L	2.191242, std 9.68e-14, n=10000, 95%CI: 2.191e+00 - 2.191e+00
% 	Lr	2.191520, std 1.44e-04, n=10000, 95%CI: 2.191e+00 - 2.192e+00
% 	Ll	2.191941, std 4.63e-05, n=10000, 95%CI: 2.192e+00 - 2.192e+00
% 	C	0.500911, std 4.49e-14, n=10000, 95%CI: 5.009e-01 - 5.009e-01
% 	Cr	0.499391, std 9.52e-06, n=10000, 95%CI: 4.994e-01 - 4.994e-01
% 	Cl	0.499369, std 2.08e-07, n=10000, 95%CI: 4.994e-01 - 4.994e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.003768, std 2.41e-02, n=10000, 95%CI: 000 - 5.163e-02
% 	SWP	0.997336, std 1.70e-02, n=10000, 95%CI: 9.635e-01 - 001
% 	SWI	-162.639683, std 7.02e+02, n=10000, 95%CI: -5.551e+02 - -4.644e+01
% 	delta	-0.917800, std 3.97e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.003042, std 1.91e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	lambda	0.999873, std 6.59e-05, n=10000, 95%CI: 9.998e-01 - 1.000e+00
% 	sigma	1.003169, std 6.81e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	omega	-0.002961, std 6.59e-05, n=10000, 95%CI: -3.106e-03 - -2.850e-03
% [*] Subject 21
% p - human fraction of ROIs significant: 0.9998
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 38.899264156818 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=90
% 	L	1.917937, std 1.51e-13, n=10000, 95%CI: 1.918e+00 - 1.918e+00
% 	Lr	1.918165, std 4.52e-05, n=10000, 95%CI: 1.918e+00 - 1.918e+00
% 	Ll	1.930608, std 8.39e-05, n=10000, 95%CI: 1.930e+00 - 1.931e+00
% 	C	0.530424, std 1.84e-14, n=10000, 95%CI: 5.304e-01 - 5.304e-01
% 	Cr	0.529381, std 4.38e-06, n=10000, 95%CI: 5.294e-01 - 5.294e-01
% 	Cl	0.524649, std 5.68e-08, n=10000, 95%CI: 5.246e-01 - 5.246e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000000, std 3.20e-05, n=10000, 95%CI: 000 - 000
% 	SWP	1.000000, std 2.26e-05, n=10000, 95%CI: 001 - 001
% 	SWI	-0.224519, std 1.36e-03, n=10000, 95%CI: -2.270e-01 - -2.217e-01
% 	delta	-0.999600, std 2.83e-02, n=10000, 95%CI: -001 - -001
% 	gamma	1.001971, std 8.29e-06, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	lambda	0.999881, std 2.35e-05, n=10000, 95%CI: 9.998e-01 - 9.999e-01
% 	sigma	1.002090, std 2.44e-05, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	omega	-0.010889, std 2.35e-05, n=10000, 95%CI: -1.094e-02 - -1.085e-02
% [*] Subject 22
% p - human fraction of ROIs significant: 0.9583
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 7.975282847881 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=48
% 	L	2.051075, std 2.78e-13, n=10000, 95%CI: 2.051e+00 - 2.051e+00
% 	Lr	2.052319, std 1.41e-03, n=10000, 95%CI: 2.050e+00 - 2.055e+00
% 	Ll	2.073286, std 3.05e-04, n=10000, 95%CI: 2.073e+00 - 2.074e+00
% 	C	0.514363, std 5.03e-14, n=10000, 95%CI: 5.144e-01 - 5.144e-01
% 	Cr	0.489293, std 1.97e-04, n=10000, 95%CI: 4.890e-01 - 4.897e-01
% 	Cl	0.488638, std 1.43e-06, n=10000, 95%CI: 4.886e-01 - 4.886e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.006456, std 1.77e-02, n=10000, 95%CI: 000 - 6.570e-02
% 	SWP	0.995435, std 1.26e-02, n=10000, 95%CI: 9.535e-01 - 001
% 	SWI	-45.118351, std 1.67e+01, n=10000, 95%CI: -8.600e+01 - -2.421e+01
% 	delta	-0.619800, std 7.85e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.051239, std 4.24e-04, n=10000, 95%CI: 1.050e+00 - 1.052e+00
% 	lambda	0.999394, std 6.88e-04, n=10000, 95%CI: 9.981e-01 - 1.001e+00
% 	sigma	1.051877, std 7.85e-04, n=10000, 95%CI: 1.050e+00 - 1.053e+00
% 	omega	-0.052040, std 6.88e-04, n=10000, 95%CI: -5.341e-02 - -5.073e-02
% [*] Subject 23
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 29.748655527830 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=76
% 	L	2.175924, std 1.82e-14, n=10000, 95%CI: 2.176e+00 - 2.176e+00
% 	Lr	2.175923, std 9.07e-06, n=10000, 95%CI: 2.176e+00 - 2.176e+00
% 	Ll	2.192968, std 1.25e-04, n=10000, 95%CI: 2.193e+00 - 2.193e+00
% 	C	0.547233, std 4.33e-14, n=10000, 95%CI: 5.472e-01 - 5.472e-01
% 	Cr	0.546546, std 4.02e-06, n=10000, 95%CI: 5.465e-01 - 5.466e-01
% 	Cl	0.540479, std 9.78e-08, n=10000, 95%CI: 5.405e-01 - 5.405e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000107, std 5.30e-04, n=10000, 95%CI: 000 - 1.584e-03
% 	SWP	0.999924, std 3.75e-04, n=10000, 95%CI: 9.989e-01 - 001
% 	SWI	-0.113281, std 7.41e-04, n=10000, 95%CI: -1.146e-01 - -1.117e-01
% 	delta	0.418800, std 9.08e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.001258, std 7.37e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000001, std 4.17e-06, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001257, std 8.52e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.012497, std 4.17e-06, n=10000, 95%CI: -1.251e-02 - -1.250e-02
% [*] Subject 24
% p - human fraction of ROIs significant: 0.9985
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 28.657501518726 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=90
% 	L	2.387372, std 2.98e-13, n=10000, 95%CI: 2.387e+00 - 2.387e+00
% 	Lr	2.389685, std 4.67e-04, n=10000, 95%CI: 2.389e+00 - 2.391e+00
% 	Ll	2.403979, std 8.97e-05, n=10000, 95%CI: 2.404e+00 - 2.404e+00
% 	C	0.427046, std 4.61e-15, n=10000, 95%CI: 4.270e-01 - 4.270e-01
% 	Cr	0.425634, std 7.33e-06, n=10000, 95%CI: 4.256e-01 - 4.256e-01
% 	Cl	0.422436, std 8.34e-08, n=10000, 95%CI: 4.224e-01 - 4.224e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	SWP	1.000000, std 00, n=10000, 95%CI: 001 - 001
% 	SWI	-0.513436, std 1.73e-02, n=10000, 95%CI: -5.475e-01 - -4.798e-01
% 	delta	-1.000000, std 00, n=10000, 95%CI: -001 - -001
% 	gamma	1.003317, std 1.73e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	lambda	0.999032, std 1.95e-04, n=10000, 95%CI: 9.987e-01 - 9.994e-01
% 	sigma	1.004289, std 1.97e-04, n=10000, 95%CI: 1.004e+00 - 1.005e+00
% 	omega	-0.009944, std 1.96e-04, n=10000, 95%CI: -1.034e-02 - -9.575e-03
% [*] Subject 25
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 36.059015274048 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=90
% 	L	2.187302, std 2.00e-13, n=10000, 95%CI: 2.187e+00 - 2.187e+00
% 	Lr	2.187302, std 6.93e-07, n=10000, 95%CI: 2.187e+00 - 2.187e+00
% 	Ll	2.204015, std 9.74e-05, n=10000, 95%CI: 2.204e+00 - 2.204e+00
% 	C	0.563785, std 3.55e-15, n=10000, 95%CI: 5.638e-01 - 5.638e-01
% 	Cr	0.563466, std 1.82e-06, n=10000, 95%CI: 5.635e-01 - 5.635e-01
% 	Cl	0.557975, std 4.47e-08, n=10000, 95%CI: 5.580e-01 - 5.580e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000002, std 4.15e-05, n=10000, 95%CI: 000 - 000
% 	SWP	0.999999, std 2.94e-05, n=10000, 95%CI: 001 - 001
% 	SWI	-0.058096, std 3.52e-04, n=10000, 95%CI: -5.874e-02 - -5.735e-02
% 	delta	-0.981000, std 1.94e-01, n=10000, 95%CI: -001 - -001
% 	gamma	1.000566, std 3.24e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000000, std 3.17e-07, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.000566, std 3.26e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.010412, std 3.28e-07, n=10000, 95%CI: -1.041e-02 - -1.041e-02
% [*] Subject 26
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 32.795680373907 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=90
% 	L	2.098083, std 2.27e-13, n=10000, 95%CI: 2.098e+00 - 2.098e+00
% 	Lr	2.096862, std 2.39e-04, n=10000, 95%CI: 2.096e+00 - 2.097e+00
% 	Ll	2.109659, std 6.24e-05, n=10000, 95%CI: 2.110e+00 - 2.110e+00
% 	C	0.497411, std 1.24e-14, n=10000, 95%CI: 4.974e-01 - 4.974e-01
% 	Cr	0.495301, std 7.94e-06, n=10000, 95%CI: 4.953e-01 - 4.953e-01
% 	Cl	0.491447, std 4.80e-08, n=10000, 95%CI: 4.914e-01 - 4.914e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.095142, std 1.68e-02, n=10000, 95%CI: 6.342e-02 - 1.295e-01
% 	SWP	0.932725, std 1.19e-02, n=10000, 95%CI: 9.084e-01 - 9.551e-01
% 	SWI	-0.495299, std 9.59e-03, n=10000, 95%CI: -5.135e-01 - -4.757e-01
% 	delta	1.000000, std 00, n=10000, 95%CI: 001 - 001
% 	gamma	1.004260, std 1.61e-05, n=10000, 95%CI: 1.004e+00 - 1.004e+00
% 	lambda	1.000583, std 1.14e-04, n=10000, 95%CI: 1.000e+00 - 1.001e+00
% 	sigma	1.003675, std 1.15e-04, n=10000, 95%CI: 1.003e+00 - 1.004e+00
% 	omega	-0.012719, std 1.14e-04, n=10000, 95%CI: -1.296e-02 - -1.251e-02
% [*] Subject 27
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 34.938577771187 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=76
% 	L	1.832233, std 5.93e-14, n=10000, 95%CI: 1.832e+00 - 1.832e+00
% 	Lr	1.832233, std 5.78e-14, n=10000, 95%CI: 1.832e+00 - 1.832e+00
% 	Ll	1.849333, std 5.14e-05, n=10000, 95%CI: 1.849e+00 - 1.849e+00
% 	C	0.670679, std 8.46e-14, n=10000, 95%CI: 6.707e-01 - 6.707e-01
% 	Cr	0.669979, std 5.02e-06, n=10000, 95%CI: 6.700e-01 - 6.700e-01
% 	Cl	0.662541, std 2.84e-08, n=10000, 95%CI: 6.625e-01 - 6.625e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000000, std 1.54e-13, n=10000, 95%CI: 000 - 5.074e-13
% 	SWP	1.000000, std 1.27e-13, n=10000, 95%CI: 1.000e+00 - 001
% 	SWI	-0.094220, std 7.39e-04, n=10000, 95%CI: -9.556e-02 - -9.265e-02
% 	delta	0.314600, std 9.49e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.001046, std 7.50e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000000, std 2.15e-15, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001046, std 7.50e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.012283, std 4.33e-08, n=10000, 95%CI: -1.228e-02 - -1.228e-02
% [*] Subject 28
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 28.014603286982 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=62
% 	L	1.762562, std 1.23e-13, n=10000, 95%CI: 1.763e+00 - 1.763e+00
% 	Lr	1.762558, std 1.17e-05, n=10000, 95%CI: 1.763e+00 - 1.763e+00
% 	Ll	1.779462, std 6.57e-05, n=10000, 95%CI: 1.779e+00 - 1.780e+00
% 	C	0.644091, std 8.66e-14, n=10000, 95%CI: 6.441e-01 - 6.441e-01
% 	Cr	0.642177, std 1.23e-05, n=10000, 95%CI: 6.422e-01 - 6.422e-01
% 	Cl	0.634267, std 6.76e-08, n=10000, 95%CI: 6.343e-01 - 6.343e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000202, std 6.89e-04, n=10000, 95%CI: 000 - 2.506e-03
% 	SWP	0.999857, std 4.87e-04, n=10000, 95%CI: 9.982e-01 - 001
% 	SWI	-0.241916, std 1.94e-03, n=10000, 95%CI: -2.454e-01 - -2.379e-01
% 	delta	0.729800, std 6.84e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.002980, std 1.92e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	lambda	1.000002, std 6.63e-06, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.002978, std 2.03e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	omega	-0.015491, std 6.63e-06, n=10000, 95%CI: -1.551e-02 - -1.549e-02
% [*] Subject 29
% p - human fraction of ROIs significant: 0.9664
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 11.642734229565 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=111
% 	L	1.790157, std 4.31e-14, n=10000, 95%CI: 1.790e+00 - 1.790e+00
% 	Lr	1.790201, std 3.85e-04, n=10000, 95%CI: 1.789e+00 - 1.791e+00
% 	Ll	1.815808, std 6.24e-05, n=10000, 95%CI: 1.816e+00 - 1.816e+00
% 	C	0.577062, std 5.50e-14, n=10000, 95%CI: 5.771e-01 - 5.771e-01
% 	Cr	0.558412, std 5.06e-05, n=10000, 95%CI: 5.583e-01 - 5.585e-01
% 	Cl	0.558388, std 1.78e-07, n=10000, 95%CI: 5.584e-01 - 5.584e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.005054, std 8.02e-03, n=10000, 95%CI: 000 - 2.709e-02
% 	SWP	0.996426, std 5.67e-03, n=10000, 95%CI: 9.808e-01 - 001
% 	SWI	-3445.440910, std 3.48e+05, n=10000, 95%CI: -5.433e+03 - 5.077e+03
% 	delta	-0.100200, std 9.95e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.033397, std 9.37e-05, n=10000, 95%CI: 1.033e+00 - 1.034e+00
% 	lambda	0.999976, std 2.15e-04, n=10000, 95%CI: 9.996e-01 - 1.000e+00
% 	sigma	1.033423, std 2.34e-04, n=10000, 95%CI: 1.033e+00 - 1.034e+00
% 	omega	-0.033417, std 2.15e-04, n=10000, 95%CI: -3.384e-02 - -3.300e-02
% [*] Subject 30
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 6.615686357021 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=19
% 	L	2.231136, std 1.80e-13, n=10000, 95%CI: 2.231e+00 - 2.231e+00
% 	Lr	2.231110, std 1.37e-04, n=10000, 95%CI: 2.231e+00 - 2.231e+00
% 	Ll	2.231136, std 1.79e-13, n=10000, 95%CI: 2.231e+00 - 2.231e+00
% 	C	0.618998, std 5.85e-14, n=10000, 95%CI: 6.190e-01 - 6.190e-01
% 	Cr	0.617522, std 9.53e-05, n=10000, 95%CI: 6.174e-01 - 6.177e-01
% 	Cl	0.617284, std 3.33e-06, n=10000, 95%CI: 6.173e-01 - 6.173e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	NaN, std NaN, n=10000, 95%CI: -005 - 001
% 	SWP	NaN, std NaN, n=10000, 95%CI: -2.536e+00 - 001
% 	SWI	NaN, std NaN, n=10000, 95%CI: -Inf - 3.703e+01
% 	delta	NaN, std NaN, n=10000, 95%CI: -003 - 001
% 	gamma	1.002390, std 1.55e-04, n=10000, 95%CI: 1.002e+00 - 1.003e+00
% 	lambda	1.000012, std 6.12e-05, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.002379, std 1.66e-04, n=10000, 95%CI: 1.002e+00 - 1.003e+00
% 	omega	-0.002788, std 6.15e-05, n=10000, 95%CI: -2.959e-03 - -2.765e-03
% [*] Subject 31
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 11.157685965300 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=26
% 	L	1.903606, std 1.91e-13, n=10000, 95%CI: 1.904e+00 - 1.904e+00
% 	Lr	1.903606, std 1.92e-13, n=10000, 95%CI: 1.904e+00 - 1.904e+00
% 	Ll	1.952813, std 3.74e-04, n=10000, 95%CI: 1.952e+00 - 1.954e+00
% 	C	0.660796, std 2.54e-14, n=10000, 95%CI: 6.608e-01 - 6.608e-01
% 	Cr	0.659431, std 4.10e-05, n=10000, 95%CI: 6.594e-01 - 6.595e-01
% 	Cl	0.636503, std 1.02e-06, n=10000, 95%CI: 6.365e-01 - 6.365e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000000, std 2.08e-14, n=10000, 95%CI: 000 - 7.192e-14
% 	SWP	1.000000, std 2.17e-14, n=10000, 95%CI: 1.000e+00 - 001
% 	SWI	-0.059557, std 1.89e-03, n=10000, 95%CI: -6.277e-02 - -5.545e-02
% 	delta	0.506400, std 8.62e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.002071, std 6.23e-05, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	lambda	1.000000, std 8.42e-16, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.002071, std 6.23e-05, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	omega	-0.038167, std 1.67e-06, n=10000, 95%CI: -3.817e-02 - -3.816e-02
% [*] Subject 32
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 9.955208182335 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=31
% 	L	2.621203, std 4.24e-13, n=10000, 95%CI: 2.621e+00 - 2.621e+00
% 	Lr	2.621145, std 1.72e-04, n=10000, 95%CI: 2.621e+00 - 2.621e+00
% 	Ll	2.621203, std 4.25e-13, n=10000, 95%CI: 2.621e+00 - 2.621e+00
% 	C	0.500479, std 7.31e-14, n=10000, 95%CI: 5.005e-01 - 5.005e-01
% 	Cr	0.500037, std 1.97e-05, n=10000, 95%CI: 5.000e-01 - 5.001e-01
% 	Cl	0.499980, std 1.78e-06, n=10000, 95%CI: 5.000e-01 - 5.000e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	NaN, std NaN, n=10000, 95%CI: -6.667e-01 - NaN
% 	SWP	NaN, std NaN, n=10000, 95%CI: 2.929e-01 - NaN
% 	SWI	NaN, std NaN, n=10000, 95%CI: -6.340e+01 - Inf
% 	delta	NaN, std NaN, n=10000, 95%CI: -003 - NaN
% 	gamma	1.000882, std 3.93e-05, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000022, std 6.55e-05, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.000860, std 7.77e-05, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.001020, std 6.56e-05, n=10000, 95%CI: -1.247e-03 - -9.902e-04
% [*] Subject 33
% p - human fraction of ROIs significant: 0.9993
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 30.689817279577 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=76
% 	L	2.026370, std 2.04e-14, n=10000, 95%CI: 2.026e+00 - 2.026e+00
% 	Lr	2.026874, std 8.66e-05, n=10000, 95%CI: 2.027e+00 - 2.027e+00
% 	Ll	2.042578, std 8.86e-05, n=10000, 95%CI: 2.042e+00 - 2.043e+00
% 	C	0.518498, std 9.50e-14, n=10000, 95%CI: 5.185e-01 - 5.185e-01
% 	Cr	0.517878, std 5.26e-06, n=10000, 95%CI: 5.179e-01 - 5.179e-01
% 	Cl	0.512574, std 6.53e-08, n=10000, 95%CI: 5.126e-01 - 5.126e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	SWP	1.000000, std 00, n=10000, 95%CI: 001 - 001
% 	SWI	-0.120646, std 1.25e-03, n=10000, 95%CI: -1.229e-01 - -1.181e-01
% 	delta	-1.000000, std 00, n=10000, 95%CI: -001 - -001
% 	gamma	1.001197, std 1.02e-05, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	0.999752, std 4.27e-05, n=10000, 95%CI: 9.997e-01 - 9.998e-01
% 	sigma	1.001446, std 4.28e-05, n=10000, 95%CI: 1.001e+00 - 1.002e+00
% 	omega	-0.011310, std 4.28e-05, n=10000, 95%CI: -1.140e-02 - -1.124e-02
% [*] Subject 34
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 23.016997545958 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=64
% 	L	2.312965, std 3.74e-13, n=10000, 95%CI: 2.313e+00 - 2.313e+00
% 	Lr	2.312964, std 7.82e-06, n=10000, 95%CI: 2.313e+00 - 2.313e+00
% 	Ll	2.335360, std 1.04e-04, n=10000, 95%CI: 2.335e+00 - 2.336e+00
% 	C	0.563212, std 8.15e-14, n=10000, 95%CI: 5.632e-01 - 5.632e-01
% 	Cr	0.562579, std 7.46e-06, n=10000, 95%CI: 5.626e-01 - 5.626e-01
% 	Cl	0.555318, std 8.31e-08, n=10000, 95%CI: 5.553e-01 - 5.553e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000046, std 3.48e-04, n=10000, 95%CI: 000 - 3.223e-04
% 	SWP	0.999968, std 2.46e-04, n=10000, 95%CI: 9.998e-01 - 001
% 	SWI	-0.087201, std 1.12e-03, n=10000, 95%CI: -8.922e-02 - -8.483e-02
% 	delta	-0.201000, std 9.80e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.001126, std 1.33e-05, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000000, std 3.38e-06, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001125, std 1.37e-05, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.014216, std 3.38e-06, n=10000, 95%CI: -1.422e-02 - -1.422e-02
% [*] Subject 35
% p - human fraction of ROIs significant: 0.9340
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 6.515794515610 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=76
% 	L	1.701250, std 8.28e-14, n=10000, 95%CI: 1.701e+00 - 1.701e+00
% 	Lr	1.696758, std 6.78e-04, n=10000, 95%CI: 1.695e+00 - 1.698e+00
% 	Ll	1.754081, std 6.34e-05, n=10000, 95%CI: 1.754e+00 - 1.754e+00
% 	C	0.608787, std 1.29e-13, n=10000, 95%CI: 6.088e-01 - 6.088e-01
% 	Cr	0.584540, std 1.93e-04, n=10000, 95%CI: 5.842e-01 - 5.849e-01
% 	Cl	0.584767, std 3.92e-08, n=10000, 95%CI: 5.848e-01 - 5.848e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.078235, std 1.09e-02, n=10000, 95%CI: 5.696e-02 - 9.925e-02
% 	SWP	0.944680, std 7.71e-03, n=10000, 95%CI: 9.298e-01 - 9.597e-01
% 	SWI	112.689760, std 5.79e+03, n=10000, 95%CI: -8.091e+02 - 9.909e+02
% 	delta	1.000000, std 00, n=10000, 95%CI: 001 - 001
% 	gamma	1.041480, std 3.45e-04, n=10000, 95%CI: 1.041e+00 - 1.042e+00
% 	lambda	1.002648, std 4.00e-04, n=10000, 95%CI: 1.002e+00 - 1.003e+00
% 	sigma	1.038730, std 5.10e-04, n=10000, 95%CI: 1.038e+00 - 1.040e+00
% 	omega	-0.043717, std 3.98e-04, n=10000, 95%CI: -4.450e-02 - -4.295e-02
% [*] Subject 36
% p - human fraction of ROIs significant: 0.9990
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 31.439570337534 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=76
% 	L	1.853524, std 1.14e-13, n=10000, 95%CI: 1.854e+00 - 1.854e+00
% 	Lr	1.853979, std 1.18e-04, n=10000, 95%CI: 1.854e+00 - 1.854e+00
% 	Ll	1.866583, std 8.05e-05, n=10000, 95%CI: 1.866e+00 - 1.867e+00
% 	C	0.570401, std 5.19e-14, n=10000, 95%CI: 5.704e-01 - 5.704e-01
% 	Cr	0.568899, std 9.38e-06, n=10000, 95%CI: 5.689e-01 - 5.689e-01
% 	Cl	0.563756, std 6.92e-08, n=10000, 95%CI: 5.638e-01 - 5.638e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000000, std 1.96e-05, n=10000, 95%CI: 000 - 000
% 	SWP	1.000000, std 1.39e-05, n=10000, 95%CI: 001 - 001
% 	SWI	-0.302618, std 3.55e-03, n=10000, 95%CI: -3.092e-01 - -2.953e-01
% 	delta	-0.999800, std 2.00e-02, n=10000, 95%CI: -001 - -001
% 	gamma	1.002640, std 1.65e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	lambda	0.999755, std 6.39e-05, n=10000, 95%CI: 9.996e-01 - 9.999e-01
% 	sigma	1.002886, std 6.45e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	omega	-0.011541, std 6.39e-05, n=10000, 95%CI: -1.168e-02 - -1.143e-02
% [*] Subject 37
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 21.770607680082 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=62
% 	L	2.060114, std 1.59e-13, n=10000, 95%CI: 2.060e+00 - 2.060e+00
% 	Lr	2.052466, std 7.07e-04, n=10000, 95%CI: 2.051e+00 - 2.054e+00
% 	Ll	2.073945, std 9.10e-05, n=10000, 95%CI: 2.074e+00 - 2.074e+00
% 	C	0.543434, std 1.32e-14, n=10000, 95%CI: 5.434e-01 - 5.434e-01
% 	Cr	0.540196, std 2.23e-05, n=10000, 95%CI: 5.402e-01 - 5.402e-01
% 	Cl	0.534719, std 9.96e-08, n=10000, 95%CI: 5.347e-01 - 5.347e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.355352, std 2.13e-02, n=10000, 95%CI: 3.133e-01 - 3.960e-01
% 	SWP	0.748728, std 1.50e-02, n=10000, 95%CI: 7.200e-01 - 7.785e-01
% 	SWI	-0.381221, std 1.28e-02, n=10000, 95%CI: -4.068e-01 - -3.566e-01
% 	delta	1.000000, std 00, n=10000, 95%CI: 001 - 001
% 	gamma	1.005996, std 4.15e-05, n=10000, 95%CI: 1.006e+00 - 1.006e+00
% 	lambda	1.003726, std 3.46e-04, n=10000, 95%CI: 1.003e+00 - 1.004e+00
% 	sigma	1.002261, std 3.43e-04, n=10000, 95%CI: 1.002e+00 - 1.003e+00
% 	omega	-0.020012, std 3.43e-04, n=10000, 95%CI: -2.070e-02 - -1.936e-02
% [*] Subject 38
% p - human fraction of ROIs significant: 0.9974
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 40.116965681314 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=104
% 	L	1.718523, std 2.84e-13, n=10000, 95%CI: 1.719e+00 - 1.719e+00
% 	Lr	1.718396, std 1.64e-04, n=10000, 95%CI: 1.718e+00 - 1.719e+00
% 	Ll	1.727613, std 2.22e-05, n=10000, 95%CI: 1.728e+00 - 1.728e+00
% 	C	0.597518, std 6.55e-14, n=10000, 95%CI: 5.975e-01 - 5.975e-01
% 	Cr	0.594835, std 1.03e-05, n=10000, 95%CI: 5.948e-01 - 5.949e-01
% 	Cl	0.592097, std 1.43e-08, n=10000, 95%CI: 5.921e-01 - 5.921e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.015708, std 1.42e-02, n=10000, 95%CI: 000 - 4.744e-02
% 	SWP	0.988893, std 1.01e-02, n=10000, 95%CI: 9.664e-01 - 001
% 	SWI	-0.966646, std 1.81e-02, n=10000, 95%CI: -1.002e+00 - -9.312e-01
% 	delta	0.560600, std 8.28e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.004510, std 1.74e-05, n=10000, 95%CI: 1.004e+00 - 1.005e+00
% 	lambda	1.000074, std 9.52e-05, n=10000, 95%CI: 9.999e-01 - 1.000e+00
% 	sigma	1.004436, std 9.59e-05, n=10000, 95%CI: 1.004e+00 - 1.005e+00
% 	omega	-0.009229, std 9.51e-05, n=10000, 95%CI: -9.419e-03 - -9.045e-03
% [*] Subject 39
% p - human fraction of ROIs significant: 0.9993
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 30.585854589939 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=90
% 	L	2.304371, std 2.30e-13, n=10000, 95%CI: 2.304e+00 - 2.304e+00
% 	Lr	2.304425, std 2.22e-04, n=10000, 95%CI: 2.304e+00 - 2.305e+00
% 	Ll	2.317450, std 7.88e-05, n=10000, 95%CI: 2.317e+00 - 2.318e+00
% 	C	0.462232, std 6.35e-14, n=10000, 95%CI: 4.622e-01 - 4.622e-01
% 	Cr	0.460834, std 6.40e-06, n=10000, 95%CI: 4.608e-01 - 4.608e-01
% 	Cl	0.457155, std 6.17e-08, n=10000, 95%CI: 4.572e-01 - 4.572e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.004911, std 8.81e-03, n=10000, 95%CI: 000 - 2.974e-02
% 	SWP	0.996527, std 6.23e-03, n=10000, 95%CI: 9.790e-01 - 001
% 	SWI	-0.381602, std 6.91e-03, n=10000, 95%CI: -3.943e-01 - -3.678e-01
% 	delta	-0.234600, std 9.72e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.003033, std 1.39e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	lambda	0.999977, std 9.65e-05, n=10000, 95%CI: 9.998e-01 - 1.000e+00
% 	sigma	1.003057, std 9.77e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	omega	-0.011082, std 9.65e-05, n=10000, 95%CI: -1.128e-02 - -1.091e-02
% [*] Subject 40
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 59.124808937311 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=153
% 	L	2.150920, std 1.30e-13, n=10000, 95%CI: 2.151e+00 - 2.151e+00
% 	Lr	2.150919, std 3.06e-06, n=10000, 95%CI: 2.151e+00 - 2.151e+00
% 	Ll	2.150920, std 1.36e-13, n=10000, 95%CI: 2.151e+00 - 2.151e+00
% 	C	0.564092, std 7.37e-14, n=10000, 95%CI: 5.641e-01 - 5.641e-01
% 	Cr	0.563476, std 1.38e-06, n=10000, 95%CI: 5.635e-01 - 5.635e-01
% 	Cl	0.563475, std 1.28e-08, n=10000, 95%CI: 5.635e-01 - 5.635e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	NaN, std NaN, n=10000, 95%CI: 000 - 001
% 	SWP	NaN, std NaN, n=10000, 95%CI: 2.929e-01 - 001
% 	SWI	NaN, std NaN, n=10000, 95%CI: -1.271e+04 - 1.727e+04
% 	delta	NaN, std NaN, n=10000, 95%CI: -001 - 001
% 	gamma	1.001093, std 2.45e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	lambda	1.000000, std 1.42e-06, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001093, std 2.85e-06, n=10000, 95%CI: 1.001e+00 - 1.001e+00
% 	omega	-0.001095, std 1.42e-06, n=10000, 95%CI: -1.100e-03 - -1.095e-03
% [*] Subject 41
% p - human fraction of ROIs significant: 0.9996
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 23.138008564711 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=67
% 	L	2.247762, std 2.68e-13, n=10000, 95%CI: 2.248e+00 - 2.248e+00
% 	Lr	2.248111, std 1.12e-04, n=10000, 95%CI: 2.248e+00 - 2.248e+00
% 	Ll	2.248549, std 4.32e-05, n=10000, 95%CI: 2.248e+00 - 2.249e+00
% 	C	0.543166, std 1.18e-14, n=10000, 95%CI: 5.432e-01 - 5.432e-01
% 	Cr	0.541666, std 1.08e-05, n=10000, 95%CI: 5.416e-01 - 5.417e-01
% 	Cl	0.541646, std 1.61e-07, n=10000, 95%CI: 5.416e-01 - 5.416e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000371, std 5.99e-03, n=10000, 95%CI: 000 - 000
% 	SWP	0.999738, std 4.23e-03, n=10000, 95%CI: 001 - 001
% 	SWI	-194.093030, std 3.62e+03, n=10000, 95%CI: -9.511e+02 - -4.653e+01
% 	delta	-0.988600, std 1.51e-01, n=10000, 95%CI: -001 - -001
% 	gamma	1.002769, std 2.00e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	lambda	0.999844, std 4.98e-05, n=10000, 95%CI: 9.998e-01 - 1.000e+00
% 	sigma	1.002925, std 5.26e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	omega	-0.002651, std 4.98e-05, n=10000, 95%CI: -2.763e-03 - -2.567e-03
% [*] Subject 42
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 30.304634630680 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=83
% 	L	2.313146, std 4.70e-13, n=10000, 95%CI: 2.313e+00 - 2.313e+00
% 	Lr	2.313133, std 2.29e-05, n=10000, 95%CI: 2.313e+00 - 2.313e+00
% 	Ll	2.313146, std 4.65e-13, n=10000, 95%CI: 2.313e+00 - 2.313e+00
% 	C	0.521813, std 2.55e-15, n=10000, 95%CI: 5.218e-01 - 5.218e-01
% 	Cr	0.520989, std 4.25e-06, n=10000, 95%CI: 5.210e-01 - 5.210e-01
% 	Cl	0.520976, std 7.93e-08, n=10000, 95%CI: 5.210e-01 - 5.210e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	NaN, std NaN, n=10000, 95%CI: -006 - 001
% 	SWP	NaN, std NaN, n=10000, 95%CI: -3.243e+00 - 001
% 	SWI	NaN, std NaN, n=10000, 95%CI: -1.530e+03 - 4.637e+02
% 	delta	NaN, std NaN, n=10000, 95%CI: -003 - 001
% 	gamma	1.001581, std 8.17e-06, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	lambda	1.000005, std 9.91e-06, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.001576, std 1.30e-05, n=10000, 95%CI: 1.002e+00 - 1.002e+00
% 	omega	-0.001611, std 9.91e-06, n=10000, 95%CI: -1.641e-03 - -1.606e-03
% [*] Subject 43
% p - human fraction of ROIs significant: 0.9997
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 30.079671561718 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=79
% 	L	2.065805, std 3.25e-13, n=10000, 95%CI: 2.066e+00 - 2.066e+00
% 	Lr	2.065447, std 1.97e-04, n=10000, 95%CI: 2.065e+00 - 2.066e+00
% 	Ll	2.066337, std 2.43e-05, n=10000, 95%CI: 2.066e+00 - 2.066e+00
% 	C	0.498696, std 5.70e-14, n=10000, 95%CI: 4.987e-01 - 4.987e-01
% 	Cr	0.496635, std 9.67e-06, n=10000, 95%CI: 4.966e-01 - 4.967e-01
% 	Cl	0.496622, std 1.03e-07, n=10000, 95%CI: 4.966e-01 - 4.966e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.373746, std 1.41e-01, n=10000, 95%CI: 1.068e-02 - 5.960e-01
% 	SWP	0.735722, std 1.00e-01, n=10000, 95%CI: 5.784e-01 - 9.923e-01
% 	SWI	39.856160, std 1.18e+04, n=10000, 95%CI: -1.084e+03 - 8.045e+02
% 	delta	0.955200, std 2.96e-01, n=10000, 95%CI: 001 - 001
% 	gamma	1.004149, std 1.95e-05, n=10000, 95%CI: 1.004e+00 - 1.004e+00
% 	lambda	1.000173, std 9.54e-05, n=10000, 95%CI: 1.000e+00 - 1.000e+00
% 	sigma	1.003975, std 9.76e-05, n=10000, 95%CI: 1.004e+00 - 1.004e+00
% 	omega	-0.004348, std 9.53e-05, n=10000, 95%CI: -4.552e-03 - -4.177e-03
% [*] Subject 44
% p - human fraction of ROIs significant: 0.9867
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 16.650530606508 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=62
% 	L	1.964147, std 1.59e-13, n=10000, 95%CI: 1.964e+00 - 1.964e+00
% 	Lr	1.985106, std 7.42e-04, n=10000, 95%CI: 1.984e+00 - 1.987e+00
% 	Ll	1.997011, std 1.07e-04, n=10000, 95%CI: 1.997e+00 - 1.997e+00
% 	C	0.512211, std 5.72e-14, n=10000, 95%CI: 5.122e-01 - 5.122e-01
% 	Cr	0.508211, std 4.72e-05, n=10000, 95%CI: 5.081e-01 - 5.083e-01
% 	Cl	0.506950, std 2.22e-07, n=10000, 95%CI: 5.069e-01 - 5.070e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	SWP	1.000000, std 00, n=10000, 95%CI: 001 - 001
% 	SWI	-8.803601, std 6.47e-01, n=10000, 95%CI: -1.010e+01 - -7.570e+00
% 	delta	-1.000000, std 00, n=10000, 95%CI: -001 - -001
% 	gamma	1.007872, std 9.36e-05, n=10000, 95%CI: 1.008e+00 - 1.008e+00
% 	lambda	0.989442, std 3.70e-04, n=10000, 95%CI: 9.887e-01 - 9.902e-01
% 	sigma	1.018626, std 3.79e-04, n=10000, 95%CI: 1.018e+00 - 1.019e+00
% 	omega	0.000292, std 3.78e-04, n=10000, 95%CI: -4.679e-04 - 1.008e-03
% [*] Subject 45
% p - human fraction of ROIs significant: 0.9996
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 41.319707363844 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=97
% 	L	1.780932, std 2.31e-14, n=10000, 95%CI: 1.781e+00 - 1.781e+00
% 	Lr	1.781007, std 5.65e-05, n=10000, 95%CI: 1.781e+00 - 1.781e+00
% 	Ll	1.781434, std 1.45e-05, n=10000, 95%CI: 1.781e+00 - 1.781e+00
% 	C	0.579171, std 1.39e-13, n=10000, 95%CI: 5.792e-01 - 5.792e-01
% 	Cr	0.576744, std 8.10e-06, n=10000, 95%CI: 5.767e-01 - 5.768e-01
% 	Cl	0.576751, std 5.23e-08, n=10000, 95%CI: 5.768e-01 - 5.768e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.005463, std 2.21e-02, n=10000, 95%CI: 000 - 8.155e-02
% 	SWP	0.996137, std 1.56e-02, n=10000, 95%CI: 9.421e-01 - 001
% 	SWI	1805.741329, std 1.10e+05, n=10000, 95%CI: -3.250e+03 - 3.525e+03
% 	delta	-0.805400, std 5.93e-01, n=10000, 95%CI: -001 - 001
% 	gamma	1.004209, std 1.41e-05, n=10000, 95%CI: 1.004e+00 - 1.004e+00
% 	lambda	0.999958, std 3.17e-05, n=10000, 95%CI: 9.999e-01 - 1.000e+00
% 	sigma	1.004251, std 3.44e-05, n=10000, 95%CI: 1.004e+00 - 1.004e+00
% 	omega	-0.004154, std 3.17e-05, n=10000, 95%CI: -4.221e-03 - -4.097e-03
% [*] Subject 46
% p - human fraction of ROIs significant: 0.5683
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 0.000000000000 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=104
% 	L	1.767345, std 3.85e-13, n=10000, 95%CI: 1.767e+00 - 1.767e+00
% 	Lr	1.810987, std 8.44e-04, n=10000, 95%CI: 1.809e+00 - 1.813e+00
% 	Ll	2.336869, std 1.08e-04, n=10000, 95%CI: 2.337e+00 - 2.337e+00
% 	C	0.554398, std 8.75e-14, n=10000, 95%CI: 5.544e-01 - 5.544e-01
% 	Cr	0.382990, std 1.07e-03, n=10000, 95%CI: 3.810e-01 - 3.851e-01
% 	Cl	0.527832, std 8.00e-07, n=10000, 95%CI: 5.278e-01 - 5.278e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	SWP	1.000000, std 00, n=10000, 95%CI: 001 - 001
% 	SWI	1.281635, std 2.79e-03, n=10000, 95%CI: 1.276e+00 - 1.287e+00
% 	delta	-1.000000, std 00, n=10000, 95%CI: -001 - -001
% 	gamma	1.447562, std 4.05e-03, n=10000, 95%CI: 1.439e+00 - 1.455e+00
% 	lambda	0.975902, std 4.55e-04, n=10000, 95%CI: 9.750e-01 - 9.768e-01
% 	sigma	1.483307, std 4.04e-03, n=10000, 95%CI: 1.475e+00 - 1.491e+00
% 	omega	-0.025637, std 4.78e-04, n=10000, 95%CI: -2.656e-02 - -2.468e-02
% [*] Subject 47
% p - human fraction of ROIs significant: 0.9145
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 4.424398332834 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=104
% 	L	1.845552, std 4.49e-13, n=10000, 95%CI: 1.846e+00 - 1.846e+00
% 	Lr	1.835374, std 7.86e-04, n=10000, 95%CI: 1.834e+00 - 1.837e+00
% 	Ll	1.975078, std 3.36e-05, n=10000, 95%CI: 1.975e+00 - 1.975e+00
% 	C	0.548021, std 6.54e-14, n=10000, 95%CI: 5.480e-01 - 5.480e-01
% 	Cr	0.520909, std 1.44e-04, n=10000, 95%CI: 5.206e-01 - 5.212e-01
% 	Cl	0.523168, std 1.32e-07, n=10000, 95%CI: 5.232e-01 - 5.232e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	0.072823, std 5.22e-03, n=10000, 95%CI: 6.233e-02 - 8.291e-02
% 	SWP	0.948506, std 3.69e-03, n=10000, 95%CI: 9.414e-01 - 9.559e-01
% 	SWI	11.169681, std 6.78e-01, n=10000, 95%CI: 1.005e+01 - 1.267e+01
% 	delta	1.000000, std 00, n=10000, 95%CI: 001 - 001
% 	gamma	1.052048, std 2.90e-04, n=10000, 95%CI: 1.051e+00 - 1.053e+00
% 	lambda	1.005546, std 4.31e-04, n=10000, 95%CI: 1.005e+00 - 1.006e+00
% 	sigma	1.046246, std 5.07e-04, n=10000, 95%CI: 1.045e+00 - 1.047e+00
% 	omega	-0.053019, std 4.26e-04, n=10000, 95%CI: -5.385e-02 - -5.217e-02
% [*] Subject 48
% p - human fraction of ROIs significant: 1.0000
% === Frequency band: 1 ===
% [*] Clustering random network..
% cluster clash: 23.689980894327 mm
% [!] === Weighted small world index calculation ===
% [*] Calculating human SWP with 10000 permutations..
% [*] Adj, nodes=69
% 	L	2.346720, std 6.48e-14, n=10000, 95%CI: 2.347e+00 - 2.347e+00
% 	Lr	2.346584, std 1.01e-04, n=10000, 95%CI: 2.346e+00 - 2.347e+00
% 	Ll	2.346720, std 6.30e-14, n=10000, 95%CI: 2.347e+00 - 2.347e+00
% 	C	0.499354, std 1.51e-14, n=10000, 95%CI: 4.994e-01 - 4.994e-01
% 	Cr	0.497781, std 9.98e-06, n=10000, 95%CI: 4.978e-01 - 4.978e-01
% 	Cl	0.497761, std 8.19e-08, n=10000, 95%CI: 4.978e-01 - 4.978e-01
% 	DC	0.000000, std 00, n=10000, 95%CI: 000 - 000
% 	DL	NaN, std NaN, n=10000, 95%CI: 000 - 001
% 	SWP	NaN, std NaN, n=10000, 95%CI: 2.929e-01 - 001
% 	SWI	NaN, std NaN, n=10000, 95%CI: -1.106e-06 - 2.534e-08
% 	delta	NaN, std NaN, n=10000, 95%CI: -001 - 001
% 	gamma	1.003160, std 2.01e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	lambda	1.000058, std 4.29e-05, n=10000, 95%CI: 001 - 1.000e+00
% 	sigma	1.003102, std 4.74e-05, n=10000, 95%CI: 1.003e+00 - 1.003e+00
% 	omega	-0.003259, std 4.29e-05, n=10000, 95%CI: -3.364e-03 - -3.201e-03