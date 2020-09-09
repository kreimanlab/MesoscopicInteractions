close all;
clear;
rng('shuffle');

%/media/klab/internal/data/coreg/fsaverage_sym/label/all_surf_ielvis_sub18.label
dir_artLp = '/media/klab/internal/data/h5_notch20/art';
dir_corLp = '/media/klab/internal/data/coreg';
dir_resLp = '/home/jerry/data/results/coh_w10';
setenv('SUBJECTS_DIR',dir_corLp);
dir_cacheLp = './cache';
dir_h5Lp = '/media/klab/internal/data/h5_notch20';


metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
%metricsp = {'pcBroadband'};

% Corresponding variable (11: CT, 12: Mag)
variable_idx = 12;
variable_idx_ct = 11;

% Coherence value for no significant interaction, NaN means don't include
% in average during dimensionality reduction.
val_no_sig = NaN;

% Final dimension of dimensionality reduction
ss_thresh_dim = 150;
        
% Load fsaverage_sym surface
trig_ignore_cache_reduce = false; %true;
trig_plot_all = false;
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_corLp,'fsaverage_sym','r','pial'));
[s_vert_inf, faces_inf] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_corLp,'fsaverage_sym','r','inflated'));
[s_vert_sphere, faces_sphere] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_corLp,'fsaverage_sym','r','sphere'));
sphere_r = max(s_vert_sphere(:));

trig_plot_mag_dist = false;

% Plot arc
trig_plot_arc = true;
tube_radius = 0.15;
tube_n_edges = 3;

trig_kmeans = false;

trig_ss = true;
n_cmap = 100;
for iM_local = 1:length(metricsp)
    
    Kmaster = 12; %8:2:48;
    for km = 1:length(Kmaster)
        
        fn_cache = [dir_cacheLp,'/xsub_out_all_',num2str(iM_local)];
        Ca = load(fn_cache);
        % build master coherence table
        D = [];
        E = [];
        E_sphere = [];
        
        n_total_pairs = 0;
        n_art_skip = 0;
        for iS = 1:length(Ca.Subjects)
            % read subject cache
            sid = Ca.Subjects{iS};
            sidint = str2num(sid(2:end));
            fprintf('[%s]\n',sid)
            Cas = load(sprintf('cache/xsub_out_%s_%i.mat',sid,iM_local));
            La = read_label('fsaverage_sym',sprintf('ielvis_%s',sid));

            if isempty(La)
                fprintf(2,'[*] Skip %s, no fsaverage_sym label\n',sid)
            else

                % Get all electrode coordinates
                for j = 1:Cas.ecog.n_bchan
                    b1c1 = Cas.ecog.bip(j,1);
                    % get fsaverage_sym coordinates of electrodes
                    %e_coords = s_vert(La(:,1) + 1,:);
                    lidx = La(:,end) == b1c1;
                    %coord1 = La(lidx,2:4);
                    coord1 = s_vert(La(lidx,1) + 1,:);
                    E = [E; [sidint, j, coord1]];
                    
                    coord1_sphere = s_vert_sphere(La(lidx,1) + 1,:);
                    E_sphere = [E_sphere; [sidint, j, coord1_sphere]];
                end

                % Get significant interaction coordinates
                art_frac_thresh = 0; % 9: keep everything, 1: keep nothing
                idx_afrac = ((Cas.n_no + Cas.n_yes) ./ Cas.n_graph) > art_frac_thresh;
                n_skipped_art = 0;
                
                for j = 1:length(Cas.ct)
                    % Consider bip pairs in subject that are significant and
                    % far enough apart
                    if ((Cas.ct(j) > Cas.ct_thresh) && (Cas.Dmats(j) > Cas.dist_thresh))
                        b1 = Cas.chan1(j) + 1;
                        b2 = Cas.chan2(j) + 1;
                        % get electrode numbers of first in bipolar pair
                        b1c1 = Cas.ecog.bip(b1,1);
                        b2c1 = Cas.ecog.bip(b2,1);
                        % get fsaverage_sym coordinates of electrodes
                        %e_coords = s_vert(La(:,1) + 1,:);
                        lidx = La(:,end) == b1c1;
                        %coord1 = La(lidx,2:4);
                        coord1 = s_vert(La(lidx,1) + 1,:);

                        lidx = La(:,end) == b2c1;
                        %coord2 = La(lidx,2:4);
                        coord2 = s_vert(La(lidx,1) + 1,:);


                        if (isempty(coord1) || isempty(coord2))
                            fprintf(2,'[!] coordinates empty for %s.\n',sid)
                        else
                            % --------------------------------------------------------------------------------
                            if (idx_afrac(j))
                                D = [D; [sidint, b1, b2, Cas.Dmats(j), coord1, coord2, Cas.ct(j), Cas.mag(j)]];
                            else
                                %fprintf('\t(!) skip %s: pair %i of %i - artifact fraction above %.3f\n',sid,j,length(Cas.ct),1-art_frac_thresh)
                                % D = [D; [sidint, b1, b2, Cas.Dmats(j), coord1, coord2, NaN, NaN]];
                                n_skipped_art = n_skipped_art + 1;
                                
%                                 % DEBUG: check original coh values
%                                 coh_tmp = h5read(sprintf('%s/%s_graph-pcBroadband.h5',dir_resLp,sid),'/R',[j 1],[1 Cas.n_graph]);
%                                 art_idx_tmp = h5read(sprintf('%s/%s_art.h5',dir_artLp,sid),'/art_idx',[j 1],[1 Cas.n_graph]);
%                                 fprintf(2,'[!] Debug percent artifacts: %.2f%%\n',100*sum((art_idx_tmp) > 0)/Cas.n_graph)
%                                 if (rand([1 1]) > 0)
%                                     return
%                                 end
                                
                            end
                        end
                    end
                end
                
                fprintf('[*] %s total skipped due to artifacts: %.2f%%\n',sid,100*(n_skipped_art/length(Cas.ct)));
                n_total_pairs = n_total_pairs + length(Cas.ct);
                n_art_skip = n_art_skip + n_skipped_art;
                
            end

        end
        sidintu = unique(D(:,1));
        n_sidintu = length(unique(D(:,1)));

        fprintf('[!] Total dataset skipped due to artifacts: %.2f%%\n',100*(n_art_skip)/n_total_pairs);
        fprintf('\t(*) Threshold for artifact frac: %.3f\n',art_frac_thresh);
    %     for i1 = 1:(Ca.n_rois-1)
    %         for i2 = (i1+1):Ca.n_rois
    %             AA = Ca.AdjAtl{i1,i2};
    %             AAs = Ca.AdjAtl_sid{i1,i2};
    %             for j = 1:length(AA)
    %                 % only consider significant interactions
    %                 if (AA(j) > 0)
    %                 end
    %             end
    %         end
    %     end


        if (trig_plot_mag_dist)
            h = figure('visible','off');

            subplot(2,1,1)
            % Plot magnitude distance for individual subjects
            cc_sub = jet(n_sidintu);
            for i = 1:n_sidintu
                Ds = D(D(:,1) == sidintu(i),:);
                dists = Ds(:,4);
                vars = Ds(:,variable_idx);
                plot(dists,vars,'.','color',cc_sub(i,:));
                hold on;
                % fit trendline
                p = polyfit(dists,vars,1);
                x = [Cas.dist_thresh max(dists)];
                plot(x,p(2) + x.*p(1),'-','color',cc_sub(i,:));
                hold on;
            end
            lz = cell(n_sidintu,1);
            for i = 1:length(lz)
                lz{i} = sprintf('m%i',sidintu(i));
            end
            %legend(lz)
            axis([Cas.dist_thresh 120 0 1]);
            set(gca,'TickDir','Out');
            set(gca,'Box','Off');
            xlabel('Distance (mm)');
            ylabel('Coherence');


            subplot(2,1,2)
            variable_idx = 11;
            % Plot magnitude distance for individual subjects
            cc_sub = jet(n_sidintu);
            for i = 1:n_sidintu
                Ds = D(D(:,1) == sidintu(i),:);
                dists = Ds(:,4);
                vars = Ds(:,variable_idx);
                plot(dists,vars,'.','color',cc_sub(i,:));
                hold on;
                % fit trendline
                p = polyfit(dists,vars,1);
                x = [Cas.dist_thresh max(dists)];
                plot(x,p(2) + x.*p(1),'-','color',cc_sub(i,:));
                hold on;
            end
            lz = cell(n_sidintu,1);
            for i = 1:length(lz)
                lz{i} = sprintf('m%i',sidintu(i));
            end
            %legend(lz)
            axis([Cas.dist_thresh 120 0 1]);
            set(gca,'TickDir','Out');
            set(gca,'Box','Off');
            xlabel('Distance (mm)');
            ylabel('Consistency across time');

        end


        % Spatial smoothing
        % collapse electrode-electrode adjacency matrix by combining electrodes
        % that are spatially close
        if (trig_ss)
            ss_thresh_mm = 10; % threshold in millimeters
            %ss_thresh_dim = 100;

            [n_E,~] = size(E);
            [n_D,~] = size(D);


            A = nan(n_E,n_E);
            A_nocov = nan(n_E,n_E);
            Act = nan(n_E,n_E);
            Act_nocov = nan(n_E,n_E);
            
            Ad = nan(n_E,n_E);
            A_xsub = false(n_E,n_E);

            fn_ca = sprintf('%s/fig_cluster2_A_%i.mat',dir_cacheLp,iM_local);
            if (exist(fn_ca,'file'))
                fprintf('[*] Found cached initial adjacency matrix..\n')
                Ca = load(fn_ca);
                A = Ca.A;
                Ad = Ca.Ad;
                Act = Ca.Act;
                Act_nocov = Ca.Act_nocov;
                A_xsub = Ca.A_xsub;
                A_nocov = Ca.A_nocov;
                fprintf('\tDone.\n')
            else
                fprintf('[*] Building initial adjacency matrix..\n')
                parfor i = 1:n_E % parfor
                    for j = 1:n_E


                        % Find D index
                        sidint_i = E(i,1);
                        b_i = E(i,2);
                        sidint_j = E(j,1);
                        b_j = E(j,2);

                        if (sidint_i == sidint_j)
                            % Same patients
                            A_xsub(i,j) = false;
                            a_idx = find( (D(:,1) == sidint_i) & ( ((D(:,2) == b_i) & (D(:,3) == b_j)) | ((D(:,2) == b_j) & (D(:,3) == b_i)) ) );
                            if (isempty(a_idx))
                                % No significant interaction
                                if (b_i == b_j)
                                    % Diagonal
                                    a = NaN;
                                    a_ct = NaN;
                                    a2 = NaN;
                                    a_ct2 = NaN;
                                else
                                    a = val_no_sig;
                                    a_ct = val_no_sig;
                                    a2 = 0;
                                    a_ct2 = 0;
                                end
                            else
                                % Significant interaction
                                a = D(a_idx,variable_idx);
                                a_ct = D(a_idx,variable_idx_ct);
                                a2 = a;
                                a_ct2 = a_ct;
                            end
                        else
                            % Cross patient comparisons are undefined
                            A_xsub(i,j) = true;
                            a = NaN;
                            a_ct = NaN;
                            a2 = NaN;
                            a_ct2 = NaN;
                        end


                        A(i,j) = a;
                        A_nocov(i,j) = a2;
                        Act(i,j) = a_ct;
                        Act_nocov(i,j) = a_ct2;

                        % Distance matrix
                        c_i = E(i,3:5);
                        c_j = E(j,3:5);
                        ad = sqrt(sum((c_i - c_j).^2));
                        if (ad == 0)
                            Ad(i,j) = NaN;
                        else
                            Ad(i,j) = ad;
                        end
                    end
                end
                save(fn_ca,'A','Ad','Act','Act_nocov','A_xsub','A_nocov','E','E_sphere');
                fprintf('\tDone.\n')

            end

            % Dimensionality reduction by iterative spatial combination
            A_sav = A;
            Ad_sav = Ad;
            A_nocov_sav = A_nocov;
            E_sav = E;
            E_sphere_sav = E_sphere;
            fprintf('[*] Starting spatial combination with threshold: %i mm\n',ss_thresh_mm);
            fn_ca2 = sprintf('%s/fig_cluster2_reduce_%i.mat',dir_cacheLp,iM_local);
            if (exist(fn_ca2,'file') && (~trig_ignore_cache_reduce))

                Ca2 = load(fn_ca2);
                A = Ca2.A;
                Ad = Ca2.Ad;
                A_nocov = Ca2.A_nocov;
                Act = Ca2.Act;
                Act_nocov = Ca2.Act_nocov;
                E = Ca2.E;
                E_sphere = Ca2.E_sphere;
                ss_thresh_mm = Ca2.ss_thresh_mm;
                ss_samesid = Ca2.ss_samesid;
                Es = Ca2.Es;
                Es2 = Ca2.Es2;
                Es_sphere = Ca2.Es_sphere;
                fprintf('\tFound cache with threshold: %i mm\n',ss_thresh_mm);
            else

                ss_samesid = [];
                
                % Build list of combined electrodes
                Es = cell(n_E,1);
                Es2 = cell(n_E,1);
                Es_sphere = cell(n_E,1);
                for ies = 1:n_E
                    n_sig = sum((~isnan(A_sav(ies,:)))&(A_sav(ies,:) ~= 0));
                    Es{ies} = E(ies,:);
                    Es2{ies} = [E(ies,:),n_sig];
                    Es_sphere{ies} = E_sphere(ies,:);
                end
                
                sat = true;
                while(sat)
                    % Find electrode pair with shortest distance
                    [ad_v1, ad_i1] = min(Ad); % row
                    [ad_v2, ad_i2] = min(ad_v1); % col
                    i_row = ad_i1(ad_i2);
                    i_col = ad_i2;

                    sidint1 = E(i_row,1);
                    b1 = E(i_row,2);
                    sidint2 = E(i_col,1);
                    b2 = E(i_col,2);
                    ad = Ad(i_row,i_col);

                    % Break condition
                    %if (ad > ss_thresh_mm)
                    if (length(A) < (ss_thresh_dim + 1))
                        fprintf('[*] While loop break.\n')
                        sat = false;
                        break
                    end

                    % Save whether elec pair is from same subject
                    if (sidint1 == sidint2)
                        is_same = 1;
                        %ss_samesid = [ss_samesid; [1, ad]];
                    else
                        is_same = 0;
                        %ss_samesid = [ss_samesid; [0, ad]];
                    end
                    ss_samesid = [ss_samesid; [is_same, ad]];

                    % overwrite E
                    
                    % compute new midpoint iteratively
                    %combE = [E(i_row,1:2),(E(i_row,3:5) + E(i_col,3:5))/2];
                    
                    % compute new midpoint from full history of electrodes
                    % Collect points
                    % -----------------------------------------------------------------------------------------------
                    Epts = [E(i_row,3:5) ; E(i_col,3:5)];
                    %Epts = [Es{i_row}(3:5); Es{i_col}(3:5)]; % full history
                    % Compute center
                    Ecenter = mean(Epts);
                    
                    % Assign to point "with the most territory"
                    % recompute distances
                    nn_territory = 4; % number of neighbors to use in territory calculation
                    Epd1 = nan(1,length(E));
                    c_i1 = Epts(1,:);
                    Epd2 = nan(1,length(E));
                    c_i2 = Epts(2,:);
                    for j = 1:length(E)
                        c_j = E(j,3:5);
                        d1 = sqrt(sum((c_i1 - c_j).^2));
                        d2 = sqrt(sum((c_i2 - c_j).^2));
                        Epd1(j) = d1;
                        Epd2(j) = d2;
                    end
                    Epd1s = sort(Epd1);
                    Epd2s = sort(Epd2);
                    area1 = sum(Epd1s(1:nn_territory));
                    area2 = sum(Epd2s(1:nn_territory));
                    
                    if (area1 == area2)
                        % in unlikely case of equal area, randomly choose
                        Ecenter = Epts(randi([1 2]),:);
                    elseif (area1 > area2)
                        Ecenter = Epts(1,:);
                    elseif (area2 > area1)
                        Ecenter = Epts(1,:);
                    end
                    
                    % Project center out to sphere surface
                    %Ecenter = sphere_r * Ecenter/norm(Ecenter);
                    % Find nearest vertex
                    [to_sphere_mm,new_i] = min(sqrt(sum((s_vert - Ecenter).^2,2)));

                    % Get coordinate at vertex
                    combE = [E(i_row,1:2), s_vert(new_i,:) ];
                    %combE = [E(i_row,1:2),(E(i_row,3:5) + E(i_col,3:5))/2];
                    
                    combE_sphere = [E(i_row,1:2), s_vert_sphere(new_i,:) ];
                    
                    E(i_row,:) = combE;
                    E_sphere(i_row,:) = combE_sphere;
                    Es{i_row} = [Es{i_row}; Es{i_col}];
                    Es2{i_row} = [Es2{i_row}; Es2{i_col}];
                    Es_sphere{i_row} = [Es_sphere{i_row}; Es_sphere{i_col}];

                    % compute combined row

                    % %Number of shared electrodes between i_row, i_col
                    % sum((~isnan(A(i_row,:))) & (~isnan(A(i_col,:))));
                   
                    combA = nanmean([A(i_row,:);A(i_col,:)]);
                    combAnc = nanmean([A_nocov(i_row,:);A_nocov(i_col,:)]);
                    combAct = nanmean([Act(i_row,:);Act(i_col,:)]);
                    combAct_nc = nanmean([Act_nocov(i_row,:);Act_nocov(i_col,:)]);

                    % recompute distances
                    combAd = nan(1,length(E));
                    c_i = combE(3:5);
                    for j = 1:length(E)
                        c_j = E(j,3:5);
                        ad2 = sqrt(sum((c_i - c_j).^2));
                        combAd(j) = ad2;
                    end
                    combAd(combAd == 0) = NaN;
                    %combAd = nanmean([Ad(i_row,:);Ad(i_col,:)]);

                    % overwrite at i_row
                    A(i_row,:) = combA;
                    A(:,i_row) = combA';
                    A_nocov(i_row,:) = combAnc;
                    A_nocov(:,i_row) = combAnc';
                    % consistency across time
                    Act(i_row,:) = combAct;
                    Act(:,i_row) = combAct';
                    Act_nocov(i_row,:) = combAct_nc;
                    Act_nocov(:,i_row) = combAct_nc';

                    % overwrite
                    Ad(i_row,:) = combAd;
                    Ad(:,i_row) = combAd';

                    % pop i_col
                    Amask = true(1,length(A));
                    Amask(i_col) = false;
                    A = A(Amask,Amask);
                    Ad = Ad(Amask,Amask);
                    A_nocov = A_nocov(Amask,Amask);
                    Act = Act(Amask,Amask);
                    Act_nocov = Act_nocov(Amask,Amask);
                    
                    % slow
%                     A(i_col,:) = [];
%                     A(:,i_col) = [];
%                     Ad(i_col,:) = [];
%                     Ad(:,i_col) = [];
                    
                    E(i_col,:) = [];
                    Es(i_col) = [];
                    Es2(i_col) = [];
                    E_sphere(i_col,:) = [];
                    Es_sphere(i_col) = [];


                    % Combine
                    fprintf('\t[%.2f mm] subject %i bchan %i, subject %i bchan %i\tA dim: %i\t%i\ttravel mm: %.2f\n',...
                        ad,sidint1,b1,sidint2,b2,length(A),is_same,to_sphere_mm);

                end
                
                
                n_loc_before = length(A);
                
                % -------------------- electrode exclude ----------------
                %rem = [5; 29];
%                 rem = [];
%                 
%                 A(rem,:) = [];
%                 A(:,rem) = [];
%                 A_nocov(rem,:) = [];
%                 A_nocov(:,rem) = [];
%                 n_loc_after = length(A);
%                 Ad(rem,:) = [];
%                 Ad(:,rem) = [];
%                 E(rem,:) = [];
%                 E_sphere(rem,:) = [];
%                 Es(rem,:) = [];
%                 Es2(rem,:) = [];
%                 Es_sphere(rem,:) = [];
                

                save(fn_ca2,'A','Ad','A_nocov','Act','Act_nocov','E','E_sphere','ss_thresh_mm','ss_samesid','Es','Es2','Es_sphere');
                fprintf('\tDone.\n')

            end

            % Plot
            fprintf('[*] Plotting\n')
            system('mkdir figures');
            system('mkdir figures/cluster2');

    %         % Full adjacency matrix
    %         h = figure;
    %         imagesc(A_sav);
    %         xlabel('Bipolar Electrode')
    %         ylabel('Bipolar Electrode')
    %         colorbar;
    %         colormap(inferno(100))
    %         print(h,'figures/cluster2/A_full','-dpng')
    %         close(h);
    %         
    %         % Reduced matrix
    %         h = figure;
    %         imagesc(A);
    %         xlabel('Bipolar Electrode')
    %         ylabel('Bipolar Electrode')
    %         colorbar;
    %         colormap(inferno(100))
    %         print(h,sprintf('figures/cluster2/A_reduced-%i',ss_thresh_mm),'-dpng')
    %         close(h);

            % Compute spatial variances of nodes
            [n_Es,~] = size(Es);
            Es_var = zeros(n_Es,4);
            Es2_var = zeros(n_Es,4);
            for iv = 1:n_Es
                Es_var(iv,1:3) = var(Es{iv}(:,3:5));
                Es_var(iv,4) = (det(cov(Es{iv}(:,3:5))))^(1/3);
                
                cond = Es2{iv}(:,6) ~= 0;
                Es2_var(iv,1:3) = var(Es{iv}(cond,3:5));
                Es2_var(iv,4) = (det(cov(Es{iv}(cond,3:5))))^(1/3);
                % hist(sqrt(Es_var(:,4)));
            end
    
            
            fprintf('[!] Adjacency missing fraction: %.2f%%\n',100*sum(sum(isnan(A)))/numel(A))
            
            fn_ca3 = sprintf('%s/fig_cluster2_fill_%i.mat',dir_cacheLp,iM_local);
            
            
            if (exist(fn_ca3,'file'))
                Ca3 = load(fn_ca3);
                A2 = Ca3.A2;
                k_fill = Ca3.k_fill;
                fprintf('\tFound cached missing data fill using %i-NN interpolation.\n',k_fill);
            else
                A2 = A;
                % Fill missing data
                k_fill = 4;
                fprintf('[*] Filling in missing data using %i-NN interpolation..\n',k_fill);
                [n_A,~] = size(A2);
                c = 0;
                for ii = 1:(n_A-1)
                    fprintf('[*] Progress: %i of %i (%.2f %%)\n',c,(n_A-1),100*c/(n_A-1));
                    c = c + 1;
                    parfor jj = (ii+1):n_A
                        %c = c + 1;
                        %c = 0;
                        if (isnan(A(ii,jj)))
                            % compute distances
                            coord_ii = E(ii,3:5);
                            coord_jj = E(jj,3:5);

                            d_ii = sqrt(sum((coord_ii - E(:,3:5)).^2,2));
                            d_jj = sqrt(sum((coord_jj - E(:,3:5)).^2,2));

                            d_ii(ii) = NaN;
                            d_jj(jj) = NaN;

    %                         % ii goes first unless jj has closer electrode
    %                         if ( min(d_jj) < min(d_i) )
    %                             i_1 = jj;
    %                             i_2 = ii;
    %                         else
    %                             i_1 = ii;
    %                             i_2 = jj;
    %                         end

                            k_fill2 = length(d_ii);
                            a_nn = zeros(k_fill2,1);
                            d_nn = zeros(k_fill2,1);
                            for kk = 1:k_fill2
                                % Pick nearest neighbors without replacement
                                best_d = Inf;
                                best_i = NaN;
                                best_j = NaN;
                                
                                maxdnn = max(d_nn);
                                % Seek best neighbor
                                for iii = 1:(n_A-1)
                                    for jjj = (iii+1):n_A
                                        if ( (iii ~= ii) && (jjj ~= jj) )
                                            curr_d = d_ii(iii) + d_jj(jjj);
                                            if ((curr_d < best_d) && (curr_d> maxdnn))
                                                best_d = curr_d;
                                                best_i = iii;
                                                best_j = jjj;
                                            end
                                        end
                                    end
                                end

                                % No replacement
    %                             d_ii(best_i) = NaN;
    %                             d_jj(best_j) = NaN;

                                a_nn(kk) = A(best_i,best_j);
                                d_nn(kk) = best_d;
                                
                                if (sum(~isnan(a_nn)) >= k_fill)
                                    break;
                                end

                                if (all(isnan(a_nn)))
                                    fprintf(2,'[WARN] All neighbors have missing data, not enough neighbors.\n');
                                end
                            end


                            % Fill in missing data
                            a_nn(isnan(a_nn)) = []; % pop nans
                            a = nanmean(a_nn(1:k_fill));
                            if (isnan(a))
                                fprintf(2,'[WARN] Increase k_fill. Value to fill is still nan\n')
                            end
                            A2(ii,jj) = a;
                            %A2(jj,ii) = a;
                        end
                    end
                end
                
                % fill in A2
                for ii = 1:(n_A-1)
                    for jj = (ii+1):n_A
                        A2(jj,ii) = A2(ii,jj);
                    end
                end
                
                fprintf('\tDone.\n');
                
                save(fn_ca3,'A2','k_fill');
            
            end
            %A2(isnan(A2)) = 0;
            
            % clustering coefficient calculation
            ci = clusteringcoef(A2);
            c_WS = nanmean(ci);
            fprintf('[!] Watts and Strogatz C_actual: %.6f\n',c_WS);
            
            [n_A2,~] = size(A2);
            G = graph(A2 > 0);
            A2_pathl = nan(size(A2));
            A2_paths = cell(size(A2));
            pathl = nan(nchoosek(n_A2,2),1);
            ka2 = 1;
            for ia2 = 1:(n_A2-1)
                for ja2 = (ia2+1):n_A2
                    path = shortestpath(G,ia2,ja2);
                    spl = length(path);
                    A2_pathl(ia2,ja2) = spl;
                    A2_pathl(ja2,ia2) = spl;
                    A2_paths{ia2,ja2} = path;
                    A2_paths{ja2,ia2} = path;
                    pathl(ka2) = spl;
                    ka2 = ka2 + 1;
                end
            end
            
            l_WS = nanmean(pathl(:));
            fprintf('[!] Watts and Strogatz L_actual: %.6f\n',l_WS);
            
            %return
            
            % Plot 300 by 300 adjacency
            h = figure('visible','off','Position',[0 0 1920 1080]);
            imagesc(A2);
            colormap(corrcmap(n_cmap));
            set(gca,'TickDir','out')
            cb = colorbar;
            xticks(1:n_A2);
            yticks(1:n_A2);
            set(cb,'TickLength',0)
            print(h,sprintf('./figures/cluster2/Adj_reduced-%i_%s',n_A2,metricsp{iM_local}),'-depsc')
            print(h,sprintf('./figures/cluster2/Adj_reduced-%i_%s',n_A2,metricsp{iM_local}),'-dpng','-r300')
            close(h);
            
            % Plot filled in values
            h = figure('visible','off');
            imagesc(isnan(A));
            colormap(corrcmap(n_cmap));
            set(gca,'TickDir','out')
            cb = colorbar;
            xticks(1:n_A2);
            yticks(1:n_A2);
            set(cb,'TickLength',0)
            frac_fill = sum(sum(isnan(A)))/numel(A);
            frac_xsub = sum(sum(A_xsub))/numel(A_xsub);
            print(h,sprintf('./figures/cluster2/Adj_reduced-%i_fFill-%i_fXsub-%i_%s',...
                n_A2,round(1000*frac_fill),round(1000*frac_xsub),metricsp{iM_local}),'-depsc')
            print(h,sprintf('./figures/cluster2/Adj_reduced-%i_fFill-%i_fXsub-%i_%s',...
                n_A2,round(1000*frac_fill),round(1000*frac_xsub),metricsp{iM_local}),'-dpng','-r300')
            close(h);


            % kmeans
            %K = 8;
            % ---------------------------------------------------------------------------------------------------------------------------------
            K = Kmaster(km); 
            Ssumd = zeros(length(K),1);
            Res = cell(length(K),4);
            A2km = A2;
            A2km(isnan(A2km)) = 0; % set diagonals to zeros;
            for ki = 1:length(K)
                k = K(ki);
                stream = RandStream('mlfg6331_64');  % Random number stream
                options = statset('UseParallel',1,'UseSubstreams',1,'Streams',stream);
                [idx,C,sumd,D] = kmeans(A2km,k,'Options',options,'MaxIter',10000,'Display','final','Replicates',18000);
                Res{ki,1} = idx;
                Res{ki,2} = C;
                Res{ki,3} = sumd;
                Res{ki,4} = D;

                % cross-val metric
                if (length(K) > 1)
%                     cv = 0;
%                     for j = 1:k
%                         cv = cv + var(D(idx==j,j));
%                     end
                    cv = sum(sumd); 
                    Ssumd(ki) = cv; %sum(sumd);
                end
            end


            % tsne

            dimout = 2;
            options = statset('MaxIter',10000,'TolFun',1e-14);
            Y = tsne(A2km,'Options',options,'Algorithm','exact','Standardize',true,'NumDimensions',dimout,'Verbose',1);
            h = figure('visible','off');
            if (dimout == 2)
                %cck = jet(k);
                cc_p = corrcmap(k+2);
                cck = cc_p(2:(end-1),:);
                %cck = corrcmap(k);
                for i = 1:length(Y(:,1))
                    plot(Y(i,1),Y(i,2),'.','MarkerSize',15,'Color',cck(idx(i),:));
                    hold on;
                end
            else
                plot3(Y(:,1),Y(:,2),Y(:,3),'black.');
            end
            xlabel('t-SNE first dimension');
            ylabel('t-SNE second dimension');
            set(gca,'Box','off')
            set(gca,'TickDir','out')
            print(h,sprintf('figures/cluster2/A_reduced-%i_kmeans-%i_d-%i_tsne-%i_%s',...
                ss_thresh_mm,k,round(sum(sumd)*1000),dimout,metricsp{iM_local}),'-dpng')
            close(h);



            % Plot brain ======================================================
            % Plot on brain
            h = figure('visible','off');
            set(h,'Position',[0 0 1920 1080])

            SMOOTH_FACE = true;
            BRAIN_MESH_ALPHA = 0.1;

            p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
            set(p,'FaceVertexCdata',0.7*[1 1 1]);
            hold on;
            ax = gca;
            ax.Clipping = 'off';
            hold all;

            % colormap params
            
            cc_p = corrcmap(K+2);
            cc = cc_p(2:(end-1),:);
            %cc = corrcmap(K); %jet(K);
            cc_kmeans = cc; %corrcmap(K); %jet(K);

            col_var_min = min(A(:));
            col_var_max = max(A(:));


            if (trig_plot_arc)
                Cb = load('brainexport/red_6all_fsaverage.mat');
            end
            
            % Plot lines
            [n_E,~] = size(E); 
            count = 1;
            for i = 1:(n_E-1)
                for j = (i+1):n_E
                    a = A(i,j);
                    if (~isnan(a))
                        
                        if (trig_plot_arc)
                            % kmeans membership color
                            col_tube = cc_kmeans(idx(i),:);
                            
                            path = Cb.Paths{count};
                            path = downsample(path,1);
                            [n_path,~] = size(path);
                            
                            % plot
                            if (idx(i) == idx(j))
                                % within-cluster
                                tube_radius1 = tube_radius*2;
                                [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius1,ones(1,n_path),tube_n_edges,[0 0 1]);
                                surf(X,Y,Z,'FaceColor',col_tube,'EdgeColor','none'); hold on;
                            else
                                % across-cluster
                                tube_radius2 = tube_radius;
                                col_tube2 = cc_kmeans(idx(j),:);
                                col_tube = mean([col_tube; col_tube2]);
                                [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius2,ones(1,n_path),tube_n_edges,[0 0 1]);
                                surf(X,Y,Z,'FaceColor',col_tube,'EdgeColor','none'); hold on;
                                colorspec = sin((1:n_path)*(0.3 * 2*pi)); %sin((1:n_path)*); %ones(1,n_path);
                            end
                            
                            %return
                        else
                        
                            % Significant lines
                            if (idx(i) == idx(j))
                                % line connects within cluster
                                plot3([E(i,3),E(j,3)],[E(i,4),E(j,4)],[E(i,5),E(j,5)],'-','Color',cc_kmeans(idx(i),:));
                                hold on;
                            else
                                % line connects different clusters
                                plot3([E(i,3),E(j,3)],[E(i,4),E(j,4)],[E(i,5),E(j,5)],'-','Color',cc_kmeans(idx(i),:));
                                hold on;
                                plot3([E(i,3),E(j,3)],[E(i,4),E(j,4)],[E(i,5),E(j,5)],':','Color',cc_kmeans(idx(j),:));
                                hold on;
                            end
                        end
                    end
                    count = count + 1;
                end
            end

            % Plot electrodes
            for i = 1:length(E)
                plot3(E(i,3),E(i,4),E(i,5),'.','MarkerSize',15,'Color',cc_kmeans(idx(i),:));
                hold all;
            end

            axis tight;
            daspect([1 1 1]);
            

            view(90,0);
            axis off;
            print(h,sprintf('figures/cluster2/A_reduced-%i_kmeans-%i_d-%i_tsne-%i_lat_%s',...
                ss_thresh_mm,k,round(sum(sumd)*1000),dimout,metricsp{iM_local}),'-dpng')

%             view(90,-90);
%             print(h,sprintf('figures/cluster2/A_reduced-%i_kmeans-%i_d-%i_tsne-%i_inf',ss_thresh_mm,k,round(sum(sumd)*1000),dimout),'-dpng')
% 
%             view(90,90);
%             print(h,sprintf('figures/cluster2/A_reduced-%i_kmeans-%i_d-%i_tsne-%i_sup',ss_thresh_mm,k,round(sum(sumd)*1000),dimout),'-dpng')
% 
%             view(0,0);
%             print(h,sprintf('figures/cluster2/A_reduced-%i_kmeans-%i_d-%i_tsne-%i_pos',ss_thresh_mm,k,round(sum(sumd)*1000),dimout),'-dpng')
% 
%             view(180,0);
%             print(h,sprintf('figures/cluster2/A_reduced-%i_kmeans-%i_d-%i_tsne-%i_ant',ss_thresh_mm,k,round(sum(sumd)*1000),dimout),'-dpng')
% 
%             view(-90,0);
%             print(h,sprintf('figures/cluster2/A_reduced-%i_kmeans-%i_d-%i_tsne-%i_med',ss_thresh_mm,k,round(sum(sumd)*1000),dimout),'-dpng')

            % End plot brain ==================================================
            close(h)
            
            
            fprintf('[!] Number of unique subjects after reduction: %i\n',length(unique(E(:,1))))
            fprintf('[!] Fraction of combinations with same subs: %.2f%%\n',mean(100*ss_samesid(:,1)))

        end
% 
% 
%         % K-means
%         if (trig_kmeans)
%             [n_E,~] = size(E);
%             [n_D,~] = size(D);
% 
%             K = [4:2:33];
%             n_inits = 20;
%             Iidx = cell(1,length(K));
%             Iidx_mindist = cell(1,length(K));
%             Sum_mindist = zeros(1,length(K));
% 
%             parfor i = 1:length(K)
%                 k = K(i);
% 
%                 var_init_prev = Inf;
%                 var_init_final = NaN;
%                 iidx_final = [];
%                 iidx_mindist_final = [];
%                 % init
%                 for j = 1:n_inits
%                     iidx = zeros(n_E,1);
%                     iidx(1:k) = 1:k;
%                     iidx = iidx(randperm(n_E));
%                     inits = find(iidx);
% 
%                     % Extras
%                     iidx_mindist = zeros(n_E,1);
% 
%                     % Coords of each init electrode
%                     coord_init = nan(k,3);
%                     for l = 1:k
%                         e = E(inits(l),:);
%                         e_sidint = e(1);
%                         e_b = e(2);
%                         % Find all interactions for this subject-bipolar set
%     %                     ee = find((D(:,1) == e_sidint) & ((D(:,2)==e_b) | (D(:,3)==e_b)));
%     %                     coords_all = D(ee,1:7);
%                         ee1 = find((D(:,1) == e_sidint) & ((D(:,3)==e_b)));
%                         ee2 = find((D(:,1) == e_sidint) & ((D(:,2)==e_b)));
%                         coords_all = [D(ee1,5:7); D(ee2,8:10)];
%                         % collapse all electrode pairs to single average
%                         if (isempty(coords_all))
%                             coord_init(l,:) = nan(1,3);
%                         else
%                             coord_init(l,:) = mean(coords_all,1);
%                         end
%                         %disp(coords_all)
%                     end
%                     %disp(coord_init)
%                     %return
% 
% 
%                     % first cluster loop
%                     % assign all points to one of k clusters
%                     for k2 = 1:n_E
%                         if (iidx(k2) == 0) % skip clustering inits
%                             % calculate distance to each init cluster
%                             cdist = nan(k,1);
% 
%                             coord_q = nan(1,3);
% 
%                             e = E(k2,:);
%                             e_sidint = e(1);
%                             e_b = e(2);
%                             % Find all interactions for this subject-bipolar set
%         %                     ee = find((D(:,1) == e_sidint) & ((D(:,2)==e_b) | (D(:,3)==e_b)));
%         %                     coords_all = D(ee,1:7);
%                             ee1 = find((D(:,1) == e_sidint) & ((D(:,3)==e_b)));
%                             ee2 = find((D(:,1) == e_sidint) & ((D(:,2)==e_b)));
%                             coords_all = [D(ee1,5:7); D(ee2,8:10)];
%                             % collapse all electrode pairs to single average
%                             if (isempty(coords_all))
%                                 coord_q(1,:) = nan(1,3);
%                             else
%                                 coord_q(1,:) = mean(coords_all,1);
%                             end
% 
%                            % cdist(l) = cluster_dist;
% 
%                             % compute distances
%                             for l = 1:k
%                                 if (isnan(coord_q(1)) && isnan(coord_init(l,1)))
%                                     % Both cluster and query have no interactions
%                                     cdist(l) = 0;
%                                 elseif (~isnan(coord_q(1)) && ~isnan(coord_init(l,1)))
%                                     % Both cluster and query have interactions
%                                     % return euclidean distances between
%                                     % geometric centers of interacting
%                                     % electrode set
%                                     cdist(l) = sqrt(sum((coord_q - coord_init(l,:)).^2));
%                                 else
%                                     % Only one has no interactions
%                                     cdist(l) = Inf;
%                                 end
%                             end
% 
%                             % Assign minimum distance cluster
%                             [mindist,cluster_id] = min(cdist);
%                             iidx(k2) = cluster_id;
%                             iidx_mindist(k2) = mindist;
% 
% 
%                         end
%                     end
% 
%                     % Calculate goodness of init
%                     var_init = 0;
%                     for k3 = 1:k
%                         var_init = var_init + mean(iidx_mindist(iidx==k3));
%                     end
%                     var_init = var_init / k;
%                     %var_init = sum(iidx_mindist);
% 
% 
%                     if (var_init < var_init_prev)
%                         var_init_final = var_init;
%                         var_init_prev = var_init;
% 
%                         iidx_final = iidx;
%                         iidx_mindist_final = iidx_mindist;
% 
%                         fprintf('[*] k = %i init %i of %i, log sum mindist: %.6f\n',k,j,n_inits,log(var_init_final));
%                     end
% 
%                 end % init loop
% 
%                 Iidx{i} = iidx_final;
%                 Iidx_mindist{i} = iidx_mindist_final;
%                 Sum_mindist(i) = var_init_final;
% 
%             end % large K loop
% 
%             % Plot variance against k
%             h = figure;
%             plot(K,Sum_mindist,'black-');
%             xlabel('k')
%             ylabel('Average cluster variance (mm)')
% 
%             [~,mIdx] = min(Sum_mindist);
% 
%             iidx = Iidx{mIdx};
%             iidx_mindist = Iidx_mindist{mIdx};
%             k = K(mIdx);
% 
% 
% 
% 
%         end
% 
% 
% 
%         if (trig_plot_all)
%             h = figure;
%             % plot surface
%             SMOOTH_FACE = true;
%             BRAIN_MESH_ALPHA = 0.1;
%             p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
%             set(p,'FaceVertexCdata',0.7*[1 1 1]);
%             ax = gca;
%             ax.Clipping = 'off';
%             hold all;
% 
%             % colormap params
%             n_cmap = 100;
%             cc = jet(n_cmap);
% 
% 
%             col_var_min = min(D(:,variable_idx)); %0.1;
%             col_var_max = max(D(:,variable_idx)); %1;
%             % Plot lines
%             for i = 1:length(D)
%                 % Plot first electrode
%                 %plot3(D(i,5),D(i,6),D(i,7),'black.')
%                 % plot second electrode
%                 %plot3(D(i,8),D(i,9),D(i,10),'black.')
%                 % plot interaction
%                 % - magnitude -------------------------------------------------
%                 col_var = D(i,variable_idx);
%                 %col_line = cc(round((col_var/col_var_max)*n_cmap),:);
%                 col_line = cc(round(( (col_var - col_var_min)/(col_var_max - col_var_min) ) *(n_cmap-1)+1),:);
% 
% 
%                 if (trig_kmeans)
%                     % find E index for this D
%                     cc_kmeans = jet(k);
% 
%                     d_sidint = D(i,1);
%                     d_e = D(i,2);
%                     e_idx = find( (E(:,1) == d_sidint) & (E(:,2) == d_e) );
% 
%                     % find cluster id for this line
%                     plot3([D(i,5),D(i,8)],[D(i,6),D(i,9)],[D(i,7),D(i,10)],'Color',cc_kmeans(iidx(e_idx),:));
%                 else
%                     plot3([D(i,5),D(i,8)],[D(i,6),D(i,9)],[D(i,7),D(i,10)],'Color',col_line);
%                 end
%         %         if (i == 8000)
%         %             break
%         %         end
%             end
% 
%             %plot electrodes
%             for i = 1:length(E)
%                 if (trig_kmeans)
%                     plot3(E(i,3),E(i,4),E(i,5),'.','Color',cc_kmeans(iidx(i),:));
%                 else
%                     plot3(E(i,3),E(i,4),E(i,5),'black.');
%                 end
%             end
% 
%             axis tight;
%             daspect([1 1 1]);
%             view(90,0);
%             if (SMOOTH_FACE)
%                 p.FaceLighting = 'gouraud';
%             end
%             colorbar
%             if (trig_kmeans)
%                 colormap(cc_kmeans)
%                 caxis([1 k])
%             else
%                 colormap(cc)
%                 caxis([col_var_min col_var_max])
%             end
%         end




        % Show numbers
        fprintf('Number of electrodes: %i\n',length(E))
        fprintf('Number of significant interactions: %i\n',length(D))
    
    end

end