% Figure W3


close all;
clear;
clc;
rng shuffle;

n_perm = 10000;
perm_alpha_sec = 20;
cp_thresh_override = 0.05; 
n_pairs_thresh = 10;% 10; % at least this many electrode pairs to be considered
n_subs_thresh = 2;% 2; % at least this many subjects to be considered
n_subs_ct_thresh = 0;% 2; % significant CTs in region pair must be from at least this many subjects
fig_size_scale = 0.75;
stride_axislabel = 1; % stride for x and y axis labels
%n_resample = 20; % number of times to resample with n_pairs_thresh pairs and n_subs_thresh subjects

% n_pairs_thresh = 10 * 4; % at least this many electrode pairs to be considered
% n_subs_thresh = 2; % at least this many subjects to be considered
% [*] Original ROI coverage: 78.0952 %
% [*] Thresholded ROI coverage: 61.4286 %

fig_fmt = '-dpng';
trig_eps = true;
trig_mag_no_cp = true;
trig_plot_mag_cov = true;
system('mkdir figures');
system('mkdir figures/T14_allatl_resampSage');



[~,host] = system('hostname');

if contains(host,'ubuntu_1604')
    dir_artL = '/nas_share/RawData/data/h5_notch20/art';
    dir_resL = '/nas_share/RawData/data/results/coh_w10';
    dir_corL = '/nas_share/RawData/data/coreg';
    dir_cacheL = './cache';
    subjects_dirL = '/mnt/cuenap_ssd/coregistration';
    dir_h5L = '/nas_share/RawData/data/h5_notch20';
else
%     dir_artL = '/media/klab/internal/data/h5_notch20/art';
%     dir_resL = '/media/klab/internal/data/results/coh_w10';
%     dir_corL = '/media/klab/internal/data/coreg';
    dir_artL = '/media/klab/KLAB101/h5_notch20/art_nosz';
    dir_resL = '/media/klab/KLAB101/results/coh_w10';
    dir_corL = '/media/klab/internal/data/coreg';
    dir_cacheL = './cache';
    subjects_dirL = '/mnt/cuenap_ssd/coregistration';
    dir_h5L = '/media/klab/KLAB101/h5_notch20';
end

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
metrics_suffix = {'0.5-125 Hz','3-8 Hz','8-12 Hz','12-30 Hz','30-100 Hz'};
%metrics = {'pcBroadband'};

for iAgeG = 1:2

    % Get atlas names
    CaAtl = load(sprintf('./cache/xsub_out_all_%i',1));
    AtlNames = CaAtl.C.AtlNames;
    for atl = 2 %1:20

        system(sprintf('mkdir figures/T14_allatl_resampSage/atl%i_%s',atl,AtlNames{atl}));

        for iM = 1 %[1 5] %1:length(metrics) % [1 5] %
            metric = metrics{iM};

            % Load human cache
            Ca_hum = load(sprintf('%s/xsub_out_all_%i_atl%i.mat',dir_cacheL,iM,atl));
            roi_labels = Ca_hum.C.AtlROIs{atl}.RH.struct_names;
            n_rois = Ca_hum.n_rois;
            %AdjAA = cell(n_rois,n_rois);


            % Vectorize matrix
            n_roi_pairs = nchoosek(n_rois,2);
            Vec_ij = zeros(n_roi_pairs,2);
            c = 1;
            for i1 = 1:(n_rois-1)
                for i2 = (i1+1):n_rois
                    Vec_ij(c,1) = i1;
                    Vec_ij(c,2) = i2;
                    c = c + 1;
                end
            end


            % Unknown areas filter
            rois = Ca_hum.rois;
            cond_empty = false(length(rois),1);
            for ce = 1:length(cond_empty)
                elem = rois{ce};
                cond_empty(ce) = (all(isspace(elem)) | isempty(elem));
            end
            cond_contains = (contains(lower(rois),'unknown')) | (contains(lower(rois),'???')) | (contains(lower(rois),'wall')) | (cond_empty);
            known_idx = (~ cond_contains);


            % Permutation parameters <-- important
            %DropoutN = [0 8]; %0:length(Ca_hum.Subjects);
            DropoutN = 0;
            n_perm = 100000;


            % age
            SubAge = [...
                21,...
                14,...
                3,...
                32,...
                46,...
                10,...
                20,...
                39,...
                47,...
                21,...
                26,...
                18,...
                9,...
                32,...
                21,...
                8,...
                10,...
                9,...
                18,...
                15,...
                12,...
                17,...
                9,...
                16,...
                9,...
                17,...
                3,...
                11,...
                19,...
                44,...
                31,...
                42,...
                18,...
                31,...
                16,...
                18,...
                18,...
                19,...
                11,...
                17,...
                31,...
                16,...
                13,...
                14,...
                10,...
                7,...
                10,...
                17,...
            ];

            % make age lookup
            n_sub = length(Ca_hum.Subjects);
            [~,sIdx_SubAge] = sort(SubAge);
            i_group = false(1,n_sub);
            i_group(sIdx_SubAge(1:round(n_sub/2))) = true;
            fprintf('[*] Total subjects: %i\n',n_sub);
            fprintf('\tCutoff median age: %.2f\n',median(SubAge));
            fprintf('\tGroup 1 - %i, mean age: %.3f, std: %.4f\n',sum(i_group),mean(SubAge(i_group)),std(SubAge(i_group)));
            fprintf('\tGroup 2 - %i, mean age: %.3f, std: %.4f\n',sum(~i_group),mean(SubAge(~i_group)),std(SubAge(~i_group)));
            [pval,~,~] = ranksum(SubAge(i_group),SubAge(~i_group));
            fprintf('\tAge ranksum p: %.8d\n',pval);


            % Drop out subjects
            i_nSub = 1;
            AdjMag_Hi = nan(n_roi_pairs,length(DropoutN));
            AdjMag_Lo = nan(n_roi_pairs,length(DropoutN));
            for nSub = DropoutN

                % Set number of resampling permutations
                if (nSub == 0)
                    n_perm_ac = 1;
                else
                    n_perm_ac = n_perm;
                end

                fprintf('[n_perm: %i] atl: %i\tmetric: %i\t#drop: %i\tof %i\n',...
                    n_perm_ac,atl,iM,nSub,max(DropoutN));

                AdjMag3 = nan(n_rois,n_rois,n_perm_ac);
                AdjMag3_1 = nan(n_rois,n_rois,n_perm_ac);
                AdjMag3_2 = nan(n_rois,n_rois,n_perm_ac);
                tic;

                for ip = 1:n_perm_ac

                    % Pick nSub subjects to drop out by randperm
                    sIdx = 1:n_sub;
                    sIdx = sIdx(randperm(n_sub));
                    sub_ex = sIdx(1:nSub);

                    % Rebuild AdjMag excluding these subjects
                    Adj = nan(n_rois,n_rois);
                    AdjMag = nan(n_rois,n_rois);
                    AdjMag_1 = nan(n_rois,n_rois);
                    AdjMag_2 = nan(n_rois,n_rois);
                    AdjMag_1std = nan(n_rois,n_rois);
                    AdjMag_2std = nan(n_rois,n_rois);
                    AdjMag4cl = nan(n_rois,n_rois);
                    AdjNpairs = nan(n_rois,n_rois);
                    AdjNpairs_sig = nan(n_rois,n_rois);
                    AdjNusubs = nan(n_rois,n_rois);
                    AdjNusubs_sig = nan(n_rois,n_rois);
                    AdjCP = nan(n_rois,n_rois);
                    AdjMagVar = nan(n_rois,n_rois);
                    AdjMagReS = nan(n_rois,n_rois);
                    AdjMagL = cell(n_rois,n_rois);
                    N_bchan = nan(n_rois,n_rois);
                    cp_thresh = cp_thresh_override;
                    Dmat = Inf(n_rois,n_rois); % set to inf to avoid removing from every instance
                    DistsAtl = [];

                    elec_count = 0;
                    V_mag = [];
                    V_age = [];
                    for i1 = 1:(n_rois)
                        for i2 = 1:n_rois %(i1):n_rois

                            cond_known = (known_idx(i1) && known_idx(i2));

                            if (cond_known)

                                AA_sub_all = Ca_hum.AdjAtl_sid{i1,i2};
                                AA_all = Ca_hum.AdjAtl{i1,i2};
                                AA_dist_all = Ca_hum.adjct_dist{i1,i2};

                                % =============================================
                                % Filter sub
    %                             i_trim = false(size(AA_sub_all));
    %                             for i3 = 1:length(i_trim)
    %                                 if (mod(elec_count,2) == 0)
    %                                     i_trim(i3) = true;
    %                                 else
    %                                     i_trim(i3) = false;
    %                                 end
    %                                 elec_count = elec_count + 1;
    %                                 %i_trim = i_trim | (AA_sub == sub_ex(i3));
    %                             end
                                i_trim = i_group(AA_sub_all);

                                AA_age = SubAge(AA_sub_all);
                                % =============================================

                                % apply filter
                                AA_sub = AA_sub_all(i_trim);
                                AA = AA_all(i_trim);
                                AA_dist = AA_dist_all(i_trim);

                                % apply complement filter
                                AA_sub_c = AA_sub_all(~i_trim);
                                AA_c = AA_all(~i_trim);
                                AA_dist_c = AA_dist_all(~i_trim);


                                n_pairs = length(AA_sub_all);
                                n_subs = length(unique(AA_sub_all));
                                n_subs_ct = length(unique(AA_sub_all(AA_all ~= 0)));
                                % ROI pair coverage condition (grey)
                                if ( (~ isempty(AA_all))  && (n_pairs >= n_pairs_thresh) && (n_subs >= n_subs_thresh))


                                    N_bchan(i1,i2) = length(AA_all);
                                    frac_cp = sum(AA_all ~= 0)/length(AA_all);
                                    AdjCP(i1,i2) = frac_cp;
                                    % ROI pair significance condition (white)
                                    if ( (frac_cp > cp_thresh) && (n_subs_ct >= n_subs_ct_thresh) )

                                        % =========================================
                                        % save
                                        AA_all2 = AA_all; %((AA_all~=0));
                                        AA_age2 = AA_age; %((AA_all~=0));
                                        V_mag = [V_mag, AA_all2'];
                                        V_age = [V_age, AA_age2];
                                        % =========================================

                                        %return;
                                        Adj(i1,i2) = 1;
                                        AdjMag(i1,i2) = mean(AA_all(AA_all ~= 0));

                                        % Separate calculations for pair resamp
                                        AdjMag_1(i1,i2) = mean(AA(AA ~= 0));
                                        AdjMag_2(i1,i2) = mean(AA_c(AA_c ~= 0));

                                        AdjMag_1std(i1,i2) = std(AA(AA ~= 0));
                                        AdjMag_2std(i1,i2) = std(AA_c(AA_c ~= 0));

                                        AdjNpairs(i1,i2) = n_pairs;
                                        AdjNpairs_sig(i1,i2) = length(AA(AA ~= 0));
                                        AdjNusubs(i1,i2) = n_subs;
                                        AdjNusubs_sig(i1,i2) = n_subs_ct;

                                        AdjMagVar(i1,i2) = std(AA_all(AA_all ~= 0));
                                        AdjMagL{i1,i2} = AA_all(AA_all~=0);
                                        DistsAtl = [DistsAtl; [mean(AA_dist_all(AA_all ~= 0)), mean(AA_all(AA_all ~= 0))]];
                                    else
                                        Adj(i1,i2) = 0;
                                        AdjMag(i1,i2) = 0;
                                        AdjMag_1(i1,i2) = 0;
                                        AdjMag_2(i1,i2) = 0;
                                        AdjMag_1std(i1,i2) = 0;
                                        AdjMag_2std(i1,i2) = 0;
                                        AdjMagVar(i1,i2) = 0;
                                        AdjNpairs(i1,i2) = 0;
                                        AdjNpairs_sig(i1,i2) = 0;
                                        AdjNusubs(i1,i2) = 0;
                                        AdjNusubs_sig(i1,i2) = 0;
                                        AdjMagL{i1,i2} = [];
                                    end
                                end
                                % AdjMag for clustering
                                AdjMag4cl(i1,i2) = mean(AA(AA~=0));
                            end
                        end
                    end

                    fprintf('[*] Correlation age vs. coherence\n');
                    [r,pval] = corr(V_mag',V_age','Type','Spearman');
                    fprintf('\tr=%.6f, p=%.6d, n=%i\n',r,pval,length(V_mag));

                    AdjMag3(:,:,ip) = AdjMag;
                    AdjMag3_1(:,:,ip) = AdjMag_1;
                    AdjMag3_2(:,:,ip) = AdjMag_2;
                end


                % 95%CI for AdjMag3
                CI_alpha = 0.05;
                AdjMag_hi = nan(n_rois,n_rois);
                AdjMag_lo = nan(n_rois,n_rois);
                for i1 = 1:n_rois
                    for i2 = (i1):n_rois
                        mags = squeeze(AdjMag3(i1,i2,:));
                        s_mags = sort(mags,'descend');
                        sIdx = (~isnan(s_mags)) & (s_mags ~= 0);
                        s_mags = s_mags(sIdx);

                        n_mags = length(s_mags);
                        idx_hi = round((CI_alpha/2) * n_mags);
                        idx_lo = round((1-(CI_alpha/2)) * n_mags);

                        % Exception for no resampling (nSub = 0)
                        if (idx_hi == 0)
                            idx_hi = 1;
                        end

                        if (isempty(s_mags))
                            hi = NaN;
                            lo = NaN;
                        else
                            hi = s_mags(idx_hi);
                            lo = s_mags(idx_lo);
                        end

                        AdjMag_hi(i1,i2) = hi;
                        AdjMag_lo(i1,i2) = lo;

    %                     if (i1 == 4) && (i2 == 10)
    %                         return
    %                     end
                    end
                end

                % Vectorize CI
                AdjMag_hi_v = nan(n_roi_pairs,1);
                AdjMag_lo_v = nan(n_roi_pairs,1);
                for c = 1:n_roi_pairs
                    i1 = Vec_ij(c,1);
                    i2 = Vec_ij(c,2);
                    AdjMag_hi_v(c) = AdjMag_hi(i1,i2);
                    AdjMag_lo_v(c) = AdjMag_lo(i1,i2);
                end

                % Store
                AdjMag_Hi(:,i_nSub) = AdjMag_hi_v;
                AdjMag_Lo(:,i_nSub) = AdjMag_lo_v;

                t_sing = toc;
                fprintf('per atl-metric ETA: %.1f mins\n',(t_sing * (length(DropoutN) - i_nSub))/(60))

    %             if (i_nSub == 3)
    %                 return
    %             end

                i_nSub = i_nSub + 1;
            end

            save(sprintf('./cache/figure_t14_allatl_resampSage_im-%i_ageg-%i',iM,iAgeG));



            %% Plot

            % --- PLOT  --------------------------------------------------
            colordef_adj;
            color_nocov = COLOR_ADJ_NOCOV;
            color_distt = color_nocov;
            color_isnan = color_nocov;
            color_not_sig0 = COLOR_ADJ_NOSIG; %0.6*[1 1 1]; %[0 0.5 0]; %0.4*[0.5 1 0.5];
            fontsz = 10;
            fsize = fontsz;

            % Calculate rows of nans
        %     Adj2 = Adj;
        %     ind_isnan = false(length(Adj2),1);
        %     for i = 1:length(ind_isnan)
        %         ind_isnan(i) = (sum(isnan(Adj2(i,:))) == length(ind_isnan));
        %     end

            % overall 36-area isnan storage
            ind_isnan_master = false(n_rois,1);
            for iii = 3%1:3

                %##################################################################
                % Make figure
                %##################################################################
                trig_ct = false;

                for iii2 = 1:6 % 1:2

                    %h = figure('visible','on');
                    h = figure('visible','off');
                    %set(h,'Position',round(fig_size_scale*[0 0 1*1080 1*1080]));
                    set(h,'PaperUnits','Inches');
                    set(h,'PaperPosition',[0 0 8.5 6.8]);

                    if (iii == 1)
                        if (trig_ct)
                            Adj_plt = AdjCT;
                        else
                            if (trig_mag_no_cp)
                                if (iAgeG == 1)
                                    Adj_plt = AdjMag_1; %AdjMag;
                                elseif (iAgeG == 2)
                                    Adj_plt = AdjMag_2; %AdjMag;
                                end
                            else
                                Adj_plt = AdjCP;
                            end
                        end
                    elseif (iii == 2)
                        Adj_plt = Adj;
                    elseif (iii == 3)
                        if (trig_mag_no_cp)
                            if (iii2 == 1)
                                if (trig_ct)
                                    Adj_plt = AdjCT;
                                else
                                    if (iAgeG == 1)
                                        Adj_plt = AdjMag_1; %AdjMag;
                                    elseif (iAgeG == 2)
                                        Adj_plt = AdjMag_2; %AdjMag;
                                    end
                                    %Adj_plt = AdjMag_1-AdjMag_2; %AdjMag;
                                end
                            elseif (iii2 == 2)
                                if (trig_ct)
                                    Adj_plt = AdjCTVar;
                                else
                                    if (iAgeG == 1)
                                        Adj_plt = AdjMag_1std; %AdjMag;
                                    elseif (iAgeG == 2)
                                        Adj_plt = AdjMag_2std; %AdjMag;
                                    end
                                    %Adj_plt = AdjMagVar;
                                end
                            elseif (iii2 == 3)
                                Adj_plt = AdjNpairs_sig;
                            elseif (iii2 == 4)
                                Adj_plt = AdjNpairs;
                            elseif (iii2 == 5)
                                Adj_plt = AdjNusubs_sig;
                            elseif (iii2 == 6)
                                Adj_plt = AdjNusubs;
                            end
                        else
                            Adj_plt = AdjCP;
                        end
                        %Adj_plt_null = AdjMagNull;
                        %Adj_plt(Adj == 0) = 0;
                        color_not_sig = color_not_sig0;
                    end
                    Adj_plt_var = AdjMagVar;
                    Adj_plt_npairs = AdjNpairs;
                    Adj_plt_npairs_sig = AdjNpairs_sig;
                    Adj_plt_nusubs = AdjNusubs;
                    Adj_plt_nusubs_sig = AdjNusubs_sig;
                    Adj_plt2_cl = AdjMag4cl;

                    rois = Ca_hum.rois;
                    Dmat_plt = Dmat;
                    rois_plt = rois;

                    % checkpoint:
                    % Adjacency matrix Adj_plt: 36 by 36
                    % Adjacency rois_plt: 36 by 1
                    % Distance matrix Dmat_plt: 


                    % Filter out unknown from adjacency to plot



                    cond_empty = false(length(rois),1);
                    for ce = 1:length(cond_empty)
                        elem = rois{ce};
                        cond_empty(ce) = (all(isspace(elem)) | isempty(elem));
                    end

                    cond_contains = (contains(lower(rois),'unknown')) | (contains(lower(rois),'???')) | (contains(lower(rois),'wall')) | (cond_empty);
                    known_idx = (~ cond_contains);
                    ind_isnan_master(~known_idx) = true;
                    Adj_plt = Adj_plt(known_idx,known_idx);
                    %Adj_plt_null = Adj_plt_null(known_idx,known_idx);
                    Adj_plt_var = Adj_plt_var(known_idx,known_idx);
                    Adj_plt2_cl = Adj_plt2_cl(known_idx,known_idx);
                    Adj_plt_npairs = Adj_plt_npairs(known_idx,known_idx);
                    Adj_plt_npairs_sig = Adj_plt_npairs_sig(known_idx,known_idx);
                    Adj_plt_nusubs = Adj_plt_nusubs(known_idx,known_idx);
                    Adj_plt_nusubs_sig = Adj_plt_nusubs_sig(known_idx,known_idx);
                    rois_plt = rois_plt(known_idx);
                    Dmat_plt = Dmat_plt(known_idx,known_idx);

                    % mag and CT
                    AM = AdjMag(known_idx,known_idx);
                    %ACT = AdjCT(known_idx,known_idx);

                    % Clean ROI labels for human readability
                    for i = 1:length(rois_plt)
                        % Single atlas exceptions
                        if (startsWith(rois_plt{i},'L_'))
                            rpt = strsplit(rois_plt{i},'L_');
                            rpt = strsplit(rpt{end},'_ROI');
                            rois_plt{i} = rpt{1};
                        end
                        if (endsWith(rois_plt{i},'_M132'))
                            rpt = strsplit(rois_plt{i},'_M132');
                            rois_plt{i} = rpt{1};
                        end
                        if (startsWith(rois_plt{i},'FVE_all_'))
                            rpt = strsplit(rois_plt{i},'FVE_all_');
                            rois_plt{i} = rpt{end};
                        elseif (startsWith(rois_plt{i},'FVE_'))
                            rpt = strsplit(rois_plt{i},'FVE_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'LVE00_'))
                            rpt = strsplit(rois_plt{i},'LVE00_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'PHT00_'))
                            rpt = strsplit(rois_plt{i},'PHT00_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'BoninBailey_'))
                            rpt = strsplit(rois_plt{i},'BoninBailey_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'FerryEtAl_00'))
                            rpt = strsplit(rois_plt{i},'FerryEtAl_00');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'UD86_'))
                            rpt = strsplit(rois_plt{i},'UD86_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'PGR91_'))
                            rpt = strsplit(rois_plt{i},'PGR91_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'LyonKaas02_'))
                            rpt = strsplit(rois_plt{i},'LyonKaas02_');
                            rois_plt{i} = rpt{end};
                        end
                        if (startsWith(rois_plt{i},'BRL87_'))
                            rpt = strsplit(rois_plt{i},'BRL87_');
                            rois_plt{i} = rpt{end};
                        end

                        % Remove underscores
                        %rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '_',' ');
                        %rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '.',' ');

                        if (atl == 2)
                            rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '_',' ');
                            rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '.',' ');
                            rois_plt{i} = convertRoiDK(rois_plt{i});
                        end
                    end



                    % Distance threshold: add nans to Adj_plt
                    %Adj_plt(Dmat_plt <= dist_thresh) = nanmean(Adj_plt(:));
                    Adj_plt2 = Adj_plt;
                    Adj_plt2_var = Adj_plt_var;
                    Adj_plt2_npairs = Adj_plt_npairs;
                    Adj_plt2_npairs_sig = Adj_plt_npairs_sig;
                    Adj_plt2_nusubs = Adj_plt_nusubs;
                    Adj_plt2_nusubs_sig = Adj_plt_nusubs_sig;
                    dist_thresh = Ca_hum.dist_thresh;
                    Adj_plt(Dmat_plt <= dist_thresh) = nan;
                    %Adj_plt_null(Dmat_plt <= dist_thresh) = nan;
                    Adj_plt2_cl(Dmat_plt <= dist_thresh) = nan;

                    % mag and Ct
                    AM(Dmat_plt <= dist_thresh) = nan;
                    %ACT(Dmat_plt <= dist_thresh) = nan;

                    %return
                    % Filter out nans nodes
                    if (iii2 == 1)
                        % only compute cov_idx for mean, not std
                        cov_idx = ~ all(isnan(Adj_plt));
                        %cov_idx = ~ all(isnan(Adj_plt) | (Adj_plt==0));
                    end

                    Adj_plt = Adj_plt(cov_idx,cov_idx);
                    %Adj_plt_null = Adj_plt_null(cov_idx,cov_idx);
                    Adj_plt2 = Adj_plt2(cov_idx,cov_idx);
                    Adj_plt2_var = Adj_plt2_var(cov_idx,cov_idx);
                    Adj_plt2_npairs = Adj_plt2_npairs(cov_idx,cov_idx);
                    Adj_plt2_npairs_sig = Adj_plt2_npairs_sig(cov_idx,cov_idx);
                    Adj_plt2_nusubs = Adj_plt2_nusubs(cov_idx,cov_idx);
                    Adj_plt2_nusubs_sig = Adj_plt2_nusubs_sig(cov_idx,cov_idx);
                    Adj_plt2_cl = Adj_plt2_cl(cov_idx,cov_idx);
                    rois_plt = rois_plt(cov_idx);
                    Dmat_plt = Dmat_plt(cov_idx,cov_idx);
                    ind_isnan_master(~cov_idx) = true;

                    % mag and ct
    %                 AM = AM(cov_idx,cov_idx); %(Dmat_plt <= dist_thresh) = nan;
    %                 ACT = ACT(cov_idx,cov_idx); %(Dmat_plt <= dist_thresh) = nan;
    %                 vAM = nan(nchoosek(length(AM),2),1);
    %                 vACT = nan(nchoosek(length(ACT),2),1);
    %                 cvam = 1;
    %                 for ivam = 1:(length(AM)-1)
    %                     for ivam2 = (ivam+1):length(AM)
    %                         vAM(cvam) = AM(ivam,ivam2);
    %                         vACT(cvam) = ACT(ivam,ivam2);
    %                         cvam = cvam + 1;
    %                     end
    %                 end
    %                 pIdx = ((~isnan(vAM)) & (vAM ~= 0)) & ((~isnan(vACT)) & (vACT ~= 0));
    %                 vAM = vAM(pIdx);
    %                 vACT = vACT(pIdx);
    %                 [r,p] = corr(vAM,vACT);
    %                 fprintf('[*] corr Coh vs. CT: %.4f, p=%.4d, n=%i\n',r,p,length(vAM));
                    %[*] corr Coh vs. CT: 0.1103, p=1.2676e-01, n=193

                    % Get number of nodes
                    [n_roi_p,~] = size(Adj_plt);

                    % Main setting for plot
                    v = Adj_plt;
                    [v_n, v_m] = size(v);
                    %map = corrcmap(100);
                    map = COLOR_ADJ_CMAP; %inferno(100);
                    if (iii == 3)
                        minv = min(v(v~=0));
                        maxv = max(v(v~=0));
                    else
                        minv = min(v(:));
                        maxv = max(v(:));
                    end
                    
                    % -----------------------------------------------------
                    % Colorscale max mins for mean coherence
                    % 	[*] minv: 0.137898698449
                    % 	[*] maxv: 0.510548055172
                    % for std coherence
                    % 	[*] minv: 0.001970777856
                    % 	[*] maxv: 0.165467841509

                    % Second age group mean coh
                    % 	[*] minv: 0.143031621973
                    % 	[*] maxv: 0.617667578161
                    % Second age group std coh
                    % 	[*] minv: 0.000529849106
                    % 	[*] maxv: 0.206156606033

                    % Remap min/max colorscale
                    if (abs(minv - 0.143031621973) <= 1e-11)
                        minv = 0.137898698449;
                        %caxis([0.137898698449 maxv]);
                    elseif (abs(minv - 0.137898698449) <= 1e-11)
                        maxv = 0.617667578161;
                        %caxis([minv 0.617667578161]);
                    elseif (abs(minv - 0.001970777856) <= 1e-11)
                        minv = 0.000529849106;
                        maxv = 0.206156606033;
                        %caxis([0.000529849106 0.206156606033]);
                    elseif (abs(minv - 0.000529849106) <= 1e-11)
                        minv = minv;
                        maxv = maxv;
                        %caxis([minv maxv]);
                    end
                    % -----------------------------------------------------
                    
                    
                    ncol = size(map,1);
                    s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
                    Im = ind2rgb(s,map);
                    for i = 1:n_roi_p
                        for j = 1:n_roi_p
                            if (isnan(Adj_plt(i,j)))
                                Im(i,j,:) = color_isnan;
                            end
                            if (Dmat_plt(i,j) <= dist_thresh)
                                Im(i,j,:) = color_distt;
                            end
                            if ((iii == 3) && (Adj_plt(i,j) == 0))
                                Im(i,j,:) = color_not_sig;
                            end
                        end
                    end

                    %Im = Im(~ind_isnan,~ind_isnan,:); % remove rows of nans
                    %rois_plt = rois_plt(~ind_isnan);
                    %Adj_plt2 = Adj_plt(~ind_isnan,~ind_isnan);
                    %rois_plt2 = rois_plt;
                    %Adj_plt2 = Adj_plt;
                    %rois_plt2 = rois_plt;




                    % --- Clustering (output: cluster_i) --------------------------------------
    %adjct_dist_cl = adjct_dist(ind_hum2mac,ind_hum2mac);



                    %if ((iM == 1) && (iii2 == 1))
                    if ((iii2 == 1))
                        % TODO: RERUN xsub_out_stats to enable adjct_dist clustering
                        %adjct_dist = Ca_hum.adjct_dist;
                        %adjct_dist_cl = adjct_dist;
                        %adjct_dist_cl = adjct_dist_cl(~ind_isnan_master,~ind_isnan_master);
                        %adjct_dist = adjct_dist(~ind_isnan)
                        n_rois_cl = length(rois_plt);
                        Y = zeros(1,nchoosek(n_rois_cl,2));
                        yc = 1;
                        for m1 = 1:(n_rois_cl-1)
                            for m2 = (m1+1):n_rois_cl

                                % ========= CLUSTER BY DISTANCE ===========================
                                %ad = adjct_dist_cl{m2,m1};
                %                 if (isempty(ad))
                %                     Y(yc) = Inf;%Inf; 600;
                %                 else
                %                     Y(yc) = mean(ad);
                %                 end
                %                 yc = yc + 1;

                                % ========= CLUSTER BY COHERENCE ==========================
                                %ad = Adj_plt(m2,m1);
                                ad = Adj_plt(m1,m2);
                                if (isnan(ad))
                                    Y(yc) = 0;
                                else
                                    Y(yc) = ad;
                                end
                                yc = yc + 1;


                            end
                        end
                        %X = Y;
                        %Y = pdist(X);
                        roi_dist = squareform(Y);
                        % average   2778.832643576550 mm
                        % centroid  2773.791496333409 mm
                        % complete  2808.894068525726 mm
                        % median    2767.246370358684 mm
                        % single    2836.305543514634 mm
                        % ward      2778.832643576550 mm
                        % weighted  2778.918494966365 mm

                        %

                        % Imputed
                        %Z = linkage(Y,'complete'); %,'centroid'

                        % Marginalized
                        %M = Adj_plt2_cl;
                        M = Adj_plt2;
                        M = (M + M.')/2;

                        fn_cluster_ca = sprintf('./cache/figure_t14_%i_atl%i_%s',iM,atl,AtlNames{atl});
                        if (exist([fn_cluster_ca,'.mat'],'file'))
                            fprintf('[*] Loading existing hier-clustering: %s\n',fn_cluster_ca)
                            Ca_clu = load(fn_cluster_ca);
                            cluster_i = Ca_clu.cluster_i;
                            Z = Ca_clu.Z;
                            roi_dist = roi_dist(cluster_i,cluster_i);
                            clash = nansum(nansum(triu(roi_dist,1) - triu(roi_dist,2)));
                            fprintf('cluster clash: %.12f mm\n',clash)
                        else

                            %Z = linkage(M,'ward',m);
                            %return
                            Z = linkage(M,'ward',@nanDist4clustering);

                            cluster_i = optimalleaforder(Z,Y); % ,'transformation','inverse'
                            roi_dist = roi_dist(cluster_i,cluster_i);
                            %roi_dist(isinf(roi_dist)) = ;
                            clash = nansum(nansum(triu(roi_dist,1) - triu(roi_dist,2)));
                            fprintf('cluster clash: %.12f mm\n',clash)
                        end
                        %save(sprintf('./cache/figure_t14_%i_atl%i_%s',iM,atl,AtlNames{atl}),'AdjMag','AdjMagVar','Adj_plt','Adj_plt2','Adj_plt_null','Adj_plt2_cl','Y','Z','M','cluster_i','rois_plt','atl');
                        %return
                        %cluster_i = optimalleaforder(Z,Y);
                    end

                    % -------------------------------------------------------------------------

                    % moved saving up
    %                 if (iii2 == 1)
    %                     save(sprintf('./cache/figure_t14_%i_atl%i_%s',iM,atl,AtlNames{atl}),'Adj_plt','Adj_plt2','Adj_plt2_cl','Y','Z','M','cluster_i','rois_plt','atl');
    %                 end

                    % bypass clustering
                    %cluster_i = 1:length(cluster_i);
                    %Cai = load('cache/fig_cluster3_cluster_i.mat');
                    %cluster_i = Cai.cluster_i;

                    % Apply clustering
                    Im = Im(cluster_i,cluster_i,:);
                    %rois_plt = rois_plt(1:length(cluster_i));
                    rois_plt = rois_plt(cluster_i);

                    % Print Figure 2A coherence
                    if (atl==2)
                        if (iii2 == 1)
                            coh_mtxt = 'mean';
                            Adj_tmp = Adj_plt2(cluster_i,cluster_i);
                            %close all;
                        elseif (iii2 == 2)
                            coh_mtxt = 'std';
                            Adj_tmp = Adj_plt2(cluster_i,cluster_i);
                            %imagesc(Adj_tmp); colormap inferno; colorbar;
                            %return
                        elseif (iii2 == 3)
                            coh_mtxt = 'nsigpair';
                            Adj_tmp = Adj_plt2(cluster_i,cluster_i);
                        elseif (iii2 == 4)
                            coh_mtxt = 'npair';
                            Adj_tmp = Adj_plt2(cluster_i,cluster_i);
                        elseif (iii2 == 5)
                            coh_mtxt = 'nsigusub';
                            Adj_tmp = Adj_plt2(cluster_i,cluster_i);
                        elseif (iii2 == 6)
                            coh_mtxt = 'nusub';
                            Adj_tmp = Adj_plt2(cluster_i,cluster_i);
                        end

                        idx_1 = strcmp(rois_plt,'Pars Opercularis');
                        idx_2 = strcmp(rois_plt,'Superior Temporal');
                        fprintf('[*] %s - %s, coh %s: %.2f\n',...
                            rois_plt{idx_1},rois_plt{idx_2},coh_mtxt,Adj_tmp(idx_1,idx_2));
                    end


                    % Manual adjustment
                    n_wrap = 0;
                    %n_wrap = round(0.2258 * n_rois_cl); % 9 - Translate the whole clustering and wrap around
                    cluster_i_manual = [(n_rois_cl - (n_wrap-1)):n_rois_cl, 1:(n_rois_cl - (n_wrap))];
                    Im = Im(cluster_i_manual,cluster_i_manual,:);
                    rois_plt = rois_plt(cluster_i_manual);

                    rois_plt_sav = rois_plt;


                    imagesc(Im); hold all;
                    %colordef_draw_cross_nocov;
                    
                    
                    %imagesc(Im, 'Parent', ax);

                    % 
                    % % --- Axis label staggering---------------------
                    % bip_labels = rois_plt;
                    % stagger_labels = true;
                    % xtick_vec = 1:length(bip_labels);
                    % ytick_vec = 1:length(bip_labels);
                    % ax1_xind = true(1,length(xtick_vec));
                    % ax1_yind = true(1,length(ytick_vec));
                    % ax2_xind = true(1,length(xtick_vec));
                    % ax2_yind = true(1,length(ytick_vec));
                    % if (stagger_labels)
                    %     odd_i = downsample(1:length(ax1_xind),2,1);
                    %     even_i = downsample(1:length(ax1_xind),2,0);
                    %     ax1_xind(even_i) = false;
                    %     ax2_xind(odd_i) = false;
                    %     ax1_yind(odd_i) = false;
                    %     ax2_yind(even_i) = false;
                    % end
                    % 
                    % % Show axis labels
                    % %ax1 = ax;
                    % ax = gca;
                    % set(ax,'xtick',xtick_vec(ax1_xind),'Ticklength',[0.001, 0.001],...
                    %     'xticklabel',bip_labels(ax1_xind),'fontsize',fsize,'TickDir','out');
                    % xtickangle(90);
                    % set(ax,'ytick',ytick_vec(ax1_yind),'Ticklength',[0.001, 0.001],...
                    %     'yticklabel',bip_labels(ax1_yind),'fontsize',fsize,'TickDir','out');
                    % 
                    % ax2 = copyobj(ax,ax.Parent);
                    % set(ax2,'XAxisLocation','top');
                    % set(ax2,'YAxisLocation','right');
                    % set(ax2,'xtick',xtick_vec(ax2_xind),'xticklabel',bip_labels(ax2_xind));
                    % set(ax2,'ytick',ytick_vec(ax2_yind),'yticklabel',bip_labels(ax2_yind));
                    % return

                    % --- END Axis label staggering---------------------

                     % colormap
                    % colormap(ax,map)
                    % colormap(ax2,map)
                    % if (minv == maxv)
                    %     caxis(ax,[(minv - 1e-6), (maxv + 1e-6)]); 
                    %     caxis(ax2,[(minv - 1e-6), (maxv + 1e-6)]); 
                    % else
                    %     caxis(ax,[minv maxv]);
                    %     caxis(ax2,[minv maxv]);
                    % end

                    set(gca,'YTickMode','manual')
                    set(gca,'XTickMode','manual')

                    % numbering
                    stride_num = 1;
                    rois_plt2 = cell(size(rois_plt));
                    for inum = 1:length(rois_plt)
                        if (mod(inum-1,stride_num) == 0)
                            n_suf = sprintf('%i',inum);
                            while (length(n_suf) < 4)
                                n_suf = [' ',n_suf];
                            end
                        else
                            n_suf = '    ';
                        end
                        rois_plt{inum} = [rois_plt{inum},'\bf{',n_suf,'}'];
                        rois_plt2{inum} = rois_plt{inum}; %['\bf{',n_suf,'}'];
                        %rois_plt{inum} = ['\bf{',rois_plt{inum},'} \t'];
                    end

                    sub_idx = stride_axislabel:stride_axislabel:length(rois_plt);
                    rty = 1:length(rois_plt);
                    yticks(rty(sub_idx));
                    yticklabels(rois_plt(sub_idx));

                    % Horizontal axis
                    rtx = 1:length(rois_plt2);
                    xticks(rtx(sub_idx));
                    xticklabels(rois_plt2(sub_idx));
                    xtickangle(90);

        %             yticks(1:length(rois_plt));
        %             yticklabels(rois_plt);
        %             xticks(1:length(rois_plt));
        %             xticklabels(rois_plt);
        %             xtickangle(90);

                    set(gca,'tickdir','out');
                    set(gca,'fontsize',fontsz*(31/n_rois_cl)^(0.5));
                    set(gca,'TickLength',[0.001, 0.001])
                    daspect([1 1 1]);
                    set(gca,'FontName','Arial');

                    % Move axis if needed
                    ax = gca;
                    tp = ax.Position;
                    ax1_xoffset = 0.05; % 0.05
                    tp(1) = tp(1) - ax1_xoffset;
                    ax.Position = tp;
                    if(ax.FontSize > 10)
                        ax.FontSize = 10;
                    end
            %         
            %         if (strcmp(metric,'pcBroadband'))
            %             title('Functional Interactions')
            %         else
            %             title(sprintf('Functional Interactions - %s',metric(3:end)))
            %         end

                    %return

                    if (iii == 1)
                        if (trig_mag_no_cp)
                            ss = 'Mag';
                        else
                            ss = 'CP';
                        end
                    elseif (iii == 2)
                        ss = 'bin';
                    elseif (iii == 3)
                        if (trig_mag_no_cp)
                            if (trig_ct)
                                ss = ['CT-bin_',AtlNames{atl}];
                            else
                                ss = ['Mag-bin_',AtlNames{atl}];
                            end
                            %ss = sprintf('%s_cov-',ss,round(1000*));
            %                 mag_min = min(AdjMag(AdjMag ~= 0 ));
            %                 mag_max = max(AdjMag(AdjMag ~= 0 ));
                        else
                            ss = 'CP-bin';
                        end
                    end
                    %set(gca,'FontName','Times New Roman');

            %         colormap(map);
            %         if (((~isnan(minv)) && (~isnan(maxv))) && (minv ~= maxv))
            %             ctick = linspace(minv,maxv,5);
            %             ctickT = cell(size(ctick));
            %             for ict = 1:length(ctickT)
            %                 ctickT{ict} = sprintf('%.2f',ctick(ict));
            %             end
            %             cb = colorbar('Ytick',ctick,'TickLabels',ctickT);
            %             set(cb,'TickLength',0);
            %             caxis([minv maxv]);
            %         else
            %             cb = colorbar;
            %             set(cb,'TickLength',0);
            %         end

                    % New colormap handling
                    %-----------------------------------------
                    colormap(map);
                    if (((~isnan(minv)) && (~isnan(maxv))) && (minv ~= maxv))
                        if (iii <= 2)
                            % Coherence values use 3 ticks
                            cytick = linspace(minv,maxv,3);
                            cytick_str = cell(1,length(cytick));
                        else
                            % Integer values use increments
                            cytick = linspace(minv,maxv,3);
                            cytick_str = cell(1,length(cytick));
                        end

                        for kk = 1:length(cytick)
                            if (iii2 == 1) % Mean of Coherence
                                cytick_str{kk} = sprintf('%.2f',cytick(kk));
                            elseif (iii2 == 2) % Std of Coherence
                                cytick_str{kk} = sprintf('%.3f',cytick(kk));
                            elseif (iii2 == 3) % Number of significant pairs
                                cytick_str{kk} = sprintf('%.0f',cytick(kk));
                            elseif (iii2 == 4) % Number of pairs
                                cytick_str{kk} = sprintf('%.0f',cytick(kk));
                            elseif (iii2 == 5) % Number of significant subjects
                                cytick_str{kk} = sprintf('%.0f',cytick(kk));
                            elseif (iii2 == 6) % Number of subjects
                                cytick_str{kk} = sprintf('%.0f',cytick(kk));
                            end
                        end

                        if ((iii2 <= 2))
                            % Coherence values
                            cb = colorbar('ytick',cytick,'yticklabel',cytick_str,'FontSize',fontsz);
                        else
                            % Integer values
                            cb = colorbar('FontSize',fontsz);
                        end
                        
                        caxis([minv maxv])
                        
                        
                        fprintf('\t[*] minv: %.12f\n',minv);
                        fprintf('\t[*] maxv: %.12f\n',maxv);
                        
                        
                        
                    else
                        cb = colorbar('ytick',minv);
                    end
                    set(cb,'TickDir','out');

                    if ((iii2 <= 2))
                        % Coherence values
                        set(cb,'TickLength',0);
                    else
                        % Integer values
                        %set(cb,'TickLength',0);
                    end

                    %--- display frequency band range ----
                    metricTxt = metric(3:end);
                    if (strcmp(metricTxt,'Broadband'))
                        if (iii2 == 1)
                            if (trig_ct)
                                ylabel(cb,sprintf('Time Consistency'));
                            else
                                ylabel(cb,sprintf('Coherence'));
                            end
                        elseif (iii2 == 2)
                            if (trig_ct)
                                ylabel(cb,sprintf('Std of Time Consistency'));
                            else
                                ylabel(cb,sprintf('Std of Coherence'));
                            end
                        elseif (iii2 == 3)
                            ylabel(cb,sprintf('Number of significant pairs'));
                        elseif (iii2 == 4)
                            ylabel(cb,sprintf('Number of pairs'));
                        elseif (iii2 == 5)
                            ylabel(cb,sprintf('Number of significant subjects'));
                        elseif (iii2 == 6)
                            ylabel(cb,sprintf('Number of subjects'));
                        end
                    else
                        if (iii2 == 1)
                            if (trig_ct)
                                ylabel(cb,sprintf('%s (%s) Time Consistency',metric(3:end),metrics_suffix{iM}));
                            else
                                ylabel(cb,sprintf('%s Coherence (%s)',metric(3:end),metrics_suffix{iM}));
                            end
                        elseif (iii2 == 2)
                            if (trig_ct)
                                ylabel(cb,sprintf('Std of %s Time Consistency (%s)',metric(3:end),metrics_suffix{iM}));
                            else
                                ylabel(cb,sprintf('Std of %s Coherence (%s)',metric(3:end),metrics_suffix{iM}));
                            end
                        elseif (iii2 == 3)
                            ylabel(cb,sprintf('Number of significant pairs'))
                        elseif (iii2 == 4)
                            ylabel(cb,sprintf('Number of pairs'))
                        elseif (iii2 == 5)
                            ylabel(cb,sprintf('Number of significant subjects'))
                        elseif (iii2 == 6)
                            ylabel(cb,sprintf('Number of subjects'))
                        end
                    end
                    %-------------------------------------

                    colormap(map);
                    cb.Location = 'eastoutside';
                    cbpos = cb.Position;
                    cbpos(1) = ax.Position(1)+ax.Position(3)+0.11; %0.86;
                    cbpos(2) = ax.Position(2)+ax.Position(4)/3; %cbpos(2) + 0.2; %2.6*cbpos(2);
                    cbpos(3) = 0.02;
                    cbpos(4) = ax.Position(4)*0.5; %0.4; %0.6*cbpos(4);
                    cb.Position = cbpos;


                    annotation('rectangle',[cb.Position(1) cb.Position(2)-0.06 cb.Position(3) 0.02],'FaceColor',color_not_sig0);
                    annotation('textbox',[cb.Position(1)+0.02 cb.Position(2)-0.06 0.2 0.02],...
                        'String','Not Significant','FitBoxToText','on','EdgeColor','none',...
                        'VerticalAlignment','middle','FontName','Arial','FontSize',9);
                    annotation('rectangle',[cb.Position(1) cb.Position(2)-0.09 cb.Position(3) 0.02],'FaceColor',color_nocov);
                    annotation('textbox',[cb.Position(1)+0.02 cb.Position(2)-0.09 0.2 0.02],...
                        'String','No Coverage','FitBoxToText','on','EdgeColor','none',...
                        'VerticalAlignment','middle','FontName','Arial','FontSize',9);
                    %-------------------------------------------------



                    % Dendrogram
                    if (iii2 == 1)

                        ax2_offset = 0; %(-1)*(ax1.Position(1)+ax1.Position(3))*0.004; %-0.016;
                        switch (atl)
                            case (1)
                                ax2_offset = -0.016;
                            case (2)
                                ax2_offset = -0.022;
                            case (3)
                                ax2_offset = 0;
                            case (4)
                                ax2_offset = -0.02;
                            case (5)
                                ax2_offset = -0.02;
                            case (6)
                                ax2_offset = -0.003;
                                ax.Position(1) = ax.Position(1)-0.01;
                            case (7)
                                ax2_offset = -0.025;
                                ax.Position(1) = ax.Position(1)-0.0001;
                            case (8)
                                ax2_offset = -0.033;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (9)
                                ax2_offset = -0.003;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (10)
                                ax2_offset = -0.019;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (11)
                                ax2_offset = -0.003;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (12)
                                ax2_offset = -0.010;
                            case (13)
                                ax2_offset = -0.009;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (15)
                                ax2_offset = -0.033;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (16)
                                ax2_offset = -0.015;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (17)
                                ax2_offset = -0.08;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (18)
                                ax2_offset = -0.035;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (19)
                                ax2_offset = -0.05;
                                ax.Position(1) = ax.Position(1)-0.00001;
                            case (20)
                                ax2_offset = -0.075;
                                ax.Position(1) = ax.Position(1)-0.00001;
                        end

                        %ax2 = axes('Position',[ax1.Position(1)+ax1.Position(3)+ax2_offset,ax1.Position(2),0.1,ax1.Position(4)]);
                        ax2 = axes('Position',[ax.Position(1)+ax.Position(3)+ax2_offset,ax.Position(2),0.1,ax.Position(4)]);
                        hD = dendrogram(Z,0,'Reorder',fliplr(cluster_i),'Orientation','right','ColorThreshold',0); %[a,b,c]


                        for ihD = 1:length(hD)
                            hD(ihD).Color = 0*[1 1 1];
                        end
                        axis tight;
                        ax2.YLim = [0.5 ax2.YLim(2)+0.5];
                        axis off;
                    end


                    % Get info
                    if ((iii2 == 1) && ((iM == 1) ||(iM == 5)) && (atl == 2))
                        Adiag = diag(Adj_plt2);
                        fprintf('[!!] %i of %i in diagonal are nans.\n',sum(isnan(Adiag)),length(Adiag) );
                        % 6 of 31 in diagonal are nan, 25 are covered
                        AdiagC = Adiag(~isnan(Adiag));
                        fprintf('[!!] %i of %i covered are significant\n',sum(AdiagC~=0),length(AdiagC~=0));
                        % 17 of 25 are significant (68%)

                        AM = Adj_plt2(cluster_i,cluster_i);
                        AMvar = Adj_plt2_var(cluster_i,cluster_i);

                        AM(AM==0) = NaN;
                        [min_vals,min_idxs] = min(AM);
                        [minv,mini] = min(min_vals);
                        minv_std = sqrt(AMvar(mini,min_idxs(mini)));
                        fprintf('[!!] minimum at: %i, %i, coherence = %.4f +- %.5f\n',mini,min_idxs(mini),minv,minv_std);
                        fprintf('\t%s, %s\n',rois_plt2{mini},rois_plt2{min_idxs(mini)})

                        AM = Adj_plt2(cluster_i,cluster_i);
                        AM(AM==0) = NaN;
                        [max_vals,max_idxs] = max(AM);
                        [maxv,maxi] = max(max_vals);
                        maxv_std = sqrt(AMvar(maxi,max_idxs(maxi)));
                        fprintf('[!!] maximum at: %i, %i, coherence = %.4f +- %.5f\n',maxi,max_idxs(maxi),maxv,maxv_std);
                        fprintf('\t%s, %s\n',rois_plt2{maxi},rois_plt2{max_idxs(maxi)})

                        % min and max regions
                        AM2 = Adj_plt2(cluster_i,cluster_i);
                        [a,b] = min(sum(~isnan(AM)));
                        fprintf('[!!] Min region: %s, total pairs sig: %.2f\n',rois_plt2{b},a);
                        am = AM2(b,:);
                        am = am(~isnan(am));
                        fprintf('\t%i of %i sig\n',sum(am~=0),length(am));

                        AM2 = Adj_plt2(cluster_i,cluster_i);
                        [a,b] = max(sum(~isnan(AM)));
                        fprintf('[!!] Max region: %s, total pairs sig: %.2f\n',rois_plt2{b},a);
                        am = AM2(b,:);
                        am = am(~isnan(am));
                        fprintf('\t%i of %i sig\n',sum(am~=0),length(am));
                        %return
                    end


                    %export_fig(sprintf('figures/T14/Adj_%s_%s_2',metric,ss),'-eps');
                    if (iii2 == 2)
                        ss = [ss,'_var'];
                    elseif (iii2 == 3)
                        ss = [ss, '_nsigpairs'];
                    elseif (iii2 == 4)
                        ss = [ss, '_npairs'];
                    elseif (iii2 == 5)
                        ss = [ss, '_nsigsubs'];
                    elseif (iii2 == 6)
                        ss = [ss, '_nsubs'];
                    end

                    %
                    print(h,sprintf('figures/T14_allatl_resampSage/atl%i_%s/Adj_%s_%s_ageg-%i',atl,AtlNames{atl},metric,ss,iAgeG),fig_fmt);
                    if (trig_eps)
                        print(h,sprintf('figures/T14_allatl_resampSage/atl%i_%s/Adj_%s_%s_ageg-%i',atl,AtlNames{atl},metric,ss,iAgeG),'-depsc');
                    end

                    % save to web figure
    %                 if (iii2 == 1)
    %                     print(h,sprintf('/home/jerry/Nextcloud2/Wangetal_PutativeInteractome1_nc/_figures_web/Figure_W3/Figure_W3_atl%i_%s_Adj_%i_%s',atl,AtlNames{atl},iM,metric(3:end)),'-dpng','-r300');
    %                     print(h,sprintf('/home/jerry/Nextcloud2/Wangetal_PutativeInteractome1_nc/_figures_web/Figure_W3/Figure_W3_atl%i_%s_Adj_%i_%s',atl,AtlNames{atl},iM,metric(3:end)),'-depsc','-r300');
    %                     
    %                     % Load atlas maps
    %                     AC = load('./cache/custom_parcellation_atlas_overlap.mat');
    %                     
    %                     %return
    %                     % Save for json export
    %                     S = struct();
    %                     S.Img = Im;
    %                     S.A = Adj_plt2(cluster_i,cluster_i);
    %                     S.Astd = sqrt(Adj_plt2_var(cluster_i,cluster_i));
    %                     S.A_npairs = Adj_plt2_npairs(cluster_i,cluster_i);
    %                     S.A_npairs_sig = Adj_plt2_npairs_sig(cluster_i,cluster_i);
    %                     S.A_nusubs = Adj_plt2_nusubs(cluster_i,cluster_i);
    %                     S.A_nusubs_sig = Adj_plt2_nusubs_sig(cluster_i,cluster_i);
    %                     S.A_min = nanmin(nanmin(S.A(S.A ~= 0)));
    %                     S.A_max = nanmax(nanmax(S.A));
    %                     S.dendro_Z = Z;
    %                     S.dendro_reorder = fliplr(cluster_i);
    %                     S.labels = rois_plt_sav;
    %                     S.atl_name = AtlNames{atl};
    %                     S.atl_chans = AC.AChans{atl}; %1:length(Adj_plt2);
    %                     S.atl_chansR = AC.AChansR{atl};
    %                     save(sprintf('./brainexport/controls_data_sub-0_freq-%i_atl-%i',iM,atl),'S');
    %                 end

                    %return
                    close(h);



                    % Dendrogram
    %                 if (iii2 == 1)
    %                     h = figure('visible','off');
    %                     set(h,'PaperUnits','Inches');
    %                     set(h,'PaperPosition',[0 0 1.5 5]);
    %     %                 cluster_i2 = 1:length(cluster_i);
    %     %                 cluster_i2 = cluster_i2(cluster_i);
    %     %                 cluster_i2 = cluster_i2(cluster_i_manual);
    %                     hD = dendrogram(Z,0,'Reorder',fliplr(cluster_i),'Orientation','right','ColorThreshold',0); %[a,b,c]
    %                     for ihD = 1:length(hD)
    %                         hD(ihD).Color = 0*[1 1 1];
    %                     end
    %                     axis off;
    %                     print(h,sprintf('figures/T14_allatl/atl%i_%s/Adj_%s_%s_dendro',atl,AtlNames{atl},metric,ss),fig_fmt);
    %                     if (trig_eps)
    %                         print(h,sprintf('figures/T14_allatl/atl%i_%s/Adj_%s_%s_dendro',atl,AtlNames{atl},metric,ss),'-depsc');
    %                     end
    %                     close(h);
    %                 end


                    % print info
                    if (iii2 == 1)
                        frac_nocov = sum(sum(isnan(Adj_plt2)))/numel(Adj_plt2);
                        frac_nosig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 == 0)))/numel(Adj_plt2);
                        frac_sig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 ~= 0)))/numel(Adj_plt2);
                        fprintf('[*] saved file: %s\n',sprintf('figures/T14_allatl_resampSage/atl%i_%s/Adj_%s_%s',atl,AtlNames{atl},metric,ss));

                        if (false)
                            fprintf('[*] no coverage: %.6f (%i of %i)\n',frac_nocov,sum(sum(isnan(Adj_plt2))),numel(Adj_plt2));
                            n_nosig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 == 0)));
                            fprintf('[*] not significant: %.6f (%i of %i)\n',frac_nosig,n_nosig,numel(Adj_plt2));
                            n_sig = sum(sum(~isnan(Adj_plt2) & (Adj_plt2 ~= 0)));
                            fprintf('[*] significant: %.6f (%i of %i)\n',frac_sig,n_sig,numel(Adj_plt2));
                            fprintf('[*] significance fraction: %.6f (%i of %i)\n',100*frac_sig/(frac_sig + frac_nosig),n_sig,(n_sig + n_nosig));
                        end

                        % calculate fraction sig
                        Adj_plt = Adj_plt2;
                        n_sig2 = 0;
                        n_all2 = 0;

                        % Vectorized Adj_plt
                        Adj_plt_vec = nan(nchoosek(length(Adj_plt),2),1);
                        c_tmp = 1;

                        for isig1 = 1:(length(Adj_plt)-1)
                            for isig2 = (isig1+1):length(Adj_plt)
                                cond_count = ~isnan(Adj_plt(isig1,isig2));
                                if (cond_count)
                                    if (Adj_plt(isig1,isig2) ~= 0)
                                        n_sig2 = n_sig2 + 1;
                                    end
                                    n_all2 = n_all2 + 1;
                                end

                                % Vectorized Adj_plt
                                Adj_plt_vec(c_tmp) = Adj_plt(isig1,isig2);
                                c_tmp = c_tmp + 1;
                            end
                        end
                        n_nosig2 = n_all2 - n_sig2;
                        fprintf('[!] Significant in : %i of %i (%.2f%%), total possible (%i choose 2): %i\n',n_sig2,(n_sig2 + n_nosig2),100*n_sig2/(n_sig2 + n_nosig2),length(Adj_plt),nchoosek(length(Adj_plt),2));

                        % Save values from mean

                        Adj_plt_mu = Adj_plt;
                        Adj_plt_mu_vec = Adj_plt_vec;
                    end

                    % GAMMA:
                    % [!] Significant in : 150 of 378 (39.68%), total possible: 465
                end
                %##################################################################
                % End mean/variance loop
                %return
            end
            %return


            %%
            % Important variables:
            %
            % AdjMag_Hi
            % AdjMag_Lo
            % Vec_ij
            % DropoutN
            %
            %clear;
            %load('./cache/figure_t14_allatl_resamp');
            load(sprintf('./cache/figure_t14_allatl_resampSage_im-%i',iM));





            % Filter
            AdjMag3_k = AdjMag3(known_idx,known_idx,1);
            AdjMag3_1_k = AdjMag3_1(known_idx,known_idx,1);
            AdjMag3_2_k = AdjMag3_2(known_idx,known_idx,1);

            % Vectorize
            n_known = sum(known_idx);
            Av_0 = zeros(nchoosek(n_known,2),1);
            Av_1 = zeros(nchoosek(n_known,2),1);
            Av_2 = zeros(nchoosek(n_known,2),1);

            c = 0;
            for i1 = 1:(n_known - 1)
                for i2 = (i1 + 1):n_known
                    Av_0(c+1) = AdjMag3_k(i1,i2);
                    Av_1(c+1) = AdjMag3_1_k(i1,i2);
                    Av_2(c+1) = AdjMag3_2_k(i1,i2);
                    c = c + 1;
                end
            end


            for iz = 1:length(AdjMag_Hi)
                i1 = Vec_ij(iz,1);
                i2 = Vec_ij(iz,2);

                if ( (~ known_idx(i1)) || (~ known_idx(i1)) )
                    AdjMag_Hi(iz) = NaN;
                    AdjMag_Lo(iz) = NaN;
                end
            end



            pct_hi = nanmean(100 * ((AdjMag_Hi(:,2:end) ./ AdjMag_Hi(:,1)) - 1));
            pct_lo = nanmean(100 * (1 - (AdjMag_Lo(:,2:end) ./ AdjMag_Lo(:,1))));
            pct_hi_s = nanstd(100 * ((AdjMag_Hi(:,2:end) ./ AdjMag_Hi(:,1)) - 1));
            pct_lo_s = nanstd(100 * (1 - (AdjMag_Lo(:,2:end) ./ AdjMag_Lo(:,1))));

            fprintf('[*] Number of permutations: %i\n',n_perm);
            for i = 1:length(pct_hi)
                fprintf('# subject drop: %i, 95%% CI Hi: %.2f%% (+- %.2f%%), Lo:%.2f%% (+- %.2f%%)\n',DropoutN(i+1),pct_hi(i),pct_hi_s(i),pct_lo(i),pct_lo_s(i));

                [p,H,stats] = ranksum(v_A,v_A1);
            end


            % Remove zeros and nans
            key = AdjMag_Hi(:,1);
            idx_sig = (~isnan(key)) & ((key) ~= 0);
            AdjMag_Hi = AdjMag_Hi(idx_sig,:);
            AdjMag_Lo = AdjMag_Lo(idx_sig,:);
            Vec_ij = Vec_ij(idx_sig,:);


            % Print information
            fprintf('[*] metric: %i, atl: %i\n',iM,atl);
            fprintf('\tBip-pair full n=%i\n',elec_count);
            fprintf('\tFull\tn_sig=%i\n',sum(Av_0>0));
            fprintf('\tSplit1\tn_sig=%i\n',sum(Av_1>0));
            fprintf('\tSplit2\tn_sig=%i\n',sum(Av_2>0));
            idx_full = (Av_0 > 0);
            diff = (Av_1(idx_full) - Av_2(idx_full));
            diffPct = 100*(Av_1(idx_full) - Av_2(idx_full))./Av_1(idx_full);

            [p,h,s] = ranksum(Av_1(idx_full),Av_2(idx_full));
            fprintf('\tRanksum\tp=%.4d\n',p);
            [h,p,s] = ttest(Av_1(idx_full),Av_2(idx_full),'tail','both');
            fprintf('\tT-test\tp=%.4d\n',p);
            fprintf('\tdiff mean: %.4f (%.2f%%)\n',nanmean(diff),nanmean(diffPct));
            fprintf('\tdiff std: %.4f (%.2f%%)\n',nanstd(diff),nanstd(diffPct));
            fprintf('\tdiff n: %i\n',sum(~isnan(diff)));


            close all;
            h = figure('Position',[0 0 400 300],'visible','off');
            histogram(diffPct,'Normalization','pdf','FaceColor',[1 1 1]*0.5);
            set(gca,'TickDir','Out');

            xlabel('Coherence Change (%, Younger - Older)','FontSize',10);
            ylabel('pdf');
            title('Dataset split by Age');

            axis tight;
            xlim(max(abs(diffPct)) * [-1 1]);
            grid off;
            box off;
            set(gcf,'renderer','painters');
            print(h,sprintf('figures/figure_t14_allatl_resampSage_metric%i_atl%i_diff',iM,atl),'-depsc');
            close(h);

            %--------------------------------------------------------



            % Sort
            key = AdjMag_Hi(:,1);
            [~,idx_sort] = sort(key,'descend');
            AdjMag_Hi = AdjMag_Hi(idx_sort,:);
            AdjMag_Lo = AdjMag_Lo(idx_sort,:);
            Vec_ij = Vec_ij(idx_sort,:);

            fprintf('[*] Number of edges compared: %i\n',length(AdjMag_Hi(:,1)));

            % X-axis
            X = 1:length(AdjMag_Hi);

            close all;
            h = figure('Position',[0 0 800 200],'visible','off');
            hold all;

            % Colormap
            try
                cmap = inferno(length(DropoutN) - 1);
                cmap = [(0.5*[1 1 1]); cmap];
            catch
                cmap = [1 1 1]*0.5;
            end
            % Plot
            DropoutNp = DropoutN(1:end);
            for iR = 1:length(DropoutNp)
                i = length(DropoutNp) - iR + 1;
                if (i == 1)
                    mk_size = 5;
                else
                    mk_size = 2;
                end
                Y = AdjMag_Hi(:,i);
                plot(X(Y~=0),Y(Y~=0),'.','Color',cmap(i,:),'MarkerSize',mk_size);
                Y = AdjMag_Lo(:,i);
                plot(X(Y~=0),Y(Y~=0),'.','Color',cmap(i,:),'MarkerSize',mk_size);
            end

            % X labels
            xlab_txt = cell(length(X),1);
            for iL = 1:length(X)
                xlab_txt{iL} = sprintf('%2i,%2i',Vec_ij(iL,1),Vec_ij(iL,2));
            end
            xticks(X);
            xticklabels(xlab_txt);
            xtickangle(90);

            % Y labels
            coh_max = nanmax(AdjMag_Hi(:));
            coh_min = nanmin(AdjMag_Lo(:));
            yt_vals = linspace(coh_min,coh_max,5);
            yt_txt = cell(size(yt_vals));
            yticks(yt_vals);
            for iY = 1:length(yt_vals)
                yt_txt{iY} = sprintf('%6.2f',yt_vals(iY));
            end
            yticklabels(yt_txt);

            % Axis font sizes
            xAX = get(gca,'XAxis');
            set(xAX,'FontSize',3);
            set(xAX,'FontName','FixedWidth'); %FixedWidth
            set(xAX,'TickLength',[0.001, 0.001]);
            yAX = get(gca,'YAxis');
            set(yAX,'TickLength',[0.005, 0.005]);

            % Plot format
            set(gca,'TickDir','Out');
            xlabel('Desikan-Killiany area pairs','FontSize',10);
            ylabel('Coherence')
            grid off;
            axis tight;

            print(h,sprintf('figures/figure_t14_allatl_resampSage_metric%i_atl%i',iM,atl),'-depsc');
            close(h);

        end

    end

end


% 3 jun 2021

% 
% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T14_allatl_resampSage’: File exists
% mkdir: cannot create directory ‘figures/T14_allatl_resampSage/atl2_Desikan-Killiany’: File exists
% [*] Total subjects: 48
% 	Cutoff median age: 17.00
% 	Group 1 - 24, mean age: 11.208, std: 3.9450
% 	Group 2 - 24, mean age: 26.917, std: 10.2317
% 	Age ranksum p: 3.68728174e-09
% [n_perm: 1] atl: 2	metric: 1	#drop: 0	of 0
% [*] Correlation age vs. coherence
% 	r=0.105252, p=000000, n=168282
% per atl-metric ETA: 0.0 mins
% [*] Loading existing hier-clustering: ./cache/figure_t14_1_atl2_Desikan-Killiany
% cluster clash: 3.118823393796 mm
% [*] Pars Opercularis - Superior Temporal, coh mean: 0.24
% 	[*] minv: 0.137898698449
% 	[*] maxv: 0.510548055172
% [!!] 9 of 31 in diagonal are nans.
% [!!] 16 of 22 covered are significant
% [!!] minimum at: 4, 4, coherence = 0.1379 +- 0.17348
% 	Inferior Parietal\bf{   4}, Inferior Parietal\bf{   4}
% [!!] maximum at: 26, 30, coherence = 0.5105 +- 0.28933
% 	Pericalcarine\bf{  26}, Precuneus\bf{  30}
% [!!] Min region: Paracentral\bf{  18}, total pairs sig: 1.00
% 	1 of 22 sig
% [!!] Max region: Fusiform\bf{   5}, total pairs sig: 23.00
% 	23 of 29 sig
% [*] saved file: figures/T14_allatl_resampSage/atl2_Desikan-Killiany/Adj_pcBroadband_Mag-bin_Desikan-Killiany
% [!] Significant in : 168 of 353 (47.59%), total possible (31 choose 2): 465
% [*] Pars Opercularis - Superior Temporal, coh std: 0.13
% 	[*] minv: 0.001970777856
% 	[*] maxv: 0.165467841509
% [*] Pars Opercularis - Superior Temporal, coh nsigpair: 16.00
% 	[*] minv: 1.000000000000
% 	[*] maxv: 344.000000000000
% [*] Pars Opercularis - Superior Temporal, coh npair: 441.00
% 	[*] minv: 10.000000000000
% 	[*] maxv: 3258.000000000000
% [*] Pars Opercularis - Superior Temporal, coh nsigusub: 10.00
% 	[*] minv: 1.000000000000
% 	[*] maxv: 28.000000000000
% [*] Pars Opercularis - Superior Temporal, coh nusub: 29.00
% 	[*] minv: 2.000000000000
% 	[*] maxv: 44.000000000000
% [*] Number of permutations: 100000
% [*] metric: 1, atl: 2
% 	Bip-pair full n=0
% 	Full	n_sig=193
% 	Split1	n_sig=168
% 	Split2	n_sig=165
% 	Ranksum	p=1.5682e-17
% 	T-test	p=5.9283e-16
% 	diff mean: -0.0945 (-54.91%)
% 	diff std: 0.1220 (61.63%)
% 	diff n: 140
% [*] Number of edges compared: 193
% mkdir: cannot create directory ‘figures/T14_allatl_resampSage/atl2_Desikan-Killiany’: File exists
% [*] Total subjects: 48
% 	Cutoff median age: 17.00
% 	Group 1 - 24, mean age: 11.208, std: 3.9450
% 	Group 2 - 24, mean age: 26.917, std: 10.2317
% 	Age ranksum p: 3.68728174e-09
% [n_perm: 1] atl: 2	metric: 1	#drop: 0	of 0
% [*] Correlation age vs. coherence
% 	r=0.105252, p=000000, n=168282
% per atl-metric ETA: 0.0 mins
% [*] Loading existing hier-clustering: ./cache/figure_t14_1_atl2_Desikan-Killiany
% cluster clash: 3.623762491379 mm
% [*] Pars Opercularis - Superior Temporal, coh mean: 0.33
% 	[*] minv: 0.143031621973
% 	[*] maxv: 0.617667578161
% [!!] 8 of 31 in diagonal are nans.
% [!!] 17 of 23 covered are significant
% [!!] minimum at: 7, 24, coherence = 0.1430 +- 0.05451
% 	Inferior Temporal\bf{   7}, Superior Parietal\bf{  24}
% [!!] maximum at: 1, 10, coherence = 0.6177 +- 0.50298
% 	Parahippocampal\bf{   1}, Pars Opercularis\bf{  10}
% [!!] Min region: Paracentral\bf{  18}, total pairs sig: 0.00
% 	0 of 21 sig
% [!!] Max region: Fusiform\bf{   5}, total pairs sig: 22.00
% 	22 of 28 sig
% [*] saved file: figures/T14_allatl_resampSage/atl2_Desikan-Killiany/Adj_pcBroadband_Mag-bin_Desikan-Killiany
% [!] Significant in : 165 of 350 (47.14%), total possible (31 choose 2): 465
% [*] Pars Opercularis - Superior Temporal, coh std: 0.13
% 	[*] minv: 0.000529849106
% 	[*] maxv: 0.206156606033
% [*] Pars Opercularis - Superior Temporal, coh nsigpair: 16.00
% 	[*] minv: 1.000000000000
% 	[*] maxv: 344.000000000000
% [*] Pars Opercularis - Superior Temporal, coh npair: 441.00
% 	[*] minv: 10.000000000000
% 	[*] maxv: 3258.000000000000
% [*] Pars Opercularis - Superior Temporal, coh nsigusub: 10.00
% 	[*] minv: 1.000000000000
% 	[*] maxv: 28.000000000000
% [*] Pars Opercularis - Superior Temporal, coh nusub: 29.00
% 	[*] minv: 2.000000000000
% 	[*] maxv: 44.000000000000
% [*] Number of permutations: 100000
% [*] metric: 1, atl: 2
% 	Bip-pair full n=0
% 	Full	n_sig=193
% 	Split1	n_sig=168
% 	Split2	n_sig=165
% 	Ranksum	p=1.5682e-17
% 	T-test	p=5.9283e-16
% 	diff mean: -0.0945 (-54.91%)
% 	diff std: 0.1220 (61.63%)
% 	diff n: 140
% [*] Number of edges compared: 193


% 2 mars 2021

% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T14_allatl_resampSage’: File exists
% mkdir: cannot create directory ‘figures/T14_allatl_resampSage/atl2_Desikan-Killiany’: File exists
% [*] Total subjects: 48
% 	Cutoff median age: 17.00
% 	Group 1 - 24, mean age: 11.208, std: 3.9450
% 	Group 2 - 24, mean age: 26.917, std: 10.2317
% 	Age ranksum p: 3.68728174e-09
% [n_perm: 1] atl: 2	metric: 1	#drop: 0	of 0
% [*] Correlation age vs. coherence
% 	r=0.105252, p=000000, n=168282
% per atl-metric ETA: 0.0 mins
% [*] Loading existing hier-clustering: ./cache/figure_t14_1_atl2_Desikan-Killiany
% cluster clash: 2.781464304392 mm
% [*] Pars Opercularis - Superior Temporal, coh mean: 0.24
% [!!] 9 of 31 in diagonal are nans.
% [!!] 16 of 22 covered are significant
% [!!] minimum at: 13, 13, coherence = 0.1379 +- 0.17348
% 	Inferior Parietal\bf{  13}, Inferior Parietal\bf{  13}
% [!!] maximum at: 26, 30, coherence = 0.5105 +- 0.28933
% 	Pericalcarine\bf{  26}, Precuneus\bf{  30}
% [!!] Min region: Paracentral\bf{  24}, total pairs sig: 1.00
% 	1 of 22 sig
% [!!] Max region: Fusiform\bf{   4}, total pairs sig: 23.00
% 	23 of 29 sig
% [*] saved file: figures/T14_allatl_resampSage/atl2_Desikan-Killiany/Adj_pcBroadband_Mag-bin_Desikan-Killiany
% [!] Significant in : 168 of 353 (47.59%), total possible (31 choose 2): 465
% [*] Pars Opercularis - Superior Temporal, coh std: 0.13
% [*] Pars Opercularis - Superior Temporal, coh nsigpair: 16.00
% [*] Pars Opercularis - Superior Temporal, coh npair: 441.00
% [*] Pars Opercularis - Superior Temporal, coh nsigusub: 10.00
% [*] Pars Opercularis - Superior Temporal, coh nusub: 29.00
% [*] Number of permutations: 100000
% [*] metric: 1, atl: 2
% 	Bip-pair full n=0
% 	Full	n_sig=193
% 	Split1	n_sig=168
% 	Split2	n_sig=165
% 	Ranksum	p=1.5682e-17
% 	T-test	p=5.9283e-16
% 	diff mean: -0.0945 (-54.91%)
% 	diff std: 0.1220 (61.63%)
% 	diff n: 140
% [*] Number of edges compared: 193
% [*] Total subjects: 48
% 	Cutoff median age: 17.00
% 	Group 1 - 24, mean age: 11.208, std: 3.9450
% 	Group 2 - 24, mean age: 26.917, std: 10.2317
% 	Age ranksum p: 3.68728174e-09
% [n_perm: 1] atl: 2	metric: 5	#drop: 0	of 0
% [*] Correlation age vs. coherence
% 	r=0.118723, p=000000, n=153018
% per atl-metric ETA: 0.0 mins
% [*] Loading existing hier-clustering: ./cache/figure_t14_5_atl2_Desikan-Killiany
% cluster clash: 2.810592244817 mm
% [*] Pars Opercularis - Superior Temporal, coh mean: 0.24
% [!!] 10 of 31 in diagonal are nans.
% [!!] 15 of 21 covered are significant
% [!!] minimum at: 12, 21, coherence = 0.1550 +- 0.36022
% 	Entorhinal\bf{  12}, Inferior Parietal\bf{  21}
% [!!] maximum at: 5, 9, coherence = 0.4921 +- 0.29170
% 	Pericalcarine\bf{   5}, Precuneus\bf{   9}
% [!!] Min region: Frontal Pole\bf{   3}, total pairs sig: 0.00
% 	0 of 16 sig
% [!!] Max region: Precentral\bf{  28}, total pairs sig: 20.00
% 	20 of 27 sig
% [*] saved file: figures/T14_allatl_resampSage/atl2_Desikan-Killiany/Adj_pcGamma_Mag-bin_Desikan-Killiany
% [!] Significant in : 159 of 353 (45.04%), total possible (31 choose 2): 465
% [*] Pars Opercularis - Superior Temporal, coh std: 0.13
% [*] Pars Opercularis - Superior Temporal, coh nsigpair: 18.00
% [*] Pars Opercularis - Superior Temporal, coh npair: 441.00
% [*] Pars Opercularis - Superior Temporal, coh nsigusub: 11.00
% [*] Pars Opercularis - Superior Temporal, coh nusub: 29.00
% [*] Number of permutations: 100000
% [*] metric: 5, atl: 2
% 	Bip-pair full n=0
% 	Full	n_sig=184
% 	Split1	n_sig=159
% 	Split2	n_sig=150
% 	Ranksum	p=1.5364e-19
% 	T-test	p=1.4560e-19
% 	diff mean: -0.1032 (-53.83%)
% 	diff std: 0.1069 (52.88%)
% 	diff n: 125
% [*] Number of edges compared: 184
% mkdir: cannot create directory ‘figures/T14_allatl_resampSage/atl2_Desikan-Killiany’: File exists

% [*] Total subjects: 48
% 	Cutoff median age: 17.00
% 	Group 1 - 24, mean age: 11.208, std: 3.9450
% 	Group 2 - 24, mean age: 26.917, std: 10.2317
% 	Age ranksum p: 3.68728174e-09
% [n_perm: 1] atl: 2	metric: 1	#drop: 0	of 0
% [*] Correlation age vs. coherence
% 	r=0.105252, p=000000, n=168282
% per atl-metric ETA: 0.0 mins
% [*] Loading existing hier-clustering: ./cache/figure_t14_1_atl2_Desikan-Killiany
% cluster clash: 4.704779810071 mm
% [*] Pars Opercularis - Superior Temporal, coh mean: 0.33
% [!!] 8 of 31 in diagonal are nans.
% [!!] 17 of 23 covered are significant
% [!!] minimum at: 1, 23, coherence = 0.1430 +- 0.05451
% 	Inferior Temporal\bf{   1}, Superior Parietal\bf{  23}
% [!!] maximum at: 16, 17, coherence = 0.6177 +- 0.50298
% 	Pars Opercularis\bf{  16}, Parahippocampal\bf{  17}
% [!!] Min region: Caudal Anterior Cingulate\bf{  18}, total pairs sig: 0.00
% 	0 of 3 sig
% [!!] Max region: Fusiform\bf{   4}, total pairs sig: 22.00
% 	22 of 28 sig
% [*] saved file: figures/T14_allatl_resampSage/atl2_Desikan-Killiany/Adj_pcBroadband_Mag-bin_Desikan-Killiany
% [!] Significant in : 165 of 350 (47.14%), total possible (31 choose 2): 465
% [*] Pars Opercularis - Superior Temporal, coh std: 0.13
% [*] Pars Opercularis - Superior Temporal, coh nsigpair: 16.00
% [*] Pars Opercularis - Superior Temporal, coh npair: 441.00
% [*] Pars Opercularis - Superior Temporal, coh nsigusub: 10.00
% [*] Pars Opercularis - Superior Temporal, coh nusub: 29.00
% [*] Number of permutations: 100000
% [*] metric: 1, atl: 2
% 	Bip-pair full n=0
% 	Full	n_sig=193
% 	Split1	n_sig=168
% 	Split2	n_sig=165
% 	Ranksum	p=1.5682e-17
% 	T-test	p=5.9283e-16
% 	diff mean: -0.0945 (-54.91%)
% 	diff std: 0.1220 (61.63%)
% 	diff n: 140
% [*] Number of edges compared: 193
% [*] Total subjects: 48
% 	Cutoff median age: 17.00
% 	Group 1 - 24, mean age: 11.208, std: 3.9450
% 	Group 2 - 24, mean age: 26.917, std: 10.2317
% 	Age ranksum p: 3.68728174e-09
% [n_perm: 1] atl: 2	metric: 5	#drop: 0	of 0
% [*] Correlation age vs. coherence
% 	r=0.118723, p=000000, n=153018
% per atl-metric ETA: 0.0 mins
% [*] Loading existing hier-clustering: ./cache/figure_t14_5_atl2_Desikan-Killiany
% cluster clash: 4.135055819556 mm
% [*] Pars Opercularis - Superior Temporal, coh mean: 0.31
% [!!] 8 of 31 in diagonal are nans.
% [!!] 17 of 23 covered are significant
% [!!] minimum at: 1, 7, coherence = 0.1582 +- 0.27858
% 	Lingual\bf{   1}, Superior Parietal\bf{   7}
% [!!] maximum at: 13, 20, coherence = 0.5692 +- 0.31978
% 	Parahippocampal\bf{  13}, Caudal Middle Frontal\bf{  20}
% [!!] Min region: Paracentral\bf{   2}, total pairs sig: 0.00
% 	0 of 17 sig
% [!!] Max region: Middle Temporal\bf{  22}, total pairs sig: 18.00
% 	18 of 31 sig
% [*] saved file: figures/T14_allatl_resampSage/atl2_Desikan-Killiany/Adj_pcGamma_Mag-bin_Desikan-Killiany
% [!] Significant in : 150 of 344 (43.60%), total possible (31 choose 2): 465
% [*] Pars Opercularis - Superior Temporal, coh std: 0.13
% [*] Pars Opercularis - Superior Temporal, coh nsigpair: 18.00
% [*] Pars Opercularis - Superior Temporal, coh npair: 441.00
% [*] Pars Opercularis - Superior Temporal, coh nsigusub: 11.00
% [*] Pars Opercularis - Superior Temporal, coh nusub: 29.00
% [*] Number of permutations: 100000
% [*] metric: 5, atl: 2
% 	Bip-pair full n=0
% 	Full	n_sig=184
% 	Split1	n_sig=159
% 	Split2	n_sig=150
% 	Ranksum	p=1.5364e-19
% 	T-test	p=1.4560e-19
% 	diff mean: -0.1032 (-53.83%)
% 	diff std: 0.1069 (52.88%)
% 	diff n: 125
% [*] Number of edges compared: 184


% 18 fevr. 2021

% [*] Total subjects: 48
% 	Cutoff median age: 17.00
% 	Group 1 - 24, mean age: 11.208, std: 3.9450
% 	Group 2 - 24, mean age: 26.917, std: 10.2317
% 	Age ranksum p: 3.68728174e-09
% [n_perm: 1] atl: 2	metric: 1	#drop: 0	of 0
% [*] Correlation age vs. coherence
% 	r=0.099571, p=5.742354e-199, n=90898
% per atl-metric ETA: 0.0 mins
% [*] Number of permutations: 100000
% [*] metric: 1, atl: 2
% 	Bip-pair full n=0
% 	Full	n_sig=193
% 	Split1	n_sig=168
% 	Split2	n_sig=165
% 	Ranksum	p=1.5682e-17
% 	T-test	p=5.9283e-16
% 	diff mean: -0.0945 (-54.91%)
% 	diff std: 0.1220 (61.63%)
% 	diff n: 140
% [*] Number of edges compared: 193
% [*] Total subjects: 48
% 	Cutoff median age: 17.00
% 	Group 1 - 24, mean age: 11.208, std: 3.9450
% 	Group 2 - 24, mean age: 26.917, std: 10.2317
% 	Age ranksum p: 3.68728174e-09
% [n_perm: 1] atl: 2	metric: 5	#drop: 0	of 0
% [*] Correlation age vs. coherence
% 	r=0.112721, p=1.536568e-233, n=83266
% per atl-metric ETA: 0.0 mins
% [*] Number of permutations: 100000
% [*] metric: 5, atl: 2
% 	Bip-pair full n=0
% 	Full	n_sig=184
% 	Split1	n_sig=159
% 	Split2	n_sig=150
% 	Ranksum	p=1.5364e-19
% 	T-test	p=1.4560e-19
% 	diff mean: -0.1032 (-53.83%)
% 	diff std: 0.1069 (52.88%)
% 	diff n: 125
% [*] Number of edges compared: 184




% 9 fevr. 2021

% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/T14_allatl_resampSage’: File exists
% mkdir: cannot create directory ‘figures/T14_allatl_resampSage/atl2_Desikan-Killiany’: File exists
% [*] Total subjects: 48
% 	Cutoff median age: 17.00
% 	Group 1 - 24, mean age: 11.208, std: 3.9450
% 	Group 2 - 24, mean age: 26.917, std: 10.2317
% 	Age ranksum p: 3.68728174e-09
% [n_perm: 1] atl: 2	metric: 1	#drop: 0	of 0
% [*] Correlation age vs. coherence
% 	r=0.099571, p=5.742354e-199, n=90898
% per atl-metric ETA: 0.0 mins
% [*] Number of permutations: 100000
% [*] metric: 1, atl: 2
% 	Bip-pair full n=0
% 	Full	n_sig=193
% 	Split1	n_sig=168
% 	Split2	n_sig=165
% 	Ranksum	p=1.5682e-17
% 	T-test	p=5.9283e-16
% 	diff mean: -0.0945 (-54.91%)
% 	diff std: 0.1220 (61.63%)
% 	diff n: 140
% [*] Number of edges compared: 193
% [*] Total subjects: 48
% 	Cutoff median age: 17.00
% 	Group 1 - 24, mean age: 11.208, std: 3.9450
% 	Group 2 - 24, mean age: 26.917, std: 10.2317
% 	Age ranksum p: 3.68728174e-09
% [n_perm: 1] atl: 2	metric: 5	#drop: 0	of 0
% [*] Correlation age vs. coherence
% 	r=0.112721, p=1.536568e-233, n=83266
% per atl-metric ETA: 0.0 mins
% [*] Number of permutations: 100000
% [*] metric: 5, atl: 2
% 	Bip-pair full n=0
% 	Full	n_sig=184
% 	Split1	n_sig=159
% 	Split2	n_sig=150
% 	Ranksum	p=1.5364e-19
% 	T-test	p=1.4560e-19
% 	diff mean: -0.1032 (-53.83%)
% 	diff std: 0.1069 (52.88%)
% 	diff n: 125
% [*] Number of edges compared: 184
