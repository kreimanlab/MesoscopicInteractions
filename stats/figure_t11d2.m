close all;
clear;

metric = 'pcBroadband';
n_perm = 10000;

% anatomy folder
%subjects_dir = '/mnt/cuenap_ssd/coregistration';
%subjects_dirL = '/home/klab/data/h5eeg/artifact_removed/test/opencl/coreg';
subjects_dirL = '/media/klab/internal/data/coreg';
setenv('SUBJECTS_DIR',subjects_dirL);

% Fast i/o definitions

dir_art = '/media/klab/KLAB101/h5_notch2/art_nosz';%'/media/klab/internal/data/h5_notch20/art';
dir_res = '/media/klab/KLAB101/results/coh_w10';
dir_cor = '/media/klab/internal/data/coreg';
dir_cacheL = './cache';

% Slow i/o definitions
dir_h5L = '/media/klab/KLAB101/h5_notch20';

metricsL = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
 % 'pcDelta'
% Patients
Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124',...
    'mSu'};

% Exclude monkey
Subjects = Subjects(1:(end-1));

Subjects = {'m00061'};
%Subjects = {'m00095'};


load('rois')
% skip unknown
% rois = rois(2:end);
% %rois = {'middletemporal'};
%rois = {'parstriangularis'};
rois = {'parsopercularis'};
%rois = rois(randperm(length(rois)));

for i_roi = 1:length(rois)

% ROI
%roi_1_const = 'lateraloccipital';
%roi_1_const = 'pericalcarine';
%roi_1_const = 'caudalanteriorcingulate';
%roi_1_const = 'middletemporal';
roi_1_const = rois{i_roi};

fig_fmt = '-depsc';
trig_eps = false;
trig_mag_not_ct = true;
trig_skip_plot = false;
trig_overwrite_cache = false;
system('mkdir figures');
system('mkdir figures/T11d2');

ct_range = [];

%for iM = 1:length(metrics)
for iM = 1
    metric = metricsL{iM};
    
    % Main loop
    for iL = 1:length(Subjects)
        sidL = Subjects{iL};
        %fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
        %fn_dist = sprintf('%s/%s_dists-%s-%i.mat',dir_res,sid,metric,n_perm);
        %fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
        fn_h5L = sprintf('%s/%s.h5',dir_h5L,sidL);
        %fn_coreg = sprintf('%s/%s/label/all_parcellation.mat',dir_cor,sid);

        % Check if files exist
        %ckf = {fn_art,fn_dist,fn_graph,fn_h5,fn_coreg};
        ckf = {fn_h5L};
        for jj = 1:length(ckf)
            if (~exist(ckf{jj},'file'))
                fprintf(2,'E> File not found: %s\n',ckf{jj});
                return
            end
        end

        % check hemisphere
        hemi = {};
        if (exist(sprintf('%s/%s/surf/rh.full.flat.patch.3d',subjects_dirL,sidL), 'file'))
            hemi = {hemi{:},'rh'};
        end
        if (exist(sprintf('%s/%s/surf/lh.full.flat.patch.3d',subjects_dirL,sidL), 'file'))
            hemi = {hemi{:},'lh'};
        end
        if (isempty(hemi))
            fprintf('W> file not found: %s\n',sprintf('%s/%s/surf/*h.full.flat.patch.3d',subjects_dirL,sidL));
        else
            for jj = 1:length(hemi)
                % Read in file
                fn_patch = sprintf('%s/%s/surf/%s.full.flat.patch.3d',subjects_dirL,sidL,hemi{jj});
                fn_pial = sprintf('%s/%s/surf/%s.pial',subjects_dirL,sidL,hemi{jj});
                fn_cache = sprintf('cache/%s_%s_flat.mat',sidL,hemi{jj});
    %             fn_patch_asc = sprintf('%s/%s/surf/%s.full.flat.patch.3d.asc',subjects_dir,sid,hemi{j});
    % 
    %             % if ascii doesn't exist, make it
    %             if (~ exist(fn_patch_asc, 'file'))
    %                 fprintf('[*] ASCII flat patch not found, converting...\n')
    %                 system(sprintf('cd %s/%s/surf; mris_convert -p %s.full.flat.patch.3d %s.full.flat.patch.3d.asc',subjects_dir,sid,hemi{j},hemi{j}));
    %             end


                % Check if patch cache exists:
                if (exist(fn_cache,'file') && (~trig_overwrite_cache))
                    load(fn_cache,'pa','tri_flat','pa','v_pial','tri_pial');
                else
                    [pa] = read_patch_rev(fn_patch);
                    %k = convhull(pa.x,pa.y);
                    [v_pial,tri_pial] = read_surf(fn_pial);
                    % Correct to be 0-indexed like patch indexes
                    %tri_pial = tri_pial + 1;
        %             
        %             % Build truth table for whether triangulation contains flat
        %             % patch vertex pa.ind(ii)
        % 
        %             ind_tri = false(size(tri_pial));
        %             for ii = 1:length(pa.ind)
        %                 ind_s = (tri_pial == pa.ind(ii));
        %                 ind_tri = ind_tri | ind_s;
        %                 if (mod(ii-1,round(length(pa.ind)/10)) == 0)
        %                     fprintf('pial2flat> %i of %i - %.2f %%\n',ii,length(pa.ind),100*ii/length(pa.ind));
        %                 end
        %             end
        %             % Only take triangulation if all three vertices are in patch
        %             ind_tri_pial2flat = ind_tri(:,1) & ind_tri(:,2) & ind_tri(:,3);
        %             tri_flat = tri_pial(ind_tri_pial2flat,:);

                    % old patch to surface registration
                    tri_flat = [];
                    tri_flat_i = 1;
                    for ii = 1:length(tri_pial)
                        tri_1 = find(pa.vno == tri_pial(ii,1),1);
                        tri_2 = find(pa.vno == tri_pial(ii,2),1);
                        tri_3 = find(pa.vno == tri_pial(ii,3),1);
                        if ((~isempty(tri_1)) && (~isempty(tri_2)) && (~isempty(tri_3)))
                            tri_flat(tri_flat_i,:) = [tri_1, tri_2, tri_3];
                            tri_flat_i = tri_flat_i + 1;
                        end
                        if (mod(ii-1,round(length(tri_pial)/10)) == 0)
                            fprintf('%s> %i of %i\n',sidL,ii,length(tri_pial));
                        end
                    end

                    %plot(pa.x,pa.y,'black.');
                    %h = trisurf(tri_flat,pa.x,pa.y,pa.z);
                    %set(h,'edgecolor','none')
                    save(fn_cache,'-v6');
                end

                % Load human cache
                Subjects_sav = Subjects;
                load(sprintf('%s/xsub_out_%s_%i.mat',dir_cacheL,sidL,iM));
                Subjects = Subjects_sav;
                
                chan_labels = h5readatt(fn_h5L,'/h5eeg/eeg','labels');
                bchan_labels = cell(ecog.n_bchan,1);
                for ii0 = 1:ecog.n_bchan
                    bchan_labels{ii0} = sprintf('%s-%s',chan_labels{ecog.bip(ii0,1)},chan_labels{ecog.bip(ii0,2)});
                end
                
                
                
                
                % Plot flat map

                % Read annotation
                %atl_table = AT.P.AtlROIs{atl_i}.RH.table;
                atl_table = C.AtlROIs{atl}.RH.table;
                annot_dir = sprintf('%s/%s/label',subjects_dirL,sidL);
                switch(C.AtlNames{atl})
                    case ('Desikan-Killiany')
                        atl_fname = sprintf('%s/%s.aparc.annot',annot_dir,hemi{jj});
                    case ('Destrieux-2009')
                        atl_fname = sprintf('%s/%s.aparc.a2009s.annot',annot_dir,hemi{jj});
                    case ('HCP-MMP1')
                        atl_fname = sprintf('%s/%s.HCP-MMP1.annot',annot_dir,hemi{jj});
                    case ('M132')
                        atl_fname = sprintf('%s/%s.MACAQUE_M132.annot',annot_dir,hemi{jj});
                end
                [v,l,c] = read_annotation(atl_fname);
                roi_id = l(pa.vno + 1);
                
                % ROI midpoints for fsaverage_sym
                % collect all points for each ROI
%                 roi_pts = cell(length(rois),1);
%                 for iv = 1:length(roi_id)
%                     roi_idx = find(roi_id(iv) == C.AtlROIs{2}.LH.table(:,end));
%                     if (~isempty(roi_idx))
%                         roi = C.AtlROIs{2}.LH.struct_names{roi_idx};
%                         if (~strcmp(roi,'unknown'))
%                             idx2 = find(strcmp(rois,roi));
%                             roi_pts{idx2} = [roi_pts{idx2}, [pa.x(iv); pa.y(iv)]];
%                         end
%                     else
%                         fprintf(2,'W: roi identifying number not found %i\n',roi_id(iv));
%                     end
%                 end
                
                % -- copy
                Ca_fn = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheL,sidL,iM);
                Ca = load(Ca_fn);
                ctL = Ca.ct;
                ct_threshL = Ca.ct_thresh;
                if (~ isfield(Ca,'roi_pts') )
                    roi_pts = cell(length(rois),1);
                    for iv = 1:length(roi_id)
                        roi_idx = find(roi_id(iv) == C.AtlROIs{2}.LH.table(:,end));
                        if (~isempty(roi_idx))
                            roi = C.AtlROIs{2}.LH.struct_names{roi_idx};
                            if (~strcmp(roi,'unknown'))
                                idx2 = find(strcmp(rois,roi));
                                roi_pts{idx2} = [roi_pts{idx2}, [pa.x(iv); pa.y(iv)]];
                            end
                        else
                            %fprintf(2,'W: roi identifying number not found %i\n',roi_id(iv));
                        end
                    end
                    save(Ca_fn,'roi_pts','-append');
                else
                    roi_pts = Ca.roi_pts;
                end
                
                % --
                
                % get midpoint
                roi_midpt = nan(length(rois),2);
                for iv = 1:length(rois)
                    if (~isempty(roi_pts{iv}))
                        roi_midpt(iv,:) = nanmean(roi_pts{iv},2)';
                    end
                end
                
                % read electrode locations
                fprintf('[!] Ignore any errors:\n')
                l = read_label(sprintf('%s/%s',subjects_dirL,sidL),'all_surf_ielvis');
                if (isempty(l))
                    l = read_label(sidL,'all_surf_ielvis');
                end
                
                % constant color
                col_unknown = 0.333*[1 1 1];
                col_src = 0.6*[1 1 1]; % 0.5
                col_surf = 0.9*[1 1 1]; % 0.85
                col_roi_hase = col_surf; %0.8*[1 1 1]; % 0.7
                col_bound = 0.4*[1 1 1];
                col_elec = 0.5*[1 1 1];
                col_elec_hi = col_elec; % 0*[1 1 1];
                col_elec_border = [0 0 0];
                border_fac = 1.2; % how much larger to draw border
                col_elec_nosig = col_roi_hase; %0.8*[1 1 1]; % 0.3
                size_bound = 3; %5
                bound_repeat = 1;
                bound_loop_rep = 1;
                size_elec = 15;
                size_elec_hi = size_elec;
                size_elec_nosig = size_elec;
                trig_show_txt_mid = false;

                
                %col_unknown = 0.333*[1 1 1];
                %col_src = 0.5*[1 1 1];
                %col_surf = 0.85*[1 1 1];
                %col_bound = 0.2*[1 1 1];
                %col_elec = 0*[1 1 1];
                %col_elec_hi = 0.8*[1 1 1];
                col_bip = 0*[1 1 1];
                %size_elec = 20;
                %size_elec_hi = 20;
                %trig_show_txt_mid = true;
                % 
                % ind_lingual = 14;
                % ind_cuneus = 6;
                % ind_inferiorparietal = 9;
                % source = ind_inferiorparietal;
                % %source = 34;

                tri_color = mean(col_surf)*ones(length(roi_id),3);
                n_colors = 100;
                cmapL = corrcmap(n_colors);

                pname = sprintf('%s_%s-all_%.0f-%.0fmm_%s',...
                    sidL,roi_1_const,min(Dmat(:)),max(Dmat(:)),metric);
                
                if (~trig_skip_plot)
                
                    % =======================================================================================

                    % Plot flat patch
                    % -------------------------------------------------------------------------------
                    h = figure;
                    fig_w = 10.5;
                    fig_h = 10.5;
                    set(h,'Position',[0 0 fig_w*100 fig_h*100])
                    set(h, 'PaperUnits', 'Inches')
                    set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
                    set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
                    a = trisurf(tri_flat,pa.x',pa.y',pa.z'); hold all;
                    set(a,'edgecolor','none');

                    % Show boundary
                    %isEven = true;
                    count_bound = 0;
                    count_thresh = 10000;
                    for i3 = 1:bound_loop_rep:length(tri_flat)
                        vB = tri_flat(i3,:); % 6.8 seconds
                        tv = roi_id(vB); % 9.6 seconds
                        if ((length(unique(tv)) > 1))
                            %if (isEven)
                            if (mod(count_bound,bound_repeat) == 0)
                                plot(mean(pa.x(vB)),mean(pa.y(vB)),'.','color',col_bound,'MarkerSize',size_bound)
                            end
                            %isEven = ~ isEven;
                            count_bound = count_bound + 1;
                            if (count_bound > count_thresh)
                                fprintf(2,'W> Boundary point threshold reached. Stopping boundary draw prematurely.\n');
                                break
                            end
                        end
                    end


                    % Plot electrodes
                    elec_xy = nan(length(l(:,end)),2);
                    for i4 = 1:length(l(:,end))
                        ind_pa = find(pa.ind == l(i4,1),1);
                        if (~isempty(ind_pa))
                            %plot(pa.x(ind_pa),pa.y(ind_pa),'.','Color',col_elec,'MarkerSize',size_elec)

%                             plot(pa.x(ind_pa),pa.y(ind_pa),'.','Color',...
%                                 col_elec_border,'MarkerSize',size_elec_hi*1.2)
%                             plot(pa.x(ind_pa),pa.y(ind_pa),'.','Color',...
%                                 col_elec_nosig,'MarkerSize',size_elec_hi)

                            % Save for future
                            elec_xy(l(i4,5),1) = pa.x(ind_pa);
                            elec_xy(l(i4,5),2) = pa.y(ind_pa);
                            % col_roi_hase
                        else
                            fprintf('[%s] skipped electrode: %i\n',sidL,l(i4,5));
                        end
                    end

                    view(90,90);
                    axis off;

    %                 roi_xy = cell(1,ecog.n_bchan);
                    for i2 = 1:length(roi_id)

                        c_i = find(roi_id(i2) == atl_table(:,end),1);

                        if (~isempty(c_i))
                            roi_name = C.AtlROIs{atl}.RH.struct_names{c_i};
                            if (strcmp(roi_name,roi_1_const))
                                tri_color(i2,:) = col_src;
                            elseif ( sum(strcmp(roi_name,C.AtlLabels{atl})) > 0)
                                %tri_color(i2,:) = atl_table(c_i,1:3)/255;
                                tri_color(i2,:) = col_roi_hase;
                            else
                                tri_color(i2,:) = col_surf;
                            end

    %                         roi_xy{c_i} = [roi_xy{c_i}; [pa.x(i2) pa.y(i2)]];
                        end
                    end

                    set(a,'FaceVertexCdata',tri_color); % color by ROI
                    %set(a,'FaceVertexCdata',tri_color*0.08 + 0.8*[1 1 1]); % color by trans ROI
                    %set(a,'FaceVertexCdata',col_surf) % color bland
                    daspect([1 1 1]); % fix aspect ratio to data
                    % -------------------------------------------------------------------------------

                    tri_color_val = zeros(length(roi_id),3);

                    % Plot 1 pair
                    elec_xy_sig = [];
                    elec_xy_sig_roi = {};
                    elec_xy_nosig = [];
                    elec_xy_nosig_roi = {};
                    elec_xy_src = [];
                    roi_1 = roi_1_const;
                    %roi_2 = 'inferiorparietal';
                    coh_t = [];
                    coh_t_sid = {};
                    countL = 1;
                    plotted_b1 = false;
                    for ii1 = 1:(ecog.n_bchan-1)
                        for ii2 = (ii1+1):ecog.n_bchan
                            %fprintf('chan1:%i ii1:%i chan2:%i ii2:%i\n',chan1(count)+1,ii1,chan2(count)+1,ii2)

                            has_roi_f = ( strcmp(atl_labels{ecog.bip(ii1,1)},roi_1) );
                            
%                             has_roi_f = ( strcmp(atl_labels{ecog.bip(ii1,1)},roi_1)  )...
%                                 & ( strcmp(atl_labels{ecog.bip(ii1,2)},roi_1)  );

                            %has_roi_r = (  strcmp(atl_labels{ecog.bip(ii2,1)},roi_1) )...
                            %        | (  strcmp(atl_labels{ecog.bip(ii2,1)},roi_1) );
    %                         has_roi_f = ( strcmp(atl_labels{ecog.bip(ii1,1)},roi_1)  )...
    %                                 | ( strcmp(atl_labels{ecog.bip(ii1,1)},roi_1)  )...
    %                                 | ( strcmp(atl_labels{ecog.bip(ii1,2)},roi_1)  )...
    %                                 | ( strcmp(atl_labels{ecog.bip(ii1,2)},roi_1)  );
    %                         has_roi_r = (  strcmp(atl_labels{ecog.bip(ii2,1)},roi_1) )...
    %                                 | (  strcmp(atl_labels{ecog.bip(ii2,2)},roi_1) )...
    %                                 | (  strcmp(atl_labels{ecog.bip(ii2,1)},roi_1) )...
    %                                 | (  strcmp(atl_labels{ecog.bip(ii2,2)},roi_1) );

                            %if (((has_roi_f || has_roi_r) && (Dmat(ii1,ii2) > dist_thresh)) && (ct(count) > ct_thresh))
                            %if (((has_roi_f || has_roi_r) && (Dmat(ii1,ii2) > dist_thresh)) && (ct(count) > ct_thresh))
                            
                            if (((has_roi_f) && (Dmat(ii1,ii2) > dist_thresh)))
                            %if (((has_roi_f) && (Dmat(ii1,ii2) > dist_thresh)) && (ct(count) > ct_thresh))

                                b1 = ecog.bip(ii1,1:2);
                                b2 = ecog.bip(ii2,1:2);

                                % Electrode marker size modulation
%                                 size_elec_hi_max = size_elec_hi * (1 + 1);
%                                 size_elec_hi_plt = size_elec_hi * (1 + 1*ct(count));

%                                 % Constant electrode size
%                                 size_elec_hi_max = size_elec_hi;
%                                 size_elec_hi_plt = size_elec_hi;
%                                 size_elec_border = size_elec_hi_plt*1.2;
% 
%                                 % set significant electrode color
%                                 elec_xy_src = [elec_xy_src; elec_xy(b1(:),:),repmat([col_elec_hi, size_elec_hi_max],2,1)];
%                                 if (ct(count) > ct_thresh)
%                                     if (trig_mag_not_ct)
%                                         col_elec_local = cmap(round(mag(count)*n_colors),:);
%                                         ct_range = [ct_range; mag(count); mag(count)]; % double for both bipolar electrodes
%                                     else
%                                         col_elec_local = cmap(round(ct(count)*n_colors),:);
%                                         ct_range = [ct_range; ct(count); ct(count)];
%                                     end
%                                     coh_t = [coh_t;coh_thresh(count)];
%                                     elec_xy_sig = [elec_xy_sig; [elec_xy(b2(:),:),repmat([col_elec_local, size_elec_hi_plt],2,1)]];
%                                     elec_xy_sig_roi = [elec_xy_sig_roi; {C.AtlLabels{2}{b2(1)}; C.AtlLabels{2}{b2(2)}}];
%                                 else
%                                     col_elec_local = col_elec_nosig;
%                                     elec_xy_nosig = [elec_xy_nosig; [elec_xy(b2(:),:),repmat([col_elec_local, size_elec_hi],2,1)]];
%                                     elec_xy_nosig_roi = [elec_xy_nosig_roi; {C.AtlLabels{2}{b2(1)}; C.AtlLabels{2}{b2(2)}}];
%                                 end

                                elec_xy_src = [elec_xy_src; elec_xy(b1(1),:), col_elec_hi];
                                if (ctL(countL) > ct_threshL)
                                    if (trig_mag_not_ct)
                                        col_elec_local = cmapL(round(mag(countL)*n_colors),:);
                                        ct_range = [ct_range; mag(countL)];
                                    else
                                        col_elec_local = cmapL(round(ctL(countL)*n_colors),:);
                                        ct_range = [ct_range; ctL(countL)];
                                    end
                                    coh_t = [coh_t;Ca.coh_thresh(countL)];
                                    coh_t_sid = [coh_t_sid; Subjects{iL}];
                                    elec_xy_sig = [elec_xy_sig; [elec_xy(b2(1),:), col_elec_local ]];
                                    elec_xy_sig_roi = [elec_xy_sig_roi; C.AtlLabels{2}{b2(1)} ];
                                else
                                    col_elec_local = col_elec_nosig;
                                    elec_xy_nosig = [elec_xy_nosig; [elec_xy(b2(1),:), col_elec_local ]];
                                    elec_xy_nosig_roi = [elec_xy_nosig_roi; C.AtlLabels{2}{b2(1)} ];
                                end
                                

%                                 if (~plotted_b1)
%                                     for jj1 = 1:length(b1)
%                                         plot(elec_xy(b1(jj1),1),elec_xy(b1(jj1),2),'.','Color',...
%                                             col_elec_border,'MarkerSize',size_elec_hi*1.2)
%                                         plot(elec_xy(b1(jj1),1),elec_xy(b1(jj1),2),'.','Color',...
%                                             col_elec_hi,'MarkerSize',size_elec_hi)
%                                     end
%                                 end
%                                 for jj2 = 1:length(b2)
%                                     plot(elec_xy(b2(jj2),1),elec_xy(b2(jj2),2),'.','Color',...
%                                         col_elec_border,'MarkerSize',size_elec_border)
%                                     plot(elec_xy(b2(jj2),1),elec_xy(b2(jj2),2),'.','Color',...
%                                         col_elec_local,'MarkerSize',size_elec_hi_plt)
%                                 end
                                % comment below to draw b1 multiple times
                                %plotted_b1 = true;

                            end
                            countL = countL + 1;
                        end
                    end
                    
                    % show coherence threshold stats
                    fprintf('Coherence threshold mean: %.3f\n',mean(coh_t));
                    fprintf('Coherence threshold stdev: %.3f\n',std(coh_t));
                    fprintf('Coherence threshold median: %.3f\n',median(coh_t));
                    fprintf('Number of significant bchan pairs: %i\n',length(coh_t));
                    
                    % Plot electrodes
                    border_fac = 1.2;
%                     elec_xy_nosig = unique(elec_xy_nosig,'rows');
%                     elec_xy_sig = unique(elec_xy_sig,'rows');
                    
%                     [elec_xy_nosig,uidx] = unique(elec_xy_nosig,'rows');
%                     elec_xy_nosig_roi = elec_xy_nosig_roi(uidx);
% 
%                     [elec_xy_sig,uidx] = unique(elec_xy_sig,'rows');
%                     elec_xy_sig_roi = elec_xy_sig_roi(uidx);
%                     
%                     elec_xy_src = unique(elec_xy_src,'rows');
                    
                    % Plot order
                    %elec_xy_all = [elec_xy_nosig; elec_xy_sig; elec_xy_src];
                    elec_xy_all = [elec_xy_nosig; elec_xy_sig];
                    if (~ isempty(elec_xy_all))
                    
                        % --- COPY PASTA ---
                        % Plot electrodes
                        %elec_xy_nosig = unique(elec_xy_nosig,'rows');
                        %elec_xy_sig = unique(elec_xy_sig,'rows');
                        %elec_xy_src = unique(elec_xy_src,'rows');
                        % Plot order
                        %elec_xy_all = [elec_xy_nosig; elec_xy_sig; elec_xy_src];
                        %elec_xy_all = [elec_xy_nosig; elec_xy_src]; % modified from copy pasta
                        elec_xy_all = [elec_xy_nosig];
                        [n_pe,~] = size(elec_xy_all);
                        for ip = 1:n_pe
                            % TODO: [x_remap, y_remap] = electrode_remap(x,y,rois,roi_1_const,roi_midpt,roi_id,c,pa)
                            x = elec_xy_all(ip,1);
                            y = elec_xy_all(ip,2);
                            [x_remap, y_remap] = electrode_remap(x,y,rois,elec_xy_nosig_roi(ip),roi_midpt,roi_id,c,pa);
                            plot(x_remap,y_remap,'.','Color',...
                                col_elec_border,'MarkerSize',size_elec_hi*border_fac)
                            plot(x_remap,y_remap,'.','Color',...
                                elec_xy_all(ip,3:5),'MarkerSize',size_elec_hi)
                        end

                        %return

                        [n_pe,~] = size(elec_xy_sig);
                        ct_abs = max(ct_range) - min(ct_range);
                        elec_xy_all = elec_xy_sig;
                        for ip = 1:n_pe
                            % TODO: [x_remap, y_remap] = electrode_remap(x,y,rois,roi_1_const,roi_midpt,roi_id,c,pa)
                            if (ct_abs == 0)
                                ct_norm = 1;
                            else
                                ct_norm = (ct_range(ip) - min(ct_range) )/ct_abs;
                            end
                            %fprintf('ct_norm: %.4f\n',ct_norm);
                            col_ct = cmapL(ceil((n_colors-1)*ct_norm)+1,:);

                            x = elec_xy_all(ip,1);
                            y = elec_xy_all(ip,2);
                            [x_remap, y_remap] = electrode_remap(x,y,rois,elec_xy_sig_roi(ip),roi_midpt,roi_id,c,pa);
                            plot(x_remap,y_remap,'.','Color',...
                                col_elec_border,'MarkerSize',size_elec_hi*border_fac)
                            plot(x_remap,y_remap,'.','Color',...
                                col_ct,'MarkerSize',size_elec_hi)
                        end
                        
                        elec_xy_all = [elec_xy_src];
                        [n_pe,~] = size(elec_xy_all);
                        for ip = 1:n_pe
                            % TODO: [x_remap, y_remap] = electrode_remap(x,y,rois,roi_1_const,roi_midpt,roi_id,c,pa)
                            x = elec_xy_all(ip,1);
                            y = elec_xy_all(ip,2);
                            [x_remap, y_remap] = electrode_remap(x,y,rois,roi_1_const,roi_midpt,roi_id,c,pa);
                            plot(x_remap, y_remap,'.','Color',...
                                col_elec_border,'MarkerSize',size_elec_hi*border_fac)
                            plot(x_remap, y_remap,'.','Color',...
                                elec_xy_all(ip,3:5),'MarkerSize',size_elec_hi)
                        end
                        % --- end COPY PASTA ---
                        
%                         % Filter electrodes to those in roi, dist thresh
%                         [n_pe,n_pe_u_cols] = size(elec_xy_all);
%                         elec_xy_u = unique(elec_xy_all(:,1:2),'rows');
%                         [n_pe_u,~] = size(elec_xy_u);
%                         elec_xy_all_c = zeros(n_pe_u,n_pe_u_cols);
% 
%                         for ip = 1:n_pe_u
%                             to_c = elec_xy_all(elec_xy_u(ip,1) == elec_xy_all(:,1),:); % matrix to collapse
%                             if (~ isempty(to_c))
%                                 % Collapsing function
%                                 cp = sum(to_c(:,6) ~= size_elec)/length(to_c(:,6));
%                                 if (cp > cp_thresh)
%                                     elec_xy_all_c(ip,:) = mean(to_c(to_c(:,6) ~= size_elec,:),1);
%                                 else
%                                     elec_xy_all_c(ip,:) = mean(to_c(to_c(:,6) == size_elec,:),1);
%                                 end
%                                 %elec_xy_all_c(ip,:) = mean(to_c,1);
%                             else
%                                 fprintf(2,'E: electrodes to collapse is empty for electrode: %i\n',ip)
%                             end
%                         end
% 
%                         elec_xy_all_c = [elec_xy_all_c; elec_xy_src];
%                         [n_pe_u,~] = size(elec_xy_all_c);
%                         for ip = 1:n_pe_u
%                             plot(elec_xy_all_c(ip,1),elec_xy_all_c(ip,2),'.','Color',...
%                                 col_elec_border,'MarkerSize',elec_xy_all_c(ip,6)*border_fac)
%                             plot(elec_xy_all_c(ip,1),elec_xy_all_c(ip,2),'.','Color',...
%                                 elec_xy_all_c(ip,3:5),'MarkerSize',elec_xy_all_c(ip,6))
%                         end

                        
                        % Plot scale
                        xrange = max(pa.x) - min(pa.x);
                        yrange = max(pa.y) - min(pa.y);
                        scale_fac = 10; % mm
                        bar_horiz_x = [min(pa.x), min(pa.x) + scale_fac];
                        bar_horiz_y = [min(pa.y), min(pa.y)];
                        bar_vert_x = [min(pa.x), min(pa.x)];
                        bar_vert_y = [min(pa.y), min(pa.y) + scale_fac];
                        plot(bar_horiz_x,bar_horiz_y,'black-','LineWidth',5);
                        text(min(pa.x)+0.5*scale_fac,min(pa.y)+0.5*scale_fac,'1 cm','FontSize',8);

%                         if (trig_mag_not_ct)
%                             magoct = 'mag';
%                         else
%                             magoct = 'ct';
%                         end
%                         
                        %set(a,'FaceVertexCdata',tri_color_val);
%                         colormap(cmap);
%                         caxis([0 1]);
%                         colorbar;
%                         colormap(cmap);
%                         caxis([min(ct_range) max(ct_range)]);
%                         colorbar;

                        colormap(cmapL);
                        if (~isempty(ct_range) && ((max(ct_range)-min(ct_range)) ~= 0) )
                            caxis([min(ct_range) max(ct_range)]);
                            cbh = colorbar('XTick', linspace(min(ct_range),max(ct_range),5));
                        else
                            colorbar;
                        end

                        %return

                        if (trig_mag_not_ct)
                            magoct = 'mag';
                        else
                            magoct = 'ct';
                        end
                        print(h,sprintf('figures/T11d2/%s_%s_%s',pname,hemi{jj},magoct),fig_fmt);
                        if (trig_eps)
                            print(h,sprintf('figures/T11d2/%s_%s_%s',pname,hemi{jj},magoct),'-depsc');
                        end
                        
                    end
                    
                    
                    close(h);
                    
                    % reverse labels
                    % -------------------------------------------------------------------------------
                end
                
            end
        end

    end

end

% Clear loop indices
clear i;
clear j;

% Print finish message
fprintf('Done.\n')

end
