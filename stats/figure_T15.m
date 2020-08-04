close all;
clear;


load('rois');
rois = rois(2:end);
% 
% pp = gcp('nocreate');
% 
% delete(pp);
% pp = parpool(2);

% ROI remapping
n_remap_steps = 1000; % set to 0 to skip, default 1000
step_size = 1/n_remap_steps;

trig_debug = false;
trig_roi_txt = true;

%for i_roi = 1:length(rois)

    
%metric = 'pcBroadband';
n_perm = 10000;

% anatomy folder
%subjects_dir = '/mnt/cuenap_ssd/coregistration';
subjects_dir = '/home/klab/data/h5eeg/artifact_removed/test/opencl/coreg';
setenv('SUBJECTS_DIR',subjects_dir);

% Fast i/o definitions
dir_art = '/media/klab/internal/data/h5_notch20/art';
dir_res = '/media/klab/internal/data/results/coh_w10';
dir_cor = '/media/klab/internal/data/coreg';
dir_cacheL = './cache';

% Slow i/o definitions
dir_h5 = '/media/klab/KLAB101/h5_notch20';

metrics = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};

% Patients
% TODO: figure out where Subjects is being loaded
SubjectsL = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044'         ,'m00047',... % ,'m00045'
    'm00048','m00049','m00052','m00053','m00055',         'm00058','m00059',... % ,'m00056'
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124',...
    'mSu'};

% SubjectsL = {'m00001','mSu'};

% Exclude monkey
SubjectsL = SubjectsL(1:(end-1));


% ROI
%roi_1_const = 'lateraloccipital';
%roi_1_const = 'inferiorparietal';
roi_1_const = '';
%roi_1_const = 'middletemporal';
%roi_1_const = rois{i_roi};

fig_fmt = '-djpeg';
trig_eps = true;
trig_mag_not_ct = true;
%trig_skip_plot = false;
trig_overwrite_cache = false;
system('mkdir figures');
system('mkdir figures/T15');

fn_cache = sprintf('cache/fsaverage_sym_flat.mat');
fn_patch = sprintf('%s/%s/surf/%s.full.flat.patch.3d',subjects_dir,'fsaverage_sym','lh');
% fn_curv = sprintf('%s/%s/surf/%s.curv',subjects_dir,'fsaverage_sym','lh');
% fn_sulc = sprintf('%s/%s/surf/%s.sulc',subjects_dir,'fsaverage_sym','lh');
fn_pial = sprintf('%s/%s/surf/%s.pial',subjects_dir,'fsaverage_sym','lh');

n_colors = 100;
cmap = corrcmap(n_colors);
            
%for iM = 1:length(metrics)
for iM = 1
    metric = metrics{iM};
    %fn_cache = sprintf('cache/%s_%s_flat.mat','fsaverage_sym',hemi{jj});
        
    T = load(fn_cache,'pa','tri_flat','pa','v_pial','tri_pial');
    pa = T.pa;
    tri_flat = T.tri_flat;
    v_pial = T.v_pial;
    tri_pial = T.tri_pial;
    
    covered_rois = {};
    elec_xy_sig = [];
    elec_xy_nosig = [];
    elec_xy_nosig_roi = {};
    elec_xy_src = [];
    % Main loop
    pname = sprintf('%s_%s-all_%s','fsaverage_sym',roi_1_const,metric);
            
    
    
    annot_dir = sprintf('%s/%s/label',subjects_dir,'fsaverage_sym');
    switch('Desikan-Killiany')
        case ('Desikan-Killiany')
            atl_fname = sprintf('%s/%s.aparc.annot',annot_dir,'lh');
        case ('Destrieux-2009')
            atl_fname = sprintf('%s/%s.aparc.a2009s.annot',annot_dir,'lh');
        case ('HCP-MMP1')
            atl_fname = sprintf('%s/%s.HCP-MMP1.annot',annot_dir,'lh');
        case ('M132')
            atl_fname = sprintf('%s/%s.MACAQUE_M132.annot',annot_dir,'lh');
    end
    [v,l2,c] = read_annotation(atl_fname);
    roi_id = l2(pa.vno + 1);
    
    Ca = load(sprintf('%s/xsub_out_%s_%i.mat',dir_cacheL,SubjectsL{1},iM));
    ecog = Ca.ecog;
    C1 = Ca.C;
    
    % ROI midpoints for fsaverage_sym
    % collect all points for each ROI
    roi_pts = cell(length(rois),1);
    for iv = 1:length(roi_id)
        roi_idx = find(roi_id(iv) == C1.AtlROIs{2}.LH.table(:,end));
        if (~isempty(roi_idx))
            roi = C1.AtlROIs{2}.LH.struct_names{roi_idx};
            if (~strcmp(roi,'unknown'))
                idx2 = find(strcmp(rois,roi));
                roi_pts{idx2} = [roi_pts{idx2}, [pa.x(iv); pa.y(iv)]];
            end
        else
            fprintf(2,'W: roi identifying number not found %i\n',roi_id(iv));
        end
    end
   
    % get midpoint
    roi_midpt = nan(length(rois),2);
    for iv = 1:length(rois)
        if (~isempty(roi_pts{iv}))
            roi_midpt(iv,:) = nanmean(roi_pts{iv},2)';
        end
    end
    
    % ParC - parallel cache
    % 1 - elec_xy_nosig 
    % 2 - elec_xy_nosig_roi
    % 3 - covered_rois
    
    ParC = cell(3,length(SubjectsL));
    
    for i = 1:length(SubjectsL)
        
        % ParC local vars
        elec_xy_nosig_parc = [];
        elec_xy_nosig_roi_parc = [];
        covered_rois_parc = [];
        
        sid = SubjectsL{i};
        fprintf('%s> Started.\n',sid)
        %sid = 'fsaverage_sym';
        fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
        fn_dist = sprintf('%s/%s_dists-%s-%i.mat',dir_res,sid,metric,n_perm);
        fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
        fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
        fn_coreg = sprintf('%s/%s/label/all_parcellation.mat',dir_cor,sid);

        % Check if files exist
        ckf = {fn_art,fn_dist,fn_graph,fn_h5,fn_coreg};
        for jj = 1:length(ckf)
            if (~exist(ckf{jj},'file'))
                fprintf(2,'E> File not found: %s\n',ckf{jj});
                %return
            end
        end
%         if (isempty(hemi))
%             fprintf('W> file not found: %s\n',sprintf('%s/%s/surf/*h.full.flat.patch.3d',subjects_dir,sid));
%         else
%     for jj = 1:length(hemi)
        % Read in file
        
%             fn_patch_asc = sprintf('%s/%s/surf/%s.full.flat.patch.3d.asc',subjects_dir,sid,hemi{j});
% 
%             % if ascii doesn't exist, make it
%             if (~ exist(fn_patch_asc, 'file'))
%                 fprintf('[*] ASCII flat patch not found, converting...\n')
%                 system(sprintf('cd %s/%s/surf; mris_convert -p %s.full.flat.patch.3d %s.full.flat.patch.3d.asc',subjects_dir,sid,hemi{j},hemi{j}));
%             end


        
        % Load human cache
        Ca = load(sprintf('%s/xsub_out_%s_%i.mat',dir_cacheL,sid,iM));
        ecog = Ca.ecog;
        C = Ca.C;
        atl_labels = Ca.atl_labels;
        atl = Ca.atl;
        Dmat = Ca.Dmat;
        mag = Ca.mag;
        ct = Ca.ct;
        ct_thresh = Ca.ct_thresh;
        cp_thresh = Ca.cp_thresh;
        dist_thresh = Ca.dist_thresh;
        
        chan_labels = h5readatt(fn_h5,'/h5eeg/eeg','labels');
        bchan_labels = cell(ecog.n_bchan,1);
        for ii0 = 1:ecog.n_bchan
            bchan_labels{ii0} = sprintf('%s-%s',chan_labels{ecog.bip(ii0,1)},chan_labels{ecog.bip(ii0,2)});
        end


        % Plot flat map

        % check hemisphere
%         hemi = {};
%         if (exist(sprintf('%s/%s/surf/rh.full.flat.patch.3d',subjects_dir,sid), 'file'))
%             hemi = {hemi{:},'rh'};
%         end
%         if (exist(sprintf('%s/%s/surf/lh.full.flat.patch.3d',subjects_dir,sid), 'file'))
%             hemi = {hemi{:},'lh'};
%         end
        
        %for jj = 1%:length(hemi)
        %jj = 1;
        % Read annotation
        %atl_table = AT.P.AtlROIs{atl_i}.RH.table;
        % only plot base patch once

        % read electrode locations
        fprintf('[!] Ignore any errors:\n')
        l = read_label(sprintf('%s/fsaverage_sym',subjects_dir),sprintf('ielvis_%s',sid));
        if (isempty(l))
            l = read_label('fsaverage_sym',sprintf('ielvis_%s',sid));
        end

        %if (i == 1)

        atl_table = C.AtlROIs{atl}.RH.table;
        

        % constant color
        col_unknown = 0.333*[1 1 1];
        col_src = 0.5*[1 1 1];
        col_surf = 0.85*[1 1 1];
        col_roi_hase = 0.8*[1 1 1];
        col_bound = 0.4*[1 1 1];
        col_elec = 0*[1 1 1];
        col_elec_hi = 0*[1 1 1];
        col_elec_border = [0 0 0];
        border_fac = 1.2; % how much larger to draw border
        col_elec_nosig = 0.3*[1 1 1];
        size_bound = 3; %5
        bound_repeat = 1;
        bound_loop_rep = 1;
        size_elec = 10;
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



    % =======================================================================================


        % --- Different for each patient ---
        % Plot electrodes
        elec_xy = nan(length(l(:,end)),2);
        for i4 = 1:length(l(:,end))
            ind_pa = find(pa.ind == l(i4,1),1);
            if (~isempty(ind_pa))
                %plot(pa.x(ind_pa),pa.y(ind_pa),'.','Color',col_elec,'MarkerSize',size_elec)

%                 plot(pa.x(ind_pa),pa.y(ind_pa),'.','Color',...
%                     col_elec_border,'MarkerSize',size_elec_hi*1.2)
%                 plot(pa.x(ind_pa),pa.y(ind_pa),'.','Color',...
%                     col_elec_nosig,'MarkerSize',size_elec_hi)

                % Save for future
                elec_xy(l(i4,5),1) = pa.x(ind_pa);
                elec_xy(l(i4,5),2) = pa.y(ind_pa);
                % col_roi_hase
            else
                fprintf('[%s] did not map electrode %i to flat map\n',sid,l(i4,5));
            end
        end


        % -------------------------------------------------------------------------------

        tri_color_val = zeros(length(roi_id),3);

        % Plot 1 pair
        %roi_1 = roi_1_const;
        %roi_2 = 'inferiorparietal';
        count = 1;
        plotted_b1 = false;
        for ii1 = 1:(ecog.n_bchan-1)
            for ii2 = (ii1+1):ecog.n_bchan
                %fprintf('chan1:%i ii1:%i chan2:%i ii2:%i\n',chan1(count)+1,ii1,chan2(count)+1,ii2)
% 
%                 has_roi_f = ( strcmp(atl_labels{ecog.bip(ii1,1)},roi_1)  )...
%                         & ( strcmp(atl_labels{ecog.bip(ii1,2)},roi_1)  );

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

%                 if (((has_roi_f) && (Dmat(ii1,ii2) > dist_thresh)))
                if ( (Dmat(ii1,ii2) > dist_thresh) )

                    b1 = ecog.bip(ii1,1:2);
                    b2 = ecog.bip(ii2,1:2);
                    
                    % Plot electrode interactivity
                    size_elec_hi_max = size_elec_hi * (1 + 1);
                    size_elec_hi_plt = size_elec_hi * (1 + 1*ct(count));
                    size_elec_border = size_elec_hi_plt*1.2;
                    
                    % set significant electrode color
%                     elec_xy_src = [elec_xy_src; elec_xy(b1(:),:),repmat([col_elec_hi, size_elec_hi_max],2,1)];
%                     if (ct(count) > ct_thresh)
%                         if (trig_mag_not_ct)
%                             col_elec_local = cmap(round(mag(count)*n_colors),:);
%                         else
%                             col_elec_local = cmap(round(ct(count)*n_colors),:);
%                         end
%                         elec_xy_sig = [elec_xy_sig; [elec_xy(b2(:),:),repmat([col_elec_local, size_elec_hi_plt],2,1)]];
%                     else
                    col_elec_local = col_elec_nosig;
                    %elec_xy_nosig_parc = [elec_xy_nosig_parc; [elec_xy(b2(:),:),repmat([col_elec_local, size_elec_hi],2,1)]];
                    elec_xy_nosig = [elec_xy_nosig; [elec_xy(b2(:),:),repmat([col_elec_local, size_elec_hi],2,1)]];
                    
%                     end
                    
                    % add rois to covered rois
                    %covered_rois_parc = unique([covered_rois_parc;C.AtlLabels{atl}]);
                    covered_rois = unique([covered_rois;C.AtlLabels{atl}]);
                    
                    % roi
%                     elec_xy_nosig_roi = [elec_xy_nosig_roi; C.AtlLabels{2}{b1(1)}];
%                     elec_xy_nosig_roi = [elec_xy_nosig_roi; C.AtlLabels{2}{b2(1)}];
                    %elec_xy_nosig_roi_parc = [elec_xy_nosig_roi_parc; [C.AtlLabels{2}{b1(1)}; C.AtlLabels{2}{b2(1)}]];
                    %elec_xy_nosig_roi = [elec_xy_nosig_roi; [C.AtlLabels{2}{b1(1)}; C.AtlLabels{2}{b2(1)}]];
                    elec_xy_nosig_roi = [elec_xy_nosig_roi; [C.AtlLabels{2}{b2(1)}; C.AtlLabels{2}{b2(2)}]];
                    
                    
%                             if (has_roi_r)
%                                 b2t = b2;
%                                 b2 = b1;
%                                 b1 = b2t;
%                             end

                    %return

%                     % Plot electrode interactivity
%                     if (~plotted_b1)
%                         % source electrode
%                         for jj1 = 1:length(b1)
%                             plot(elec_xy(b1(jj1),1),elec_xy(b1(jj1),2),'.','Color',...
%                                 col_elec_border,'MarkerSize',size_elec_hi*border_fac)
%                             plot(elec_xy(b1(jj1),1),elec_xy(b1(jj1),2),'.','Color',...
%                                 col_elec_hi,'MarkerSize',size_elec_hi)
%                         end
%                     end
%                     % real electrode
%                     for jj2 = 1:length(b2)
%                         plot(elec_xy(b2(jj2),1),elec_xy(b2(jj2),2),'.','Color',...
%                             col_elec_border,'MarkerSize',size_elec_hi*border_fac)
%                         plot(elec_xy(b2(jj2),1),elec_xy(b2(jj2),2),'.','Color',...
%                             col_elec_local,'MarkerSize',size_elec_hi)
% 
%                     end
                    % comment below to draw b1 multiple times
                    %plotted_b1 = true;


                end
                count = count + 1;
            end
        end




        %return
        % =========== on last subject ====================================
        if (i == length(SubjectsL))
            
            
            % Plot flat patch
            % -------------------------------------------------------------------------------
            h = figure;
            fig_w = 10.5;
            fig_h = 10.5;
            set(h,'Position',round([0 0 fig_w*100 fig_h*100]))
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
                vB = tri_flat(i3,:);
                tv = roi_id(vB);
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

            view(90,90);
            axis off;
            daspect([1 1 1]); % fix aspect ratio to data
            
            % --- ROI colors --------------------------------------------
            % PAR
            atl_tble = atl_table(:,end);
            ele_rois = covered_rois;
            for i2 = 1:length(roi_id)
                c_i = find(roi_id(i2) == atl_tble,1);
                if (~isempty(c_i))
%                     roi_name = C.AtlROIs{atl}.RH.struct_names{c_i};
%                     if ( strcmp(roi_1_const,roi_name) )
%                         tri_color(i2,:) = col_src;
%                     elseif ( ~ isempty(find(strcmp(roi_name,ele_rois),1)) )
                    tri_color(i2,:) = col_roi_hase;
%                     end
                end
            end
            set(a,'FaceVertexCdata',tri_color); % color by ROI
            %set(a,'FaceVertexCdata',tri_color*0.08 + 0.8*[1 1 1]); % color by trans ROI
            %set(a,'FaceVertexCdata',col_surf) % color bland
            
            
            
            
            % Plot electrodes
            [elec_xy_nosig,uidx] = unique(elec_xy_nosig,'rows');
            elec_xy_nosig_roi = elec_xy_nosig_roi(uidx);
            elec_xy_sig = unique(elec_xy_sig,'rows');
            elec_xy_src = unique(elec_xy_src,'rows');
            
            % Plot order
            %elec_xy_all = [elec_xy_nosig; elec_xy_sig; elec_xy_src];
            %elec_xy_all = [elec_xy_nosig; elec_xy_sig];
            elec_xy_all = elec_xy_nosig;
            elec_xy_all_roi = elec_xy_nosig_roi;
            if (isempty(elec_xy_all))
                %break;
            end
            
            
            
%             % Filter electrodes to those in roi, dist thresh
%             [n_pe,n_pe_u_cols] = size(elec_xy_all);
%             elec_xy_u = unique(elec_xy_all(:,1:2),'rows');
%             [n_pe_u,~] = size(elec_xy_u);
%             elec_xy_all_c = zeros(n_pe_u,n_pe_u_cols);
%             %elec_xy_all_c_roi = cell(n_pe_u,1);
% 
%             for ip = 1:n_pe_u
%                 to_c = elec_xy_all(elec_xy_u(ip,1) == elec_xy_all(:,1),:); % matrix to collapse
%                 if (~ isempty(to_c))
%                     % Collapsing function
%                     cp = sum(to_c(:,6) ~= size_elec)/length(to_c(:,6));
%                     if (cp > cp_thresh)
%                         elec_xy_all_c(ip,:) = mean(to_c(to_c(:,6) ~= size_elec,:),1);
%                     else
%                         elec_xy_all_c(ip,:) = mean(to_c(to_c(:,6) == size_elec,:),1);
%                     end
%                     
%                     
%                     %elec_xy_all_c_roi{ip} = ;
%                     %elec_xy_all_c(ip,:) = mean(to_c,1);
%                 else
%                     fprintf(2,'W: electrodes to collapse is empty for electrode: %i, setting to NaN.\n',ip)
%                     elec_xy_all_c(ip,:) = nan(1,n_pe_u_cols);
%                 end
%             end

            % Sort electrodes by significance
            %elec_xy_all_c = [elec_xy_all_c; elec_xy_src];
            %elec_xy_all_c_sav = elec_xy_all_c;
            elec_xy_all_c = elec_xy_all;
            elec_xy_all_roi_c = elec_xy_all_roi;
%             skey = all(diff(elec_xy_all_c(:,3:5)) == 0,2);
%             [~,sidx] = sort(skey,'descend');
%             elec_xy_all_c = elec_xy_all_c(sidx,:);
            [n_pe_u,~] = size(elec_xy_all_c);
            for ip = 1:n_pe_u
                
                roi_1_const = elec_xy_all_roi_c{ip};
                
                if (~isnan(elec_xy_all_c(ip,1)))
                    
                    %remap electrode back into roi
                    x_remap = elec_xy_all_c(ip,1);
                    y_remap = elec_xy_all_c(ip,2);
                    
                    % Get roi at current point
                    diff = (pa.x - x_remap).^2 + (pa.y - y_remap).^2;
                    [~,min_idx] = min(diff);
                    roi_id(min_idx);
                    elec_roi = c.struct_names{c.table(:,5) == roi_id(min_idx)};
                    elec_xy_roi_midpt = roi_midpt(strcmp(rois,roi_1_const),:);
                    
                    if (~isempty(elec_xy_roi_midpt))
                        % debug
    %                     if (trig_debug)
    %                         fprintf('roi: %s, midpt: %.3f %.3f\n',elec_roi,elec_xy_roi_midpt(1),elec_xy_roi_midpt(2))
    %                         plot(elec_xy_roi_midpt(1),elec_xy_roi_midpt(2),'red+');
    %                         hold on;
    %                     end
    
                        x_remap0 = x_remap;
                        y_remap0 = y_remap;
                        if (trig_debug)
                            plot(x_remap0,y_remap0,'blueo');
                            text(x_remap0,y_remap0,roi_1_const,'fontsize',5,'color',0.2*[1 1 1]);
                        end

                        t_tmp = 0;
                        remap_break_prob = 0.25;
                        if (~ strcmp(roi_1_const,elec_roi))
                            for ir = 1:n_remap_steps
                                t_tmp = t_tmp + step_size;
                                x_remap = x_remap*(1-t_tmp) + t_tmp * elec_xy_roi_midpt(1);
                                y_remap = y_remap*(1-t_tmp) + t_tmp * elec_xy_roi_midpt(2);

                                diff = (pa.x - x_remap).^2 + (pa.y - y_remap).^2;
                                [~,min_idx] = min(diff);
                                roi_id(min_idx);
                                elec_roi = c.struct_names{c.table(:,5) == roi_id(min_idx)};
                                if strcmp(roi_1_const,elec_roi)
                                    to_break = ( rand([1,1]) < remap_break_prob );
                                    if (to_break)
                                        break;
                                    end
                                end

                            end
                        end
                        if (~strcmp(elec_roi,'insula'))
                            plot(x_remap,y_remap,'.','Color',...
                                col_elec_border,'MarkerSize',elec_xy_all_c(ip,6)*border_fac)
                            plot(x_remap,y_remap,'.','Color',...
                                elec_xy_all_c(ip,3:5),'MarkerSize',elec_xy_all_c(ip,6))
                        end
                        if (trig_debug)
                            plot(x_remap,y_remap,'blue.');
                            plot([x_remap0 x_remap],[y_remap0 y_remap],'blue--')
                            %text(x_remap,y_remap,elec_roi);
                        end
                    end
%                     plot(elec_xy_all_c(ip,1),elec_xy_all_c(ip,2),'.','Color',...
%                         col_elec_border,'MarkerSize',elec_xy_all_c(ip,6)*border_fac)
%                     plot(elec_xy_all_c(ip,1),elec_xy_all_c(ip,2),'.','Color',...
%                         elec_xy_all_c(ip,3:5),'MarkerSize',elec_xy_all_c(ip,6))
                end
            end
            
%             % Plot source electrodes separately
%             elec_xy_all_c = [elec_xy_src];
%             [n_pe_u,~] = size(elec_xy_all_c);
%             elec_xy_roi_midpt = nanmean(elec_xy_all_c(:,1:2),1);
%             n_remap_steps = 1000;
%             step_size = 1/n_remap_steps;
%             
%             for ip = 1:n_pe_u
%                 if (~isnan(elec_xy_all_c(ip,1)))
%                     
%                     % map electrode onto roi
%                     diff = (pa.x - elec_xy_all_c(ip,1)).^2 + (pa.y - elec_xy_all_c(ip,2)).^2;
%                     [~,min_idx] = min(diff);
%                     roi_id(min_idx);
%                     elec_roi = c.struct_names{c.table(:,5) == roi_id(min_idx)};
%                     
%                     %remap electrode back into roi
%                     x_remap = elec_xy_all_c(ip,1);
%                     y_remap = elec_xy_all_c(ip,2);
%                     t_tmp = 0;
%                     if (~ strcmp(roi_1_const,elec_roi))
%                         for ir = 1:n_remap_steps
%                             t_tmp = t_tmp + step_size;
%                             x_remap = x_remap*(1-t_tmp) + step_size * elec_xy_roi_midpt(1);
%                             y_remap = y_remap*(1-t_tmp) + step_size * elec_xy_roi_midpt(2);
%                             
%                             
%                             diff = (pa.x - x_remap).^2 + (pa.y - y_remap).^2;
%                             [~,min_idx] = min(diff);
%                             roi_id(min_idx);
%                             elec_roi = c.struct_names{c.table(:,5) == roi_id(min_idx)};
%                             if strcmp(roi_1_const,elec_roi)
%                                 break;
%                             end
%                             
%                         end
%                     end
%                     
%                     plot(x_remap,y_remap,'.','Color',...
%                         col_elec_border,'MarkerSize',elec_xy_all_c(ip,6)*border_fac)
%                     plot(x_remap,y_remap,'.','Color',...
%                         elec_xy_all_c(ip,3:5),'MarkerSize',elec_xy_all_c(ip,6))
%                 end
%             end



            if (trig_roi_txt)
                for iv = 1:length(rois)
                    text(roi_midpt(iv,1),roi_midpt(iv,2),rois{iv},...
                        'HorizontalAlignment','center','fontsize',6);%,'BackgroundColor',col_roi_hase);
                end
            end


            %return
%             [n_pe,~] = size(elec_xy_all);
%             for ip = 1:n_pe
%                 plot(elec_xy_all(ip,1),elec_xy_all(ip,2),'.','Color',...
%                     col_elec_border,'MarkerSize',size_elec_hi*border_fac)
%                 plot(elec_xy_all(ip,1),elec_xy_all(ip,2),'.','Color',...
%                     elec_xy_all(ip,3:5),'MarkerSize',size_elec_hi)
%             end
            
            
            % - colormap
            %colormap(cmap);
            %caxis([0 1]);
            %colorbar;
            
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
            %plot(bar_vert_x,bar_vert_y,'black-','LineWidth',5);
            
            %return
            
            % Save
            if (trig_mag_not_ct)
                magoct = 'mag';
            else
                magoct = 'ct';
            end
            print(h,sprintf('figures/T15/%s_%s_%s',pname,'sym',magoct),fig_fmt);
            if (trig_eps)
                print(h,sprintf('figures/T15/%s_%s_%s',pname,'sym',magoct),'-depsc');
            end
            
            close(h);
        end
        % reverse labels
        % -------------------------------------------------------------------------------
        %end
        %end

    end

end

% Clear loop indices
% clear i;
% clear j;

% Print finish message
fprintf('Done.\n')


