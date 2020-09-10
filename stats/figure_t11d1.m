close all;
clear;

metricM = 'pcBroadband';
n_perm = 10000;

% anatomy folder
%subjects_dir = '/mnt/cuenap_ssd/coregistration';
subjects_dir = '/media/klab/internal/data/coreg';
%subjects_dir = '/home/klab/data/h5eeg/artifact_removed/test/opencl/coreg';
setenv('SUBJECTS_DIR',subjects_dir);

% Fast i/o definitions
dir_artM = '/media/klab/KLAB101/h5_notch2/art_nosz';%'/media/klab/internal/data/h5_notch20/art';
dir_resM = '/media/klab/KLAB101/results/coh_w10';
dir_corM = '/media/klab/internal/data/coreg';
dir_cacheLM = './cache';

% Slow i/o definitions
dir_h5M = '/media/klab/KLAB101/h5_notch20';

metricsM = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};

% Patients
SubjectsM = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48',...
    'mSu'};

% Exclude monkey
SubjectsM = SubjectsM(1:(end-1));

% Shuffle
%SubjectsM = SubjectsM(randperm(length(SubjectsM)));

SubjectsM = {'sub41'};
%r_samp_const = 106865001; % set to nan to pick randomly
%r_samp_const = 90545001; % strong t1, t9 good but next to art
r_samp_const = 77837501; % good t9, magnitude of t1 is only 0.339
%r_samp_const = 22547501; % bad t1, but strong t9 1 hr
%r_samp_const = 90835001; % good t9, but not as good t1, high amplitude
%r_samp_const = 91430001; % t9 missing right half
bchan1_const = 21;
bchan2_const = 46;
roi_1 = 'parsopercularis'; % source roi
roi_2 = 'superiortemporal'; % doesn't matter, plots against all
roi_1_const = roi_1;

elec_z = 1; % z coordinate of electrode points
fig_fmt = '-djpeg';
trig_eps = true;
trig_draw_lines = true;
trig_mag_no_ct = true;
system('mkdir figures');
system('mkdir figures/T11d1');


%for iMM = 1:length(metricsM)
for iMM = 1
    metricM = metricsM{iMM};
    
    % Main loop
    for i = 1:length(SubjectsM)
        sidM = SubjectsM{i};
        %fn_artM = sprintf('%s/%s_art.h5',dir_artM,sidM);
        %fn_distM = sprintf('%s/%s_dists-%s-%i.mat',dir_resM,sidM,metricM,n_perm);
        %fn_graphM = sprintf('%s/%s_graph-%s.h5',dir_resM,sidM,metricM);
        fn_h5M = sprintf('%s/%s.h5',dir_h5M,sidM);
        %fn_coregM = sprintf('%s/%s/label/all_parcellation.mat',dir_corM,sidM);

        % Check if files exist
        ckf = {fn_h5M};
        for jj = 1:length(ckf)
            if (~exist(ckf{jj},'file'))
                fprintf(2,'E> File not found: %s\n',ckf{jj});
                return
            end
        end
        

        % check hemisphere
        hemi = {};
        if (exist(sprintf('%s/%s/surf/rh.full.flat.patch.3d',subjects_dir,sidM), 'file'))
            hemi = {hemi{:},'rh'};
        elseif (exist(sprintf('%s/%s/surf/lh.full.flat.patch.3d',subjects_dir,sidM), 'file'))
            hemi = {hemi{:},'lh'};
        end
        if (isempty(hemi))
            fprintf('W> file not found: %s\n',sprintf('%s/%s/surf/*h.full.flat.patch.3d',subjects_dir,sidM));
        else
            for jj = 1:length(hemi)
                % Read in file
                fn_patchM = sprintf('%s/%s/surf/%s.full.flat.patch.3d',subjects_dir,sidM,hemi{jj});
                fn_pialM = sprintf('%s/%s/surf/%s.pial',subjects_dir,sidM,hemi{jj});
                fn_cacheM = sprintf('cache/%s_%s_flat.mat',sidM,hemi{jj});
    %             fn_patch_asc = sprintf('%s/%s/surf/%s.full.flat.patch.3d.asc',subjects_dir,sid,hemi{j});
    % 
    %             % if ascii doesn't exist, make it
    %             if (~ exist(fn_patch_asc, 'file'))
    %                 fprintf('[*] ASCII flat patch not found, converting...\n')
    %                 system(sprintf('cd %s/%s/surf; mris_convert -p %s.full.flat.patch.3d %s.full.flat.patch.3d.asc',subjects_dir,sid,hemi{j},hemi{j}));
    %             end

                % Check if patch cache exists:
                if (exist(fn_cacheM,'file'))
                    load(fn_cacheM,'pa','tri_flat');
                else
                    [pa] = read_patch_rev(fn_patchM);
                    %k = convhull(pa.x,pa.y);
                    [v_pial,tri_pial] = read_surf(fn_pialM);
                    % Correct to be 0-indexed like patch indexes
                    tri_pial = tri_pial - 1;
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
                            fprintf('%s> %i of %i\n',sidM,ii,length(tri_pial));
                        end
                    end

                    %plot(pa.x,pa.y,'black.');
                    %h = trisurf(tri_flat,pa.x,pa.y,pa.z);
                    %set(h,'edgecolor','none')
                    save(fn_cacheM,'-v6');
                end

                % Load human cache
                load(sprintf('%s/xsub_out_%s_%i.mat',dir_cacheLM,sidM,iMM));
               
                chan_labels = h5readatt(fn_h5M,'/h5eeg/eeg','labels');
                bchan_labels = cell(ecog.n_bchan,1);
                for ii0 = 1:ecog.n_bchan
                    bchan_labels{ii0} = sprintf('%s-%s',chan_labels{ecog.bip(ii0,1)},chan_labels{ecog.bip(ii0,2)});
                end
                              
                
                % Load graphs
%                 art_idx = h5read(fn_artM,'/art_idx');
%                 art_idxB = (art_idx == 1);
%                 frac_art = h5readatt(fn_artM,'/art_idx','frac_art');
%                 w_art = h5readatt(fn_artM,'/art_idx','w');
                % trim last time sample of graph
%                 R = h5read(fn_graphM,'/R',[1 1],size(art_idx));
%                 w = double(h5read(fn_graphM,'/w'));
%                 [n_comb,n_graph] = size(R);
                
                % Plot flat map

                % Read annotation
                %atl_table = AT.P.AtlROIs{atl_i}.RH.table;
                atl_table = C.AtlROIs{atl}.RH.table;
                annot_dir = sprintf('%s/%s/label',subjects_dir,sidM);
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
                tri_color = zeros(length(roi_id),3);
                
                % ROI midpoints for fsaverage_sym
                % collect all points for each ROI
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
                
                % read electrode locations
                fprintf('[!] Ignore any errors:\n')
                l = read_label(sprintf('%s/%s',subjects_dir,sidM),'all_surf_ielvis');
                if (isempty(l))
                    l = read_label(sidM,'all_surf_ielvis');
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
                
%                 % constant color
%                 col_unknown = 0.333*[1 1 1];
%                 col_src = 0.5*[1 1 1];
%                 col_surf = 0.85*[1 1 1];
%                 col_roi_hase = 0.7*[1 1 1];
%                 col_bound = 0.4*[1 1 1];
%                 col_elec = 0*[1 1 1];
%                 col_elec_hi = 0.3*[1 1 1];
%                 col_elec_border = [0 0 0];
%                 col_elec_nosig = 0.4*[1 1 1];
%                 size_elec = 15;
%                 size_elec_hi = 15;
                size_elec_border = size_elec_hi*border_fac;
%                 trig_show_txt_mid = false;
                % 
                % ind_lingual = 14;
                % ind_cuneus = 6;
                % ind_inferiorparietal = 9;
                % source = ind_inferiorparietal;
                % %source = 34;

                n_colors = 100;
                cmap = corrcmap(n_colors);
                ct_range = [];
                line_x = [];
                line_y = [];
                
                
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

                % Save plot electrodes
                elec_xy = nan(length(l(:,end)),2);
                for i4 = 1:length(l(:,end))
                    ind_pa = find(pa.ind == l(i4,1),1);
                    if (~isempty(ind_pa))
                        %plot(pa.x(ind_pa),pa.y(ind_pa),'.','Color',col_elec,'MarkerSize',size_elec)
                        % Save for future
                        elec_xy(l(i4,5),1) = pa.x(ind_pa);
                        elec_xy(l(i4,5),2) = pa.y(ind_pa);
                    else
                        fprintf('[%s] skipped electrode: %i\n',sidM,l(i4,5));
                    end
                end

                view(90,90);
                axis off;

                roi_xy = cell(1,ecog.n_bchan);
                for i2 = 1:length(roi_id)

                    c_i = find(roi_id(i2) == atl_table(:,end),1);

                    if (~isempty(c_i))
                        roi_name = C.AtlROIs{atl}.RH.struct_names{c_i};
                        tri_color(i2,:) = atl_table(c_i,1:3)/255;
                        roi_xy{c_i} = [roi_xy{c_i}; [pa.x(i2) pa.y(i2)]];
                    end
                end


                % Show boundary
                for i3 = 1:length(tri_flat)
                    vB = tri_flat(i3,:); % 6.8 seconds
                    tv = roi_id(vB); % 9.6 seconds
                    if (length(unique(tv)) > 1)
                        plot(mean(pa.x(vB)),mean(pa.y(vB)),'.','color',col_bound,'MarkerSize',5)
                    end
                end

                %set(a,'FaceVertexCdata',tri_color); % color by ROI
                %set(a,'FaceVertexCdata',tri_color*0.08 + 0.8*[1 1 1]); % color by trans ROI
                %set(a,'FaceVertexCdata',col_surf) % color bland
                daspect([1 1 1]); % fix aspect ratio to data
                % -------------------------------------------------------------------------------
                tri_color_val = mean(col_surf)*ones(length(roi_id),3);
                tri_color_val_l = cell(length(roi_id),1);
                
                % Plot 1 pair
                %roi_1 = 'parsorbitalis';
                %roi_2 = 'inferiorparietal';
                coh_t = [];
                countM = 1;
                for ii1M = 1:(ecog.n_bchan-1)
                    for ii2M = (ii1M+1):ecog.n_bchan
                        %fprintf('chan1:%i ii1:%i chan2:%i ii2:%i\n',chan1(count)+1,ii1,chan2(count)+1,ii2)

                        has_roi_f = ( strcmp(atl_labels{ecog.bip(ii1M,1)},roi_1) & strcmp(atl_labels{ecog.bip(ii2M,1)},roi_2) )...
                                | ( strcmp(atl_labels{ecog.bip(ii1M,1)},roi_1) & strcmp(atl_labels{ecog.bip(ii2M,2)},roi_2) )...
                                | ( strcmp(atl_labels{ecog.bip(ii1M,2)},roi_1) & strcmp(atl_labels{ecog.bip(ii2M,1)},roi_2) )...
                                | ( strcmp(atl_labels{ecog.bip(ii1M,2)},roi_1) & strcmp(atl_labels{ecog.bip(ii2M,2)},roi_2) );
                        has_roi_r = ( strcmp(atl_labels{ecog.bip(ii1M,1)},roi_2) & strcmp(atl_labels{ecog.bip(ii2M,1)},roi_1) )...
                                | ( strcmp(atl_labels{ecog.bip(ii1M,1)},roi_2) & strcmp(atl_labels{ecog.bip(ii2M,2)},roi_1) )...
                                | ( strcmp(atl_labels{ecog.bip(ii1M,2)},roi_2) & strcmp(atl_labels{ecog.bip(ii2M,1)},roi_1) )...
                                | ( strcmp(atl_labels{ecog.bip(ii1M,2)},roi_2) & strcmp(atl_labels{ecog.bip(ii2M,2)},roi_1) );

                        %has_chans = ( (ii1M == bchan1_const) && (ii2M == bchan2_const) );
                        has_chans = ( (ii1M == bchan1_const) );

                        %if (((has_roi_f || has_roi_r) && (Dmat(ii1,ii2) > dist_thresh)) && (ct(count) > ct_thresh))
                        %if (((has_chans) && (Dmat(ii1M,ii2M) > dist_thresh)) && (ct(countM) > ct_thresh))
                        if (((has_chans) && (Dmat(ii1M,ii2M) > dist_thresh)))
%                             p_val_ct_2 = 0.01/n_comb;
%                             ct_thresh_perm = (binoinv([p_val_ct_2 (1-p_val_ct_2)],n_graph,p_val/bonf_seconds)/n_graph);
                              ct_thresh_perm = [0 0];
%                             pname = sprintf('%s__%i_%s_%s__%i_%s_%s__%imm_ct%i_mag%i',...
%                                 sidM,ii1M,atl_labels{ecog.bip(ii1M,1)},bchan_labels{ii1M},ii2M,atl_labels{ecog.bip(ii2M,1)},...
%                                 bchan_labels{ii2M},round(Dmat(ii1M,ii2M)),round(1000*ct(count)),round(1000*mag(count)));
                            pname = sprintf('%s__%i_%s_%s__%i_%s_%s__%imm_ct%i_ctt%i-%i_mag%i_coh%it%i',...
                                sidM,ii1M,atl_labels{ecog.bip(ii1M,1)},bchan_labels{ii1M},ii2M,atl_labels{ecog.bip(ii2M,1)},...
                                bchan_labels{ii2M},round(Dmat(ii1M,ii2M)),round(1000*ct(countM)),...
                                ceil(1000*ct_thresh_perm(1)),ceil(1000*ct_thresh_perm(2)),round(1000*mag(countM)),...
                                iMM,round(1000*coh_thresh(countM)));
                
                            b1 = ecog.bip(ii1M,1:2);
                            b2 = ecog.bip(ii2M,1:2);
                            if (has_roi_r)
                                b2t = b2;
                                b2 = b1;
                                b1 = b2t;
                            end
                            
                            % Plot electrode interactivity
                            
                            % show electrodes involved
                            %return
                            if (trig_draw_lines && ( ct(countM) > ct_thresh) )
                                line_x = [line_x; [elec_xy(b1(1),1), elec_xy(b2(1),1)]];
                                line_y = [line_y; [elec_xy(b1(1),2), elec_xy(b2(1),2)]];
                                if (trig_mag_no_ct)
                                    ct_range = [ct_range; mag(countM)];
                                    %col_line = cmap(round(n_colors*mag(countM)),:);
                                else
                                    ct_range = [ct_range; ct(countM)];
                                    %col_line = cmap(round(n_colors*ct(countM)),:);
                                end
                                coh_t = [coh_t;coh_thresh(countM)];
                                
                                %plot(line_x,line_y,'-','Color',col_line,'LineWidth',5);
                            end
                            
                            for jj1 = 1 %1:length(b1)
                                % size_elec_border
                                plot3(elec_xy(b1(jj1),1),elec_xy(b1(jj1),2),elec_z,'.','Color',col_elec_border,'MarkerSize',size_elec_border)
                                plot3(elec_xy(b1(jj1),1),elec_xy(b1(jj1),2),elec_z,'.','Color',col_elec_hi,'MarkerSize',size_elec_hi)
                            end
                            for jj2 = 1 %1:length(b2)
                                plot3(elec_xy(b1(jj2),1),elec_xy(b1(jj2),2),elec_z,'.','Color',col_elec_border,'MarkerSize',size_elec_border)
                                plot3(elec_xy(b2(jj2),1),elec_xy(b2(jj2),2),elec_z,'.','Color',col_elec_hi,'MarkerSize',size_elec_hi)
                            end
                            
                            midpt_b2 = (elec_xy(b2(1),:) + elec_xy(b2(2),:))/2;
                            
                            roi_b1c1 = atl_labels{b1c1(countM)};
                            roi_b1c2 = atl_labels{b1c2(countM)};
                            roi_b2c1 = atl_labels{b2c1(countM)};
                            roi_b2c2 = atl_labels{b2c2(countM)};
                            
                            % color roi
                            for i2 = 1:length(roi_id)
                                c_i = find(roi_id(i2) == atl_table(:,end),1);
                                if (~isempty(c_i))
                                    roi_name = C.AtlROIs{atl}.RH.struct_names{c_i};
                                    % Color ROIs
                                    if (strcmp(roi_name,roi_b1c1) || strcmp(roi_name,roi_b1c2))
                                        % Source ROI color
                                        tri_color_val(i2,:) = col_src;
%                                     elseif(strcmp(roi_name,roi_b2c1) || strcmp(roi_name,roi_b2c2))
%                                         % Target ROI color
%                                         if (sum(tri_color_val(i2,:)==0) == 3)
%                                             % no pre-existing ROI color
%                                             tri_color_val(i2,:) = cmap(round(ct(countM)*n_colors),:);
%                                         else
%                                             % remap roi color
%                                             tri_color_val_l{i2} = [tri_color_val_l{i2}, ct(countM)];
%                                             tri_color_val(i2,:) = cmap(round(mean(tri_color_val_l{i2})*n_colors),:);
%                                         end
                                    else
                                        tri_color_val(i2,:) = col_surf;
                                    end
                                end
                            end
                            
                            if (trig_show_txt_mid)
                                text(midpt_b2(1),midpt_b2(2),sprintf('%.3f',mag(countM)),'HorizontalAlignment','center','FontSize',6);
                            end
                            
                        end
                        countM = countM + 1;
                    end
                end
                
                % show coherence threshold stats
                fprintf('Coherence threshold mean: %.3f\n',mean(coh_t));
                fprintf('Coherence threshold stdev: %.3f\n',std(coh_t));
                fprintf('Coherence threshold median: %.3f\n',median(coh_t));
                fprintf('Number of significant bchan pairs: %i\n',length(coh_t));

                
                % Plot lines
                ct_abs = max(ct_range) - min(ct_range);
                for ipl = 1:length(ct_range)
%                     line_x = [line_x; [elec_xy(b1(1),1), elec_xy(b2(1),1)]];
%                     line_y = [line_y; [elec_xy(b1(1),2), elec_xy(b2(1),2)]];
%                     if (trig_mag_no_ct)
%                         ct_range = [ct_range; mag(countM)];
%                         %col_line = cmap(round(n_colors*mag(countM)),:);
%                     else
%                         ct_range = [ct_range; ct(countM)];
%                         %col_line = cmap(round(n_colors*ct(countM)),:);
%                     end
                    if (ct_abs == 0)
                        ct_norm = 1;
                    else
                        ct_norm = (ct_range(ipl) - min(ct_range) )/ct_abs;
                    end
                    %ct_norm = (ct_range(ipl) - min(ct_range) )/ct_abs;
                    %fprintf('ct_norm: %.4f\n',ct_norm);
                    col_line = cmap(ceil((n_colors-1)*ct_norm)+1,:);
                    plot(line_x(ipl,:),line_y(ipl,:),'-','Color',col_line,'LineWidth',3);
                end
                
                set(a,'FaceVertexCdata',tri_color_val);
%                 colormap(cmap);
%                 %caxis([0 1]);
%                 caxis([min(ct_range),max(ct_range)]);
%                 colorbar;

                colormap(cmap);
                %if (~isempty(ct_range))
                if (~isempty(ct_range) && ((max(ct_range)-min(ct_range)) ~= 0) )
                    caxis([min(ct_range) max(ct_range)]);
                    cbh = colorbar('XTick', linspace(min(ct_range),max(ct_range),5));
                else
                    colorbar;
                end

                xrange = max(pa.x) - min(pa.x);
                yrange = max(pa.y) - min(pa.y);
                scale_fac = 10; % mm
                bar_horiz_x = [min(pa.x), min(pa.x) + scale_fac];
                bar_horiz_y = [min(pa.y), min(pa.y)];
                bar_vert_x = [min(pa.x), min(pa.x)];
                bar_vert_y = [min(pa.y), min(pa.y) + scale_fac];
                plot(bar_horiz_x,bar_horiz_y,'black-','LineWidth',5);
                text(min(pa.x)+0.5*scale_fac,min(pa.y)+0.5*scale_fac,'1 cm','FontSize',8);

                %return
                if (trig_mag_no_ct)
                    magoct = 'mag';
                else
                    magoct = 'ct';
                end
                print(h,sprintf('figures/T11d1/%s_%s',pname,magoct),fig_fmt);
                if (trig_eps)
                    print(h,sprintf('figures/T11d1/%s_%s',pname,magoct),'-depsc');
                end
                close(h);
                % reverse labels
                
                % -------------------------------------------------------------------------------
                
                
            end
        end
        
%         if ( i == length(SubjectsM) )
%         end

    end

end

% Clear loop indices
clear i;
clear j;

% Print finish message
fprintf('Done.\n')





function patch = read_patch_rev(fname)
% Original Author: Bruce Fischl
% Edited: Jerry Wang, July 31 2018

fid = fopen(fname,'r');
if (fid == -1)
   error('could not open file %s', fname) ;
end

ver = fread(fid, 1, 'int', 0, 'b');
if (ver ~= -1)
   error('incorrect version # %d (not -1) found in file',ver) ;
end

patch.npts = fread(fid, 1, 'int', 0, 'b') ;

for i=1:patch.npts
    ind = fread(fid, 1, 'int', 0, 'b') ;
    if (ind < 0)
       ind = -ind - 1 ;
    else
       ind = ind - 1 ;
    end
    patch.ind(i) = ind ;
    patch.x(i) = fread(fid, 1, 'float', 0, 'b') ;
    patch.y(i) = fread(fid, 1, 'float', 0, 'b') ;
    patch.z(i) = fread(fid, 1, 'float', 0, 'b') ;
    patch.vno(i) = ind ;
end

fclose(fid);
end