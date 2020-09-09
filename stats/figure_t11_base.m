close all;
clear;

metric = 'pcBroadband';
n_perm = 10000;

% anatomy folder
subjects_dir = '/mnt/cuenap_ssd/coregistration';

% Fast i/o definitions
dir_art = '/media/klab/internal/data/h5_notch20/art';
dir_res = '/media/klab/internal/data/results/coh_w10';
dir_cor = '/media/klab/internal/data/coreg';
dir_cacheL = './cache';

% Slow i/o definitions
dir_h5 = '/media/klab/KLAB101/h5_notch20';

metrics = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};

% Patients
Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48',...
    'mSu'};

% Exclude monkey
Subjects = Subjects(1:(end-1));


for iM = 1:length(metrics)
    metric = metrics{iM};
    
    % Main loop
    for i = 1:length(Subjects)
        sid = Subjects{i};
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
                return
            end
        end

        % check hemisphere
        hemi = {};
        if (exist(sprintf('%s/%s/surf/rh.full.flat.patch.3d',subjects_dir,sid), 'file'))
            hemi = {hemi{:},'rh'};
        elseif (exist(sprintf('%s/%s/surf/lh.full.flat.patch.3d',subjects_dir,sid), 'file'))
            hemi = {hemi{:},'lh'};
        end
        if (isempty(hemi))
            fprintf('W> file not found: %s\n',sprintf('%s/%s/surf/*h.full.flat.patch.3d',subjects_dir,sid));
        else
            for jj = 1:length(hemi)
                % Read in file
                fn_patch = sprintf('%s/%s/surf/%s.full.flat.patch.3d',subjects_dir,sid,hemi{jj});
                fn_pial = sprintf('%s/%s/surf/%s.pial',subjects_dir,sid,hemi{jj});
                fn_cache = sprintf('cache/%s_%s_flat.mat',sid,hemi{jj});
    %             fn_patch_asc = sprintf('%s/%s/surf/%s.full.flat.patch.3d.asc',subjects_dir,sid,hemi{j});
    % 
    %             % if ascii doesn't exist, make it
    %             if (~ exist(fn_patch_asc, 'file'))
    %                 fprintf('[*] ASCII flat patch not found, converting...\n')
    %                 system(sprintf('cd %s/%s/surf; mris_convert -p %s.full.flat.patch.3d %s.full.flat.patch.3d.asc',subjects_dir,sid,hemi{j},hemi{j}));
    %             end


                % Check if patch cache exists:
                if (exist(fn_cache,'file'))
                    load(fn_cache,'pa','tri_flat');
                else
                    [pa] = read_patch_rev(fn_patch);
                    %k = convhull(pa.x,pa.y);
                    [v_pial,tri_pial] = read_surf(fn_pial);
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
                            fprintf('%s> %i of %i\n',sid,ii,length(tri_pial));
                        end
                    end

                    %plot(pa.x,pa.y,'black.');
                    %h = trisurf(tri_flat,pa.x,pa.y,pa.z);
                    %set(h,'edgecolor','none')
                    save(fn_cache,'-v6');
                end

                % Load human cache
                load(sprintf('%s/xsub_out_%s_%i.mat',dir_cacheL,sid,iM));
               
                
                
                
                
                
                % Plot flat map

                % Read annotation
                %atl_table = AT.P.AtlROIs{atl_i}.RH.table;
                atl_table = C.AtlROIs{atl}.RH.table;
                annot_dir = sprintf('%s/%s/label',subjects_dir,sid);
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
                
                % read electrode locations
                l = read_label(sprintf('%s/%s',subjects_dir,sid),'all_surf_ielvis');
                
                % constant color
                col_unknown = 0.333*[1 1 1];
                col_source = [0.2 1 0.2];
                col_surf = 0.8*[1 1 1];
                col_bound = 0.2*[1 1 1];
                col_elec = 0*[1 1 1];
                size_elec = 15;
                % 
                % ind_lingual = 14;
                % ind_cuneus = 6;
                % ind_inferiorparietal = 9;
                % source = ind_inferiorparietal;
                % %source = 34;

                n_colors = 100;
                cc = corrcmap(n_colors);

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
                
                % Plot electrodes
                elec_xy = nan(length(l(:,end)),2);
                for i4 = 1:length(l(:,end))
                    ind_pa = find(pa.ind == l(i4,1),1);
                    if (~isempty(ind_pa))
                        plot(pa.x(ind_pa),pa.y(ind_pa),'.','Color',col_elec,'MarkerSize',size_elec)
                        % Save for future
                        elec_xy(l(i4,5),1) = pa.x(ind_pa);
                        elec_xy(l(i4,5),2) = pa.y(ind_pa);
                    else
                        fprintf('[%s] skipped electrode: %i\n',sid,l(i4,5));
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
                        plot(mean(pa.x(vB)),mean(pa.y(vB)),'.','color',col_bound)
                    end
                end

                %set(a,'FaceVertexCdata',tri_color); % color by ROI
                set(a,'FaceVertexCdata',tri_color*0.08 + 0.8*[1 1 1]); % color by trans ROI
                %set(a,'FaceVertexCdata',col_surf) % color bland
                daspect([1 1 1]);
                % -------------------------------------------------------------------------------
                
                return
            end
        end

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
