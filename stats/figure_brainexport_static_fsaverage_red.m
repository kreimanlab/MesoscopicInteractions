close all;
clear;

subjects_dir = '/media/jerry/internal/data/coreg';
dir_proj = '/home/jerry/Downloads/blend4web_ce/projects/brainview';
dir_assets = sprintf('%s/assets',dir_proj);

SubjectsLoc = { ...
   'fsaverage_sym','sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
   'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
   'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',... % 'sub23',
   'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',... % 'sub30',
   'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
   'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
hz_txt = {'0-125 Hz','3-8 Hz','8-12 Hz','12-30 Hz','30-100 Hz'};

path_ds_factor = 1;
scale_mm_per_obj = 30;
flag_export_obj = true;
flag_nowrite = false; % False to overwrite existing files
flag_plot = false;
flag_colorbar = true;
flag_list_sig = false; % false maps average location to *all* member electrodes regardless of significance

% Pial surface downsampling
FS_SURF_DSAMPLING_F = 0.1; % fraction from 0 to 1 of faces to keep
            
% clear temporary script
fn_sh = sprintf('%s/blender/tmp.sh',dir_proj);
%if (~flag_nowrite)
system(sprintf('rm -f %s',fn_sh));
%end

% deprecated
cluster_i = load('cache/fig_cluster3_cluster_i.mat');
cluster_i = cluster_i.cluster_i;
%cluster_i = 1:length(cluster_i); %does not work
% end deprecated

CaT14 = load(sprintf('./cache/figure_t14_%i_150',1));
cluster_i2 = zeros(size(cluster_i));
for i = 1:length(cluster_i)
    idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(i))));
    cluster_i2(i) = idx;
end
cluster_i = cluster_i2; % cluster_i2: converts Es indices to final adjacency matrix indices

% Override
%cluster_i = 1:length(cluster_i);

%for isub = 1:length(SubjectsLoc) %2:25 %length(SubjectsLoc)
for isub = 1:length(SubjectsLoc)

    for metrici = 1:length(metrics)

        % Frequency band nomenclature
        metric_str = metrics{metrici};
        metric_str = metric_str(3:end);
        if (strcmp(metric_str,'Broadband'))
            coh_str = sprintf('Coherence (%s)',hz_txt{metrici});
        else
            coh_str = sprintf('%s coherence (%s)',metric_str,hz_txt{metrici});
        end
        
        dir_cache = './cache';
        sid_const = SubjectsLoc{isub}; %'fsaverage_sym';
        if (~strcmp(sid_const,'fsaverage_sym'))
            sid_const_int = find(strcmp(SubjectsLoc(2:end),sid_const),1); %str2double(sid_const(2:end));
        else
            sid_const_int = 0;
        end

        mkdir('brainexport');
        dir_be = sprintf('brainexport/%s_%i_%i',sid_const,isub - 1,metrici);
        mkdir(dir_be);
        
        if (~flag_nowrite)
            system(sprintf('rm -rf %s/*',dir_be));
        end


        % Calculate maximum and minimum colors
        %fn_colorca = 'cache/figure_brainexport_static_fsaverage.mat';


        %load('figure_brainexport_6all_fsaverage_red');

        % Load paths
        CaRloc = load(sprintf('%s/fig_cluster2_reduce_%i_new.mat','cache',metrici));
        CaT14 = load(sprintf('./cache/figure_t14_%i_150',metrici));
        if (strcmp(sid_const,'fsaverage_sym'))
            % Load from figure_brainexport_6all_fsaverge_red.m
            load(sprintf('brainexport/red_6all_fsaverage_%i.mat',metrici));
            
            A = CaT14.Adj_plt2;
            A = A(cluster_i,cluster_i);
            A(A==0) = NaN;
            
%             A = CaRloc.A;
%             
%             %return
% 
%             % Remove electrodes in medial wall
%             A(rem,:) = [];
%             A(:,rem) = [];
            
            hemi = 'r';

            % Copy subjects list to int
            SubjectsInt = zeros(size(Subjects));
            for itm = 1:length(Subjects)
                SubjectsInt(itm) = str2double(Subjects{itm}(2:end));
            end

            % Load surface
            [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,'pial'));
            FV = triangulation(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3));
        else
            
            load(sprintf('brainexport/%s_6all_%i.mat',sid_const,metrici));
            % Load surface
            [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,'pial'));
            FV = triangulation(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3));
            center = mean(s_vert,1);
            l = read_label(sprintf('%s/%s',subjects_dir,sid_const),'all_surf_ielvis');
            if (isempty(l))
                l = read_label(sid_const,'all_surf_ielvis');
            end
            [~,sIdx] = sort(l(:,5),'ascend');
            l = l(sIdx,:);
            l = [l(:,1), l(:,5), l(:,2:4)];
            Ca = load(sprintf('%s/xsub_out_%s_%i.mat',dir_cache,sid_const,metrici));
            %if (false)            
            if (strcmp(sid_const,'sub26'))
                elec_hemi_idx = strcmp(Ca.C.EleHemi,upper(hemi));
                E = []; %zeros(sum(elec_hemi_idx),5);
                ic = 1;
                bchan_ishemi = false(Ca.ecog.n_bchan,1);
                for i = 1:(Ca.ecog.n_bchan)
                    b1c1 = Ca.ecog.bip(i,1);
                    b1c1_ishemi = elec_hemi_idx(b1c1);
                    if (b1c1_ishemi)
                        ctmp = vf_pia.Points(l(b1c1,1)+1,:);
                        E(ic,:) = l(b1c1,:);
                        E(ic,3:5) = ctmp;
                        E(ic,2) = ic;
                        ic = ic + 1;
                    end
                    bchan_ishemi(i) = b1c1_ishemi;
                end
                A = Ca.AdjMag;
                n_A_full = length(A);
                A_ind = 1:n_A_full;
                A(A == 0) = NaN;
                A = A(bchan_ishemi,bchan_ishemi);
                A_ind = A_ind(bchan_ishemi);
                
                % Index Paths by 2-channels
                Paths2 = cell(n_A_full,n_A_full);
                nap = 1;
                for iap = 1:(n_A_full-1)
                    for jap = (iap+1):n_A_full
                        Paths2{iap,jap} = Paths{nap};
                        Paths2{jap,iap} = Paths{nap};
                        nap = nap + 1;
                    end
                end
            else
                E = zeros(Ca.ecog.n_bchan,5);
                for i = 1:(Ca.ecog.n_bchan)
                    b1c1 = Ca.ecog.bip(i,1);
                    ctmp = vf_pia.Points(l(b1c1,1)+1,:);
                    E(i,:) = l(b1c1,:);
                    E(i,3:5) = ctmp;
                    E(i,2) = i;
                    %return
                end
                A = Ca.AdjMag;
                A(A == 0) = NaN;
            end
            
            % Build list of locations for subject
            sid_const_int_tmp = str2double(sid_const(2:end));
            Eb2avg = nan(Ca.ecog.n_bchan,1);
            ico = 0;
            for iet = 1:length(CaRloc.Es)
                est = CaRloc.Es{iet};
                est2 = CaRloc.Es2{iet};
                sidints = est(:,1);
                for ieu = 1:length(sidints)
                    if ( (est(ieu,1)==sid_const_int_tmp) ) % (est2(ieu,6) > 0) &&
                        Eb2avg(est(ieu,2)) = iet;
                        ico = ico + 1;
                    end
                end
            end
            
            if (ico ~= Ca.ecog.n_bchan)
                fprintf(2,'[!] WARN: %i bchans mapped to fsaverage_sym in subject with %i bchans\n',ico,Ca.ecog.n_bchan);
            end
            
        end



        Eraw = E;
        E(:,3) = E(:,3) - center(1);
        E(:,4) = E(:,4) - center(2);
        E(:,5) = E(:,5) - center(3);

    %     Ad = CaR.Ad;
    %     Ad(rem,:) = [];
    %     Ad(:,rem) = [];
        n_A = length(A);

        % colormap function
        n_col = 2^8;
        %cmap = corrcmap(n_col);
        cmap = inferno(n_col);
        
        % Load distorted colormap
        %CaCM = load('./cache/corrcmap_test-cmap.mat');
        %n_col = CaCM.n_col;
        %cmap = CaCM.cmap;

        
        mag_min = nanmin(A(:));
        mag_max = nanmax(A(:));
        mag2col = @(x) cmap(round(((x - mag_min)/(mag_max - mag_min))*(n_col-1) + 1),:);


        tube_radius = 0.15;
        tube_n_edges = 3; % number of vertices of tube cross-section
        
        % Apply rotation matrix according to hemisphere
        if (strcmp(hemi,'r'))
            rotm = rotationVectorToMatrix([pi/2, 0, 0]);
            rotm_hemi = rotationVectorToMatrix([0, 0, 0]);
            rotm_hemi_elec = rotationVectorToMatrix([0, 0, 0]);
        else
            rotm = rotationVectorToMatrix([pi/2, 0, 0]);
            rotm_hemi = rotationVectorToMatrix([0, pi, 0]);
            rotm_hemi_elec = rotationVectorToMatrix([0, 0, pi]);
            rotm = rotm_hemi * rotm;
        end
        n = 1;
        Paths = Paths(Paths_ind);
        NamesPair = {};
        for i = 1:(n_A - 1)
            for j = (i+1):n_A
                mag = A(i,j);
                if (~isnan(mag))

                    if (strcmp(sid_const,'sub26'))
                        path = Paths2{A_ind(i),A_ind(j)};
                    else
                        path = Paths{n};
                    end

                    path = downsample(path,path_ds_factor);

                    path(:,1) = path(:,1) - center(1);
                    path(:,2) = path(:,2) - center(2);
                    path(:,3) = path(:,3) - center(3);


                    [n_path,~] = size(path);

                    try
                        col_tube = mag2col(mag);
                    catch
                        col_tube = cmap(end,:);
                    end
                    
                    if (~ flag_nowrite)
                        [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
                        pa = surf2patch(X,Y,Z,'triangles');
                    end
                    
                    if (flag_plot)
                        surf(X,Y,Z,'FaceColor',col_tube,'EdgeColor','none'); hold on;
                    end

                    if (strcmp(sid_const,'fsaverage_sym'))
                        chan1 = cluster_i(i); %Count(Count(:,3) == n,1);
                        chan2 = cluster_i(j); %Count(Count(:,3) == n,2);
                    elseif (strcmp(sid_const,'sub26'))
                        chan1 = A_ind(i);
                        chan2 = A_ind(j);
                    else
                        chan1 = i;%Ca.ecog.bip(i,1);
                        chan2 = j;%Ca.ecog.bip(j,1);
                    end
                    

                    if ((flag_export_obj) && (~ flag_nowrite))
                        clear obj;
                        s_vert_obj = pa.vertices;
                        s_vert_obj = (rotm * (s_vert_obj'))';
                        %s_vert_obj(:,1) = pa.vertices(:,1);
                        %s_vert_obj(:,2) = pa.vertices(:,3);
                        %s_vert_obj(:,3) = pa.vertices(:,2);
                        
                        % stable naming scheme
                        %mtl_name = sprintf('pair_%i_%i_%i_mag-%i',n,chan1,chan2,round(mag*1000));
                        
                        % new naming scheme
                        rgb_tube = round(col_tube*255);
                        mtl_name = sprintf('pair_%i_%i_%i_mag-%i_rgb-%i-%i-%i',n,chan1,chan2,...
                            round(mag*1000),rgb_tube(1),rgb_tube(2),rgb_tube(3));
                        
                        obj.vertices = s_vert_obj / scale_mm_per_obj;
                        obj.objects(1).type='usemtl'; % use material
                        obj.objects(1).data=mtl_name;
                        obj.objects(2).type='s'; % smooth shading
                        obj.objects(2).data='off'; % 1 for on, "off" for off
        %                 obj.objects(3).type='o';
        %                 obj.objects(3).data=mtl_name;
        %                 obj.objects(4).type='g';
        %                 obj.objects(4).data='pair';
                        obj.objects(3).type='f';
                        obj.objects(3).data.vertices = pa.faces;
                        obj.material(1).type='newmtl';
                        obj.material(1).data=mtl_name;
                        obj.material(2).type='Ka'; % ambient color
                        obj.material(2).data=col_tube;
                        obj.material(3).type='Kd'; % diffuse color
                        obj.material(3).data=col_tube;
                        obj.material(4).type='Ks'; % specular color
                        obj.material(4).data=col_tube; %col_tube_spec;
                        obj.material(5).type='illum';
                        obj.material(5).data=2;
                        obj.material(6).type='Ns'; % 0 to 1000, specular exponent
                        obj.material(6).data=0;
                        currdir = pwd;
                        fn_obj = sprintf('%s/%i_%s.obj',dir_be,sid_const_int,mtl_name);
                        fprintf('[*] Saving pair to: %s\n',fn_obj);
                        if (~ flag_nowrite)
                            write_wobj(obj,fn_obj);
                        end
                        NamesPair = [NamesPair; {mtl_name}];
                    end

                end

                n = n + 1;
            end
        end


        %return

        % Read projected electrodes
        % label is zero-indexed, checked by freeview
        % l = read_label(sprintf('%s/%s',subjects_dir,'fsaverage_sym'),sprintf('ielvis_%s',sid));
        % if (isempty(l))
        %     l = read_label(sprintf('%s','fsaverage_sym'),sprintf('ielvis_%s',sid));
        % end
        % [~,sIdx] = sort(l(:,end));
        % l = l(sIdx,:);

        [n_locations,~] = size(E);

        % Plot electrodes
        n_elec_pts = 4;
        r_elec = 2; % (mm)
        col_elec = 0.5*[1 1 1];
        col_elec_spec = [1 1 1];
        %for i = 1:(ecog.n_bchan)
        %%
        NamesElec = {};
        DescElec = {};
        TitleElec = {};
        CoordElec = {};
        
%         if (strcmp(sid_const,'sub26'))
%             return
%             n_locations = n_A;
%         end
        for i = 1:n_locations
            [x,y,z] = sphere(n_elec_pts);
            %coord = ecog.bip(i,4:6);
            %coord = l(l(:,end) == ecog.bip(i,1),2:4);
            %coord = s_vert(l(l(:,end) == ecog.bip(i,1),1)+1,:);
            coord = E(i,3:5);
            x = x + coord(1);
            y = y + coord(2);
            z = z + coord(3);

            if (flag_plot)
                mesh(x,y,z,'EdgeColor','none','FaceColor',col_elec); hold on;
            end

            pa = surf2patch(x,y,z,'triangles');

            % List subject membership for location
    %         if (strcmp(sid_const,'fsaverage_sym'))
    %             %strcmp(SubjectsInt,usubs)
    %             es = CaR.Es{i};
    %             usubs = sort(unique(es(:,1)));
    %             mstr = 'sub-';
    %             for iml = 1:length(usubs)
    %                 if (iml ~= 1)
    %                     mstr = [mstr,'-'];
    %                 end
    %                 mstr = [mstr,sprintf('%i', find(SubjectsInt == usubs(iml)) )];
    %             end
    %         else
    %             mstr = sprintf('sub-%i',find(strcmp(Subjects,sid_const)));
    %         end
    
            if (strcmp(sid_const,'fsaverage_sym'))
                mstr = 'sub-0';
            else
                % use 1-48 naming instead of m-numbers
                mstr = sprintf('sub-%i',find(strcmp(Subjects,sid_const)));
            end

            % Export to obj file
             if ((flag_export_obj) && (~ flag_nowrite))
                clear obj;
                s_vert_obj = pa.vertices;
                s_vert_obj = (rotm * (s_vert_obj'))';
                %s_vert_obj(:,1) = pa.vertices(:,1);
                %s_vert_obj(:,2) = pa.vertices(:,3);
                %s_vert_obj(:,3) = pa.vertices(:,2);
                %mtl_name = sprintf('electrode%i',ecog.bip(i,1));
                
                if (strcmp(sid_const,'fsaverage_sym'))
                    elec_name_int = cluster_i(i);
                else
                    elec_name_int = i; %E(i,2); %Ca.ecog.bip(E(i,2),1);
                end
                
                mtl_name = sprintf('elec-%i_%s',elec_name_int,mstr);
                obj.vertices = s_vert_obj / scale_mm_per_obj;
                obj.objects(1).type='usemtl'; % use material
                obj.objects(1).data=mtl_name;
                obj.objects(2).type='s'; % smooth shading
                obj.objects(2).data='off'; % 1 for on, "off" for off
        %         obj.objects(3).type='o';
        %         obj.objects(3).data=mtl_name;
                obj.objects(3).type='f';
                obj.objects(3).data.vertices = pa.faces;

                obj.material(1).type='newmtl';
                obj.material(1).data=mtl_name;
                obj.material(2).type='Ka'; % ambient color
                obj.material(2).data=col_elec;
                obj.material(3).type='Kd'; % diffuse color
                obj.material(3).data=col_elec;
                obj.material(4).type='Ks'; % specular color
                obj.material(4).data=col_elec; %col_elec_spec;
                obj.material(5).type='illum';
                obj.material(5).data=2;
                obj.material(6).type='Ns'; % 0 to 1000, specular exponent
                obj.material(6).data=0;
                currdir = pwd;
                fn_obj = sprintf('%s/%i_%s.obj',dir_be,sid_const_int,mtl_name);
                fprintf('[*] Saving electrode to: %s\n',fn_obj);
                if (~ flag_nowrite)
                    write_wobj(obj,fn_obj);
                end
                NamesElec = [NamesElec; {sprintf('%i_%s',sid_const_int,mtl_name)}];


                % Write electrode annotation description
                desc = '';
                titl = sprintf('%i',elec_name_int);
                if (strcmp(sid_const,'fsaverage_sym'))
                    %Location 
                    desc = sprintf('%sLocation: %i<br>',desc,elec_name_int);
                else
                    %Electrode 
                    desc = sprintf('%sElectrode: %i<br>',desc,elec_name_int);
                end
                
                % Average coherence for node
                mag_mean = nanmean(A(i,:));
                mag_std = nanstd(A(i,:));
                mag_nsig = sum(~isnan(A(i,:)));
                node_nsig = sum(A(i,:) > 0);
                
                coh_str2 = sprintf('Coherence');
                if (mag_nsig > 0)
%                     desc = sprintf('%s%s: %.2f (\x03c3=%.2f, n=%i)',...
%                         desc,coh_str2,mag_mean,mag_std,mag_nsig);
                    desc = sprintf('%sTotal interactions: %i',desc,node_nsig);
                else
                    %desc = sprintf('%s%s: Not available',desc,coh_str2);
                    desc = sprintf('%sTotal interactions: N/A',desc);
                end

                % get roi;
                coord = E(i,3:5);
                coord_raw = Eraw(i,3:5);
                %coord_scaled = coord / scale_mm_per_obj; %((rotm * (coord'))') / scale_mm_per_obj;
                % works, as of April 29, 2020:
                coord_scaled = ((rotm_hemi_elec * (coord'))') / scale_mm_per_obj;

                %desc = sprintf('%s<br>',desc);
%                 atlas = {'Brodmann','Desikan-Killiany','Glasser-Van Essen','Markov-Kennedy (Monkey)'};
%                 atlases = {'PALS_B12_Brodmann','aparc','HCP-MMP1','MACAQUE_M132'};
                
                % Updated May 6, 2020
                atlas = {'Desikan','Brodmann','Markov(monkey)'};
                atlases = {'aparc','PALS_B12_Brodmann','MACAQUE_M132'};
                for iatl = 1:length(atlases)
                    atlname = atlases{iatl};
                    fn_annot = sprintf('%s/%s/label/%sh.%s.annot',subjects_dir,sid_const,hemi,atlname);
                    [vertices,label,ctab] = read_annotation(fn_annot);
                    xIdx = ~ (all(ctab.table==0,2));
                    roi_table = ctab.table(xIdx,:);
                    roi_names = ctab.struct_names(xIdx,:);
                    roi_id = label(FV.nearestNeighbor(coord_raw));
                    try
                        roi = roi_names{roi_table(:,end)==roi_id};
                    catch
                        roi = 'Unknown';
                    end
                    roi = strsplit(roi,'_ROI');
                    roi = strsplit(roi{1},'_M132');
                    roi = strsplit(roi{1},'Brodmann.');
                    roi = strsplit(roi{end},sprintf('%s_',upper(hemi)));
                    roi = roi{end};
                    %roi = replace(roi,'.',' ');
                    %roi = replace(roi,'_',' ');
                    if (strcmp(atlname,'aparc'))
                        roi = convertRoiDK(roi);
                    end
                    desc = sprintf('%s<br>Area %s: %s',desc,atlas{iatl},roi);
                end

                
                if (strcmp(sid_const,'fsaverage_sym'))
                    es = CaRloc.Es{i};
                    if (flag_list_sig) % true
                        % remove nonsignificant individual electrodes from list
                        es2 = CaRloc.Es2{i};
                        sigidx = (es2(:,end) ~= 0);
                        es = es(sigidx,:);
                    end
                    
                    n_el = length(es(:,1));
                    [eusub,~,cu] = unique(es(:,1));
                    desc = sprintf('%s<br>Total electrodes: %i<br>',desc,n_el);

                    for ius = 1:length(eusub)
                        % convert subject number
                        sid_tmp = sprintf('%i',eusub(ius));
                        while (length(sid_tmp) < 5)
                            sid_tmp = ['0',sid_tmp];
                        end
                        sid_tmp = ['m',sid_tmp];
                        sid_tmp_int = find(strcmp(SubjectsLoc(2:end),sid_tmp),1);
                        CaT = load(sprintf('%s/xsub_out_%s_%i.mat',dir_cache,sid_tmp,metrici));
                        
                        % number of electrode in unique subject number
                        if (~isempty(sid_tmp_int))
                            bchans = es(cu == ius,2);
                            btxt = 'elec=';
                            for ib = 1:length(bchans)
                                b1c1t = bchans(ib); %CaT.ecog.bip(bchans(ib),1);
                                if (ib == 1)
                                    btxt = sprintf('%s%i',btxt,b1c1t);
                                else
                                    btxt = sprintf('%s,%i',btxt,b1c1t);
                                end
                            end
                            
                            % Construct links to individual subjects
                            n_elu = sum(cu == ius);
                            desc = sprintf('%s<a href="brainview.html?subject=%i&freq=1&%s" class="anchorLink">    Subject %i</a> # electrodes: %i',...
                                desc,sid_tmp_int,btxt,sid_tmp_int,n_elu); %eusub(ius));
                            %if (ius ~= length(eusub))
                            desc = sprintf('%s<br>',desc);
                            %end
                        else
                            fprintf(2,'[!] Warn: could not find patient %s\n',sid_tmp);
                        end
                    end
                else
                    % i - bchan number
                    b1 = Eb2avg(i);
                    
                    % Add 150 custom parcellation text
                    desc = sprintf('%s<br>Area Custom: %i',desc,cluster_i(b1));
                    
                    % Refer back to 150
                    %btxt = sprintf('elec=%i',b1);
                    btxt = sprintf('elec=%i',cluster_i(b1));
                    desc = sprintf('%s<br><a href="brainview.html?subject=%i&freq=1&%s" class="anchorLink">Explore on average brain</a>',desc,0,btxt);
                end

                %desc = sprintf('<div id="elec_anchor" style="overflow:auto; height:200px;">%s</div>',desc);
                desc = sprintf('<div class="elec_anchor">%s</div>',desc);

                DescElec = [DescElec; {desc}];
                TitleElec = [TitleElec; {titl}];
                CoordElec = [CoordElec; {coord_scaled}];

            end
        end

        % Export pair and elec names to json
        %currdir = pwd;
        %addpath('jsonlab-1.5');
        %savejson('NamesElec',NamesElec,sprintf('%s/brainexport/%s/NamesElec.json',currdir,sid_const));
        %savejson('NamesPair',NamesPair,sprintf('%s/brainexport/%s/NamesPair.json',currdir,sid_const));
        if (true) %(~flag_nowrite)
            s = struct();
            s.elec = NamesElec;
            s.pair = NamesPair;
            if (strcmp(sid_const,'fsaverage_sym'))
                s.A = isnan(CaRloc.A_nocov);
                s.n_days = NaN;
                s.n_chan = length(CaRloc.A);
            else
%                 % Convert bipolar adjacency matrix to electrode
%                 Asub = nan(Ca.ecog.n_chan,Ca.ecog.n_chan);
%                 for itp = 1:Ca.ecog.n_chan
%                     for itq = 1:Ca.ecog.n_chan
%                         b1 = Ca.ecog.bip(:,1) == itp;
%                         b2 = Ca.ecog.bip(:,1) == itq;
%                         if ((sum(b1) ~= 0) && (sum(b2) ~= 0))
%                             Asub(itp,itq) = Ca.AdjMag(b1,b2);
%                         end
%                     end
%                 end
                s.A = isnan(Ca.AdjMag); %isnan(Asub);
                s.n_days = Ca.n_graph*Ca.w/(3600*24);
                s.n_chan = Ca.ecog.n_chan;
                s.n_bchan = Ca.ecog.n_bchan;
            end
            S = jsonencode(s);
            %S = replace(S,'],','],\n\t');
            %S = sprintf('export default {\n%s\n}',S);
            S_fn = sprintf('%s/data_%i_%i.json',dir_be,sid_const_int,metrici);
            sf = fopen(S_fn,'w');
            fprintf(sf,S);
            fclose(sf);

            %if (sid_const_int == 0)
            system(sprintf('cp %s %s/data_%i_%i.json',S_fn,dir_assets,sid_const_int,metrici));
            %end
            
        end

        col_pial = 0.6*[1 1 1];
        col_pial_spec = [1 1 1];
        alpha_pial = 0.07; % was 0.2
        %subjects_dir = '/media/jerry/internal/data/coreg';
        SURFACE_TYPE = 'pial';

        % Save surface object
        currdir = pwd;
        [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,SURFACE_TYPE));
        s_vert(:,1) = s_vert(:,1) - center(1);
        s_vert(:,2) = s_vert(:,2) - center(2);
        s_vert(:,3) = s_vert(:,3) - center(3);

         if ((flag_export_obj) && (~ flag_nowrite))
            fn_obj = sprintf('%s/%i_%s.obj',dir_be,sid_const_int,SURFACE_TYPE);
            fprintf('[*] Saving brain to: %s\n',fn_obj);
            clear obj;
            s_vert_obj = s_vert;
            s_vert_obj = (rotm * (s_vert_obj'))';
            %s_vert_obj(:,1) = s_vert(:,1);
            %s_vert_obj(:,2) = s_vert(:,3);
            %s_vert_obj(:,3) = s_vert(:,2);
            
            % Pial surface downsampling
            %FS_SURF_DSAMPLING_F = 0.1;
            nFV = reducepatch((faces+1),s_vert_obj,FS_SURF_DSAMPLING_F);
            
            obj.vertices = nFV.vertices / scale_mm_per_obj; %s_vert_obj / scale_mm_per_obj;
            obj.objects(1).type='usemtl'; % use material
            obj.objects(1).data='pial';
            obj.objects(2).type='s'; % smooth shading
            obj.objects(2).data='1'; % 1 for on, "off" for off
            obj.objects(3).type='f';
            obj.objects(3).data.vertices = nFV.faces; %faces + 1;
            obj.material(1).type='newmtl';
            obj.material(1).data='pial';
            obj.material(2).type='Ka'; % ambient color
            obj.material(2).data=col_pial;
            obj.material(3).type='Kd'; % diffuse color
            obj.material(3).data=col_pial;
            obj.material(4).type='Ks'; % specular color
            obj.material(4).data=col_pial_spec;
            obj.material(5).type='illum';
            obj.material(5).data=9;
            obj.material(6).type='Ns'; % 0 to 1000, specular exponent
            obj.material(6).data=5;
            obj.material(7).type='d';
            obj.material(7).data=alpha_pial;
            obj.material(7).type='Tr';
            obj.material(7).data=(1-alpha_pial);
            if (~flag_nowrite)
                write_wobj(obj,fn_obj);
            end

        end



        % Plot surface
        if (flag_plot)
            p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
                'EdgeColor','none','FaceColor',col_pial,'FaceAlpha',alpha_pial); 
            daspect([1 1 1]);
            p.AmbientStrength = 0.3 ;
            p.DiffuseStrength = 0.4 ;
            p.SpecularStrength = 0;
            p.SpecularExponent = 1;
            p.BackFaceLighting = 'lit';
            cam_elev = 0;
            camlight(-135,cam_elev);
            camlight(45,cam_elev);
            camlight(-225,cam_elev);
            camlight(-45,cam_elev);

            view(90,0);
            axis off;
            set(gcf,'Position',[0 0 1920 1080]);
            colormap(cmap);
            cb = colorbar;
            set(cb,'Location','east');
            set(cb,'TickLength',0);
            set(cb,'Ticks',linspace(mag_min,mag_max,2));
            %cb.Label.String = 'Coherence';
            caxis([mag_min mag_max]);
            print(gcf,'figures/brainexport_static_fsaverage_red','-depsc');
            print(gcf,'figures/brainexport_static_fsaverage_red','-dpng','-r300');
        end

        %%
        % Print colorbar
        if (flag_colorbar && (~ flag_nowrite) && (mag_min ~= mag_max) && (~isnan(mag_min)) && (~isnan(mag_max)))
            fprintf('[*] Exporting colorbar..\n');
            FONT_SIZE = 16;
            N_TICKS = 4;
            h = figure('Position',3*[0 0 40 300],'visible','off');
            
            % Extend colorbar to shorter margin, but will cause figure to
            % be drawn regardless of 'visible','off'
            %[ha, pos] = tight_subplot(1,1,0,0,0);
            %axes(ha);
            
            set(gcf, 'Color', 'None');
            %h.Color(4) = 0;
            %set(gca,'Color',zeros(1,4));
            set(gca, 'Color', 'None');
            colormap(cmap);
            cb = colorbar;
            set(cb,'Location','west');
            set(cb,'TickLength',0);
            xt = linspace(mag_min,mag_max,N_TICKS);
            xts = cell(1,length(xt));
            for i = 1:length(xt)
                xts{i} = sprintf('%.2f',xt(i));
            end
            set(cb,'Ticks',xt);
            set(cb,'TickLabels',xts);
            cb.Label.String = coh_str; %'Coherence';
            caxis([mag_min mag_max]);
            cb.Color = 0.8*[1 1 1];
            cb.FontSize = FONT_SIZE;
            %h.Children(1).Color = 'none';
            axis off;
            fn_cb = sprintf('%s/cb_%i_%i',dir_assets,sid_const_int,metrici);

            im = frame2im(getframe(h));
            %im_sav = im;
            cutoff = 0; %round(mean(cb.Color) * (2^8) * 0.1);
            im_alpha = all(im<=cutoff,3) ~= 1;
            idx_crop = ~all(im_alpha == 0);

            % crop from left
            state = false;
            idx_crop2 = false(size(idx_crop));
            for icrop = 1:length(idx_crop)
                if (idx_crop(icrop))
                    state = true;
                end
                idx_crop2(icrop) = state;
            end
            
            imwrite(im(:,idx_crop2,:),[fn_cb,'.png'],'png','Alpha',double(im_alpha(:,idx_crop2)));
            %print(h,fn_cb,'-dsvg','-painters');
            %addpath('altmany-export_fig-d570645');
            %export_fig(fn_cb,'-png','-painters','-transparent');
            %addpath('plot2svg_20120915');
            %plot2svg([fn_cb,'.svg'], h); 
            %saveas(h,[fn_cb,'.png']);
            close(h);
        end

        %return
        %%
        % Import .obj into blender file
        if (flag_export_obj)
            % Make a copy from template .blend
            system(sprintf('cp %s/blender/brainview_template.blend %s/blender/brainview_%i_%i.blend',...
                dir_proj,dir_proj,sid_const_int,metrici));

            % Append to temporary script to be executed in terminal
            f_sh = fopen(fn_sh,'a');
            fprintf(f_sh,'blender --background %s/blender/brainview_%i_%i.blend --threads 6 --python %s/blender/batch_import_arg.py -- %i %s_%i\n',...
                dir_proj,sid_const_int,metrici, dir_proj,metrici,sid_const,sid_const_int);
            fprintf(f_sh,'blender --background %s/blender/brainview_%i_%i.blend --threads 6 --python %s/blender/batch_selectable.py --\n',...
                dir_proj,sid_const_int,metrici, dir_proj);
            fclose(f_sh);

            % Annotation coordinates
            fn_annot = sprintf('%s/annot_%i_%i.mat',dir_be,sid_const_int,metrici);
            fprintf('[*] Saving: %s\n',fn_annot);
            save(fn_annot,'E','NamesElec','TitleElec','DescElec','CoordElec');
        end

        %return

        fprintf('[*] Metric %i done.\n',metrici)

    end

end

fprintf('[!] All Done.\n')
