close all;
clear;

dir_cache = './cache';
subjects_dir = '/media/klab/internal/data/coreg';
%sid_const = 'fsaverage_sym';
scale_mm_per_obj = 30;
flag_export_obj = true;
flag_plot = false;
tube_radius = 1;

mkdir('figures');
mkdir(sprintf('figures/T16'));

%metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
imetrics = [1 5];
%imetrics = 1;

flag_cone = false;

trig_show_subject_num = false;
trig_show_colorbar = trig_show_subject_num;
trig_color_rois = true;
trig_show_esm_elecs = false;


try
    Esm = load('./cache/figure_t19');
catch
    fprintf('[*] Could not find ./cache/figure_t19.mat, did not load.\n')
end

for imm = 1:length(imetrics)
    iM = imetrics(imm);
    fprintf('[*] metric %i\n',iM);
    
    for trig_show_path = [0] %[0 1]
    
        % Load fsaverage_sym paths
        fn_paths = sprintf('brainexport/red_6all_%s_%i.mat','fsaverage',iM);
        fprintf('[*] Loading %s ..\n',fn_paths)
        CaP = load(fn_paths);
        Paths = CaP.Paths;
        Paths_ind = CaP.Paths_ind;
        Paths = Paths(Paths_ind);
        [n_loc,~] = size(CaP.E);

        PathsA = zeros(n_loc,n_loc);
        count = 1;
        for i = 1:(n_loc-1)
            for j = (i+1):n_loc
                PathsA(i,j) = count;
                PathsA(j,i) = count;
                count = count + 1;
            end
        end


        % original subject
        sid = 'm00005';
        sid_i = str2double(sid(2:end));
        b1 = 35;
        b2 = 84;
        CaR = load(sprintf('%s/fig_cluster2_reduce_%i_new',dir_cache,iM));
        CaA = load(sprintf('%s/xsub_out_all_%i',dir_cache,iM));

        % save vars
        sid_const = sid;
        b1_const = b1;
        b2_const = b2;
        
        % convert to integer subject number
        sid_int = find(strcmp(CaA.Subjects,sid));
        Ca = load(sprintf('%s/xsub_out_%s_%i',dir_cache,sid,iM));

        % Find ROIs bipolar electrodes map onto
        rois = Ca.C.AtlLabels{2};
        b1c1 = Ca.ecog.bip(b1,:);
        b2c1 = Ca.ecog.bip(b2,:);
        b1_roi = rois{b1c1};
        b2_roi = rois{b2c1};
        
        
        % Search all subs
        E_reduce = [];
        E_reduce_sid = {};
        E_reduce_a = [];
        E_reduce_sid_a = {};
        mag_atl_l = CaA.AdjAtl{strcmp(Ca.rois,b1_roi),strcmp(Ca.rois,b2_roi)};
        mag_atl_l_sid = CaA.AdjAtl_sid{strcmp(Ca.rois,b1_roi),strcmp(Ca.rois,b2_roi)};
        mag_atl_l_n = CaA.AdjAtlN{strcmp(Ca.rois,b1_roi),strcmp(Ca.rois,b2_roi)};
        
        usubs = unique(mag_atl_l_sid);
        for isub = 1:length(usubs)
            s1 = usubs(isub);
            sid_t = CaA.Subjects{s1};
            Ca_t = load(sprintf('%s/xsub_out_%s_%i',dir_cache,sid_t,iM));
            
            
            % find magnitudes
            mag_seek = mag_atl_l(mag_atl_l_sid==s1);
            %mag_seek = mag_seek(mag_seek ~= 0);
            for k = 1:length(mag_seek)
                if (mag_seek(k) == 0)
                    mag = 0;
                    b1 = 1;
                    b2 = 1;
                    hemi_isr = strcmp(Ca_t.hemi_list{1},'R');
                else
                    sIdx = find(Ca_t.mag == mag_seek(k));
                    if (isempty(sIdx))
                        fprintf(2,'[W] coherence value not found in AdjAtl.\n');
                    end
                    b1 = double(Ca_t.chan1(sIdx) + 1);
                    b2 = double(Ca_t.chan2(sIdx) + 1);
                    b1c1 = Ca_t.ecog.bip(b1,1);
                    b1c1_roi = Ca_t.C.EleHemi{b1c1};
                    % xsub_out_stats guarantees same hemisphere
%                     b2c1 = Ca_t.ecog.bip(b1,1);
%                     b2c1_roi = Ca_t.C.EleHemi{b1c1};
                    hemi_isr = strcmp(b1c1_roi,'R');
                    mag = mag_seek(k);
                    
                    E_reduce = [E_reduce; [s1,b1,b2,mag,Ca_t.ecog.bip(b1,4:6),hemi_isr]];
                    E_reduce_sid = [E_reduce_sid; {sid_t}];
                end
                
                E_reduce_a = [E_reduce_a; [s1,b1,b2,mag,Ca_t.ecog.bip(b1,4:6),hemi_isr]];
                E_reduce_sid_a = [E_reduce_sid_a; {sid_t}];
                
                
            end
            
            
            
        end
        
        
%         E_reduce_a = [];
%         E_reduce_sid_a = {};
%         for isub = 1:length(Ca.Subjects)
%             sid_t = Ca.Subjects{isub};
%             s1 = str2double(sid_t(2:end));
%             Ca_t = load(sprintf('%s/xsub_out_%s_%i',dir_cache,sid_t,iM));
%             count = 1;
%             rois_t = Ca_t.C.AtlLabels{2};
% %             hemi_list = Ca_t.hemi_list;
% %             isSameHemi = (length(unique(hemi_list)) == 1);
%             for i = 1:(Ca_t.ecog.n_bchan-1)
%                 b1c1 = Ca_t.ecog.bip(i,1);
%                 hemi_t = Ca_t.C.EleHemi{b1c1};
%                 hemi_isr = strcmp(hemi_t,'R');
%                 b1c1_roi_t = rois_t{b1c1};
%                 
%                 for j = (i+1):Ca_t.ecog.n_bchan
%                     b2c1 = Ca_t.ecog.bip(j,1);
%                     b2c1_roi_t = rois_t{b2c1};
%                     mag = Ca_t.AdjMag(count);
%                     cond_dist =  (Ca_t.Dmats(count) > Ca_t.dist_thresh);
%                     %cond_issig = ((mag > 0) && (~isnan(mag))) && cond_dist;
%                     cond_roipair1 = (strcmp(b1c1_roi_t,b1_roi) && strcmp(b2c1_roi_t,b2_roi));
%                     cond_roipair2 = (strcmp(b2c1_roi_t,b1_roi) && strcmp(b1c1_roi_t,b2_roi));
%                     cond_roipair = cond_roipair1 || cond_roipair2;
%                     
%                     if ((((mag > 0) && (~isnan(mag))) && cond_dist) && cond_roipair)
%                         E_reduce = [E_reduce; [s1,i,j,mag,Ca_t.ecog.bip(i,4:6),hemi_isr]];
%                         E_reduce_sid = [E_reduce_sid; {sid_t}];
%                     end
%                     
%                     if (((~isnan(mag)) && cond_dist) && cond_roipair)
%                         E_reduce_a = [E_reduce_a; [s1,i,j,mag,Ca_t.ecog.bip(i,4:6),hemi_isr]];
%                         E_reduce_sid_a = [E_reduce_sid_a; {sid_t}];
%                     end
%                     
%                     count = count + 1;
%                 end
%             end
%         end
        
        
        mag_atl = nanmean(mag_atl_l(mag_atl_l>0));
        mag_atl_cp = sum((mag_atl_l > 0) & (~isnan(mag_atl_l)))/numel(mag_atl_l);
        fprintf('[*] %s - %s: %.6f, consistency across pairs: %.6f (n=%i)\n',...
            b1_roi,b2_roi,mag_atl,mag_atl_cp,numel(mag_atl_l));
        
        %return
        %%
%         % Find locations bipolar electrodes map onto
%         l1 = NaN;
%         l2 = NaN;
%         n_Es = length(CaR.Es);
%         for i = 1:n_Es
%             ce = CaR.Es{i};
%             [n_ce,~] = size(ce);
%             for j = 1:n_ce
%                 subj = ce(j,1);
%                 subj_b = ce(j,2);
%                 if ((subj == sid_i) && (subj_b == b1))
%                     l1 = i;
%                 elseif ((subj == sid_i) && (subj_b == b2))
%                     l2 = i;
%                 end
%             end
%         end

%         % Find all interactions between location pair
%         e1 = CaR.Es{l1};
%         [n_e1,~] = size(e1);
%         [~,sIdx] = sort(e1(:,1));
%         e1 = e1(sIdx,:);
% 
%         e2 = CaR.Es{l2};
%         [n_e2,~] = size(e2);
%         [~,sIdx] = sort(e2(:,1));
%         e2 = e2(sIdx,:);
% 
%         E_reduce = [];
%         E_reduce_sid = {};
%         for i = 1:n_e1
%             for j = 1:n_e2
%                 % is there a significant interaction between i and j?
%                 s1 = e1(i,1);
%                 b1 = e1(i,2);
%                 sid1 = sprintf('%i',s1);
%                 while (length(sid1) < 5)
%                     sid1 = ['0',sid1];
%                 end
%                 sid1 = ['m',sid1];
%                 %Ca1 = load(sprintf('%s/xsub_out_%s_%i',dir_cache,sid1,iM));
% 
%                 s2 = e2(j,1);
%                 b2 = e2(j,2);
%                 sid2 = sprintf('%i',s2);
%                 while (length(sid2) < 5)
%                     sid2 = ['0',sid2];
%                 end
%                 sid2 = ['m',sid2];
%                 %Ca2 = load(sprintf('%s/xsub_out_%s_%i',dir_cache,sid2,iM));
% 
%                 if (s1 == s2)
%                     % If same subject
%                     Ca1 = load(sprintf('%s/xsub_out_%s_%i',dir_cache,sid1,iM));
%                     mag = Ca1.AdjMag(b1,b2);
%                     dis = Ca1.Dmat(b1,b2);
% 
%                     % remember hemisphere
%                     hemi_isr = 0;
%                     if (strcmp(Ca1.hemi_list{1},'R'))
%                         hemi_isr = 1;
%                     end
% 
%                     % check for significance
%                     %if (((~isnan(mag)) && (mag > 0)) && (dis > Ca1.dist_thresh))
%                     if (((~isnan(mag))) && (dis > Ca1.dist_thresh))
%                         fprintf('[!] %s: bchan %i, %i\n',sid1,b1,b2);
%                         E_reduce = [E_reduce; [s1,b1,b2,mag,e1(i,3:5),hemi_isr]];
%                         E_reduce_sid = [E_reduce_sid; {sid1}];
%                     end
%                 end
%                 %return
%             end
%         end

        % Plot
        [usub,uIdx] = unique(E_reduce(:,1));
        usub_str = E_reduce_sid(uIdx);
        usub_isr = E_reduce(uIdx,end);
        n_usub = length(usub);

        % colorscale calculations
        %mag_all = [E_reduce(:,4); CaR.A(l1,l2)];
        %fprintf('[*] average roi pair coherence: %.4f\n',mag_atl)
        fprintf('[*] component bipolar pair coherence mean: %.4f\n',mean(E_reduce(:,4)));
        fprintf('[*] component bipolar pair coherence median: %.4f\n',median(E_reduce(:,4)));
        fprintf('[*] component bipolar pair coherence stdev: %.4f\n',std(E_reduce(:,4)));
        fprintf('[*] component bipolar pair coherence min: %.4f\n',min(E_reduce(:,4)));
        fprintf('[*] component bipolar pair coherence max: %.4f\n',max(E_reduce(:,4)));
        fprintf('[*] component bipolar pair coherence n: %.4f\n',numel(E_reduce(:,4)));
        
        
        % Sort out which significant subjects to plot
        s_d1 = 4;
        s_d2 = 3;
        n_sub_toplot = s_d1 * s_d2;
        plot_sid = unique(E_reduce(:,1));
        if (n_sub_toplot <= length(plot_sid))
            plot_sid = plot_sid(1:n_sub_toplot);
        end
        E_reduce_mag = [];
        for ips = 1:n_sub_toplot
            sIdx = (E_reduce(:,1) == plot_sid(ips));
            E_reduce_mag = [E_reduce_mag; E_reduce(sIdx,4)];
            if (ips == length(plot_sid))
                break;
            end
        end
        
        mag_all = [E_reduce_mag; mag_atl]; %[E_reduce(:,4); mag_atl];
        n_col = 2^8;
        %cmap = corrcmap(n_col);
        cmap = inferno(n_col);
        %mag_all = Ca.mag(pass_dist & pass_ct);
        mag_all_s = sort(mag_all);
        %pcent = 0.01;
        pcent = 0;
        mag_all_s = mag_all_s(mag_all_s > 0);
        mag_max = mag_all_s(round((1-pcent)*length(mag_all_s)));
        mag_min = mag_all_s(1);
        iif = @(varargin) varargin{3-(varargin{1}>0)};
        mag2col = @(x) iif(x < mag_max, cmap(round(((x - mag_min)/(mag_max - mag_min))*(n_col-1) + 1),:), cmap(n_col,:));


        %h = figure('visible','on','Units','pixels','Position',[0 0 1920 1080]);
        h = figure('visible','off'); set(h,'PaperUnits','inches'); set(h,'PaperPosition',[0 0 8 8]);
        [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,'fsaverage_sym','r','pial'));

        % Show average plot
        %const_elec_surface = 0.7*[1 1 1];
        const_col_surface = 0.9*[1 1 1];
        const_alpha_surface = 0.8; % 0.4
        const_elec_surface = 0*[1 1 1];
        const_elec_low = 0.45*[1 1 1]; % 0.4

        %[ha, pos] = tight_subplot(1,n_usub,[.01 .01],[.01 .01],[.01 .01]);
        
        [ha, pos] = tight_subplot(s_d1,s_d2,[.01 .01],[.01 .01],[.01 .01]);

%         % Plot the location-location interaction on the average brain
%         axes(ha(1));
%         p = trisurf(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
%             'EdgeColor','none','FaceColor',const_col_surface,'FaceAlpha',const_alpha_surface);
%         hold all;
%         brainlight;
%         view(90,0)
%         axis off
% 
%         % plot locatons
%         coord = CaR.E(:,3:5);
%         [n_chan,~] = size(coord);
%         [Xe,Ye,Ze] = sphere(20);
%         elec_pa = cell(1,n_chan);
%         const_elec_radius = 2;
%         for j = 1:n_chan
%             X = Xe*const_elec_radius + coord(j,1);
%             Y = Ye*const_elec_radius + coord(j,2);
%             Z = Ze*const_elec_radius + coord(j,3);
%             elec_pa{j} = surf2patch(X,Y,Z);
%         end
%         for i = 1:n_chan
%             col_e = const_elec_low;
% %             if ((i == l1) || (i == l2))
% %                 col_e = const_elec_surface;
% %             end
%             vert = elec_pa{i}.vertices;
%             q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
%                 'EdgeColor','none','FaceColor',col_e);
%             q.SpecularStrength = 0;
%             q.SpecularExponent = 1;
%         end
% 
% 



%         % show line
%         tube_n_edges = 3;
%         %tube_radius = 0.4;
%         col_path = mag2col(CaR.A(l1,l2));%[0 0 1];
%         fprintf('Average coherence: %.6f\n',CaR.A(l1,l2));
% 
%         pIdx = PathsA(l1,l2);
%         path = Paths{pIdx};
%         if (trig_show_path ~= 1)
%             path = [path(1,:); path(end,:)];
%         end
%         [n_path,~] = size(path);
%         [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
%         fv = surf2patch(X,Y,Z);
%         p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
%             'EdgeColor','none','FaceColor',col_path); hold on;
%         p2.SpecularStrength = 0;
%         p2.SpecularExponent = 1;


        %axes(ha(1));

        

        if (imm == 1)
            save('./cache_figure_t16_dk','usub','usub_str');
            %return
        end


        % Plot the individual subjects
        n_sig_lines = 0;
        n_sig_lines_sid = {};
        for i = (1+1):(n_usub+1)
            sid = usub_str{i-1};
            E_reduce_s = E_reduce(usub(i-1)==E_reduce(:,1),:);
            %E_reduce_s_a = E_reduce_a(usub(i-1)==E_reduce_a(:,1),:);

            if ((i-1) > (s_d1*s_d2))
                break;
            end

            % subaxis
            axes(ha(i-1));
            

            if (usub_isr(i-1) == 1)
                hemi = 'r';
                view_ang = 90;
            elseif (usub_isr(i-1) == 0)
                hemi = 'l';
                view_ang = -90;
            end

            % read subject surface 
            [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,'pial'));

            % Highlight areas
            if (trig_color_rois)
                [v_a, l, ct] = read_annotation(sprintf('%s/%s/label/%sh.aparc.annot',subjects_dir,sid,hemi));
                [n_vert,~] = size(s_vert);
                const_col_surface2 = repmat(const_col_surface,[n_vert 1]);
                
                for ihi = 1:ct.numEntries
                    cond_to_color = strcmp(ct.struct_names{ihi},b1_roi) | strcmp(ct.struct_names{ihi},b2_roi);
                    if (cond_to_color)
                        roi_id = ct.table(ihi,end);
                        roi_col = ct.table(ihi,1:3) / 255;
                        n_vtx = sum(l==roi_id);
                        const_col_surface2(l==roi_id,:) = repmat(roi_col,[n_vtx 1]);
                    end
                end
                
                p = trisurf(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
                    'EdgeColor','none','FaceVertexCData',const_col_surface2,'FaceAlpha',const_alpha_surface);
                %return
            else
                p = trisurf(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
                    'EdgeColor','none','FaceColor',const_col_surface,'FaceAlpha',const_alpha_surface);
            end
            
            
            hold all;
            brainlight;
            view(view_ang,0)
            axis off








            % Load subject paths
            fn_paths = sprintf('brainexport/%s_6all_%i.mat',sid,iM);
            fprintf('[*] Loading %s ..\n',fn_paths)
            CaP = load(fn_paths);
            Paths = CaP.Paths;
            Paths_ind = CaP.Paths_ind;
            Paths = Paths(Paths_ind);
            %[n_loc,~] = size(CaP.E);
            n_loc = CaP.ecog.n_bchan;

            PathsA = zeros(n_loc,n_loc);
            count = 1;
            for i2 = 1:(n_loc-1)
                for j = (i2+1):n_loc
                    PathsA(i2,j) = count;
                    PathsA(j,i2) = count;
                    count = count + 1;
                end
            end


            % plot interactions

            % show line
            tube_n_edges = 3;
            %tube_radius = 0.4;

            [n_int,~] = size(E_reduce_s);

            for k = 1:n_int
                if (E_reduce_s(k,4) > 0)
                    col_path = mag2col(E_reduce_s(k,4));
                    fprintf('[%i] %s coherence: %.6f\n',i,sid,E_reduce_s(k,4));
                    l1b = E_reduce_s(k,2);
                    l2b = E_reduce_s(k,3);
                    pIdx = PathsA(l1b,l2b);
                    path = Paths{pIdx};
                    if (trig_show_path ~= 1)
                        path = [path(1,:); path(end,:)];
                    end
                    [n_path,~] = size(path);
                    [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
                    fv = surf2patch(X,Y,Z);
                    p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
                        'EdgeColor','none','FaceColor',col_path); hold on;
                    p2.SpecularStrength = 0;
                    p2.SpecularExponent = 1;
                    n_sig_lines = n_sig_lines + 1;
                    n_sig_lines_sid = [n_sig_lines_sid; {sid}];
                end
            end
            bchans = unique(E_reduce_s(:,2:3));



            % Plot electrodes
            % Read subject electrodes
            l = read_label(sprintf('%s/%s',subjects_dir,sid),sprintf('all_surf_ielvis'));
            if (isempty(l))
                l = read_label(sprintf('%s',sid),sprintf('all_surf_ielvis'));
            end
            [~,sIdx] = sort(l(:,end));
            l = l(sIdx,:); % zero-indexed
            [n_chan,~] = size(l);
            [Xe,Ye,Ze] = sphere(20);
            elec_pa = cell(1,n_chan);
            const_elec_radius = 2;
            elec_ras = nan(n_chan,3);
            for j = 1:n_chan
                X = Xe*const_elec_radius + s_vert(l(j,1)+1,1);
                Y = Ye*const_elec_radius + s_vert(l(j,1)+1,2);
                Z = Ze*const_elec_radius + s_vert(l(j,1)+1,3);
                elec_pa{j} = surf2patch(X,Y,Z);
                elec_ras(j,1:3) = s_vert(l(j,1)+1,1:3);
            end
            
            roisP = CaP.Ca.C.AtlLabels{2};
            for k = 1:n_chan
                vert = elec_pa{k}.vertices;
                col_e = const_elec_low;
                
                % highlight electrodes in areas
                is_elec_hi = (strcmp(roisP{k},b1_roi) || strcmp(roisP{k},b2_roi));
                if (is_elec_hi)
                    col_hiroi = ct.table(strcmp(ct.struct_names,roisP{k}),1:3)/255;
                    col_e = const_elec_surface * 0.75 + col_hiroi * 0.25;
                end

                b1 = find(k == CaP.ecog.bip(:,1));
                if ((~isempty(b1)) && any(b1 == bchans))
                    col_e = const_elec_surface;
                end
                
                % Highlight original Subject 3 electrodes
                if (strcmp(CaP.sid,sid_const) && (( CaP.ecog.bip(b1_const,1) == k ) || ( CaP.ecog.bip(b2_const,1) == k )) )
                    col_e = 0.95*[1 1 1];
                end
                
                % Putative language cortex electrodes
                if (trig_show_esm_elecs)
                    cond_isesm = false;
                    n_esm = length(Esm.esm_chans);
                    sid_int = find(strcmp(CaA.Subjects,sid));
                    for i_esm = 1:n_esm
                        esm_is_el = (Esm.esm_chans(i_esm) == k);
                        esm_is_sub = (Esm.esm_subs(i_esm) == (sid_int));
                        cond_isesm = cond_isesm | (esm_is_el & esm_is_sub);
                    end
                    if (cond_isesm)
                        col_e = [1 0 0];
                    end
                end
                
                % Show bipolar electrodes only
                is_bip = (sum(CaP.ecog.bip(:,1) == k) > 0);
                if (is_bip)
                    if (flag_cone) %((i == b1c1) || (i == b2c1))
                        %CaTmp = load([subjects_dir,'/',sid,'/label/all_parcellation.mat']);
                        %ec = CaTmp.EleCoords;
                        cone_radius = [1.8 0.1];
                        cone_n = 32;
                        %cone_col = [1 1 1]*0.15;
                        cone_col = col_e;
                        cone_rref = 1.1;
                        
                        %bc = find(CaP.Ca.ecog.bip(:,1)==k);
                        %bc1 = CaP.Ca.ecog.bip(bc,1);
                        %bc2 = CaP.Ca.ecog.bip(bc,2);
                        %coord_1 = ec(ec(:,5)==bc1,2:4);
                        %coord_2 = ec(ec(:,5)==bc2,2:4);
                        
                        ecog_local = H5eeg(['/media/jerry/untitled/h5_notch20/',sid,'.h5']);
                        cone_bip = ecog_local.bip;
                        coord_1 = elec_ras(cone_bip(cone_bip(:,1)==k,1),:); %cone_bip(cone_bip(:,1)==k,4:6);
                        coord_2 = elec_ras(cone_bip(cone_bip(:,1)==k,2),:); %cone_bip(cone_bip(:,1)==k,7:9);
                        [ConeH,EndPlate1,EndPlate2] = Cone(coord_1,coord_2,cone_radius,cone_n,cone_col,1,0);
                        ConeH.SpecularStrength = 0;
                        ConeH.SpecularExponent = 1;
                        % ref marker
                        [X,Y,Z] = sphere(cone_n);
                        q = surf(coord_2(1)+X*cone_rref,coord_2(2)+Y*cone_rref,coord_2(3)+Z*cone_rref,'EdgeColor','none','FaceColor',cone_col);
                        q.SpecularStrength = 0;
                        q.SpecularExponent = 1;
                    end
                    
                    q = trisurf(elec_pa{k}.faces,vert(:,1),vert(:,2),vert(:,3),...
                        'EdgeColor','none','FaceColor',col_e);
                    q.SpecularStrength = 0;
                    q.SpecularExponent = 1;
                end
            end
            

            % show patient number
            if (trig_show_subject_num)
                sid_int = find(strcmp(CaA.Subjects,sid));
                text(min(s_vert(:,1)),min(s_vert(:,2)),min(s_vert(:,3)),...
                    sprintf('%i',sid_int),'FontSize',10);
            end


            if ((i == 2) && (trig_show_colorbar))
                % Show colorbar
                colormap(cmap);
                cb = colorbar;
                color_pts = linspace(mag_min,mag_max,5);
                color_str = cell(1,length(color_pts));
                for i2 = 1:length(color_pts)
                    color_str{i2} = sprintf('%.3f',color_pts(i2));
                end
                caxis([mag_min, mag_max])
                set(cb,'TickLength',0);
                set(cb,'Location','south');
                set(cb,'Ticks',color_pts,'TickLabels',color_str);
            end
        end

    %     axes(ha(i+1));
    %     axis off;


        % Save figure
        print(h,sprintf('figures/T16/figure_t16_%i_curve-%i_DK_%i_cone-%i',iM,trig_show_path,trig_show_subject_num,flag_cone),'-depsc','-r400');
        print(h,sprintf('figures/T16/figure_t16_%i_curve-%i_DK_%i_cone-%i',iM,trig_show_path,trig_show_subject_num,flag_cone),'-dpng','-r600');

        close(h);

        
        %%
        mag_atl_l_sig = mag_atl_l(mag_atl_l > 0);
        n_sig_lines2 = sum(mag_atl_l > 0);
        n_sig_lines_sid2 = mag_atl_l_sid(mag_atl_l > 0); % length(unique(
        n_sig_lines_sid2_nosig = mag_atl_l_sid((mag_atl_l == 0)); % not significant (ignore nocov)
        [n_Er,~] = size(E_reduce_a);
        n_usub = length(unique(E_reduce_sid_a));
        fprintf('Number of pairs: %i\n',n_Er);
        fprintf('Number of subjects: %i\n',n_usub);
        fprintf('Consistency across pairs: %.4f (%i of %i)\n',...
            n_sig_lines2/n_Er,n_sig_lines2,n_Er);
        n_usub_sig = length(unique(n_sig_lines_sid2));
        fprintf('Consistency across subjects: %.4f (%i of %i)\n',...
            n_usub_sig/n_usub,n_usub_sig,n_usub);
        
        % Left/Right
        % DEBUG: CaA.AdjAtl_sid is publication subject IDs
%         a = reshape(CaA.AdjAtl_sid,1,[]);
%         a_all = [];
%         for ia = 1:length(a)
%             at = a{ia};
%             a_all = [a_all, at];
%         end
%         sid_a_all = unique(a_all);
        
        usubLR = unique(n_sig_lines_sid2);
        is_L = false(1,length(usubLR));
        is_R = false(1,length(usubLR));
        pairs_isL = nan(1,length(n_sig_lines_sid2));
        for ius = 1:length(usubLR)
            sid = CaA.Subjects{usubLR(ius)};
            CaTmp = load(sprintf('./cache/xsub_out_%s_1.mat',sid));
            
            % CaA.Subjects(usubLR) - does not include m00049, can assume
            % all electrodes are either L or R
            is_L(ius) = (upper(CaTmp.C.EleHemi{1}) == 'L');
            is_R(ius) = (upper(CaTmp.C.EleHemi{1}) == 'R');
            
            % count pairs L/R
            pairs_isL(n_sig_lines_sid2 == usubLR(ius)) = (upper(CaTmp.C.EleHemi{1}) == 'L');
        end
        
        fprintf('Of %i significant subjects, %i had left coverage, %i had right (%.2f%%)\n',...
            numel(is_L),sum(is_L),sum(is_R),100*sum(is_L)/numel(is_L));
        fprintf('Of %i significant pairs, %i had left coverage, %i had right (%.2f %%)\n',...
            numel(pairs_isL),sum(pairs_isL),sum(~pairs_isL),100*sum(pairs_isL)/numel(pairs_isL))
        
        % L/R coherence average differences
        coh_sig_L_mean = mean(mag_atl_l_sig(pairs_isL==1));
        fprintf('coh_sig_L_mean: %.6f\n',coh_sig_L_mean);
        coh_sig_L_len = length(mag_atl_l_sig(pairs_isL==1));
        fprintf('coh_sig_L_len: %i\n',coh_sig_L_len);
        coh_sig_L_std = std(mag_atl_l_sig(pairs_isL==1));
        fprintf('coh_sig_L_std: %.6f\n',coh_sig_L_std);
        
        coh_sig_R_mean = mean(mag_atl_l_sig(pairs_isL==0));
        fprintf('coh_sig_R_mean: %.6f\n',coh_sig_R_mean);
        coh_sig_R_len = length(mag_atl_l_sig(pairs_isL==0));
        fprintf('coh_sig_R_len: %i\n',coh_sig_R_len);
        coh_sig_R_std = std(mag_atl_l_sig(pairs_isL==0));
        fprintf('coh_sig_R_std: %.6f\n',coh_sig_R_std);
        
        [p,h,stats] = ranksum(mag_atl_l_sig(pairs_isL==0), mag_atl_l_sig(pairs_isL==1));
        fprintf('Ranksum p: %.6f, zval: %.6f\n',p,stats.zval)
        %return
        
        % === Same thing, but for unsignificant subs and pairs ===========
        %usubLR_nosig = unique(n_sig_lines_sid2_nosig);
        usubLR = unique(n_sig_lines_sid2_nosig);
        is_L = false(1,length(usubLR));
        is_R = false(1,length(usubLR));
        pairs_isL = nan(1,length(n_sig_lines_sid2_nosig));
        for ius = 1:length(usubLR)
            sid = CaA.Subjects{usubLR(ius)};
            CaTmp = load(sprintf('./cache/xsub_out_%s_1.mat',sid));
            
            % CaA.Subjects(usubLR) - does not include m00049, can assume
            % all electrodes are either L or R
            is_L(ius) = (upper(CaTmp.C.EleHemi{1}) == 'L');
            is_R(ius) = (upper(CaTmp.C.EleHemi{1}) == 'R');
            
            % count pairs L/R
            pairs_isL(n_sig_lines_sid2_nosig == usubLR(ius)) = (upper(CaTmp.C.EleHemi{1}) == 'L');
        end
        
        fprintf('Of %i not significant subjects, %i had left coverage, %i had right (%.2f%%)\n',...
            numel(is_L),sum(is_L),sum(is_R),100*sum(is_L)/numel(is_L));
        fprintf('Of %i not significant pairs, %i had left coverage, %i had right (%.2f %%)\n',...
            numel(pairs_isL),sum(pairs_isL),sum(~pairs_isL),100*sum(pairs_isL)/numel(pairs_isL))
        
        
    end
    
end




return

% == Broadband ==
% Number of pairs: 441
% Number of subjects: 29
% Consistency across pairs: 0.1111 (49 of 441)
% Consistency across subjects: 0.3448 (10 of 29)
% Of 10 significant subjects, 5 had left coverage, 5 had right (50.00%)
% Of 49 significant pairs, 30 had left coverage, 19 had right (61.22 %)
% coh_sig_L_mean: 0.331574
% coh_sig_L_len: 30
% coh_sig_L_std: 0.145328
% coh_sig_R_mean: 0.247583
% coh_sig_R_len: 19
% coh_sig_R_std: 0.096179
% Ranksum p: 0.030402, zval: -2.164814
% Of 29 not significant subjects, 16 had left coverage, 13 had right (55.17%)
% Of 392 not significant pairs, 251 had left coverage, 141 had right (64.03 %)

% == Gamma ==
% Number of pairs: 441
% Number of subjects: 29
% Consistency across pairs: 0.1156 (51 of 441)
% Consistency across subjects: 0.3793 (11 of 29)
% Of 11 significant subjects, 5 had left coverage, 6 had right (45.45%)
% Of 51 significant pairs, 28 had left coverage, 23 had right (54.90 %)
% coh_sig_L_mean: 0.312207
% coh_sig_L_len: 28
% coh_sig_L_std: 0.145659
% coh_sig_R_mean: 0.244975
% coh_sig_R_len: 23
% coh_sig_R_std: 0.101044
% Ranksum p: 0.113961, zval: -1.580638
% Of 29 not significant subjects, 16 had left coverage, 13 had right (55.17%)
% Of 390 not significant pairs, 253 had left coverage, 137 had right (64.87 %)





% 
% for imm = 1:length(imetrics)
%     iM = imetrics(imm);
%     % Read subject surface
%     sid = 'm00005';
%     fn_paths = sprintf('brainexport/%s_6all_%i.mat',sid,iM);
%     fprintf('[*] Loading %s ..\n',fn_paths)
%     load(fn_paths);
%     Paths = Paths(Paths_ind);
%     fn_cache = sprintf('%s/xsub_out_%s_%i.mat',dir_cache,sid,iM);
%     Ca = load(fn_cache);
%     currdir = pwd;
%     surface_type = 'pial';
%     hemi = lower(Ca.roi_b1c1_hemi);
%     [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,surface_type));
%     faces = faces + 1;
% 
%     % Read subject electrodes
%     l = read_label(sprintf('%s/%s',subjects_dir,sid),sprintf('all_surf_ielvis'));
%     if (isempty(l))
%         l = read_label(sprintf('%s',sid),sprintf('all_surf_ielvis'));
%     end
%     [~,sIdx] = sort(l(:,end));
%     l = l(sIdx,:); % zero-indexed
%     [n_chan,~] = size(l);
%     [Xe,Ye,Ze] = sphere(20);
%     elec_pa = cell(1,n_chan);
%     const_elec_radius = 2;
%     for i = 1:n_chan
%         X = Xe*const_elec_radius + s_vert(l(i,1)+1,1);
%         Y = Ye*const_elec_radius + s_vert(l(i,1)+1,2);
%         Z = Ze*const_elec_radius + s_vert(l(i,1)+1,3);
%         elec_pa{i} = surf2patch(X,Y,Z);
%     end
% 
%     %colormap
%     n_col = 2^8;
%     %cmap = corrcmap(n_col);
%     cmap = inferno(n_col);
%     pass_dist = (Ca.Dmats > Ca.dist_thresh);
%     pass_ct = (Ca.ct > Ca.ct_thresh);
%     mag_all = Ca.mag(pass_dist & pass_ct);
%     mag_all_s = sort(mag_all);
%     pcent = 0.01;
% %     mag_max = max(mag_all);
% %     mag_min = min(mag_all);
%     mag_max = mag_all_s(round((1-pcent)*length(mag_all_s)));
%     mag_min = mag_all_s(1);
%     iif = @(varargin) varargin{3-(varargin{1}>0)};
%     %mag2col = @(x) cmap(round(((x - mag_min)/(mag_max - mag_min))*(n_col-1) + 1),:);
%     mag2col = @(x) iif(x < mag_max, cmap(round(((x - mag_min)/(mag_max - mag_min))*(n_col-1) + 1),:), cmap(n_col,:));
% 
%     % color constants
%     const_col_surface = 0.85*[1 1 1];
%     const_alpha_surface = 0.4;
%     const_elec_surface = 0*[1 1 1];
%     const_elec_low = 0.4*[1 1 1];
% 
% 
%     h = figure('visible','off','Units','pixels','Position',[0 0 1920 1080]);
%     [ha, pos] = tight_subplot(2,3,[.01 .01],[.01 .01],[.01 .01]);
% 
% 
% 
%     
%     
%     
%     
%     
%     
%     %subplot(2,3,4)
%     axes(ha(4));
%     hold all;
%     p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
%         'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
%     for i = 1:n_chan
%         vert = elec_pa{i}.vertices;
%     %     if ((i == b1c1) || (i == b2c1))
%     %         col_elec = const_elec_surface;
%     %     else
%     %         col_elec = const_elec_low;
%     %     end
%         q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
%             'EdgeColor','none','FaceColor',const_elec_surface);
%         q.SpecularStrength = 0;
%         q.SpecularExponent = 1;
%     end
%     % --- Plot all interactions to both bchan ------------------------
%     b1 = 35;
%     b2 = 84;
%     tube_n_edges = 3;
%     %tube_radius = 1.5; %0.4;
%     col_path = [0 0 1];
%     count = 1;
%     mag_range = [];
%     for i = 1:(ecog.n_bchan - 1)
%         for j = (i+1):ecog.n_bchan
%             cond_fwd = true; %((i == b1) || (j == b2));
%             cond_rev = true; %((i == b2) || (j == b1));
%             if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
%                 path = Paths{count};
%                 [n_path,~] = size(path);
%                 if (Ca.mag(count) < mag_max)
%                     col_path = mag2col(Ca.mag(count));
%                 else
%                     col_path = cmap(end,:);
%                 end
%                 mag_range = [mag_range; Ca.mag(count)];
%                 [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
%                 fv = surf2patch(X,Y,Z);
%                 p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
%                     'EdgeColor','none','FaceColor',col_path); hold on;
%                 p2.SpecularStrength = 0;
%                 p2.SpecularExponent = 1;
%                 %break;
%             end
%             count = count + 1;
%         end
%     end
%     fprintf('[!] Mag range for subplot 4: interaction to both: %.9f - %.9f\n',...
%         min(mag_range),max(mag_range));
%     brainlight;
%     view(90,0)
%     axis off
%     
%     
%     
%     
%     %subplot(2,3,1)
%     axes(ha(1));
% 
%     hold all;
%     p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
%         'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
%     b1 = 35;
%     b2 = 84;
%     b1c1 = ecog.bip(b1,1);
%     b2c1 = ecog.bip(b2,1);
%     for i = 1:n_chan
%         vert = elec_pa{i}.vertices;
%         if ((i == b1c1) || (i == b2c1))
%             col_elec = const_elec_surface;
%         else
%             col_elec = const_elec_low;
%         end
%         q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
%             'EdgeColor','none','FaceColor',col_elec);
%         q.SpecularStrength = 0;
%         q.SpecularExponent = 1;
%     end
%     % --- Plot one interaction ------------------------------------------------
%     tube_n_edges = 9;
%     %tube_radius = 1.5;
%     %col_path = [0 0 1];
%     count = 1;
%     for i = 1:(ecog.n_bchan - 1)
%         for j = (i+1):ecog.n_bchan
%             cond_fwd = ((i == b1) && (j == b2));
%             cond_rev = ((i == b2) && (j == b1));
%             if ( cond_fwd || cond_rev )
%                 path = Paths{count};
%                 [n_path,~] = size(path);
%                 [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
%                 fv = surf2patch(X,Y,Z);
%                 %col_path = mag2col(Ca.mag(count));
%                 if (Ca.mag(count) < mag_max)
%                     col_path = mag2col(Ca.mag(count));
%                 else
%                     col_path = cmap(end,:);
%                 end
%                 p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
%                     'EdgeColor','none','FaceColor',col_path); hold on;
%                 p2.SpecularStrength = 0;
%                 p2.SpecularExponent = 1;
%                 break;
%             end
%             count = count + 1;
%         end
%     end
%     brainlight;
%     view(90,0)
%     axis off
% 
% 
%     %subplot(2,3,2)
%     axes(ha(2));
%     hold all;
%     p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
%         'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
%     for i = 1:n_chan
%         vert = elec_pa{i}.vertices;
%         if ((i == b1c1) )
%             col_elec = const_elec_surface;
%         else
%             col_elec = const_elec_low;
%         end
%         q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
%             'EdgeColor','none','FaceColor',col_elec);
%         q.SpecularStrength = 0;
%         q.SpecularExponent = 1;
%     end
%     % --- Plot all interactions to one bchan ----------------------------------
%     b1 = 35;
%     b2 = 84;
%     tube_n_edges = 9;
%     %tube_radius = 1.5;
%     col_path = [0 0 1];
%     count = 1;
%     for i = 1:(ecog.n_bchan - 1)
%         for j = (i+1):ecog.n_bchan
%             cond_fwd = (i == b1); %((i == b1) && (j == b2));
%             cond_rev = (j == b1); %((i == b2) && (j == b1));
%             if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
%                 path = Paths{count};
%                 [n_path,~] = size(path);
%                 %col_path = mag2col(Ca.mag(count));
%                 if (Ca.mag(count) < mag_max)
%                     col_path = mag2col(Ca.mag(count));
%                 else
%                     col_path = cmap(end,:);
%                 end
%                 [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
%                 fv = surf2patch(X,Y,Z);
%                 p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
%                     'EdgeColor','none','FaceColor',col_path); hold on;
%                 p2.SpecularStrength = 0;
%                 p2.SpecularExponent = 1;
%                 %break;
%             end
%             count = count + 1;
%         end
%     end
%     brainlight;
%     view(90,0)
%     axis off
% 
% 
%     %subplot(2,3,3)
%     axes(ha(3));
%     hold all;
%     p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
%         'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
%     for i = 1:n_chan
%         vert = elec_pa{i}.vertices;
%         if ( (i == b2c1))
%             col_elec = const_elec_surface;
%         else
%             col_elec = const_elec_low;
%         end
%         q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
%             'EdgeColor','none','FaceColor',col_elec);
%         q.SpecularStrength = 0;
%         q.SpecularExponent = 1;
%     end
%     % --- Plot all interactions to the other one bchan ------------------------
%     b1 = 35;
%     b2 = 84;
%     tube_n_edges = 9;
%     %tube_radius = 1.5;
%     col_path = [0 0 1];
%     count = 1;
%     for i = 1:(ecog.n_bchan - 1)
%         for j = (i+1):ecog.n_bchan
%             cond_fwd = (i == b2); %((i == b1) && (j == b2));
%             cond_rev = (j == b2); %((i == b2) && (j == b1));
%             if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
%                 path = Paths{count};
%                 [n_path,~] = size(path);
%                 %col_path = mag2col(Ca.mag(count));
%                 if (Ca.mag(count) < mag_max)
%                     col_path = mag2col(Ca.mag(count));
%                 else
%                     col_path = cmap(end,:);
%                 end
%                 [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
%                 fv = surf2patch(X,Y,Z);
%                 p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
%                     'EdgeColor','none','FaceColor',col_path); hold on;
%                 p2.SpecularStrength = 0;
%                 p2.SpecularExponent = 1;
%                 %break;
%             end
%             count = count + 1;
%         end
%     end
%     brainlight;
%     view(90,0)
%     axis off
% 
% 
% 
% 
% 
% 
% 
%     
% 
% 
% 
%     %subplot(2,3,5)
%     axes(ha(5));
%     hold all;
%     p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
%         'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
%     b1c1 = ecog.bip(b1,1);
%     b2c1 = ecog.bip(b2,1);
%     roi_labels = Ca.C.AtlLabels{2};
%     b1_roi = roi_labels{b1c1};
%     b2_roi = roi_labels{b2c1};
%     for i = 1:n_chan
%         vert = elec_pa{i}.vertices;
%         b1_roi_t = roi_labels{i};
%         if (((i == b1c1) ) || strcmp(b1_roi,b1_roi_t))
%             col_elec = const_elec_surface;
%         else
%             col_elec = const_elec_low;
%         end
%         q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
%             'EdgeColor','none','FaceColor',col_elec);
%         q.SpecularStrength = 0;
%         q.SpecularExponent = 1;
%     end
%     % hold all;
%     % p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
%     %     'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
%     % for i = 1:n_chan
%     %     vert = elec_pa{i}.vertices;
%     %     if ((i == b1c1) || (i == b2c1))
%     %         col_elec = const_elec_surface;
%     %     else
%     %         col_elec = const_elec_low;
%     %     end
%     %     q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
%     %         'EdgeColor','none','FaceColor',col_elec);
%     %     q.SpecularStrength = 0;
%     %     q.SpecularExponent = 1;
%     % end
%     % --- Plot all interactions in bchan 1 roi ------------------------
%     b1 = 35;
%     b2 = 84;
%     roi_labels = Ca.C.AtlLabels{2};
%     b1c1 = ecog.bip(b1,1);
%     b2c1 = ecog.bip(b2,1);
%     b1_roi = roi_labels{b1c1};
%     b2_roi = roi_labels{b2c1};
%     tube_n_edges = 3;
%     %tube_radius = 0.8;
%     col_path = [0 0 1];
%     count = 1;
%     for i = 1:(ecog.n_bchan - 1)
%         for j = (i+1):ecog.n_bchan
%             b1c1 = ecog.bip(i,1);
%             b2c1 = ecog.bip(j,1);
%             b1_roi_t = roi_labels{b1c1};
%             b2_roi_t = roi_labels{b2c1};
% 
%             cond_fwd = (strcmp(b1_roi,b1_roi_t)); %|| strcmp(b1_roi,b1_roi_t)); %((i == b1) || (j == b2));
%             cond_rev = (strcmp(b1_roi,b2_roi_t)); %|| strcmp(b1_roi,b2_roi_t)); %((i == b2) || (j == b1));
%             if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
%                 path = Paths{count};
%                 [n_path,~] = size(path);
%                 %col_path = mag2col(Ca.mag(count));
%                 if (Ca.mag(count) < mag_max)
%                     col_path = mag2col(Ca.mag(count));
%                 else
%                     col_path = cmap(end,:);
%                 end
%                 [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
%                 fv = surf2patch(X,Y,Z);
%                 p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
%                     'EdgeColor','none','FaceColor',col_path); hold on;
%                 p2.SpecularStrength = 0;
%                 p2.SpecularExponent = 1;
%                 %break;
%             end
%             count = count + 1;
%         end
%     end
%     brainlight;
%     view(90,0)
%     axis off
% 
%     % Show colorbar
%     colormap(cmap);
%     cb = colorbar;
%     color_pts = linspace(mag_min,mag_max,5);
%     color_str = cell(1,length(color_pts));
%     for i = 1:length(color_pts)
%         color_str{i} = sprintf('%.2f',color_pts(i));
%     end
%     caxis([mag_min, mag_max])
%     set(cb,'TickLength',0);
%     set(cb,'Location','south');
%     set(cb,'Ticks',color_pts,'TickLabels',color_str);
% 
% 
%     %subplot(2,3,6)
%     axes(ha(6));
%     hold all;
%     p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
%         'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
%     b1c1 = ecog.bip(b1,1);
%     b2c1 = ecog.bip(b2,1);
%     roi_labels = Ca.C.AtlLabels{2};
%     b1_roi = roi_labels{b1c1};
%     b2_roi = roi_labels{b2c1};
%     for i = 1:n_chan
%         vert = elec_pa{i}.vertices;
%         b1_roi_t = roi_labels{i};
%         if (( (i == b2c1)) || strcmp(b2_roi,b1_roi_t))
%             col_elec = const_elec_surface;
%         else
%             col_elec = const_elec_low;
%         end
%         q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
%             'EdgeColor','none','FaceColor',col_elec);
%         q.SpecularStrength = 0;
%         q.SpecularExponent = 1;
%     end
%     % --- Plot all interactions in bchan 1 roi ------------------------
%     b1 = 35;
%     b2 = 84;
%     roi_labels = Ca.C.AtlLabels{2};
%     b1c1 = ecog.bip(b1,1);
%     b2c1 = ecog.bip(b2,1);
%     b1_roi = roi_labels{b1c1};
%     b2_roi = roi_labels{b2c1};
%     tube_n_edges = 3;
%     %tube_radius = 0.8;
%     col_path = [0 0 1];
%     count = 1;
%     for i = 1:(ecog.n_bchan - 1)
%         for j = (i+1):ecog.n_bchan
%             b1c1 = ecog.bip(i,1);
%             b2c1 = ecog.bip(j,1);
%             b1_roi_t = roi_labels{b1c1};
%             b2_roi_t = roi_labels{b2c1};
% 
%             cond_fwd = (strcmp(b2_roi,b1_roi_t)); %|| strcmp(b1_roi,b1_roi_t)); %((i == b1) || (j == b2));
%             cond_rev = (strcmp(b2_roi,b2_roi_t)); %|| strcmp(b1_roi,b2_roi_t)); %((i == b2) || (j == b1));
% 
%             if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
%                 path = Paths{count};
%                 [n_path,~] = size(path);
%                 %col_path = mag2col(Ca.mag(count));
%                 if (Ca.mag(count) < mag_max)
%                     col_path = mag2col(Ca.mag(count));
%                 else
%                     col_path = cmap(end,:);
%                 end
%                 [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
%                 fv = surf2patch(X,Y,Z);
%                 p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
%                     'EdgeColor','none','FaceColor',col_path); hold on;
%                 p2.SpecularStrength = 0;
%                 p2.SpecularExponent = 1;
%                 %break;
%             end
%             count = count + 1;
%         end
%     end
%     brainlight;
%     view(90,0)
%     axis off
% 
%     fn_fig = sprintf('figures/T16/figure_t16_metric-%i_DK_%i',iM,trig_show_subject_num);
%     print(h,fn_fig,'-depsc');
%     print(h,fn_fig,'-dpng','-r300');
%     
%     close(h);
%     
%     fprintf('[*] rois: %s, %s\n',b1_roi,b2_rois);
% end


% fprintf('[*] All Done.\n');

% figure_t16_dk
% Warning: Directory already exists. 
% > In figure_t16_dk (line 12) 
% Warning: Directory already exists. 
% > In figure_t16_dk (line 13) 
% [*] metric 1
% [*] Loading brainexport/red_6all_fsaverage_1.mat ..
% [*] superiortemporal - parsopercularis: 0.299006, consistency across pairs: 0.111111 (n=441)
% [*] component bipolar pair coherence mean: 0.2990
% [*] component bipolar pair coherence median: 0.2454
% [*] component bipolar pair coherence stdev: 0.1339
% [*] component bipolar pair coherence min: 0.1475
% [*] component bipolar pair coherence max: 0.6771
% [*] Loading brainexport/m00001_6all_1.mat ..
% [2] m00001 coherence: 0.389735
% [2] m00001 coherence: 0.291535
% [2] m00001 coherence: 0.312130
% [2] m00001 coherence: 0.238188
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00001/label/all_surf_ielvis.label
% [*] Loading brainexport/m00003_6all_1.mat ..
% [3] m00003 coherence: 0.174183
% [3] m00003 coherence: 0.188159
% [3] m00003 coherence: 0.187294
% [3] m00003 coherence: 0.191479
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00003/label/all_surf_ielvis.label
% [*] Loading brainexport/m00005_6all_1.mat ..
% [4] m00005 coherence: 0.200970
% [4] m00005 coherence: 0.226948
% [4] m00005 coherence: 0.219356
% [4] m00005 coherence: 0.173139
% [4] m00005 coherence: 0.221297
% [4] m00005 coherence: 0.224905
% [4] m00005 coherence: 0.153515
% [4] m00005 coherence: 0.147452
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00005/label/all_surf_ielvis.label
% [*] Loading brainexport/m00022_6all_1.mat ..
% [5] m00022 coherence: 0.175796
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00022/label/all_surf_ielvis.label
% [*] Loading brainexport/m00025_6all_1.mat ..
% [6] m00025 coherence: 0.173750
% [6] m00025 coherence: 0.182594
% [6] m00025 coherence: 0.188519
% [6] m00025 coherence: 0.306254
% [6] m00025 coherence: 0.546720
% [6] m00025 coherence: 0.321457
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00025/label/all_surf_ielvis.label
% [*] Loading brainexport/m00026_6all_1.mat ..
% [7] m00026 coherence: 0.366665
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00026/label/all_surf_ielvis.label
% [*] Loading brainexport/m00028_6all_1.mat ..
% [8] m00028 coherence: 0.677078
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00028/label/all_surf_ielvis.label
% [*] Loading brainexport/m00045_6all_1.mat ..
% [9] m00045 coherence: 0.303511
% [9] m00045 coherence: 0.268442
% [9] m00045 coherence: 0.302794
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00045/label/all_surf_ielvis.label
% [*] Loading brainexport/m00061_6all_1.mat ..
% [10] m00061 coherence: 0.469501
% [10] m00061 coherence: 0.504897
% [10] m00061 coherence: 0.491773
% [10] m00061 coherence: 0.607767
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00061/label/all_surf_ielvis.label
% [*] Loading brainexport/m00084_6all_1.mat ..
% [11] m00084 coherence: 0.369821
% [11] m00084 coherence: 0.208354
% [11] m00084 coherence: 0.200404
% [11] m00084 coherence: 0.199672
% [11] m00084 coherence: 0.233960
% [11] m00084 coherence: 0.437001
% [11] m00084 coherence: 0.335854
% [11] m00084 coherence: 0.478321
% [11] m00084 coherence: 0.529557
% [11] m00084 coherence: 0.524800
% [11] m00084 coherence: 0.191681
% [11] m00084 coherence: 0.192209
% [11] m00084 coherence: 0.192287
% [11] m00084 coherence: 0.270909
% [11] m00084 coherence: 0.334145
% [11] m00084 coherence: 0.245414
% [11] m00084 coherence: 0.279114
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00084/label/all_surf_ielvis.label
% Number of pairs: 441
% Number of subjects: 29
% Consistency across pairs: 0.1111 (49 of 441)
% Consistency across subjects: 0.3448 (10 of 29)
% Of 10 significant subjects, 5 had left coverage, 5 had right (50.00%)
% Of 49 significant pairs, 30 had left coverage, 19 had right (61.22 %)
% [*] metric 5
% [*] Loading brainexport/red_6all_fsaverage_5.mat ..
% [*] superiortemporal - parsopercularis: 0.281887, consistency across pairs: 0.115646 (n=441)
% [*] component bipolar pair coherence mean: 0.2819
% [*] component bipolar pair coherence median: 0.2165
% [*] component bipolar pair coherence stdev: 0.1307
% [*] component bipolar pair coherence min: 0.1590
% [*] component bipolar pair coherence max: 0.7708
% [*] Loading brainexport/m00001_6all_5.mat ..
% [2] m00001 coherence: 0.370165
% [2] m00001 coherence: 0.275318
% [2] m00001 coherence: 0.282752
% [2] m00001 coherence: 0.187649
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00001/label/all_surf_ielvis.label
% [*] Loading brainexport/m00003_6all_5.mat ..
% [3] m00003 coherence: 0.199180
% [3] m00003 coherence: 0.202934
% [3] m00003 coherence: 0.197890
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00003/label/all_surf_ielvis.label
% [*] Loading brainexport/m00005_6all_5.mat ..
% [4] m00005 coherence: 0.222912
% [4] m00005 coherence: 0.216523
% [4] m00005 coherence: 0.173722
% [4] m00005 coherence: 0.206240
% [4] m00005 coherence: 0.200862
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00005/label/all_surf_ielvis.label
% [*] Loading brainexport/m00022_6all_5.mat ..
% [5] m00022 coherence: 0.211647
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00022/label/all_surf_ielvis.label
% [*] Loading brainexport/m00025_6all_5.mat ..
% [6] m00025 coherence: 0.187313
% [6] m00025 coherence: 0.201972
% [6] m00025 coherence: 0.205284
% [6] m00025 coherence: 0.334502
% [6] m00025 coherence: 0.618738
% [6] m00025 coherence: 0.357758
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00025/label/all_surf_ielvis.label
% [*] Loading brainexport/m00026_6all_5.mat ..
% [7] m00026 coherence: 0.309830
% [7] m00026 coherence: 0.305443
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00026/label/all_surf_ielvis.label
% [*] Loading brainexport/m00028_6all_5.mat ..
% [8] m00028 coherence: 0.770846
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00028/label/all_surf_ielvis.label
% [*] Loading brainexport/m00045_6all_5.mat ..
% [9] m00045 coherence: 0.217663
% [9] m00045 coherence: 0.185935
% [9] m00045 coherence: 0.230464
% [9] m00045 coherence: 0.213993
% [9] m00045 coherence: 0.181784
% [9] m00045 coherence: 0.172569
% [9] m00045 coherence: 0.170303
% [9] m00045 coherence: 0.350004
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00045/label/all_surf_ielvis.label
% [*] Loading brainexport/m00047_6all_5.mat ..
% [10] m00047 coherence: 0.158963
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00047/label/all_surf_ielvis.label
% [*] Loading brainexport/m00061_6all_5.mat ..
% [11] m00061 coherence: 0.444970
% [11] m00061 coherence: 0.480437
% [11] m00061 coherence: 0.462720
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00061/label/all_surf_ielvis.label
% [*] Loading brainexport/m00084_6all_5.mat ..
% [12] m00084 coherence: 0.317456
% [12] m00084 coherence: 0.184946
% [12] m00084 coherence: 0.168706
% [12] m00084 coherence: 0.170113
% [12] m00084 coherence: 0.192149
% [12] m00084 coherence: 0.244291
% [12] m00084 coherence: 0.356397
% [12] m00084 coherence: 0.395376
% [12] m00084 coherence: 0.281756
% [12] m00084 coherence: 0.453921
% [12] m00084 coherence: 0.514736
% [12] m00084 coherence: 0.507080
% [12] m00084 coherence: 0.208059
% [12] m00084 coherence: 0.268381
% [12] m00084 coherence: 0.186495
% [12] m00084 coherence: 0.213882
% [12] m00084 coherence: 0.203184
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00084/label/all_surf_ielvis.label
% Number of pairs: 441
% Number of subjects: 29
% Consistency across pairs: 0.1156 (51 of 441)
% Consistency across subjects: 0.3793 (11 of 29)
% Of 11 significant subjects, 5 had left coverage, 6 had right (45.45%)
% Of 51 significant pairs, 28 had left coverage, 23 had right (54.90 %)
