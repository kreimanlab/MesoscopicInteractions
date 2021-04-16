
close all;
clear;

dir_cache = './cache';
%sid_const = 'fsaverage_sym';
scale_mm_per_obj = 30;
flag_export_obj = true;
flag_plot = false;
tube_radius = 1;

mkdir('figures');
mkdir(sprintf('figures/T11_new'));

%metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
imetrics = [1 5]; % 1:5; 

flag_cone = true;

for imm = 1:length(imetrics)
    iM = imetrics(imm);
    
    for trig_show_path = [0] %[0 1]
    
        % Read subject surface
        sid = 'm00005';
        fn_paths = sprintf('brainexport/%s_6all_%i.mat',sid,iM);
        fprintf('[*] Loading %s ..\n',fn_paths)
        load(fn_paths);
        Paths = Paths(Paths_ind);
        fn_cache = sprintf('%s/xsub_out_%s_%i.mat',dir_cache,sid,iM);
        Ca = load(fn_cache);
        currdir = pwd;
        surface_type = 'pial';
        hemi = lower(Ca.roi_b1c1_hemi);
        [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,surface_type));
        faces = faces + 1;

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
        for i = 1:n_chan
            X = Xe*const_elec_radius + s_vert(l(i,1)+1,1);
            Y = Ye*const_elec_radius + s_vert(l(i,1)+1,2);
            Z = Ze*const_elec_radius + s_vert(l(i,1)+1,3);
            elec_pa{i} = surf2patch(X,Y,Z);
        end

        %colormap
        n_col = 2^8;
        %cmap = corrcmap(n_col);
        cmap = inferno(n_col);
        pass_dist = (Ca.Dmats > Ca.dist_thresh);
        pass_ct = (Ca.ct > Ca.ct_thresh);
        mag_all = Ca.mag(pass_dist & pass_ct);
        mag_all_s = sort(mag_all);
        pcent = 0.01;
    %     mag_max = max(mag_all);
    %     mag_min = min(mag_all);
        mag_max = mag_all_s(round((1-pcent)*length(mag_all_s)));
        mag_min = mag_all_s(1);
        iif = @(varargin) varargin{3-(varargin{1}>0)};
        %mag2col = @(x) cmap(round(((x - mag_min)/(mag_max - mag_min))*(n_col-1) + 1),:);
        mag2col = @(x) iif(x < mag_max, cmap(round(((x - mag_min)/(mag_max - mag_min))*(n_col-1) + 1),:), cmap(n_col,:));

        % color constants
        %const_col_surface = 0.85*[1 1 1];
        %const_alpha_surface = 0.4;
        %const_elec_surface = 0*[1 1 1];
        %const_elec_low = 0.4*[1 1 1];
        
        % Copy pasted from figure_t16_dk.m on 19Nov2019
        const_col_surface = 0.9*[1 1 1];
        const_alpha_surface = 0.8; % 0.4
        const_elec_surface = 0*[1 1 1];
        const_elec_low = 0.45*[1 1 1]; % 0.4


        h = figure('visible','off','Units','pixels','Position',[0 0 1920 1080]);
        [ha, pos] = tight_subplot(2,3,[.01 .01],[.01 .01],[.01 .01]);









        %subplot(2,3,4)
        axes(ha(4));
        hold all;
        p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
            'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
        for i = 1:n_chan
            vert = elec_pa{i}.vertices;
        %     if ((i == b1c1) || (i == b2c1))
        %         col_elec = const_elec_surface;
        %     else
        %         col_elec = const_elec_low;
        %     end
        
            is_bip = (sum(ecog.bip(:,1) == i) > 0);
            if (is_bip && flag_cone) %((i == b1c1) || (i == b2c1))
                cone_radius = [1.8 0.1];
                cone_n = 32;
                %cone_col = [1 1 1]*0.15;
                cone_col = const_elec_surface;
                cone_rref = 1.1;
                coord_1 = ecog.bip(ecog.bip(:,1)==i,4:6);
                coord_2 = ecog.bip(ecog.bip(:,1)==i,7:9);
                [ConeH,EndPlate1,EndPlate2] = Cone(coord_1,coord_2,cone_radius,cone_n,cone_col,1,0);
                ConeH.SpecularStrength = 0;
                ConeH.SpecularExponent = 1;
                % ref marker
                [X,Y,Z] = sphere(cone_n);
                q = surf(coord_2(1)+X*cone_rref,coord_2(2)+Y*cone_rref,coord_2(3)+Z*cone_rref,'EdgeColor','none','FaceColor',cone_col);
                q.SpecularStrength = 0;
                q.SpecularExponent = 1;
            end
        
            q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
                'EdgeColor','none','FaceColor',const_elec_surface);
            q.SpecularStrength = 0;
            q.SpecularExponent = 1;
        end
        % --- Plot all interactions to both bchan ------------------------
        b1 = 35;
        b2 = 84;
        tube_n_edges = 3;
        %tube_radius = 1.5; % 0.4;
        col_path = [0 0 1];
        count = 1;
        mag_range = [];
        for i = 1:(ecog.n_bchan - 1)
            for j = (i+1):ecog.n_bchan
                cond_fwd = true; %((i == b1) || (j == b2));
                cond_rev = true; %((i == b2) || (j == b1));
                if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
                    path = Paths{count};
                    if (trig_show_path ~= 1)
                        path = [path(1,:); path(end,:)];
                    end
                    [n_path,~] = size(path);
                    if (Ca.mag(count) < mag_max)
                        col_path = mag2col(Ca.mag(count));
                    else
                        col_path = cmap(end,:);
                    end
                    mag_range = [mag_range; Ca.mag(count)];
                    [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
                    fv = surf2patch(X,Y,Z);
                    p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
                        'EdgeColor','none','FaceColor',col_path); hold on;
                    p2.SpecularStrength = 0;
                    p2.SpecularExponent = 1;
                    %break;
                end
                count = count + 1;
            end
        end
        fprintf('[!] Mag range for subplot 4: interaction to both: %.9f - %.9f\n',...
            min(mag_range),max(mag_range));
        brainlight;
        view(90,0)
        axis off




        %subplot(2,3,1)
        axes(ha(1));

        hold all;
        p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
            'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
        b1 = 35;
        b2 = 84;
        b1c1 = ecog.bip(b1,1);
        b2c1 = ecog.bip(b2,1);
        b1c2 = ecog.bip(b1,2);
        b2c2 = ecog.bip(b2,2);
        col_elec_bip = [0 1 0]*0.4;
        for i = 1:n_chan
            vert = elec_pa{i}.vertices;
            
            % Show bipolar electrodes only
            is_bip = (sum(ecog.bip(:,1) == i) > 0);
            if (is_bip)
                
                if ((i == b1c1) || (i == b2c1))
                %if ((i == b1c1) || (i == b2c1) || (i == b1c2) || (i == b2c2))
                    col_elec = const_elec_surface;
                else
                    col_elec = const_elec_low;
                end
                
                % show bipolar ground indicator
                if (flag_cone) %((i == b1c1) || (i == b2c1))
                    %cone_radius = [1.8 0.1];
                    %cone_n = 32;
                    %cone_col = [1 1 1]*0.15;
                    cone_col = col_elec;
                    %cone_rref = 1.1;
                    coord_1 = ecog.bip(ecog.bip(:,1)==i,4:6);
                    coord_2 = ecog.bip(ecog.bip(:,1)==i,7:9);
                    [ConeH,EndPlate1,EndPlate2] = Cone(coord_1,coord_2,cone_radius,cone_n,cone_col,1,0);
                    ConeH.SpecularStrength = 0;
                    ConeH.SpecularExponent = 1;
                    % ref marker
                    [X,Y,Z] = sphere(cone_n);
                    q = surf(coord_2(1)+X*cone_rref,coord_2(2)+Y*cone_rref,coord_2(3)+Z*cone_rref,'EdgeColor','none','FaceColor',cone_col);
                    q.SpecularStrength = 0;
                    q.SpecularExponent = 1;
                end
                
                % plot first electrode
                q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
                'EdgeColor','none','FaceColor',col_elec);
                q.SpecularStrength = 0;
                q.SpecularExponent = 1;
               
                
            else
%                 % plot second electrode
%                 if ((i == b1c2) || (i == b2c2))
%                     %i_2 = ecog.bip(b1,2);
%                     col_elec_bip = [0 1 0]*0.4;
%                     q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
%                     'EdgeColor','none','FaceColor',col_elec_bip);
%                     q.SpecularStrength = 0;
%                     q.SpecularExponent = 1;
%                 end
            end
            
        end
        % --- Plot one interaction ------------------------------------------------
        tube_n_edges = 9;
        %tube_radius = 1.5;
        %col_path = [0 0 1];
        count = 1;
        for i = 1:(ecog.n_bchan - 1)
            for j = (i+1):ecog.n_bchan
                cond_fwd = ((i == b1) && (j == b2));
                cond_rev = ((i == b2) && (j == b1));
                if ( cond_fwd || cond_rev )
                    path = Paths{count};
                    if (trig_show_path ~= 1)
                        path = [path(1,:); path(end,:)];
                    end
                    [n_path,~] = size(path);
                    [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
                    fv = surf2patch(X,Y,Z);
                    %col_path = mag2col(Ca.mag(count));
                    if (Ca.mag(count) < mag_max)
                        col_path = mag2col(Ca.mag(count));
                        fprintf('\t[*] Subplot 2,3,1, mag: %.4f\n',Ca.mag(count))
                        fprintf('\t[*] Subplot 2,3,1, mag_std: %.4f\n',Ca.mag_std(count))
                        fprintf('\t[*] Subplot 2,3,1, mag_null: %.4f\n',Ca.mag_null(count))
                        fprintf('\t[*] Subplot 2,3,1, mag_null_std: %.4f\n',Ca.mag_null_std(count))
                        fprintf('\t[*] Subplot 2,3,1, ct: %.4f\n',Ca.ct(count))
                        fprintf('\t[*] Subplot 2,3,1, ct_null: %.4f\n',Ca.ct_null(count))
                    else
                        col_path = cmap(end,:);
                    end
                    p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
                        'EdgeColor','none','FaceColor',col_path); hold on;
                    p2.SpecularStrength = 0;
                    p2.SpecularExponent = 1;
                    break;
                end
                count = count + 1;
            end
        end
        brainlight;
        view(90,0)
        axis off


        %subplot(2,3,2)
        axes(ha(2));
        hold all;
        p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
            'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
        for i = 1:n_chan
            vert = elec_pa{i}.vertices;
            % Show bipolar electrodes only
            is_bip = (sum(ecog.bip(:,1) == i) > 0);
            if (is_bip)
                
                if ((i == b1c1) )
                %if ((i == b1c1) || (i == b1c2))
                    col_elec = const_elec_surface;
                else
                    col_elec = const_elec_low;
                end
                
                % show bipolar ground indicator
                if (flag_cone)%((i == b1c1) )
                    cone_col = col_elec;
                    coord_1 = ecog.bip(ecog.bip(:,1)==i,4:6);
                    coord_2 = ecog.bip(ecog.bip(:,1)==i,7:9);
                    [ConeH,EndPlate1,EndPlate2] = Cone(coord_1,coord_2,cone_radius,cone_n,cone_col,1,0);
                    ConeH.SpecularStrength = 0;
                    ConeH.SpecularExponent = 1;
                    % ref marker
                    [X,Y,Z] = sphere(cone_n);
                    q = surf(coord_2(1)+X*cone_rref,coord_2(2)+Y*cone_rref,coord_2(3)+Z*cone_rref,'EdgeColor','none','FaceColor',cone_col);
                    q.SpecularStrength = 0;
                    q.SpecularExponent = 1;
                end
                
                q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
                'EdgeColor','none','FaceColor',col_elec);
                q.SpecularStrength = 0;
                q.SpecularExponent = 1;
            else
%                 % plot second electrode
%                 if ((i == b1c2))
%                     %i_2 = ecog.bip(b1,2);
%                     %col_elec_bip = [0 1 0]*0.4;
%                     q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
%                     'EdgeColor','none','FaceColor',col_elec_bip);
%                     q.SpecularStrength = 0;
%                     q.SpecularExponent = 1;
%                 end
            end
            
            
        end
        % --- Plot all interactions to one bchan ----------------------------------
        b1 = 35;
        b2 = 84;
        tube_n_edges = 9;
        %tube_radius = 1.5;
        col_path = [0 0 1];
        count = 1;
        ecount = 1;
        bcount = 1;
        for i = 1:(ecog.n_bchan - 1)
            for j = (i+1):ecog.n_bchan
                cond_fwd = (i == b1); %((i == b1) && (j == b2));
                cond_rev = (j == b1); %((i == b2) && (j == b1));
                
                if (pass_dist(count) && (cond_fwd || cond_rev))
                    bcount = bcount + 1;
                end
                
                if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
                    path = Paths{count};
                    if (trig_show_path ~= 1)
                        path = [path(1,:); path(end,:)];
                    end
                    [n_path,~] = size(path);
                    %col_path = mag2col(Ca.mag(count));
                    if (Ca.mag(count) < mag_max)
                        col_path = mag2col(Ca.mag(count));
                    else
                        col_path = cmap(end,:);
                    end
                    [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
                    fv = surf2patch(X,Y,Z);
                    p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
                        'EdgeColor','none','FaceColor',col_path); hold on;
                    p2.SpecularStrength = 0;
                    p2.SpecularExponent = 1;
                    %break;
                    ecount = ecount + 1;
                end
                count = count + 1;
            end
        end
        brainlight;
        view(90,0)
        axis off
        fprintf('[*] m00003: # edges: %i, # bip: %i (%.2f%%)\n',ecount,bcount,100*(ecount/bcount));
        


        %subplot(2,3,3)
        axes(ha(3));
        hold all;
        p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
            'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
        for i = 1:n_chan
            vert = elec_pa{i}.vertices;
            % Show bipolar electrodes only
            is_bip = (sum(ecog.bip(:,1) == i) > 0);
            if (is_bip)
                
                if ( (i == b2c1))
                %if ((i == b2c1) || (i == b2c2))
                    col_elec = const_elec_surface;
                else
                    col_elec = const_elec_low;
                end
                
                % show bipolar ground indicator
                if (flag_cone) %((i == b1c1) )
                    cone_col = col_elec;
                    coord_1 = ecog.bip(ecog.bip(:,1)==i,4:6);
                    coord_2 = ecog.bip(ecog.bip(:,1)==i,7:9);
                    [ConeH,EndPlate1,EndPlate2] = Cone(coord_1,coord_2,cone_radius,cone_n,cone_col,1,0);
                    ConeH.SpecularStrength = 0;
                    ConeH.SpecularExponent = 1;
                    % ref marker
                    [X,Y,Z] = sphere(cone_n);
                    q = surf(coord_2(1)+X*cone_rref,coord_2(2)+Y*cone_rref,coord_2(3)+Z*cone_rref,'EdgeColor','none','FaceColor',cone_col);
                    q.SpecularStrength = 0;
                    q.SpecularExponent = 1;
                end
                
                q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
                'EdgeColor','none','FaceColor',col_elec);
                q.SpecularStrength = 0;
                q.SpecularExponent = 1;
            else
%                 % plot second electrode
%                 if ((i == b2c2))
%                     %i_2 = ecog.bip(b1,2);
%                     %col_elec_bip = [0 1 0]*0.4;
%                     q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
%                     'EdgeColor','none','FaceColor',col_elec_bip);
%                     q.SpecularStrength = 0;
%                     q.SpecularExponent = 1;
%                 end
            end
            
        end
        % --- Plot all interactions to the other one bchan ------------------------
        b1 = 35;
        b2 = 84;
        tube_n_edges = 9;
        %tube_radius = 1.5;
        col_path = [0 0 1];
        count = 1;
        ecount = 1;
        bcount = 1;
        for i = 1:(ecog.n_bchan - 1)
            for j = (i+1):ecog.n_bchan
                cond_fwd = (i == b2); %((i == b1) && (j == b2));
                cond_rev = (j == b2); %((i == b2) && (j == b1));
                
                if (pass_dist(count) && (cond_fwd || cond_rev))
                    bcount = bcount + 1;
                end
                
                if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
                    path = Paths{count};
                    if (trig_show_path ~= 1)
                        path = [path(1,:); path(end,:)];
                    end
                    [n_path,~] = size(path);
                    %col_path = mag2col(Ca.mag(count));
                    if (Ca.mag(count) < mag_max)
                        col_path = mag2col(Ca.mag(count));
                    else
                        col_path = cmap(end,:);
                    end
                    [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
                    fv = surf2patch(X,Y,Z);
                    p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
                        'EdgeColor','none','FaceColor',col_path); hold on;
                    p2.SpecularStrength = 0;
                    p2.SpecularExponent = 1;
                    %break;
                    ecount = ecount + 1;
                end
                count = count + 1;
            end
        end
        brainlight;
        view(90,0)
        axis off
        %fprintf('[*] m00005: num of edges: %i\n',ecount);
        fprintf('[*] m00005: # edges: %i, # bip: %i (%.2f%%)\n',ecount,bcount,100*(ecount/bcount));
        




        %subplot(2,3,5)
        axes(ha(5));
        hold all;
        p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
            'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
        b1c1 = ecog.bip(b1,1);
        b2c1 = ecog.bip(b2,1);
        roi_labels = Ca.C.AtlLabels{2};
        b1_roi = roi_labels{b1c1};
        b2_roi = roi_labels{b2c1};
        for i = 1:n_chan
            vert = elec_pa{i}.vertices;
            b1_roi_t = roi_labels{i};
            
            % Show bipolar electrodes only
            is_bip = (sum(ecog.bip(:,1) == i) > 0);
            if (is_bip)    
                if (((i == b1c1) ) || strcmp(b1_roi,b1_roi_t))
                    col_elec = const_elec_surface;
                else
                    col_elec = const_elec_low;
                end
                
                if (flag_cone) %((i == b1c1) )
                    cone_col = col_elec;
                    coord_1 = ecog.bip(ecog.bip(:,1)==i,4:6);
                    coord_2 = ecog.bip(ecog.bip(:,1)==i,7:9);
                    [ConeH,EndPlate1,EndPlate2] = Cone(coord_1,coord_2,cone_radius,cone_n,cone_col,1,0);
                    ConeH.SpecularStrength = 0;
                    ConeH.SpecularExponent = 1;
                    % ref marker
                    [X,Y,Z] = sphere(cone_n);
                    q = surf(coord_2(1)+X*cone_rref,coord_2(2)+Y*cone_rref,coord_2(3)+Z*cone_rref,'EdgeColor','none','FaceColor',cone_col);
                    q.SpecularStrength = 0;
                    q.SpecularExponent = 1;
                end
                
                q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
                'EdgeColor','none','FaceColor',col_elec);
                q.SpecularStrength = 0;
                q.SpecularExponent = 1;
            end
            
        end
        % hold all;
        % p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
        %     'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
        % for i = 1:n_chan
        %     vert = elec_pa{i}.vertices;
        %     if ((i == b1c1) || (i == b2c1))
        %         col_elec = const_elec_surface;
        %     else
        %         col_elec = const_elec_low;
        %     end
        %     q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
        %         'EdgeColor','none','FaceColor',col_elec);
        %     q.SpecularStrength = 0;
        %     q.SpecularExponent = 1;
        % end
        % --- Plot all interactions in bchan 1 roi ------------------------
        b1 = 35;
        b2 = 84;
        roi_labels = Ca.C.AtlLabels{2};
        b1c1 = ecog.bip(b1,1);
        b2c1 = ecog.bip(b2,1);
        b1_roi = roi_labels{b1c1};
        b2_roi = roi_labels{b2c1};
        tube_n_edges = 3;
        %tube_radius = 1.5; %0.8;
        col_path = [0 0 1];
        count = 1;
        ecount = 1;
        for i = 1:(ecog.n_bchan - 1)
            for j = (i+1):ecog.n_bchan
                b1c1 = ecog.bip(i,1);
                b2c1 = ecog.bip(j,1);
                b1_roi_t = roi_labels{b1c1};
                b2_roi_t = roi_labels{b2c1};

                cond_fwd = (strcmp(b1_roi,b1_roi_t)); %|| strcmp(b1_roi,b1_roi_t)); %((i == b1) || (j == b2));
                cond_rev = (strcmp(b1_roi,b2_roi_t)); %|| strcmp(b1_roi,b2_roi_t)); %((i == b2) || (j == b1));
                if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
                    path = Paths{count};
                    if (trig_show_path ~= 1)
                        path = [path(1,:); path(end,:)];
                    end
                    [n_path,~] = size(path);
                    %col_path = mag2col(Ca.mag(count));
                    if (Ca.mag(count) < mag_max)
                        col_path = mag2col(Ca.mag(count));
                    else
                        col_path = cmap(end,:);
                    end
                    [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
                    fv = surf2patch(X,Y,Z);
                    p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
                        'EdgeColor','none','FaceColor',col_path); hold on;
                    p2.SpecularStrength = 0;
                    p2.SpecularExponent = 1;
                    %break;
                    ecount = ecount + 1;
                end
                count = count + 1;
            end
        end
        brainlight;
        view(90,0)
        axis off
        fprintf('[*] m00019: num of edges: %i\n',ecount);

        % Show colorbar
        colormap(cmap);
        cb = colorbar;
        color_pts = linspace(mag_min,mag_max,5);
        color_str = cell(1,length(color_pts));
        for i = 1:length(color_pts)
            color_str{i} = sprintf('%.2f',color_pts(i));
        end
        caxis([mag_min, mag_max])
        set(cb,'TickLength',0);
        set(cb,'Location','south');
        set(cb,'Ticks',color_pts,'TickLabels',color_str);




        %subplot(2,3,6)
        axes(ha(6));
        hold all;
        p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
            'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
        b1c1 = ecog.bip(b1,1);
        b2c1 = ecog.bip(b2,1);
        roi_labels = Ca.C.AtlLabels{2};
        b1_roi = roi_labels{b1c1};
        b2_roi = roi_labels{b2c1};
        for i = 1:n_chan
            vert = elec_pa{i}.vertices;
            b1_roi_t = roi_labels{i};
                
            % Show bipolar electrodes only
            is_bip = (sum(ecog.bip(:,1) == i) > 0);
            if (is_bip)
                if (( (i == b2c1)) || strcmp(b2_roi,b1_roi_t))
                    col_elec = const_elec_surface;
                else
                    col_elec = const_elec_low;
                end
                
                if (flag_cone) %((i == b1c1) )
                    cone_col = col_elec;
                    coord_1 = ecog.bip(ecog.bip(:,1)==i,4:6);
                    coord_2 = ecog.bip(ecog.bip(:,1)==i,7:9);
                    [ConeH,EndPlate1,EndPlate2] = Cone(coord_1,coord_2,cone_radius,cone_n,cone_col,1,0);
                    ConeH.SpecularStrength = 0;
                    ConeH.SpecularExponent = 1;
                    % ref marker
                    [X,Y,Z] = sphere(cone_n);
                    q = surf(coord_2(1)+X*cone_rref,coord_2(2)+Y*cone_rref,coord_2(3)+Z*cone_rref,'EdgeColor','none','FaceColor',cone_col);
                    q.SpecularStrength = 0;
                    q.SpecularExponent = 1;
                end
                
                q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
                'EdgeColor','none','FaceColor',col_elec);
                q.SpecularStrength = 0;
                q.SpecularExponent = 1;
            end
            
        end
        % --- Plot all interactions in bchan 1 roi ------------------------
        b1 = 35;
        b2 = 84;
        roi_labels = Ca.C.AtlLabels{2};
        b1c1 = ecog.bip(b1,1);
        b2c1 = ecog.bip(b2,1);
        b1_roi = roi_labels{b1c1};
        b2_roi = roi_labels{b2c1};
        tube_n_edges = 3;
        %tube_radius = 1.5; %0.8;
        col_path = [0 0 1];
        count = 1;
        for i = 1:(ecog.n_bchan - 1)
            for j = (i+1):ecog.n_bchan
                b1c1 = ecog.bip(i,1);
                b2c1 = ecog.bip(j,1);
                b1_roi_t = roi_labels{b1c1};
                b2_roi_t = roi_labels{b2c1};

                cond_fwd = (strcmp(b2_roi,b1_roi_t)); %|| strcmp(b1_roi,b1_roi_t)); %((i == b1) || (j == b2));
                cond_rev = (strcmp(b2_roi,b2_roi_t)); %|| strcmp(b1_roi,b2_roi_t)); %((i == b2) || (j == b1));

                if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
                    path = Paths{count};
                    if (trig_show_path ~= 1)
                        path = [path(1,:); path(end,:)];
                    end
                    [n_path,~] = size(path);
                    %col_path = mag2col(Ca.mag(count));
                    if (Ca.mag(count) < mag_max)
                        col_path = mag2col(Ca.mag(count));
                    else
                        col_path = cmap(end,:);
                    end
                    [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
                    fv = surf2patch(X,Y,Z);
                    p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
                        'EdgeColor','none','FaceColor',col_path); hold on;
                    p2.SpecularStrength = 0;
                    p2.SpecularExponent = 1;
                    %break;
                    ecount = ecount + 1;
                end
                count = count + 1;
            end
        end
        brainlight;
        view(90,0)
        axis off
        fprintf('[*] m00021: num of edges: %i\n',ecount);

        fn_fig = sprintf('figures/T11_new/figure_t11_metric-%i_curve-%i_cone-%i',iM,trig_show_path,flag_cone);
        print(h,fn_fig,'-depsc');
        print(h,fn_fig,'-dpng','-r300');
        close(h);
        
        fprintf('[*] Saved to: %s\n',fn_fig);
        
        
        % --- Show details ----------------------------------
        b1 = 35;
        b2 = 84;
        
        n_cov = 0;
        n_all = 0;
        n_either = 0;
        n_both = 0;
        n_other = 0;
        n_nosig = 0;
        n_xor = 0;
        n_bchan_cov = [];
        
        count = 1;
        for i = 1:(ecog.n_bchan - 1)
            for j = (i+1):ecog.n_bchan
                cond_fwd = (i == b1);
                cond_rev = (j == b1);
                cond_fwd_2 = (i == b2);
                cond_rev_2 = (j == b2);
                
                if (pass_dist(count))
                    n_cov = n_cov + 1;
                    
                    if (pass_ct(count))
                        if ((cond_fwd || cond_rev) && (cond_fwd_2 || cond_rev_2))
                            n_both = n_both + 1;
                        elseif (cond_fwd || cond_rev) || (cond_fwd_2 || cond_rev_2)
                            n_either = n_either + 1;
                        else
                            n_other = n_other + 1;
                        end
                        
                        if (xor((cond_fwd || cond_rev),(cond_fwd_2 || cond_rev_2)))
                            n_xor = n_xor + 1;
                        end
                    else
                        n_nosig = n_nosig + 1;
                    end
                end
                
                if (((cond_fwd || cond_rev) || (cond_fwd_2 || cond_rev_2))) && (pass_dist(count))
                    n_bchan_cov = [n_bchan_cov; i; j];
                end
                
%                 if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
%                     
%                 end
                count = count + 1;
                
                n_all = n_all + 1;
            end
        end
        
        fprintf('[*] Details:\n');
        fprintf('\t# pairs total: %i\n',n_all);
        fprintf('\t# pairs covered: %i\n',n_cov);
        fprintf('\t\t# pairs sig either: %i\n',n_either);
        fprintf('\t\t# pairs sig both: %i\n',n_both);
        fprintf('\t\t# pairs sig other: %i\n',n_other);
        fprintf('\t\t# pairs not sig: %i\n',n_nosig);
        fprintf('\t# pairs exclusively to one or another: %i\n',n_xor);
        fprintf('\t# bipolar channels covered by either: %i\n',length(unique(n_bchan_cov)));
        fprintf('\t# bipolar channels: %i\n',Ca.ecog.n_bchan);
        
        %%
        fprintf('=== newest as of 04/24/2019 ===\n');
        sig1 = Ca.AdjIsSig(b1,:);
        sig2 = Ca.AdjIsSig(b2,:);
        sIdx = (isnan(sig1) | isnan(sig2)) | ((Ca.AdjIsDistOk(b1,:) == 0) | (Ca.AdjIsDistOk(b2,:) == 0));
        sig1 = sig1(~sIdx);
        sig2 = sig2(~sIdx);
        fprintf('[*] n bchan cov: %i\n',length(sig1));
        fprintf('[*] n sig both: %i\n',sum(sig1 & sig2));
        fprintf('[*] n sig either: %i\n',sum(sig1 | sig2));
        fprintf('[*] n sig no both: %i\n',sum((~ sig1) & (~ sig2)));
        fprintf('[*] n sig xor: %i\n',sum(xor(sig1,sig2)));
        fprintf('=== end comment ===\n');
    end
end


fprintf('[*] All Done.\n');


% April 5, 2020

% figure_t11_new
% Warning: Directory already exists. 
% > In figure_t11_new (line 12) 
% Warning: Directory already exists. 
% > In figure_t11_new (line 13) 
% [*] Loading brainexport/m00005_6all_1.mat ..
% [!] Mag range for subplot 4: interaction to both: 0.142108157 - 0.272080839
% 	[*] Subplot 2,3,1, mag: 0.2213
% 	[*] Subplot 2,3,1, mag_std: 0.0492
% 	[*] Subplot 2,3,1, mag_null: 0.1980
% 	[*] Subplot 2,3,1, mag_null_std: 0.0492
% 	[*] Subplot 2,3,1, ct: 0.4052
% 	[*] Subplot 2,3,1, ct_null: 0.0018
% [*] m00003: # edges: 20, # bip: 75 (26.67%)
% [*] m00005: # edges: 23, # bip: 60 (38.33%)
% [*] m00019: num of edges: 67
% [*] m00021: num of edges: 138
% [*] Saved to: figures/T11_new/figure_t11_metric-1_curve-0
% [*] Details:
% 	# pairs total: 4095
% 	# pairs covered: 3193
% 		# pairs sig either: 39
% 		# pairs sig both: 1
% 		# pairs sig other: 517
% 		# pairs not sig: 2636
% 	# pairs exclusively to one or another: 39
% 	# bipolar channels covered by either: 89
% 	# bipolar channels: 91
% === newest as of 04/24/2019 ===
% [*] n bchan cov: 44
% [*] n sig both: 8
% [*] n sig either: 18
% [*] n sig no both: 26
% [*] n sig xor: 10
% === end comment ===
% [*] Loading brainexport/m00005_6all_5.mat ..
% [!] Mag range for subplot 4: interaction to both: 0.155799359 - 0.322029501
% 	[*] Subplot 2,3,1, mag: 0.2062
% 	[*] Subplot 2,3,1, mag_std: 0.0551
% 	[*] Subplot 2,3,1, mag_null: 0.1635
% 	[*] Subplot 2,3,1, mag_null_std: 0.0551
% 	[*] Subplot 2,3,1, ct: 0.4707
% 	[*] Subplot 2,3,1, ct_null: 0.0017
% [*] m00003: # edges: 20, # bip: 75 (26.67%)
% [*] m00005: # edges: 20, # bip: 60 (33.33%)
% [*] m00019: num of edges: 29
% [*] m00021: num of edges: 80
% [*] Saved to: figures/T11_new/figure_t11_metric-5_curve-0
% [*] Details:
% 	# pairs total: 4095
% 	# pairs covered: 3193
% 		# pairs sig either: 36
% 		# pairs sig both: 1
% 		# pairs sig other: 243
% 		# pairs not sig: 2913
% 	# pairs exclusively to one or another: 36
% 	# bipolar channels covered by either: 89
% 	# bipolar channels: 91
% === newest as of 04/24/2019 ===
% [*] n bchan cov: 44
% [*] n sig both: 5
% [*] n sig either: 16
% [*] n sig no both: 28
% [*] n sig xor: 11
% === end comment ===
% [*] All Done.