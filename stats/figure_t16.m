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

    
for imm = 1:length(imetrics)
    iM = imetrics(imm);
    
    for trig_show_path = 0 %[0 1]
    
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

        % convert to integer subject number
        sid_int = find(strcmp(CaA.Subjects,sid));
        Ca = load(sprintf('%s/xsub_out_%s_%i',dir_cache,sid,iM));

        % Find locations bipolar electrodes map onto
        l1 = NaN;
        l2 = NaN;
        n_Es = length(CaR.Es);
        for i = 1:n_Es
            ce = CaR.Es{i};
            [n_ce,~] = size(ce);
            for j = 1:n_ce
                subj = ce(j,1);
                subj_b = ce(j,2);
                if ((subj == sid_i) && (subj_b == b1))
                    l1 = i;
                elseif ((subj == sid_i) && (subj_b == b2))
                    l2 = i;
                end
            end
        end

        % Find all interactions between location pair
        e1 = CaR.Es{l1};
        [n_e1,~] = size(e1);
        [~,sIdx] = sort(e1(:,1));
        e1 = e1(sIdx,:);

        e2 = CaR.Es{l2};
        [n_e2,~] = size(e2);
        [~,sIdx] = sort(e2(:,1));
        e2 = e2(sIdx,:);

        E_reduce = [];
        E_reduce_sid = {};
        for i = 1:n_e1
            for j = 1:n_e2
                % is there a significant interaction between i and j?
                s1 = e1(i,1);
                b1 = e1(i,2);
                sid1 = sprintf('%i',s1);
                while (length(sid1) < 5)
                    sid1 = ['0',sid1];
                end
                sid1 = ['m',sid1];
                %Ca1 = load(sprintf('%s/xsub_out_%s_%i',dir_cache,sid1,iM));

                s2 = e2(j,1);
                b2 = e2(j,2);
                sid2 = sprintf('%i',s2);
                while (length(sid2) < 5)
                    sid2 = ['0',sid2];
                end
                sid2 = ['m',sid2];
                %Ca2 = load(sprintf('%s/xsub_out_%s_%i',dir_cache,sid2,iM));

                if (s1 == s2)
                    % If same subject
                    Ca1 = load(sprintf('%s/xsub_out_%s_%i',dir_cache,sid1,iM));
                    mag = Ca1.AdjMag(b1,b2);
                    dis = Ca1.Dmat(b1,b2);

                    % remember hemisphere
                    hemi_isr = 0;
                    if (strcmp(Ca1.hemi_list{1},'R'))
                        hemi_isr = 1;
                    end

                    % check for significance
                    %if (((~isnan(mag)) && (mag > 0)) && (dis > Ca1.dist_thresh))
                    if (((~isnan(mag))) && (dis > Ca1.dist_thresh))
                        fprintf('[!] %s: bchan %i, %i\n',sid1,b1,b2);
                        E_reduce = [E_reduce; [s1,b1,b2,mag,e1(i,3:5),hemi_isr]];
                        E_reduce_sid = [E_reduce_sid; {sid1}];
                    end
                end
                %return
            end
        end

        % Plot
        [usub,uIdx] = unique(E_reduce(:,1));
        usub_str = E_reduce_sid(uIdx);
        usub_isr = E_reduce(uIdx,end);
        n_usub = length(usub);

        % colorscale calculations
        mag_all = [E_reduce(:,4); CaR.A(l1,l2)];
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


        h = figure('visible','on','Units','pixels','Position',[0 0 1920 1080]);
        %h = figure; set(h,'PaperUnits','inches'); set(h,'PaperPosition',[0 0 8 4]);
        [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,'fsaverage_sym','r','pial'));

        % Show average plot
        %const_elec_surface = 0.7*[1 1 1];
        const_col_surface = 0.85*[1 1 1];
        const_alpha_surface = 0.4;
        const_elec_surface = 0*[1 1 1];
        const_elec_low = 0.75*[1 1 1]; % 0.4

        %[ha, pos] = tight_subplot(1,n_usub,[.01 .01],[.01 .01],[.01 .01]);
        s_d1 = 2;
        s_d2 = 3;
        [ha, pos] = tight_subplot(s_d1,s_d2,[.01 .01],[.01 .01],[.01 .01]);

        % Plot the location-location interaction on the average brain
        axes(ha(1));
        p = trisurf(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
            'EdgeColor','none','FaceColor',const_col_surface,'FaceAlpha',const_alpha_surface);
        hold all;
        brainlight;
        view(90,0)
        axis off

        % plot locatons
        coord = CaR.E(:,3:5);
        [n_chan,~] = size(coord);
        [Xe,Ye,Ze] = sphere(20);
        elec_pa = cell(1,n_chan);
        const_elec_radius = 2;
        for j = 1:n_chan
            X = Xe*const_elec_radius + coord(j,1);
            Y = Ye*const_elec_radius + coord(j,2);
            Z = Ze*const_elec_radius + coord(j,3);
            elec_pa{j} = surf2patch(X,Y,Z);
        end
        for i = 1:n_chan
            col_e = const_elec_low;
            if ((i == l1) || (i == l2))
                col_e = const_elec_surface;
            end
            vert = elec_pa{i}.vertices;
            q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
                'EdgeColor','none','FaceColor',col_e);
            q.SpecularStrength = 0;
            q.SpecularExponent = 1;
        end





        % show line
        tube_n_edges = 3;
        %tube_radius = 0.4;
        col_path = mag2col(CaR.A(l1,l2));%[0 0 1];
        fprintf('Average coherence: %.6f\n',CaR.A(l1,l2));

        pIdx = PathsA(l1,l2);
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



        % Plot the individual subjects
        n_sig_lines = 0;
        n_sig_lines_sid = {};
        for i = (1+1):(n_usub+1)
            sid = usub_str{i-1};
            E_reduce_s = E_reduce(usub(i-1)==E_reduce(:,1),:);

            if (i > (s_d1*s_d2))
                break;
            end

            % subaxis
            axes(ha(i));

            if (usub_isr(i-1) == 1)
                hemi = 'r';
                view_ang = 90;
            elseif (usub_isr(i-1) == 0)
                hemi = 'l';
                view_ang = -90;
            end

            % read subject surface 
            [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,'pial'));

            p = trisurf(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
                    'EdgeColor','none','FaceColor',const_col_surface,'FaceAlpha',const_alpha_surface);
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
            for j = 1:n_chan
                X = Xe*const_elec_radius + s_vert(l(j,1)+1,1);
                Y = Ye*const_elec_radius + s_vert(l(j,1)+1,2);
                Z = Ze*const_elec_radius + s_vert(l(j,1)+1,3);
                elec_pa{j} = surf2patch(X,Y,Z);
            end
            for k = 1:n_chan
                vert = elec_pa{k}.vertices;
                col_e = const_elec_low;

                b1 = find(k == CaP.ecog.bip(:,1));
                if ((~isempty(b1)) && any(b1 == bchans))
                    col_e = const_elec_surface;
                end
                q = trisurf(elec_pa{k}.faces,vert(:,1),vert(:,2),vert(:,3),...
                    'EdgeColor','none','FaceColor',col_e);
                q.SpecularStrength = 0;
                q.SpecularExponent = 1;
            end

            % show patient number
%             sid_int = find(strcmp(CaA.Subjects,sid));
%             text(min(s_vert(:,1)),min(s_vert(:,2)),min(s_vert(:,3)),...
%                 sprintf('%i',sid_int),'FontSize',10);

        end

    %     axes(ha(i+1));
    %     axis off;


        % Save figure
        print(h,sprintf('figures/T16/figure_t16_%i_curve-%i',iM,trig_show_path),'-depsc','-r400');

        close(h);

        [n_Er,~] = size(E_reduce);
        fprintf('Number of pairs: %i\n',n_Er);
        fprintf('Number of subjects: %i\n',length(unique(E_reduce_sid)));
        fprintf('Consistency across pairs: %.4f\n',n_sig_lines/n_Er);
        fprintf('Consistency across subjects: %.4f\n',length(unique(n_sig_lines_sid))/length(unique(E_reduce_sid)));


    end
    
end














return

for imm = 1:length(imetrics)
    iM = imetrics(imm);
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
    const_col_surface = 0.85*[1 1 1];
    const_alpha_surface = 0.4;
    const_elec_surface = 0*[1 1 1];
    const_elec_low = 0.4*[1 1 1];


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
        q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
            'EdgeColor','none','FaceColor',const_elec_surface);
        q.SpecularStrength = 0;
        q.SpecularExponent = 1;
    end
    % --- Plot all interactions to both bchan ------------------------
    b1 = 35;
    b2 = 84;
    tube_n_edges = 3;
    %tube_radius = 1.5; %0.4;
    col_path = [0 0 1];
    count = 1;
    mag_range = [];
    for i = 1:(ecog.n_bchan - 1)
        for j = (i+1):ecog.n_bchan
            cond_fwd = true; %((i == b1) || (j == b2));
            cond_rev = true; %((i == b2) || (j == b1));
            if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
                path = Paths{count};
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
    for i = 1:n_chan
        vert = elec_pa{i}.vertices;
        if ((i == b1c1) || (i == b2c1))
            col_elec = const_elec_surface;
        else
            col_elec = const_elec_low;
        end
        q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
            'EdgeColor','none','FaceColor',col_elec);
        q.SpecularStrength = 0;
        q.SpecularExponent = 1;
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
                [n_path,~] = size(path);
                [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
                fv = surf2patch(X,Y,Z);
                %col_path = mag2col(Ca.mag(count));
                if (Ca.mag(count) < mag_max)
                    col_path = mag2col(Ca.mag(count));
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
        if ((i == b1c1) )
            col_elec = const_elec_surface;
        else
            col_elec = const_elec_low;
        end
        q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
            'EdgeColor','none','FaceColor',col_elec);
        q.SpecularStrength = 0;
        q.SpecularExponent = 1;
    end
    % --- Plot all interactions to one bchan ----------------------------------
    b1 = 35;
    b2 = 84;
    tube_n_edges = 9;
    %tube_radius = 1.5;
    col_path = [0 0 1];
    count = 1;
    for i = 1:(ecog.n_bchan - 1)
        for j = (i+1):ecog.n_bchan
            cond_fwd = (i == b1); %((i == b1) && (j == b2));
            cond_rev = (j == b1); %((i == b2) && (j == b1));
            if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
                path = Paths{count};
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
            end
            count = count + 1;
        end
    end
    brainlight;
    view(90,0)
    axis off


    %subplot(2,3,3)
    axes(ha(3));
    hold all;
    p = trisurf(faces,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
        'FaceAlpha',const_alpha_surface,'EdgeColor','none','FaceColor',const_col_surface);
    for i = 1:n_chan
        vert = elec_pa{i}.vertices;
        if ( (i == b2c1))
            col_elec = const_elec_surface;
        else
            col_elec = const_elec_low;
        end
        q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
            'EdgeColor','none','FaceColor',col_elec);
        q.SpecularStrength = 0;
        q.SpecularExponent = 1;
    end
    % --- Plot all interactions to the other one bchan ------------------------
    b1 = 35;
    b2 = 84;
    tube_n_edges = 9;
    %tube_radius = 1.5;
    col_path = [0 0 1];
    count = 1;
    for i = 1:(ecog.n_bchan - 1)
        for j = (i+1):ecog.n_bchan
            cond_fwd = (i == b2); %((i == b1) && (j == b2));
            cond_rev = (j == b2); %((i == b2) && (j == b1));
            if ( (pass_ct(count) && pass_dist(count)) && (cond_fwd || cond_rev) )
                path = Paths{count};
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
            end
            count = count + 1;
        end
    end
    brainlight;
    view(90,0)
    axis off







    



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
        if (((i == b1c1) ) || strcmp(b1_roi,b1_roi_t))
            col_elec = const_elec_surface;
        else
            col_elec = const_elec_low;
        end
        q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
            'EdgeColor','none','FaceColor',col_elec);
        q.SpecularStrength = 0;
        q.SpecularExponent = 1;
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
    %tube_radius = 0.8;
    col_path = [0 0 1];
    count = 1;
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
            end
            count = count + 1;
        end
    end
    brainlight;
    view(90,0)
    axis off

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
        if (( (i == b2c1)) || strcmp(b2_roi,b1_roi_t))
            col_elec = const_elec_surface;
        else
            col_elec = const_elec_low;
        end
        q = trisurf(elec_pa{i}.faces,vert(:,1),vert(:,2),vert(:,3),...
            'EdgeColor','none','FaceColor',col_elec);
        q.SpecularStrength = 0;
        q.SpecularExponent = 1;
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
    %tube_radius = 0.8;
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
            end
            count = count + 1;
        end
    end
    brainlight;
    view(90,0)
    axis off

    fn_fig = sprintf('figures/T16/figure_t16_metric-%i',iM);
    print(h,fn_fig,'-depsc');
    print(h,fn_fig,'-dpng','-r300');
    
    close(h);
end


fprintf('[*] All Done.\n');

% 
% >>figure_t16
% Warning: Directory already exists. 
% > In figure_t16 (line 12) 
% Warning: Directory already exists. 
% > In figure_t16 (line 13) 
% [*] Loading brainexport/red_6all_fsaverage_1.mat ..
% [!] m00003: bchan 41, 14
% [!] m00005: bchan 35, 84
% [!] m00005: bchan 35, 91
% [!] m00006: bchan 13, 28
% [!] m00025: bchan 39, 6
% [!] m00025: bchan 39, 7
% [!] m00032: bchan 49, 55
% [!] m00039: bchan 36, 14
% [!] m00043: bchan 71, 68
% [!] m00043: bchan 71, 74
% [!] m00045: bchan 46, 20
% [!] m00045: bchan 46, 21
% [!] m00045: bchan 47, 20
% [!] m00045: bchan 47, 21
% [!] m00061: bchan 42, 21
% [!] m00071: bchan 77, 20
% [!] m00084: bchan 44, 22
% [!] m00095: bchan 34, 21
% [!] m00096: bchan 21, 79
% [!] m00096: bchan 21, 80
% [!] m00097: bchan 28, 35
% [!] m00097: bchan 56, 35
% Average coherence: 0.223101
% [*] Loading brainexport/m00003_6all_1.mat ..
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00003/label/all_surf_ielvis.label
% [*] Loading brainexport/m00005_6all_1.mat ..
% [3] m00005 coherence: 0.221297
% [3] m00005 coherence: 0.224905
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00005/label/all_surf_ielvis.label
% [*] Loading brainexport/m00006_6all_1.mat ..
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00006/label/all_surf_ielvis.label
% [*] Loading brainexport/m00025_6all_1.mat ..
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00025/label/all_surf_ielvis.label
% [*] Loading brainexport/m00032_6all_1.mat ..
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00032/label/all_surf_ielvis.label
% Number of pairs: 22
% Number of subjects: 14
% Consistency across pairs: 0.0909
% Consistency across subjects: 0.0714
% [*] Loading brainexport/red_6all_fsaverage_5.mat ..
% [!] m00003: bchan 41, 14
% [!] m00005: bchan 35, 84
% [!] m00005: bchan 35, 91
% [!] m00006: bchan 13, 28
% [!] m00025: bchan 39, 6
% [!] m00025: bchan 39, 7
% [!] m00032: bchan 49, 55
% [!] m00039: bchan 36, 14
% [!] m00043: bchan 71, 68
% [!] m00043: bchan 71, 74
% [!] m00045: bchan 46, 20
% [!] m00045: bchan 46, 21
% [!] m00045: bchan 47, 20
% [!] m00045: bchan 47, 21
% [!] m00061: bchan 42, 21
% [!] m00071: bchan 77, 20
% [!] m00084: bchan 44, 22
% [!] m00095: bchan 34, 21
% [!] m00096: bchan 21, 79
% [!] m00096: bchan 21, 80
% [!] m00097: bchan 28, 35
% [!] m00097: bchan 56, 35
% Average coherence: 0.203551
% [*] Loading brainexport/m00003_6all_5.mat ..
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00003/label/all_surf_ielvis.label
% [*] Loading brainexport/m00005_6all_5.mat ..
% [3] m00005 coherence: 0.206240
% [3] m00005 coherence: 0.200862
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00005/label/all_surf_ielvis.label
% [*] Loading brainexport/m00006_6all_5.mat ..
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00006/label/all_surf_ielvis.label
% [*] Loading brainexport/m00025_6all_5.mat ..
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00025/label/all_surf_ielvis.label
% [*] Loading brainexport/m00032_6all_5.mat ..
% ERROR: could not open /media/klab/internal/data/coreg//media/klab/internal/data/coreg/m00032/label/all_surf_ielvis.label
% Number of pairs: 22
% Number of subjects: 14
% Consistency across pairs: 0.0909
% Consistency across subjects: 0.0714

