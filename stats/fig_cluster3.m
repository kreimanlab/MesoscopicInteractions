% fig_cluster3
% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/cluster3’: File exists
% [*] pcBroadband
% Percentage no coverage: 52.76 %
% Percentage no significance: 20.55 %
% Percentage interacting: 26.69 %s
% Fraction significant ignoring nocov: 56.50 %
% cluster clash: 0.000000000000
% [*] pcTheta
% Percentage no coverage: 52.76 %
% Percentage no significance: 32.66 %
% Percentage interacting: 14.59 %
% Fraction significant ignoring nocov: 30.87 %
% [*] pcAlpha
% Percentage no coverage: 52.76 %
% Percentage no significance: 22.61 %
% Percentage interacting: 24.64 %
% Fraction significant ignoring nocov: 52.14 %
% [*] pcBeta
% Percentage no coverage: 52.76 %
% Percentage no significance: 30.92 %
% Percentage interacting: 16.32 %
% Fraction significant ignoring nocov: 34.54 %
% [*] pcGamma
% Percentage no coverage: 52.76 %
% Percentage no significance: 23.11 %
% Percentage interacting: 24.14 %
% Fraction significant ignoring nocov: 51.09 %
% [!] Done.

% == April 4 ==
% fig_cluster3
% mkdir: cannot create directory ‘figures’: File exists
% mkdir: cannot create directory ‘figures/cluster3’: File exists
% [*] pcBroadband
% Percentage no coverage: 35.02 %
% Percentage no significance: 37.74 %
% Percentage interacting: 27.24 %
% Fraction significant ignoring nocov: 41.93 %
% cluster clash: 0.000000000000
% [*] pcTheta
% Percentage no coverage: 35.02 %
% Percentage no significance: 52.80 %
% Percentage interacting: 12.18 %
% Fraction significant ignoring nocov: 18.75 %
% [*] pcAlpha
% Percentage no coverage: 35.02 %
% Percentage no significance: 40.21 %
% Percentage interacting: 24.77 %
% Fraction significant ignoring nocov: 38.12 %
% [*] pcBeta
% Percentage no coverage: 35.02 %
% Percentage no significance: 50.92 %
% Percentage interacting: 14.06 %
% Fraction significant ignoring nocov: 21.63 %
% [*] pcGamma
% Percentage no coverage: 35.02 %
% Percentage no significance: 41.52 %
% Percentage interacting: 23.46 %
% Fraction significant ignoring nocov: 36.10 %
% [!] Done.


close all;
clear;
rng('shuffle');

dir_art = '/media/klab/internal/data/h5_notch20/art';
dir_cor = '/media/klab/internal/data/coreg';
dir_res = '/home/jerry/data/results/coh_w10';
setenv('SUBJECTS_DIR',dir_cor);
dir_cache = './cache';
dir_h5 = '/media/klab/internal/data/h5_notch20';

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

% Read fsaverage_sym
fpos = 1.5*[0 0 9 7];
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_cor,'fsaverage_sym','r','pial'));
col_surf = 0.95*[1 1 1];
col_elec = 0*[1 1 1];
r_elec = (2.3); % radius (mm)
brain_alpha = 0.7;
fsize = 10;

system('mkdir figures');
system('mkdir figures/cluster3');

for iM = 1:length(metrics)
    fprintf('[*] %s\n',metrics{iM});
    
    %fn_cache = [dir_cache,'/xsub_out_all_',num2str(iM)];
    %Ca = load(fn_cache);
    fn_ca2 = sprintf('%s/fig_cluster2_reduce_%i_new.mat',dir_cache,iM);
    Ca2 = load(fn_ca2);
    
    
    
    % Plot adjacency matrix ===============================================
    
    Ca = load(sprintf('%s/xsub_out_all_%i.mat',dir_cache,iM));
    CaT14 = load(sprintf('./cache/figure_t14_%i_150',iM));
    
    % Set adjacency matrix to plot
    Adj_plt = Ca2.A;
    Adj_plt(Ca2.A_nocov == 0) = 0;
    
    % Print stats
    f_nocov = sum(sum(isnan(Adj_plt)))/numel(Adj_plt);
    f_nosig = sum(sum((Adj_plt == 0)))/numel(Adj_plt);
    f_pos = 1 - (f_nocov + f_nosig);
    fprintf('Percentage no coverage: %.2f %% (%i of %i)\n',f_nocov*100,sum(sum(isnan(Adj_plt))),numel(Adj_plt));
    fprintf('Percentage no significance: %.2f %%\n',f_nosig*100);
    fprintf('Percentage interacting: %.2f %%\n',f_pos*100);
    fprintf('Fraction significant ignoring nocov: %.2f %%\n',(f_pos/(f_pos + f_nosig))*100);
    
    % Roi list
    rois_plt = {};
    rois_plt_val = [];
    for i = 1:length(Ca2.Es)
        rois_plt = [rois_plt; {sprintf('%i',i)}];
        rois_plt_val = [rois_plt_val; i];
    end
    
    % ---- clustering --------------------
    if (iM == 1)
        %adjct_dist_cl = adjct_dist(ind_hum2mac,ind_hum2mac);

        % TODO: RERUN xsub_out_stats to enable adjct_dist clustering
%         adjct_dist_cl = adjct_dist;
%         adjct_dist_cl = adjct_dist_cl(~ind_isnan,~ind_isnan);
        n_rois_cl = length(rois_plt);
        Y = zeros(1,nchoosek(n_rois_cl,2));
        yc = 1;
        for m1 = 1:(n_rois_cl-1)
            for m2 = (m1+1):n_rois_cl

                % ========= CLUSTER BY COHERENCE ==========================
                ad = Adj_plt(m2,m1);
                if (isnan(ad))
                    Y(yc) = 0;
                else
                    Y(yc) = ad;
                end
                yc = yc + 1;
            end
        end
        roi_dist = squareform(Y);
        % average   1.588800512865
        % centroid  0.846690439765
        % complete  1.066014886514
        % median    0.455547431746
        % single    5.997111741497
        % ward      0.920860236074
        % weighted  0.412189366310

        Z = linkage(Y,'weighted'); %,'average' 'median'
        cluster_i = optimalleaforder(Z,Y); % ,'transformation','inverse'
        roi_dist = roi_dist(cluster_i,cluster_i);
        clash = nansum(nansum(triu(roi_dist,1) - triu(roi_dist,2)));
        fprintf('cluster clash: %.12f\n',clash)
        save('cache/fig_cluster3_cluster_i','cluster_i');
    end
    % ---- clustering --------------------
    
    % Apply cluster
    rois_plt = rois_plt(cluster_i);
    rois_plt_val = rois_plt_val(cluster_i);
    Adj_plt = Adj_plt(cluster_i,cluster_i);
    
    % Figure
    %fpos = figure_position_inches; %[0 0 6 8];
    h = figure('PaperUnits','inches','PaperPosition',fpos,'Position',fpos*100,'visible','off');
    %hold all;
    %fig_size_scale = 0.75;
    %set(h,'Position',round(fig_size_scale*[0 0 0.95*1080 0.8*1080]))
    
    v = Adj_plt;
    [v_n, v_m] = size(v);
    
    % Colors
    %map = corrcmap(100);
    %map = inferno(100);
%     color_isnan = 0.3*[1 1 1];
%     color_not_sig = 0.2*[1 1 1];
%     color_distt = color_isnan;
%     
    colordef_adj;
    color_nocov = COLOR_ADJ_NOCOV;
    color_distt = color_nocov;
    color_isnan = color_nocov;
    color_not_sig = COLOR_ADJ_NOSIG;
    map = COLOR_ADJ_CMAP;
    
    iii = 3;
    if (iii == 3)
        minv = nanmin(v(v~=0));
        maxv = nanmax(v(v~=0));
    else
        minv = nanmin(v(:));
        maxv = nanmax(v(:));
    end
    ncol = size(map,1);
    s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
    Im = ind2rgb(s,map);
    [n_roi_p,~] = size(Adj_plt);
    for i = 1:n_roi_p
        for j = 1:n_roi_p
            if (isnan(Adj_plt(i,j)))
                Im(i,j,:) = color_isnan;
            end
%             if (Dmat_plt(i,j) <= dist_thresh)
%                 Im(i,j,:) = color_distt;
%             end
            if ((iii == 3) && (Adj_plt(i,j) == 0))
                Im(i,j,:) = color_not_sig;
            end
        end
    end
    imagesc(Im);
    set(gca,'YTickMode','manual')
    set(gca,'XTickMode','manual')
    
    % Choose rois to label
    rois_plt2 = {};
    rois_plt_val2 = [];
    isodd = false;
    for i = 0:5:length(rois_plt)
        if (isodd)
            rois_plt2 = [rois_plt2; {''}];
        else
            rois_plt2 = [rois_plt2; {sprintf('%i',i)}]; %[rois_plt2; {rois_plt{i}}];
        end
        rois_plt_val2 = [rois_plt_val2; i];
        isodd = (~ isodd);
    end
    yticks(rois_plt_val2);
    yticklabels(rois_plt2);
    xticks(rois_plt_val2);
    xticklabels(rois_plt2);
    
    xtickangle(90);
    set(gca,'tickdir','out');
    set(gca,'fontsize',fsize);
    set(gca,'TickLength',[0.006, 0.001])
    axis tight;
    colormap(map);
    if (((~isnan(minv)) && (~isnan(maxv))) && (minv ~= maxv))
        ctick = linspace(minv,maxv,5);
        ctickT = cell(size(ctick));
        for ict = 1:length(ctickT)
            ctickT{ict} = sprintf('%.2f',ctick(ict));
        end
        cb = colorbar('Ytick',ctick,'TickLabels',ctickT);
        set(cb,'TickLength',0);
        caxis([minv maxv]);
    else
        cb = colorbar;
        set(cb,'TickLength',0);
    end
    
    % Save figure
    print(h,sprintf('figures/cluster3/fsaverage_sym_reduced_labels_%i_adj',iM),'-depsc','-r400');
    print(h,sprintf('figures/cluster3/fsaverage_sym_reduced_labels_%i_adj',iM),'-dpng','-r400');
    close(h);
    
    
    
    
    
    % Labels ========================================================
    
    % Figure
    %fpos = 1.1*[0 0 5 7]; %[0 0 6 8];
    h = figure('PaperUnits','inches','PaperPosition',fpos,'Position',fpos*100,'visible','off');
    hold all;
    
    % Offsets - medial view
    c1_coeff = (-1);
    c1_off = (-1) * (min(s_vert(:,1)*(-1)) - min(s_vert(:,1)));
    c2_off = (-1) * (min(s_vert(:,2)*(-1)) - min(s_vert(:,2)));
    c2_coeff = (-1);
    c3_off = (-1) * 130; % 210;
    
    % Offsets - bottom view
    d3_off = (-1) * 250; % 125;
    d1_coeff = (-1);
    
    % threshold for front back separation
    c3_thresh = 22;
    
    % threshold for bottom separation
    c1_thresh = -15;
    
    % Plot electrodes
    [x,y,z] = sphere(8);
    x = x * r_elec;
    y = y * r_elec;
    z = z * r_elec;
    
    % text properties
    valign = 'bottom';
    halign = 'center';
    col_tback = [0.9*[1 1 1], 0.7];
    col_tedge = 0.2*[1 1 1];
    tmargin = 0.1;
    
    % Bank
    b = [];
    bx = 11; % Second coordinate
    bx_step = 0.5;
    bbias = 0.5; % guide line slope bias (0 to 1, 0 is slope = 1)
    zlift = 60;
    
    for i = 1:length(Ca2.Es)
        
        
        % Label ROI numbers according to CaT14
        idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(i))));
        %text_val = sprintf('%i',idx);
        text_val = CaT14.rois_plt{idx};
        % old way:
        %text_val = sprintf('%i',cluster_i(i)); % sprintf('%s',dec2rom(i));
        
        elec_x = Ca2.E(i,3);
        elec_y = Ca2.E(i,4);
        elec_z = Ca2.E(i,5);
        
        % Lateral view
        e1 = elec_x;
        e2 = elec_y;
        e3 = elec_z;
        p1 = x + e1;
        p2 = y + e2;
        p3 = z + e3;
        if ((elec_x > c3_thresh) && (elec_z > c1_thresh))
        %if (true)
%             p = surf(p1,p2,p3,'FaceColor',col_elec,'EdgeColor','none');
%             p.AmbientStrength = 0.5 ;
%             p.DiffuseStrength = 0.2 ;
%             p.SpecularStrength = 0;
%             p.SpecularExponent = 1;

            % check if in bank
            e2s = e2;
            e3s = e3;
            if (isempty(b))
                sat = true;
            else
                sat = min(sqrt(sum(([e2 e3] - b).^2,2))) > bx;
            end
            while (~sat)
                e2 = e2 + bx_step*sign(elec_y)*(1 - bbias*i/length(Ca2.Es));
                e3 = e3 + bx_step*sign(elec_z)*(1 + bbias*i/length(Ca2.Es));
                sat = min(sqrt(sum(([e2 e3] - b).^2,2))) > bx;
            end
            
%             if (e2 - e2s > 0)
%                 halign = 'left';
%             else
%                 halign = 'right';
%             end
%             if (e3 - e3s > 0)
%                 valign = 'bottom';
%             else
%                 valign = 'top';
%             end

            text(e1+zlift,e2,e3,text_val,'FontSize',fsize,'BackgroundColor',col_tback,...
                'HorizontalAlignment',halign,'VerticalAlignment',valign,...
                'Margin',tmargin,'EdgeColor',col_tedge);
            hold on;
            plot3([e1 e1],[e2s e2],[e3s e3],'black.-');
            hold on;
            b = [b; [e2 e3]];
        end
        
        % Medial view
        e1 = elec_x*c1_coeff + c1_off;
        e2 = elec_y*c2_coeff + c2_off;
        e3 = elec_z + c3_off;
        p1 = x + e1;
        p2 = y + e2;
        p3 = z + e3;
        if ((elec_x <= c3_thresh) && (elec_z > c1_thresh))
        %if (true)
%             p = surf(p1,p2,p3,'FaceColor',col_elec,'EdgeColor','none');
%             p.AmbientStrength = 0.5 ;
%             p.DiffuseStrength = 0.2 ;
%             p.SpecularStrength = 0;
%             p.SpecularExponent = 1;
            
            % check if in bank
            bbias = 0.2;
            e2s = e2;
            e3s = e3;
            if (isempty(b))
                sat = true;
            else
                sat = min(sqrt(sum(([e2 e3] - b).^2,2))) > bx;
            end
            while (~sat)
%                 e2 = e2 + bx_step*sign(c2_coeff*elec_y);
%                 e3 = e3 + bx_step*sign(elec_z);
                e2 = e2 + bx_step*sign(elec_y*c2_coeff)*(1 - bbias*i/length(Ca2.Es));
                e3 = e3 + bx_step*sign(elec_z)*(1 + bbias*i/length(Ca2.Es));
                sat = min(sqrt(sum(([e2 e3] - b).^2,2))) > bx;
            end
            
            text(e1+zlift,e2,e3,text_val,'FontSize',fsize,'BackgroundColor',col_tback,...
                'HorizontalAlignment',halign,'VerticalAlignment',valign,...
                'Margin',tmargin,'EdgeColor',col_tedge);
            hold on;
            plot3([e1 e1],[e2s e2],[e3s e3],'black.-');
            hold on;
            b = [b; [e2 e3]];
        end
        
        % Bottom view
        e1 = elec_z*d1_coeff;
        e2 = elec_y;
        e3 = elec_x + d3_off;
        p1 = x + e1;
        p2 = y + e2;
        p3 = z + e3;
        if ((elec_z <= c1_thresh))
        %if (true)
%             p = surf(p1,p2,p3,'FaceColor',col_elec,'EdgeColor','none');
%             p.AmbientStrength = 0.5 ;
%             p.DiffuseStrength = 0.2 ;
%             p.SpecularStrength = 0;
%             p.SpecularExponent = 1;

            % check if in bank
            bbias = 0.02;
            bx = 9.5;
            
            e2s = e2;
            e3s = e3;
            if (isempty(b))
                sat = true;
            else
                sat = min(sqrt(sum(([e2 e3] - b).^2,2))) > bx;
            end
            if ((i == 108) || (i == 124) || (i == 107))
                elec_z = (-1)*elec_z;
            end
            while (~sat)
%                 e2 = e2 + bx_step*sign(elec_y);
%                 e3 = e3 + bx_step*sign(elec_z*d1_coeff); % elec_z*d1_coeff
                e2 = e2 + bx_step*sign(elec_y)*(1 - bbias*i/length(Ca2.Es));
                e3 = e3 + bx_step*sign(elec_z*d1_coeff)*(1 + bbias*i/length(Ca2.Es));
                sat = min(sqrt(sum(([e2 e3] - b).^2,2))) > bx;
            end
            
            text(e1+zlift,e2,e3,text_val,'FontSize',fsize,'BackgroundColor',col_tback,...
                'HorizontalAlignment',halign,'VerticalAlignment',valign,...
                'Margin',tmargin,'EdgeColor',col_tedge);
            hold on;
            plot3([e1 e1],[e2s e2],[e3s e3],'black.-');
            hold on;
            b = [b; [e2 e3]];
        end
    end
    
    CaC = load('cache/fig_cluster2_reduce_new_annot_vertex_color2.mat');
    vcolor = CaC.vertex_color2_bw;
    
    
    % Plot surface
    surf_diffuse = 0.3;
    surf_ambient = 0.5;
    
    p = trisurf(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
        'FaceColor',col_surf,'EdgeColor','none','FaceAlpha',brain_alpha);
    p.AmbientStrength = surf_ambient;
    p.DiffuseStrength = surf_diffuse;
    p.SpecularStrength = 0;
    p.SpecularExponent = 1;
    p.BackFaceLighting = 'lit';
    p.FaceLighting = 'gouraud';
    set(p,'FaceVertexCData',vcolor);
    shading('interp');
    hold on;
    
    % Plot surface, medial
    p = trisurf(faces+1,s_vert(:,1)*c1_coeff + c1_off,...
        s_vert(:,2)*c2_coeff + c2_off,... % align horizontal
        s_vert(:,3) + c3_off,...
        'FaceColor',col_surf,'EdgeColor','none','FaceAlpha',brain_alpha);
    daspect([1 1 1]);
    p.AmbientStrength = surf_ambient;
    p.DiffuseStrength = surf_diffuse;
    p.SpecularStrength = 0;
    p.SpecularExponent = 1;
    p.BackFaceLighting = 'lit';
    p.FaceLighting = 'gouraud';
    set(p,'FaceVertexCData',vcolor);
    shading('interp');
    hold on;
    
    % Plot surface, bottom
    p = trisurf(faces+1,s_vert(:,3) * d1_coeff,...
        s_vert(:,2),... % align horizontal
        s_vert(:,1) + d3_off,...
        'FaceColor',col_surf,'EdgeColor','none','FaceAlpha',brain_alpha);
    daspect([1 1 1]);
    p.AmbientStrength = surf_ambient;
    p.DiffuseStrength = surf_diffuse;
    p.SpecularStrength = 0;
    p.SpecularExponent = 1;
    p.BackFaceLighting = 'lit';
    p.FaceLighting = 'gouraud';
    set(p,'FaceVertexCData',vcolor);
    shading('interp');
    hold on;
    
    % Camera
    cam_elev = 0;
    camlight(-135,cam_elev);
    camlight(45,cam_elev);
    camlight(-225,cam_elev);
    camlight(-45,cam_elev);
    view(90,0);
    axis off;
    
    % Expand space
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];

    % Save figure
    print(h,sprintf('figures/cluster3/fsaverage_sym_reduced_labels_%i',iM),'-depsc','-r400');
    print(h,sprintf('figures/cluster3/fsaverage_sym_reduced_labels_%i',iM),'-dpng','-r400');
    close(h);
    
    
    % n_pairs vs mag
    Npairs = nan(length(Ca2.Es),1);
    Nusubs = nan(length(Ca2.Es),1);
    for i = 1:length(Ca2.Es)
        [n_pairs,~] = size(Ca2.Es{i});
        n_usubs = length(unique(Ca2.Es{i}(:,1)));
        Npairs(i) = n_pairs;
        Nusubs(i) = n_usubs;
    end
    % n_subs vs mag
    
    
end

fprintf('[!] Done.\n');


function ans = dec2rom(z)
    d = [ 1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1];
    c =  {'M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I'};
    [];
    for ii = 1:numel(d)
        if z >= d(ii)    
            ans = [ans,repmat(c{ii},1,fix(z/d(ii)))];
            z = rem(z,d(ii));
        end
    end
end