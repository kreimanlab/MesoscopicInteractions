close all;
clear;

% Colors
BRAIN_MESH_ALPHA = 1;
COLOR_MEDIAL_WALL = 0.1*[1 1 1];
COLOR_NOSIG = 0.5*[1 1 1];
COLOR_BLACK = 0.3*[1 1 1];
SURFACE_TYPE = 'pial';
paper_width = 8.5;
paper_height = 6.5;


% load surf
dir_cacheLp = './cache';
dir_corLp = '/media/klab/internal/data/coreg';
metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};
metrics_suffix = {'0.5-125 Hz','3-8 Hz','8-12 Hz','12-30 Hz','30-100 Hz'};
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_corLp,'fsaverage_sym','r',SURFACE_TYPE));



%i = 1;
%iM = i;
for iM = 1 %[1 5]

    metric = metricsp{iM};
    
    
    % plot 150 boundaries on DK areas

    h = figure('visible','off');
    set(h,'PaperUnits','inches');
    set(h,'PaperPosition',[0 0 paper_width paper_height]);
    axis off;
    hold all;

    %subplot(2,1,1);
    [ha, pos] = tight_subplot(2,2,[0.01 0.01],8*[0.01 0.01],4*[0.01 0.01]);
    axes(ha(1));
    p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
        'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
    hold all;



    %--------------------------------------------------------------------------
    %
    % original subject
    sid = 'm00005';
    sid_i = str2double(sid(2:end));
    b1 = 35;
    b2 = 84;
    %CaR = load(sprintf('%s/fig_cluster2_reduce_%i_new',dir_cache,iM));
    CaA = load(sprintf('%s/xsub_out_all_%i',dir_cacheLp,iM));

    % convert to integer subject number
    sid_int = find(strcmp(CaA.Subjects,sid));
    Ca = load(sprintf('%s/xsub_out_%s_%i',dir_cacheLp,sid,iM));

    % Find ROIs bipolar electrodes map onto
    rois = Ca.C.AtlLabels{2};
    b1c1 = Ca.ecog.bip(b1,:);
    b2c1 = Ca.ecog.bip(b2,:);
    b1_roi = rois{b1c1};
    b2_roi = rois{b2c1};
    %--------------------------------------------------------------------------
    
    
    

    % Load adjacency matrix
    CaT14 = load(sprintf('%s/figure_t14_%i_150.mat',dir_cacheLp,iM));

    % region highlight by number

    % Color DK region
    const_col_surface = 0.9 * ones(1,3);
    [v_a, l, ct] = read_annotation(sprintf('%s/%s/label/%sh.aparc.annot',dir_corLp,'fsaverage_sym','r'));
    [n_vert,~] = size(s_vert);
    const_col_surface2 = repmat(const_col_surface,[n_vert 1]);
    vtx_r1 = false(length(l),1);
    vtx_r2 = false(length(l),1);
    for ihi = 1:ct.numEntries
        roi_id = ct.table(ihi,end);
        cond_r1 = strcmp(ct.struct_names{ihi},b1_roi);
        cond_r2 = strcmp(ct.struct_names{ihi},b2_roi);
        if (cond_r1)
            vtx_r1(l==roi_id) = true;
        end
        if (cond_r2)
            vtx_r2(l==roi_id) = true;
        end
        cond_to_color = cond_r1 | cond_r2;
        if (cond_to_color)
            roi_col = ct.table(ihi,1:3) / 255;
            n_vtx = sum(l==roi_id);
            const_col_surface2(l==roi_id,:) = repmat(roi_col,[n_vtx 1]);
        end
    end
    col_r1 = ct.table(strcmp(ct.struct_names,b1_roi),1:3) / 255;
    col_r2 = ct.table(strcmp(ct.struct_names,b2_roi),1:3) / 255;
    vcolor = const_col_surface2;
    %vcolor = 0.8*ones(length(s_vert),3);
    %roi_cl = Res{k,1};

    % deprecated
    cluster_i = load('cache/fig_cluster3_cluster_i.mat');
    cluster_i = cluster_i.cluster_i;
    % end deprecated
%     CaT14 = load(sprintf('./cache/figure_t14_%i_150',iM));
%     cluster_i2 = zeros(size(cluster_i));
%     for i = 1:length(cluster_i)
%         idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(i))));
%         cluster_i2(i) = idx;
%     end
%     cluster_i = cluster_i2; 

    
    CaC2 = load('cache/fig_cluster2_reduce_annot');

    % 150 region - superior temporal overlap
    tint1_factor = 0.1;
    tint1_col = [0 0.4 0];
    tint1_border_col = [0 0.2 0];
    % 150 region - pars op overlap
    tint2_factor = 0.1;
    tint2_col = [0.4 0 0];
    tint2_border_col = [0.2 0 0];

    % First loop to find overlap areas
    %   label each vertex in s_vert with one of 150 areas
    vtx_region = nan(length(s_vert),1);
    for j = 1:length(s_vert)
        r = CaC2.region(j);
        idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(r))));
        vtx_region(j) = idx;
    end
    fprintf('[*] First loop complete.\n')

    r1_overlap_150 = zeros(length(CaT14.rois_plt),1);
    for r1 = 1:length(r1_overlap_150)
        %f_overlap = mean(vtx_r1(vtx_r1) & (vtx_region(vtx_r1) == r1));
        f_overlap = mean(vtx_r1(vtx_region == r1));
        r1_overlap_150(r1) = f_overlap;
    end
    %plot(r1_overlap_150,'black.')

    r2_overlap_150 = zeros(length(CaT14.rois_plt),1);
    for r2 = 1:length(r2_overlap_150)
        f_overlap = mean(vtx_r2(vtx_region == r2));
        r2_overlap_150(r2) = f_overlap;
    end



    thresh_overlap = 0.3; % Default 0.3
    region_const1 = find(r1_overlap_150 > thresh_overlap);
    region_const2 = find(r2_overlap_150 > thresh_overlap);
    coord1 = cell(length(region_const1),1);
    coord2 = cell(length(region_const2),1);
    
    
    % region_const1 = [ 40, 138, 64, 1, 84]; % <- region to highlight
    % region_const2 = [113, 112];

    % Second loop to calculate colors
    for j = 1:length(s_vert)

        % Color 150 region
        idx = vtx_region(j);
        %r = CaC2.region(j);
        %idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(r))));
        cond_region1 = any(idx == region_const1);
        cond_region2 = any(idx == region_const2);
        cond_tocolor = (cond_region1) || (cond_region2);

        if (cond_tocolor)

            % Check if vertex is a border
            face = faces((any((faces+1) == j,2)),:) + 1;
            vtx_span = unique(face(:));
            reg150_span = CaC2.region(vtx_span);
            cond_border = length(unique(reg150_span)) > 1;


            if (cond_border)
                t1 = 1;
                t2 = 1;
                c1 = tint1_border_col;
                c2 = tint2_border_col;
            else
                t1 = tint1_factor;
                t2 = tint2_factor;
                c1 = tint1_col;
                c2 = tint2_col;
            end
            % Region 1
            if (cond_region1)
                vcolor(j,:) = (1 - t1) * vcolor(j,:) + t1 * c1;
                coord1{(idx == region_const1)} = [coord1{(idx == region_const1)}; s_vert(j,:)];
            end
            % Region 2
            if (cond_region2)
                vcolor(j,:) = (1 - t2) * vcolor(j,:) + t2 * c2;
                coord2{(idx == region_const2)} = [coord2{(idx == region_const2)}; s_vert(j,:)];
            end

        end


    %     if ((trig_blot_nosig) && (cl_adj_nosig(idx)))
    %         col_kmeans = COLOR_NOSIG;
    %     else
    %         col_kmeans = cc_k(roi_cl(idx),:);
    %     end
    % 
    %     vcolor(j,:) = col_kmeans;
    %     if (CaC1.region_medial(j) == 1)
    %         vcolor(j,:) = COLOR_MEDIAL_WALL;
    %     end
    end

    
    fprintf('[*] Second loop complete.\n')


    %%

    set(p,'FaceVertexCData',vcolor);
    brainlight;
    view(90,0);
    shading(gca,'interp');
    axis off;
    
    
    
    % plot original electrodes
    l = [];
    if (isempty(l))
        l = read_label(sprintf('%s','fsaverage_sym'),sprintf('ielvis_%s',sid));
        if (isempty(l))
            fprintf('[*] Reading from fallback location..\n');
            l = read_label(sprintf('%s','/media/jerry/internal/data/coreg/fsaverage_sym'),sprintf('ielvis_%s',sid));
        end
    end
    e1 = Ca.ecog.bip(b1,1);
    coord1 = l(l(:,end)==e1,2:4);
    e2 = Ca.ecog.bip(b2,1);
    coord2 = l(l(:,end)==e2,2:4);

    [x,y,z] = sphere(64);
    x = x * 1.2;
    y = y * 1.2;
    z = z * 1.2;
    const_elec_surface2 = 0.000001 * [1 1 1];
%     q1 = surf(x+coord1(1),y+coord1(2),z+coord1(3),'FaceColor',const_elec_surface2,'EdgeColor','none');
%     q1.SpecularStrength = 0;
%     q2 = surf(x+coord2(1),y+coord2(2),z+coord2(3),'FaceColor',const_elec_surface2,'EdgeColor','none');
%     q2.SpecularStrength = 0;
    plot3(coord1(1)+20,coord1(2),coord1(3),'.','Color',const_elec_surface2,'MarkerSize',14);
    plot3(coord2(1)+20,coord2(2),coord2(3),'.','Color',const_elec_surface2,'MarkerSize',14);
    
    
%     e1g = Ca.ecog.bip(b1,2);
%     coord1g = l(l(:,end)==e1g,2:4);
%     e2g = Ca.ecog.bip(b2,2);
%     coord2g = l(l(:,end)==e2g,2:4);
%     plot3(coord1g(1)+0,coord1g(2),coord1g(3),'.','Color',const_elec_surface2,'MarkerSize',14);
%     plot3(coord2g(1)+0,coord2g(2),coord2g(3),'.','Color',const_elec_surface2,'MarkerSize',14);
%     
    
    %return
    

    % 150-area labels
    % r_offset = 30;
    % for i = 1:length(coord1)
    %     tx = sprintf('%i',region_const1(i));
    %     co = mean(coord1{i},1);
    %     text(co(1)+r_offset,co(2),co(3),tx);
    % end
    % for i = 1:length(coord2)
    %     tx = sprintf('%i',region_const2(i));
    %     co = mean(coord2{i},1);
    %     text(co(1)+r_offset,co(2),co(3),tx);
    % end



    legendCol(1,1:3) = col_r1;
    legendCol(2,1:3) = col_r2;
    legendStr = {'Superior Temporal','Pars Opercularis'};
    for il = 1:2
        box_width = 0.016;
        box_vert_spacer = 6*0.16 * box_width;
        horiz_offset = 0;
        horiz_text_offset = -1 * (box_width) - 0.001;
        horiz_pos = pos{1}(1) + pos{1}(3) + horiz_offset;
        vert_offset =  - il * (box_width + box_vert_spacer);
        vert_pos = pos{1}(2) + pos{1}(4) + vert_offset;
        annotation('rectangle',[horiz_pos + horiz_text_offset, vert_pos, box_width, 0.9*box_width*(paper_width/paper_height)],...
            'FaceColor',legendCol(il,:),'EdgeColor','black','LineWidth',1);
        annotation('textbox',[horiz_pos,vert_pos-0.002, 8*box_width, 2*box_width*(paper_width/paper_height)],...
            'String',legendStr{il},'EdgeColor','none',...
            'FontSize',6,'FontName','Helvetica','VerticalAlignment','middle')
    end


    axes(ha(2));
    p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
        'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
    set(p,'FaceVertexCData',vcolor);
    brainlight;
    view(-90,0);
    shading(gca,'interp');
    axis off;


    axes(ha(3));
    p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
        'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
    set(p,'FaceVertexCData',vcolor);
    brainlight;
    view(90,-90);
    shading(gca,'interp');
    axis off;

    axes(ha(4));
    axis off;

    mkdir('./figures/T18');
    print(h,sprintf('figures/T18/Overlap_%i',iM),'-depsc','-r400');
    print(h,sprintf('figures/T18/Overlap_%i',iM),'-dpng','-r400');
    close(h);



    %%

    subregion = [region_const1; region_const2];
    sub_rois = cell(length(subregion),1);
    for i = 1:length(sub_rois)
        sub_rois{i} = sprintf('%i',subregion(i));
    end
    Adj_sub = CaT14.Adj_plt2(subregion,subregion);


    % plot
    h = figure('visible','off');
    %set(h,'Position',round(fig_size_scale*[0 0 1*1080 1*1080]));
    set(h,'PaperUnits','Inches');
    set(h,'PaperPosition',[0 0 8.5 0.6*6.8]);

    rois_plt = sub_rois;
    fontsz = 10;
    colordef_adj;
    color_nocov = COLOR_ADJ_NOCOV;
    color_distt = color_nocov;
    color_isnan = color_nocov;
    color_not_sig = COLOR_ADJ_NOSIG;
    v = Adj_sub;
    [v_n, v_m] = size(v);
    %map = corrcmap(100);
    map = COLOR_ADJ_CMAP; %inferno(100);
    minv = min(v(v~=0));
    maxv = max(v(v~=0));
    ncol = size(map,1);
    s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
    Im = ind2rgb(s,map);
    n_roi_p = length(v);
    for i = 1:n_roi_p
        for j = 1:n_roi_p
            if (isnan(v(i,j)))
                Im(i,j,:) = color_isnan;
            end
    %         if (Dmat_plt(i,j) <= dist_thresh)
    %             Im(i,j,:) = color_distt;
    %         end
            if ( (v(i,j) == 0) )
                Im(i,j,:) = color_not_sig;
            end
        end
    end
    imagesc(Im);
    set(gca,'YTickMode','manual')
    set(gca,'XTickMode','manual')

    stride_axislabel = 1;
    sub_idx = stride_axislabel:stride_axislabel:length(rois_plt);
    rty = 1:length(rois_plt);
    yticks(rty(sub_idx));
    yticklabels(rois_plt(sub_idx));
    rtx = 1:length(rois_plt);
    xticks(rtx(sub_idx));
    xticklabels(rois_plt(sub_idx));
    xtickangle(90);

    %             yticks(1:length(rois_plt));
    %             yticklabels(rois_plt);
    %             xticks(1:length(rois_plt));
    %             xticklabels(rois_plt);
    %             xtickangle(90);

    set(gca,'tickdir','out');
    set(gca,'fontsize',fontsz); %fontsz*(31/n_rois_cl)^(0.7)
    set(gca,'TickLength',[0.001, 0.001])
    daspect([1 1 1]);
    set(gca,'FontName','Arial');

    % Move axis if needed
    ax = gca;
    tp = ax.Position;
    tp(1) = tp(1) - 0.05;
    ax.Position = tp;

    %-----------------------------------------
    colormap(map);
    if (((~isnan(minv)) && (~isnan(maxv))) && (minv ~= maxv))
        cytick = linspace(minv,maxv,3);
        cytick_str = cell(1,length(cytick));
        for kk = 1:length(cytick)
            cytick_str{kk} = sprintf('%.2f',cytick(kk));
        end
        cb = colorbar('ytick',cytick,'yticklabel',cytick_str,'FontSize',fontsz);
        caxis([minv maxv]);
    else
        cb = colorbar('ytick',minv);
    end
    set(cb,'TickDir','out');
    set(cb,'TickLength',0);

    %--- display frequency band range ----
    iii2 = 1;
    metricTxt = metric(3:end);
    if (strcmp(metricTxt,'Broadband'))
        if (iii2 == 1)
            ylabel(cb,sprintf('Coherence'));
        elseif (iii2 == 2)
            ylabel(cb,sprintf('Variance of Coherence'));
        end
    else
        if (iii2 == 1)
            ylabel(cb,sprintf('%s Coherence (%s)',metric(3:end),metrics_suffix{iM}));
        elseif (iii2 == 2)
            ylabel(cb,sprintf('Variance of %s Coherence (%s)',metric(3:end),metrics_suffix{iM}));
        end
    end
    %-------------------------------------

    colormap(map);
    cb.Location = 'eastoutside';
    cbpos = cb.Position;
    cbpos(1) = 0.82;
    cbpos(2) = cbpos(2) + 0.2; %2.6*cbpos(2);
    cbpos(3) = 0.02;
    cbpos(4) = 0.4; %0.6*cbpos(4);
    cb.Position = cbpos;


    annotation('rectangle',[cb.Position(1) cb.Position(2)-0.06 cb.Position(3) 0.02],'FaceColor',color_not_sig);
    annotation('textbox',[cb.Position(1)+0.02 cb.Position(2)-0.06 0.2 0.02],'String','Not Significant','FitBoxToText','on','EdgeColor','none','VerticalAlignment','middle');
    annotation('rectangle',[cb.Position(1) cb.Position(2)-0.09 cb.Position(3) 0.02],'FaceColor',color_nocov);
    annotation('textbox',[cb.Position(1)+0.02 cb.Position(2)-0.09 0.2 0.02],'String','No Coverage','FitBoxToText','on','EdgeColor','none','VerticalAlignment','middle');
    %-------------------------------------------------
    
    
    % Show details
    % cross-area
    A_bn = Adj_sub((end-2):end,1:(end-3));
    fprintf('[!] between ST, PO: %i of %i cov\n',sum(~isnan(A_bn(:))),numel(A_bn));
    A_bns = A_bn(:);
    A_bns = A_bns(~isnan(A_bns));
    fprintf('[!] between ST, PO: %i sig of %i cov\n',sum(A_bns~=0),numel(A_bns));
    A_s = A_bns(A_bns~=0);
    fprintf('\tCoherence mean: %.2f, std: %.2f, min: %.2f, max: %.2f\n',mean(A_s),std(A_s),min(A_s),max(A_s))
    
    % Superior temporal within-area
    A_st = Adj_sub(1:(end-3),1:(end-3));
    A_l = A_st;
    al = [];
    for l1 = 1:(length(A_l)-1)
        for l2 = (l1+1):length(A_l)
            al = [al; A_l(l1,l2)];
        end
    end
    
    A_po = Adj_sub((end-2):end,(end-2):end);
    
    
    
    print(h,sprintf('figures/T18/Overlap_%i_adj',iM),'-depsc','-r400');
    print(h,sprintf('figures/T18/Overlap_%i_adj',iM),'-dpng','-r400');
    close(h);
    
    save(sprintf('./cache/figure_t18_%i',iM));

end



% 23 janvier
% [*] First loop complete.
% [*] Second loop complete.
% ERROR: could not open /fsaverage_sym/label/ielvis_m00005.label
% [*] Reading from fallback location..
% Warning: Directory already exists. 
% > In figure_t18 (line 322) 
% [!] between ST, PO: 25 of 30 cov
% [!] between ST, PO: 18 sig of 25 cov
% 	Coherence mean: 0.30, std: 0.10, min: 0.19, max: 0.49
% [*] First loop complete.
% [*] Second loop complete.
% ERROR: could not open /fsaverage_sym/label/ielvis_m00005.label
% [*] Reading from fallback location..
% Warning: Directory already exists. 
% > In figure_t18 (line 322) 
% [!] between ST, PO: 25 of 30 cov
% [!] between ST, PO: 19 sig of 25 cov
% 	Coherence mean: 0.29, std: 0.08, min: 0.16, max: 0.46
%     



% 10 janvier 2020
% figure_t18
% [*] First loop complete.
% [*] Second loop complete.
% Warning: Directory already exists. 
% > In figure_t18 (line 304) 
% [!] between ST, PO: 25 of 30 cov
% [!] between ST, PO: 18 sig of 25 cov
% 	Coherence mean: 0.30, std: 0.10, min: 0.19, max: 0.49