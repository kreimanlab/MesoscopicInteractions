% One interaction, 10 seconds

close all;
clear;

dir_artL = '/media/jerry/KLAB101/h5_notch20/art_nosz';
dir_resL = '/media/jerry/KLAB101/results/coh_w10';
dir_corL = '/media/jerry/internal/data/coreg';
dir_cacheL = './cache';
dir_h5 = '/media/jerry/KLAB101/h5_notch20';

% anatomy folder
subjects_dir = dir_corL;
setenv('SUBJECTS_DIR',subjects_dir);

%metrics = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};
metrics = {'pcBroadband'};

% Patients
SubjectsL = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',... % 'sub23',
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',... % 'sub30',
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48',...
    'mSu'};

% Exclude monkey
SubjectsL = SubjectsL(1:(end-1));

% Plot options
plot_border_dk = true;
plot_border_line_not_pt = true;
plot_remap = false;
n_remap = 30; % number of steps
remap_thresh = 0; % 0 - terminate immediately, 1 - go to midpoint
n_cmap = 100; % number of colormap divisions

% Colors
COLOR_SURF = 0.9*[1 1 1];
COLOR_BORDER = 0.75*[1 1 1];
COLOR_ROI_TXT = 0.1*[1 1 1];
COLOR_ELEC = 0.5*[1 1 1];
COLOR_ELEC_BORDER = 0.1*[1 1 1];
SIZE_ELEC = 30;
FACTOR_ELEC = 0.8; % fraction of border size to make electrode
COH_LINE_WIDTH = 3;
WIDTH_BORDER = 1; % line width of region borders
SIZE_BORDER = 1; % marker size of region borders
SIGMA_NOMAP_MM = 3;

% Plot specific channel and time
r_samp_const = 72056321;
bchan1_const = 35; % 13; % sub1 13
bchan2_const = 84; % 58; % sub1 58
sid = 'sub3';
trig_coh_line_not_pt = true;

% r_samp_const = 22755001;
% bchan1_const = 13; % 13; % sub1 13
% bchan2_const = 58; % 58; % sub1 58
% sid = 'sub1';
% trig_coh_line_not_pt = true;


% ROI
load('rois');
def_rois_short;
              

fig_fmt = '-dpng';
trig_eps = true;
trig_mag_not_ct = true;
%trig_skip_plot = false;
trig_overwrite_cache = true;
trig_covered_rois = false;
system('mkdir figures');
system('mkdir figures/T11d4');

fn_cache = sprintf('cache/fsaverage_sym_flat.mat');
fn_cache_dk = sprintf('cache/fsaverage_sym_flat_dk.mat');
fn_cache_rtxt = sprintf('cache/fsaverage_sym_flat_rtxt.mat');
fn_patch = sprintf('%s/%s/surf/%s.full.flat.patch.3d',subjects_dir,'fsaverage_sym','lh');
fn_pial = sprintf('%s/%s/surf/%s.pial',subjects_dir,'fsaverage_sym','lh');


for iM = 1:length(metrics)
    metric = metrics{iM};
    Ca_fn = sprintf('%s/xsub_out_%s_%i.mat',dir_cacheL,sid,iM);
    Ca = load(Ca_fn);
    
    % load flat surface
    T = load(fn_cache,'pa','tri_flat','pa','v_pial','tri_pial');
    pa = T.pa;
    tri_flat = T.tri_flat;
    v_pial = T.v_pial;
    tri_pial = T.tri_pial;
    
    % check if roi boundary exists
    if (exist(fn_cache_dk,'file')) %  && false
        fprintf('[*] Found cached Desikan-Killiany borders.\n')
        Ca_dk = load(fn_cache_dk);
        bound_dk = Ca_dk.bound_dk;
        bound_dk_lines = Ca_dk.bound_dk_lines;
        Annot = Ca_dk.Annot;
    else
        % read desikan-killiany labels
        fprintf('[*] Computing Desikan-Killiany borders..\n')
        annot_dir = sprintf('%s/%s/label',subjects_dir,'fsaverage_sym');
        atl_fname = sprintf('%s/%s.aparc.annot',annot_dir,'lh');
        [v,l2,c] = read_annotation(atl_fname);
        roi_id = l2(pa.vno + 1);
        Annot.v = v;
        Annot.l = l2;
        Annot.c = c;
        Annot.roi_id = roi_id;
        
        count_bound = 0;
        count_thresh = 10000;
        bound_loop_rep = 1;
        bound_repeat = 1;
        bound_dk = [];
        bound_dk_lines = [];
        bound_dk_Lines = {};
        for i3 = 1:bound_loop_rep:length(tri_flat)
            vB = tri_flat(i3,:);
            tv = roi_id(vB);
            if ((length(unique(tv)) > 1))
                %if (isEven)
                if (mod(count_bound,bound_repeat) == 0)
                    mid_x = mean(pa.x(vB));
                    mid_y = mean(pa.y(vB));
                    bound_dk = [bound_dk; [mid_x,mid_y]];
                    
                    % lines for border betwen 2 regions
                    if (length(unique(tv)) == 2)
                        % find the single
                        %[a,~] = unique(tv);
                        b = sum(tv == tv');
                        sing = (tv == tv(b==1));
                        nsing = find(~sing);
                        mid_xs1 = mean([pa.x(vB(sing)),pa.x(vB(nsing(1)))]);
                        mid_ys1 = mean([pa.y(vB(sing)),pa.y(vB(nsing(1)))]);
                        mid_xs2 = mean([pa.x(vB(sing)),pa.x(vB(nsing(2)))]);
                        mid_ys2 = mean([pa.y(vB(sing)),pa.y(vB(nsing(2)))]);
                        bound_dk_lines = [bound_dk_lines; [mid_xs1, mid_ys1, mid_x, mid_y]];
                        bound_dk_lines = [bound_dk_lines; [mid_xs2, mid_ys2, mid_x, mid_y]];
                    % lines for border betwen 3 regions
                    elseif (length(unique(tv)) == 3)
                        % compute 3 midpoints
                        mid_x12 = mean([pa.x(vB(1)),pa.x(vB(2))]);
                        mid_y12 = mean([pa.y(vB(1)),pa.y(vB(2))]);
                        mid_x23 = mean([pa.x(vB(2)),pa.x(vB(3))]);
                        mid_y23 = mean([pa.y(vB(2)),pa.y(vB(3))]);
                        mid_x13 = mean([pa.x(vB(1)),pa.x(vB(3))]);
                        mid_y13 = mean([pa.y(vB(1)),pa.y(vB(3))]);
                        bound_dk_lines = [bound_dk_lines; [mid_x12, mid_y12, mid_x, mid_y]];
                        bound_dk_lines = [bound_dk_lines; [mid_x23, mid_y23, mid_x, mid_y]];
                        bound_dk_lines = [bound_dk_lines; [mid_x13, mid_y13, mid_x, mid_y]];
                    end
                    %return
                    %plot(mean(pa.x(vB)),mean(pa.y(vB)),'.','color',col_bound,'MarkerSize',size_bound)
                end
                %isEven = ~ isEven;
                count_bound = count_bound + 1;
                if (count_bound > count_thresh)
                    fprintf(2,'W> Boundary point threshold reached. Stopping boundary draw prematurely.\n');
                    break
                end
            end
        end
%         % try to connect the lines
%         bound_dk_allpts = [bound_dk_lines(:,1:2); bound_dk_lines(:,3:4)];
%         [~,sIdx] = sort(bound_dk_allpts(:,i_dim));
%         bound_dk_allpts_s = bound_dk_allpts(sIdx,:);
%         Lines = {};
%         [n_pts,~] = size(bound_dk_allpts_s);
%         % pop existing points
%         sat = true;
%         last_pt = bound_dk_allpts_s(1,:);
%         count = 2;
%         pt_thresh_mm = 1e-3;
%         while (sat)
%             new_pt = bound_dk_allpts_s(count,:);
%             pt_dist = sqrt(sum(( new_pt - last_pt ).^2,2));
%             if (pt_dist < pt_thresh_mm)
%                 bound_dk_allpts_s(count,:) = [NaN NaN];
%             end
%             count = count + 1;
% 
%             last_pt = new_pt;
%             % break on last 
%             if (count == n_pts)
%                 break
%             end
%         end
%         % pop nans
%         bound_dk_allpts_s = bound_dk_allpts_s(~any(isnan(bound_dk_allpts_s),2),:);
%         for ipl = 1:(n_pts-1)
%             % check for line breaks
%             pt_dist = sqrt(sum((bound_dk_allpts_s(ipl,:) - bound_dk_allpts_s(ipl+1,:)).^2,2));
%             break_cond = 1;
%             if (break_cond)
%             end
%         end
%         line2line_error_mm = 1e-3;
%         [n_bound_dk_lines,~] = size(bound_dk_lines);
%         longl = [bound_dk_lines(1,1:2); bound_dk_lines(1,3:4)];
%         for ibl = 2:n_bound_dk_lines
%             candidates = [bound_dk_lines((ibl):end,1:2); bound_dk_lines((ibl):end,3:4)];
%             endpt = longl(end,:);
%             cdist = sqrt(sum((endpt - candidates).^2,2));
%             if (min(cdist) < line2line_error_mm)
%                 longl
%             else
%             end
%         end
        
%         % connect boundary points
%         [n_bound_dk,~] = size(bound_dk);
%         degree = zeros(n_bound_dk,1);
%         for i3 = 1:n_bound_dk
%             dist_dki = sqrt(sum((bound_dk - bound_dk(i3,:)).^2,2));
%             dist_dki(i3) = NaN;
%             dist_dki_s = sort(dist_dki);
%             [~,i3_pair] = min(dist_dki);
% %             in_dict = any( (bdict == i3_pair) | (bdict == i3) );
% %             if ( ~in_dict )
% %                 bdict = [bdict; i3];
% %                 bdict = [bdict; i3_pair];
% %                 bdict = unique(bdict);
% %                 bound_edge = [bound_dk(i3,:) bound_dk(i3_pair,:)];
% %             end
%             
%         end
        %return
        save(fn_cache_dk,'bound_dk','bound_dk_lines','Annot');
    end
    
    
    
    % plot surface
    %h = figure('visible','off');
    fprintf('[*] Plotting fsaverage_sym flat patch..\n')
    h = figure;
    fig_w = 10.5;
    fig_h = 10.5;
    set(h,'Position',round([0 0 fig_w*100 fig_h*100]))
    set(h, 'PaperUnits', 'Inches')
    set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
    set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
    a = trisurf(tri_flat,pa.x',pa.y',pa.z'); hold all;
    
    % plot border
    fprintf('[*] Plotting Desikan-Killiany borders..\n')
%     stream = RandStream('mlfg6331_64');  % Random number stream
%     options = statset('UseParallel',1,'UseSubstreams',1,...
%         'Streams',stream);
%     k = length(rois)*8;
%     [idx,C,sumd,D] = kmeans(bound_dk,k,'Options',options,'MaxIter',10000,'Replicates',4*12);
%     for ik = 1:k
%         plot(bound_dk(idx==ik,1),bound_dk(idx==ik,2),'-','color',COLOR_BORDER);
%     end
    
    
    
    if (plot_border_dk)
        if (plot_border_line_not_pt)
            fprintf('[*] Plotting DK borders as lines..\n')
            % plot line format [mid_x12, mid_y12, mid_x, mid_y]
            [n_blines,~] = size(bound_dk_lines);
            tic
            for il = 1:n_blines
                x = bound_dk_lines(il,[1 3]);
                y = bound_dk_lines(il,[2 4]);
                %plot(x,y,'-','color',COLOR_BORDER,'LineWidth',WIDTH_BORDER);
                patch(x,y,'-','edgecolor',COLOR_BORDER,'LineWidth',WIDTH_BORDER);
            end
            toc
        else
            fprintf('[*] Plotting DK borders as dots..\n')
            plot(bound_dk(:,1),bound_dk(:,2),'.','color',COLOR_BORDER,'MarkerSize',SIZE_BORDER);
        end
        
        % Show roi names
        roi_table = Annot.c.table(~all(Annot.c.table == 0,2),:);
        [n_rois,~] = size(roi_table);
        if (exist(fn_cache_rtxt,'file'))
            fprintf('[*] Found cached ROI names.\n')
            Ca_rtxt = load(fn_cache_rtxt);
            Roitext = Ca_rtxt.Roitext;
        else
            fprintf('[*] Computing ROI names..\n')
            % compute roi midpoints
            Roitext = {};
            i_roi2 = 1;
            for i_roi = 1:n_rois
                roi_idx = ( Annot.roi_id == roi_table(i_roi,end) );
                x = pa.x(roi_idx);
                y = pa.y(roi_idx);
                z = pa.z(roi_idx);
                if (contains(rois_short{i_roi},{'ICN','PCL','PCU','RAC','FRP','LIN','MOF','UNK'}))
                    switch(rois_short{i_roi})
                        case 'ICN'
                            k = 3;
                        case 'PCL'
                            k = 2;
                        case 'PCU'
                            k = 2;
                        case 'RAC'
                            k = 2;
                        case 'FRP'
                            k = 2;
                        case 'LIN'
                            k = 2;
                        case 'MOF'
                            k = 2;
                        case 'UNK'
                            k = 2;
                    end

                    [idx,C,sumd,D] = kmeans([x', y'],k,'Replicates',12*12);
                    for ikm = 1:k
                        ci = (idx==ikm);
                        Roitext{i_roi2,1} = [mean(x(ci)),mean(y(ci)),mean(z(ci))];
                        Roitext{i_roi2,2} = rois_short{i_roi};
                        i_roi2 = i_roi2 + 1;
    %                     text(mean(x(idx==ikm)),mean(y(idx==ikm)),mean(z(idx==ikm)),sprintf('%s',rois_short{i_roi}),...
    %                         'HorizontalAlignment','center','VerticalAlignment','middle','color',COLOR_ROI_TXT)
                    end
    %                 elseif (contains(rois_short{i_roi},{'PCL'}))
    %                     k = 2;
    %                     [idx,C,sumd,D] = kmeans([x', y'],k,'Replicates',12*12);
    %                     for ikm = 1:k
    %                         text(mean(x(idx==ikm)),mean(y(idx==ikm)),mean(z(idx==ikm)),sprintf('%s',rois_short{i_roi}),...
    %                             'HorizontalAlignment','center','VerticalAlignment','middle','color',COLOR_ROI_TXT)
    %                     end
                else
                    Roitext{i_roi2,1} = [mean(x),mean(y),mean(z)];
                    Roitext{i_roi2,2} = rois_short{i_roi};
                    i_roi2 = i_roi2 + 1;
    %                 text(mean(x),mean(y),mean(z),sprintf('%s',rois_short{i_roi}),...
    %                     'HorizontalAlignment','center','VerticalAlignment','middle','color',COLOR_ROI_TXT)
                end
            end
            save(fn_cache_rtxt,'Roitext');
        end
        
        % show text
        fprintf('[*] Plotting ROI names..\n')
        [n_rois,~] = size(Roitext);
        for i_roi = 1:n_rois
            xyz = Roitext{i_roi,1};
            txt = Roitext{i_roi,2};
            text(xyz(1),xyz(2),xyz(3),sprintf('%s',txt),...
                'HorizontalAlignment','center',...
                'VerticalAlignment','middle',...
                'color',COLOR_ROI_TXT)
        end
        
        % plot scale bar
        x = 100;
        y = -140;
        w = 20;
        scalex = [x,x-w,x-w,x,x];
        scaley = [y,y,y+w,y+w,y];
        COLOR_SCALE_BAR = 0*[1 1 1];
        plot(scalex,scaley,'-','color',COLOR_SCALE_BAR)
        text(x-(w/2),y+(w/2),sprintf('%imm',w),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'color',COLOR_SCALE_BAR)
    end
    
    
    
    %return
%     if (plot_border_dk)
%         for ipb = 1:length(bound_dk)
%             plot(bound_dk(ipb,1),bound_dk(ipb,2),'.','color',COLOR_BORDER);
%         end
%     end
    
    % settings
    set(a,'edgecolor','none');
    tri_color = COLOR_SURF;
    set(a,'FaceVertexCdata',tri_color);
    axis off;
    daspect([1 1 1]);
    view(90,90);
    
    
    % Load coherences
    fn_graph = sprintf('%s/%s_graph-%s.h5',dir_resL,sid,metric);
    fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
    ecog = H5eeg(fn_h5);
    w_coh = double(h5read(fn_graph,'/w'));
    w_idx = ceil(r_samp_const/(w_coh*round(ecog.fs)));
    [n_bchan,~] = size(ecog.bip);
    n_comb = nchoosek(n_bchan,2);
    R = h5read(fn_graph,'/R',[1 w_idx],[n_comb 1]);
    chan1 = Ca.chan1 + 1;
    chan2 = Ca.chan2 + 1;
    
    % Load bipolar electrode coordinates
    l = read_label('fsaverage_sym',sprintf('ielvis_%s',sid));
    [n_lchan,~] = size(l);
    if (n_lchan ~= ecog.n_chan)
        fprintf(2,'E: number of channels in %s.label does not match .h5 file\n',sprintf('ielvis_%s',sidL))
        return
    end
    
    % Map labeled electrode to patch
    dist_remap = [];
    bchan_coords = [];
    for i = 1:n_bchan
        % Find out where electrode belongs
        b1c1 = ecog.bip(i,1);
        ind_pa = find(pa.ind == l(b1c1,1),1);
        chan = (b1c1 == l(:,5));
        roi_belong = Ca.atl_labels{chan};
        roi_idx = ( Annot.roi_id == roi_table(strcmp(rois,roi_belong),end) );
        % calculate midpoint of roi where electrode belongs
        x = pa.x(roi_idx);
        y = pa.y(roi_idx);
        xy_belong = [mean(x) mean(y)];
            
        % Make sure it's in the right roi
        if (isempty(ind_pa))
            fprintf('[%s] Electrode %i of %i did not map to patch\n',sid,i,n_bchan);
            x_plt = xy_belong(1) + normrnd(0,SIGMA_NOMAP_MM);
            y_plt = xy_belong(2) + normrnd(0,SIGMA_NOMAP_MM);
            fprintf('\tPut it near region center: %.3f, %.3f mm\n',x_plt,y_plt);
        else
            % calculate position and roi of electrode current
            roi_actual = rois{roi_table(:,end) == Annot.roi_id(ind_pa)};
            xy_actual = [pa.x(ind_pa) pa.y(ind_pa)];
            
            % remap electrode
            remap_coef = linspace(0,1,n_remap);
            if (plot_remap)
                plot(xy_actual(1),xy_actual(2),'redx'); % show remap
            end
            for ii = 1:n_remap
                x_tmp = xy_actual(1)*(1 - remap_coef(ii)) + xy_belong(1)*(remap_coef(ii));
                y_tmp = xy_actual(2)*(1 - remap_coef(ii)) + xy_belong(2)*(remap_coef(ii));
                % find nearest vertex
                [~,ind_pa2] = min(sqrt((pa.x - x_tmp).^2 + (pa.y - y_tmp).^2));
                roi_actual = rois{roi_table(:,end) == Annot.roi_id(ind_pa2)};
                if (strcmp(roi_actual,roi_belong) && (rand >= remap_thresh) )
                    break
                end
            end
            x_plt = x_tmp; %pa.x(ind_pa);
            y_plt = y_tmp; %pa.y(ind_pa);
        end
        
        if (plot_remap)
            plot([xy_actual(1) x_tmp],[xy_actual(2) y_tmp],'red--'); % show remap
        end
        plot(x_plt,y_plt,'.','MarkerSize',SIZE_ELEC,'color',COLOR_ELEC_BORDER); 
        plot(x_plt,y_plt,'.','MarkerSize',SIZE_ELEC*0.9,'color',COLOR_ELEC);
        dist_remap = [dist_remap; sqrt((xy_actual(1) - x_tmp)^2 + (xy_actual(2) - y_tmp)^2)];
        bchan_coords = [bchan_coords; [x_plt y_plt i b1c1]];
    end
    
    % Select only significant coherences
    sig_i = (Ca.Dmats > Ca.dist_thresh) & (Ca.mag > Ca.coh_thresh);
    mag = double(R(sig_i));
    chan1 = double(chan1(sig_i));
    chan2 = double(chan2(sig_i));
    cc = corrcmap(n_cmap);
    cmax = max(mag);
    cmin = min(mag);
    
    % Plot coherence
    for j = 1:length(chan1)
        b1 = chan1(j);
        b2 = chan2(j);
        coh = mag(j);
        if (((b1 == bchan1_const) && (b2 == bchan2_const)) || ((b2 == bchan1_const) && (b1 == bchan2_const)))
            if (trig_coh_line_not_pt)
                x_plt = [bchan_coords(b1,1) bchan_coords(b2,1)];
                y_plt = [bchan_coords(b1,2) bchan_coords(b2,2)];
                %col_i = round(((coh - cmin) / (cmax - cmin)) * n_cmap);
                col_i = round(((coh - cmin) / (cmax - cmin)) * (n_cmap - 1) + 1);
                patch(x_plt,y_plt,'-','EdgeColor',cc(col_i,:),'LineWidth',COH_LINE_WIDTH)
            else
            end
        end
    end
    colormap(cc);
    cb = colorbar;
    set(cb,'TickLength',0);
    caxis([cmin cmax]);
    set(cb,'Ticks',[cmin cmax]);
    
    system('mkdir figures/T011d2');
    oname = sprintf('figures/T011d2/fsaverage_sym_remap_avgMM-%i_medianMM-%i',...
        round(mean(dist_remap)),round(median(dist_remap)));
    print(h,oname,'-dpng','-r300');
    print(h,oname,'-depsc');
    
    close(h);
end
