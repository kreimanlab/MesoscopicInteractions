close all;
clear;

metric = 's';
n_sub = 51; % 51
atl_i = 2;

% constants
min_n_sub = 2;
min_n_epair_is_interact = 2; % within each patient

% read fsaverage_sym lh flattened
% [pa] = read_patch_rev('coreg/lh.full.flat.patch.3d');
% k = convhull(pa.x,pa.y);
% [v_pial,tri_pial] = read_surf('coreg/lh.pial');
% 
% plot(pa.x,pa.y,'black.');

load('coreg/fsaverage_sym_flat');


%tri = delaunay(pa.x,pa.y);
%trisurf(tri,pa.x,pa.y,pa.z);
%plot3(pa.x,pa.y,pa.z,'black.','markersize',0.1)
%view(0,90);

%n_sub = 20;

for i = 1:n_sub
    xfname = sprintf('./xsub-%s-%i',metric,i);
    fprintf('[%i/%i] Loading %s.\n',i,n_sub,xfname)
    load(xfname);
    fprintf('\tBuilding xsub matrices.\n')
    
    
    if (i == 1)
        adjct_xsub = cell(size(AT.AdjCT{atl_i}));
        adjct_sub = zeros(size(AT.AdjCT{atl_i}));
    end
    
    % For loop through roi pairs
    for j = 1:length(adjct_xsub)
        for k = 1:length(adjct_xsub)
            adjcts = AT.AdjCT{atl_i}{j,k};
            
            % if there is coverage for this roi pair
            if (~isempty(adjcts))
            
                adjcts(adjcts < AdjCT_thresh) = nan;

                % If at least n electrode pairs is significant between rois
                if (sum(isnan(adjcts)) >= min_n_epair_is_interact)
                    roi_is_interact = 1;
                else
                    roi_is_interact = 0;
                end
                adjct_xsub{j,k} = [adjct_xsub{j,k}, roi_is_interact];
                adjct_sub(j,k) = adjct_sub(j,k) + 1;

            end
            
        end
    end
    
    fprintf('\tdone.\n')
    
%     if (i == 5)
%         return
%     end
end

% consistency across subjects

adjcs = zeros(size(AT.AdjCT{atl_i}));
for j = 1:length(adjct_xsub)
    for k = 1:length(adjct_xsub)
        adjcts = adjct_xsub{j,k};

        % if there is coverage for this roi pair
        if (length(adjcts) >= min_n_sub)

            %adjct_xsub{j,k} = [adjct_xsub{j,k}, roi_is_interact];
            adjcs(j,k) = mean(adjcts)*n_sub;

        end

    end
end

% Get ROI colors
atl_table = AT.P.AtlROIs{atl_i}.RH.table;
atl_fname = 'coreg/lh.aparc.annot';
[v,l,c] = read_annotation(atl_fname);
roi_id = l(pa.vno + 1);

tri_color = zeros(length(roi_id),3);
col_unknown = 0.333*[1 1 1];

ind_lingual = 14;
ind_cuneus = 6;
ind_inferiorparietal = 9;
source = ind_inferiorparietal;
%source = 34;

n_colors = 100;
cc = corrcmap(n_colors+1);
adjcs_source = round(n_colors*(adjcs(source,:)/n_sub));

% Plot flat map
h = figure;
fig_w = 10.5;
fig_h = 10.5;
set(h,'Position',[0 0 fig_w*100 fig_h*100])
set(h, 'PaperUnits', 'Inches')
set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
a = trisurf(tri_flat,pa.x',pa.y',pa.z');
set(a,'edgecolor','none');
view(-90,90);
axis off;

% Colors
col_source = [0.2 1 0.2];

roi_xy = cell(1,length(adjct_xsub));
for i = 1:length(roi_id)
    
    c_i = find(roi_id(i) == atl_table(:,end),1);

    if (~isempty(c_i))
        %disp(c_i)
        roi_name = AT.P.AtlROIs{atl_i}.RH.struct_names{c_i};

        if (c_i == source)
            tri_color(i,:) = col_source;
        else

            tri_color(i,:) = cc(adjcs_source(c_i)+1,:);
        end


        %roi_xy{c_i} = [pa.x(), pa.y()];
    end
%     if (mod(i,round(length(roi_id)/10)) == 0)
%         text(pa.x(i),pa.y(i),roi_name)
%     end
%     if (~isempty(c_i))
%         tri_color(i,:) = atl_table(c_i,1:3)/255;
%     else
%         tri_color(i,:) = col_unknown;
%     end
end

%atl_tint_factor = 0.9;
%tri_color = tri_color + ([1 1 1] - tri_color) * atl_tint_factor;
set(a,'FaceVertexCdata',tri_color);
colormap(cc);
caxis([min(adjcs(source,:)) max(adjcs(source,:))])
cb = colorbar;

save(sprintf('xsub_atl-%i',atl_i),'-v7.3');