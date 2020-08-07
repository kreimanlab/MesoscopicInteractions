close all;

if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        resultsDir = '/media/klab/44/data/results';
        h5Dir = '/media/klab/44/h5';
    else
        resultsDir = '/mnt/cuenap2/data/results';
        h5Dir = '/mnt/cuenap2/scripts/synth/out';
    end
end

if (~exist('A','var'))
    A = Analysis(resultsDir,h5Dir);
end

% Read xsub file
load('xsub-s_Destrieux-2009.mat');
%load('xsub-sP_Destrieux-2009.mat');
%load('xsub-scBroadband_Destrieux-2009.mat');
%load('xsub-s_Desikan-Killiany.mat');
%load('xsub-sP_Desikan-Killiany.mat');
%load('result_xsub-scBroadband_Desikan-Killiany.mat');
%load('xsub-sP_HCP-MMP1.mat');

% alias


% Get ROI colors
AdjCS(AdjCS == -1) = NaN;
adjcs = AdjCS_mag;

% Read annotation
atl_table = AT.P.AtlROIs{atl_i}.RH.table;
annot_dir = './coreg';
switch(AT.P.AtlNames{atl_i})
    case ('Desikan-Killiany')
        atl_fname = [annot_dir,'/lh.aparc.annot'];
    case ('Destrieux-2009')
        atl_fname = [annot_dir,'/lh.aparc.a2009s.annot'];
    case ('HCP-MMP1')
        atl_fname = [annot_dir,'/lh.HCP-MMP1.annot'];
    case ('M132')
        atl_fname = [annot_dir,'/lh.MACAQUE_M132.annot'];
end
[v,l,c] = read_annotation(atl_fname);
roi_id = l(pa.vno + 1);
tri_color = zeros(length(roi_id),3);

% constant color
col_unknown = 0.333*[1 1 1];
col_source = [0.2 1 0.2];
% 
% ind_lingual = 14;
% ind_cuneus = 6;
% ind_inferiorparietal = 9;
% source = ind_inferiorparietal;
% %source = 34;

n_colors = 100;
cc = corrcmap(n_colors);
%adjcs_source = round(n_colors*(adjcs(source,:)));

%%
% Plot flat map
h = figure;
fig_w = 10.5;
fig_h = 10.5;
set(h,'Position',[0 0 fig_w*100 fig_h*100])
set(h, 'PaperUnits', 'Inches')
set(h, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
set(h, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
a = trisurf(tri_flat,pa.x',pa.y',pa.z'); hold all;
set(a,'edgecolor','none');
view(90,90);
axis off;

roi_xy = cell(1,length(adjct_xsub));
for i = 1:length(roi_id)
    
    c_i = find(roi_id(i) == atl_table(:,end),1);

    if (~isempty(c_i))
        %disp(c_i)
        roi_name = AT.P.AtlROIs{atl_i}.RH.struct_names{c_i};

%         if (c_i == source)
%             tri_color(i,:) = col_source;
%         else
            %tri_color(i,:) = cc(adjcs_source(c_i)+1,:);
        tri_color(i,:) = atl_table(c_i,1:3)/255;
%         end
        roi_xy{c_i} = [roi_xy{c_i}; [pa.x(i) pa.y(i)]];


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


% Show boundary
for i = 1:length(tri_flat)
    vB = tri_flat(i,:);
    tv = roi_id(vB);
    if (length(unique(tv)) > 1)
        plot(mean(pa.x(vB)),mean(pa.y(vB)),'black.')
    end
end
% for i = 1:length(roi_xy)
%     if (~isempty(roi_xy{i}))
%         color_boundary = [0 0 0];
%         k = boundary(roi_xy{i}(:,1),roi_xy{i}(:,2));
%         plot(roi_xy{i}(k,1),roi_xy{i}(k,2),'-','color',color_boundary); hold on;
%     end
% end
%atl_tint_factor = 0.9;
%tri_color = tri_color + ([1 1 1] - tri_color) * atl_tint_factor;
set(a,'FaceVertexCdata',tri_color);
%colormap(cc);
%caxis([min(adjcs(source,:)) max(adjcs(source,:))])
%cb = colorbar;

%save(sprintf('xsub_atl-%i',atl_i),'-v7.3');

