close all;
clear;

% Only show electrodes for subjects found significant for fig T16
trig_plot_t16_only = true;


Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124',...
    'mSu'};


sp_sid = {"m00001",...
    "m00001",...
    "m00001",...
    "m00001",...
    "m00001",...
    "m00006",...
    "m00006",...
    "m00032",...
    "m00032",...
    "m00044",...
    "m00071",...
    "m00071",...
    "m00071",...
    "m00071",...
    "m00071",...
    "m00075",...
    "m00075",...
    "m00075",...
    "m00084",...
    "m00084",...
    "m00084",...
    "m00084",...
    "m00084"};

sp_elec = {"LF7",...
    "LT23",...
    "LT15",...
    "PT7",...
    "LT31",...
    "RLT23",...
    "RLT31",...
    "LSF7",...
    "LSF15",...
    "AT15",...
    "TA16",...
    "TT6",...
    "TT14",...
    "TT7",...
    "TT15",...
    "AT8",...
    "HD30",...
    "HD21",...
    "LLT12",...
    "LLT20",...
    "LLT28",...
    "LLT54",...
    "LLT55"};

mkdir('figures');
mkdir(sprintf('figures/T19'));

trig_col_elec_roi = true;

% First, plot average brain
subjects_dir = '/media/klab/internal/data/coreg';
h = figure('visible','off'); set(h,'PaperUnits','inches'); set(h,'PaperPosition',[0 0 8 8]);
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,'fsaverage_sym','r','pial'));
%const_col_surface = [1 1 1];
const_col_surface = 0.9 * ones(1,3);
const_alpha_surface = 1;
const_elec_surface = 0*[1 1 1];
const_elec_surface2 = 0*[1 0 0];
const_elec_low = 0.75*[1 1 1]; % 0.4
p = trisurf(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
    'EdgeColor','none','FaceVertexCData',const_col_surface,'FaceAlpha',const_alpha_surface);

brainlight;
% daspect([1 1 1]);
% p.AmbientStrength = 0.5 ;
% p.DiffuseStrength = 0.2 ;
% p.SpecularStrength = 0;
% p.SpecularExponent = 20;
% p.BackFaceLighting = 'lit';
% p.FaceLighting = 'gouraud';
% cam_elev = 0;
% camlight(-135,cam_elev);
% camlight(45,cam_elev);
% camlight(-225,cam_elev);
% camlight(-45,cam_elev);

view(90,0)
hold all;

% Plot speech-responsive electrodes according to stimulation
n_sp = length(sp_elec);
iM = 1;
legend_rois = {};
legend_rois_col = {};

esm_sids = {};
esm_subs = [];
esm_bchans = [];
esm_chans = [];

for i = 1:n_sp
    sid = sp_sid{i};
    elec = sp_elec{i};
    
    % Read label
    l = read_label(sprintf('%s/%s',subjects_dir,'fsaverage_sym'),sprintf('ielvis_%s',sid));
    %l = [];
    if (isempty(l))
        l = read_label(sprintf('%s','fsaverage_sym'),sprintf('ielvis_%s',sid));
    end
    
    % Find electrode index
    Ca = load(sprintf('cache/xsub_out_%s_%i_atl2.mat',sid,iM));
    e1 = find(strcmp(Ca.C.EleLabels,elec));
    if (isempty(e1))
        fprintf(2,'[%i/%i] Electrode %s, subject %s not found.\n',i,n_sp,elec,sid);
    else
        fprintf('[%i/%i] Found: electrode %s, subject %s at index %i.\n',i,n_sp,elec,sid,e1);
        [x,y,z] = sphere(64);
        x = x * 2.1;
        y = y * 2.1;
        z = z * 2.1;
        coord = l(l(:,end)==e1,2:4);
        
        % remap coord
        [remap_dist,r_idx] = min(sqrt(sum((s_vert - [(1)*coord(1),coord(2),coord(3)]).^2,2)));
        coord = s_vert(r_idx,:);
        fprintf('\tremap (mm): %.12f\n',remap_dist);
        
        % show subject number
        n_sub = find(strcmp(Subjects,sid));
        n_bchan = find(Ca.ecog.bip(:,1)==e1);
        if (isempty(n_bchan))
            n_bchan = find(Ca.ecog.bip(:,2)==e1);
        end
        fprintf('\tsubject num: %i, bipolar num: %i, raw elec num: %i\n',n_sub,n_bchan,e1);
        esm_sids = [esm_sids; {sid}];
        esm_subs = [esm_subs; n_sub];
        esm_bchans = [esm_bchans; n_bchan];
        esm_chans = [esm_chans; e1];
        
        % show roi
        e1_roi = Ca.C.AtlLabels{2}{e1};
        fprintf('\tRegion: %s\n',e1_roi);
        e1_roi_col = Ca.C.AtlROIs{2}.LH.table(strcmp(Ca.C.AtlROIs{2}.LH.struct_names,e1_roi),1:3)/255;
        if (trig_col_elec_roi)
            col_ep = e1_roi_col;
            % add roi to legend
            if (~ any(strcmp(legend_rois,e1_roi)))
                legend_rois = [legend_rois; {e1_roi}];
                legend_rois_col = [legend_rois_col; {e1_roi_col}];
            end
        
            % Tint
            tint_fac = 0.15;
            col_ep = col_ep * (1 - tint_fac) + [1 1 1] * (tint_fac);
        else
            col_ep = 0.2*[1 1 1];
        end
        
        
        
%         if (strcmp(e1_roi,'superiortemporal'))
%             col_ep = [0 1 0]*0.8;
%         else
%             col_ep = [0 0 1]*0.9;
%         end

        if (trig_plot_t16_only)
            CT = load('./cache_figure_t16_dk.mat');
            
            cond_show_elec = any(strcmp(CT.usub_str,sid));
        else
            cond_show_elec = true;
        end

        if (cond_show_elec)
            q = surf(x+coord(1),y+coord(2),z+coord(3),'FaceColor',col_ep,'EdgeColor','none');
            q.SpecularStrength = 0;
        end
        %plot3(coord(1),coord(2),coord(3),'black.')
        %return
    end
    %return
end

save('./cache/figure_t19','esm_sids','esm_subs','esm_bchans','esm_chans');


% Plot m00005 example electrodes
sid = 'm00005';
roi_1 = 'superiortemporal';
roi_2 = 'parsopercularis';
bchan1_const = 35;
bchan2_const = 84;
iM = 1;
Ca = load(sprintf('cache/xsub_out_%s_%i_atl2.mat',sid,iM));

% Read label
l = read_label(sprintf('%s/%s',subjects_dir,'fsaverage_sym'),sprintf('ielvis_%s',sid));
%l = [];
if (isempty(l))
    l = read_label(sprintf('%s','fsaverage_sym'),sprintf('ielvis_%s',sid));
end

% Plot
e1 = Ca.ecog.bip(bchan1_const,1);
e2 = Ca.ecog.bip(bchan2_const,1);
coord1 = l(l(:,end)==e1,2:4);
coord2 = l(l(:,end)==e2,2:4);

% remap coord
[remap_dist,r_idx] = min(sqrt(sum((s_vert - [(1)*coord1(1),coord1(2),coord1(3)]).^2,2)));
coord1 = s_vert(r_idx,:);
fprintf('[*] Chan 1, remap (mm): %.12f\n',remap_dist);
[remap_dist,r_idx] = min(sqrt(sum((s_vert - [(1)*coord2(1),coord2(2),coord2(3)]).^2,2)));
coord2 = s_vert(r_idx,:);
fprintf('[*] Chan 1, remap (mm): %.12f\n',remap_dist);



%--------------------------------------------------------------------------
% Plot electrode locations from Figure 2A example

% %q = surf(x+coord1(1),y+coord1(2),z+coord1(3),'FaceColor',const_elec_surface2,'EdgeColor','none');
% %q.SpecularStrength = 0;
% plot3(coord1(1)+20,coord1(2),coord1(3),'.','Color',const_elec_surface2,'MarkerSize',28);
% %q = surf(x+coord2(1),y+coord2(2),z+coord2(3),'FaceColor',const_elec_surface2,'EdgeColor','none');
% %q.SpecularStrength = 0;
% plot3(coord2(1)+20,coord2(2),coord2(3),'.','Color',const_elec_surface2,'MarkerSize',28);
%--------------------------------------------------------------------------



% Recolor surface with ROI highlighted
% [v_a, l, ct] = read_annotation(sprintf('%s/%s/label/%sh.aparc.annot',subjects_dir,sid,hemi));
% [n_vert,~] = size(s_vert);
% const_col_surface2 = repmat(const_col_surface,[n_vert 1]);
% 
% for ihi = 1:ct.numEntries
%     cond_to_color = strcmp(ct.struct_names{ihi},b1_roi) | strcmp(ct.struct_names{ihi},b2_roi);
%     if (cond_to_color)
%         roi_id = ct.table(ihi,end);
%         roi_col = ct.table(ihi,1:3) / 255;
%         n_vtx = sum(l==roi_id);
%         const_col_surface2(l==roi_id,:) = repmat(roi_col,[n_vtx 1]);
%     end
% end

%return

axis off
print(h,'figures/T19/figure_t19','-depsc');
print(h,'figures/T19/figure_t19','-dpng','-r600');
close(h);



% roi legend
if (trig_col_elec_roi)
    h = figure('visible','off'); set(h,'PaperUnits','inches'); set(h,'PaperPosition',[0 0 8 8]);
    axis off;
    idx_rearrange = [3 2 1 4];
    legend_rois = legend_rois(idx_rearrange);
    legend_rois_col = legend_rois_col(idx_rearrange);
    ax = gca;
    x_offset = 0.7;
    for j = 1:length(legend_rois)
        annotation('textbox',[x_offset ax.Position(2)+j*0.02 0.2 0.02],...
            'String',sprintf('%s',convertRoiDK(legend_rois{j})),'FitBoxToText','on','EdgeColor','none','VerticalAlignment','middle');
        annotation('rectangle',[x_offset-0.02 ax.Position(2)+j*0.02 0.02 0.02],...
            'FaceColor',legend_rois_col{j});
    end
    print(h,'figures/T19/figure_t19_legend','-depsc');
    print(h,'figures/T19/figure_t19_legend','-dpng','-r600');
    close(h);
end

