%   Plot ielvis electrode localization
%
%   Input:
%       sid - freesurfer subject string
%       chan_num - index of the channels to plot in original EEG order
%       atlas - string containing which atlas to highlight
%           'mmp' - multimodal parcellation
%           'd' - Destrieux Atlas (a2009s)
%           'dk' - Desikan-Killiany Atlas
%  - Markov-Kennedy 2014 Macaque Atlas
%           'fve' - Felleman-Van Essen 1991 Macaque Atlas
%           'fvea' - Felleman-V          'mk'an Essen all Macaque Atlas
%           <other> - no atlas coloring
%
%   Output:
%       h - figure handle
%
%   Example:
%       brainplot_one('sub1',[1 2]) - plot blank pial surface with electrodes
%       brainplot_one('sub1',[1 2],'mmp') - plot surface and electrodes colored with MMP atlas
%

function [ h ] = brainplot_one( sid, chan_num, atlas, roi_str )
% Jiarui Wang :: jwang04@g.harvard.edu
% Last edit: 2018-02-21 11:12


WIN_W = 1440; % Figure height in pixels
WIN_H = 900; % Figure width in pixels
MOVE_CAM = true; % Whether to compute the view using geometric center
ROUND_CAM_ANGLE = true; % Whether to round view angles to cardinal orientations
BRAIN_MESH_ALPHA = 1; % Surface transparancy between 0,1
COLOR_UNKNOWN = [1 1 1]; % RGB color of unknown ROIs
SMOOTH_FACE = true; % Use geometric smoothing when rendering surface
BLUR_ROI = true; % Blur face colors between ROI boundaries
ATLAS_SATURATION = 1; % Multiplier between 0,1 for atlas color saturation
ATLAS_TINT = 1*[1 1 1]; % Tint color
ATLAS_TINT_FACTOR = 0.55; % Tint factor: 0 for no tint, 1 for full
ROI_TXT_DARK_FACTOR = 1; % Dark factor: 0 for black text, 1 for color
ELEC_FONT_SIZE = 0.01;
ELEC_TXT_ZOFFSET = 1; % superior axis offset
ELEC_TXT_AOFFSET = 0; % anterior axis offset
ELEC_TXT_COLOR = 0*[1 1 1];
ELEC_HORIZ_ALIGN = 'center';
ELEC_RADIUS = 4;
SHOW_ROI_TXT = false;
COLOR_ELEC_WITH_ROI = false;
ROI_FONT_SIZE = 10;
%ATLAS_SHOW_ROI = true; % Whether to show ROI name
SURFACE_TYPE = 'pial'; % Which surface to plot (e.g. pial, white, inflated)


% Make sure $SUBJECTS_DIR is set
% [~,subjects_dir] = system('echo $SUBJECTS_DIR');
% subjects_dir = strip(subjects_dir);
%subjects_dir = '/mnt/cuenap_ssd/coregistration';
%subjects_dir = '/media/klab/internal/data/coreg';
subjects_dir = '../data/coregistration';
setenv('SUBJECTS_DIR',subjects_dir);
%fprintf('subjects_dir: %s\n',subjects_dir);

% Read label files
e_labels = read_label(sid,'all_surf_ielvis');
if (isempty(e_labels))
    fprintf('Trying: %s/%s/all_surf_ielvis.label\n',subjects_dir,sid)
    e_labels = read_label([subjects_dir,'/',sid],'all_surf_ielvis');
end

% Read electrode names
ElecNames = textscan(fopen(sprintf('%s/%s/elec_recon/%s.electrodeNames',...
    subjects_dir,sid,sid),'r'),'%s %s %s','HeaderLines',2);
ENames = ElecNames{1};
EHemi = ElecNames{3};
Hemi = unique(EHemi);

% Read ielvis pial coordinates
%Ecoords = textscan(fopen(sprintf('%s/%s/elec_recon/%s.PIAL',...
%    subjects_dir,sid,sid),'r'),'%f %f %f','HeaderLines',2);
%Ecoords = [Ecoords{1},Ecoords{2},Ecoords{3}];

% Parcellation file name
if (nargin == 2)
    atlas = '';
end
readAtlas = true;
aname = '';
switch atlas
    case 'mmp'
        aname = 'HCP-MMP1';
    case 'd'
        aname = 'aparc.a2009s';
    case 'dk'
        aname = 'aparc';
    case 'mk'
        aname = 'MACAQUE_M132';
    case 'fve'
        aname = 'MACAQUE_FVE91';
    case 'fvea'
        aname = 'MACAQUE_FVEall';
    case 'pht'
        aname = 'MACAQUE_PHT00';
    case 'landmark'
        aname = 'MACAQUE_LANDMARK';
    case 'lve'
        aname = 'MACAQUE_LVE00';
    case 'bb'
        aname = 'MACAQUE_BoninBailey';
    case 'fea'
        aname = 'MACAQUE_FerryEtAl';
    case 'ud'
        aname = 'MACAQUE_UD86';
    case 'brl'
        aname = 'MACAQUE_BRL87';
    case 'sp'
        aname = 'MACAQUE_SP78';
    case 'lk'
        aname = 'MACAQUE_LyonKaas';
    case 'pgr'
        aname = 'MACAQUE_PGR91';
    case 'b'
        aname = 'PALS_B12_Brodmann';
    case '7n'
        aname = 'Yeo2011_7Networks_N1000';
    case '17n'
        aname = 'Yeo2011_17Networks_N1000';
    otherwise
        readAtlas = false;
end

%h = figure;
%set(h,'Position',[0 0 WIN_W WIN_H])

% Load freesurfer files
e_labels_sav = e_labels;

% Index electrodes to plot
[n_el,~] = size(e_labels_sav);
ploti = false(n_el,1);
[~,sIdx] = sort(e_labels(:,end));
ploti(sIdx(chan_num)) = true;
ploti_sav = ploti;
% sort to label order
%ploti = ploti(sIdx);

ENames_sav = ENames;

% Predetermine how many hemispheres need to be lit
n_hemi_for_light = 0;%length(Hemi);
for hIdx = 1:length(Hemi)
    ploti = ploti_sav(strcmp(EHemi,Hemi{hIdx}));
    if (sum(ploti) ~= 0)
        n_hemi_for_light = n_hemi_for_light + 1;
    end
end

for hIdx = 1:length(Hemi)
    hemi = lower(Hemi{hIdx});

    %disp(size(ploti))
    %disp(size(strcmp(EHemi,Hemi{hIdx})))
    % Condition plot index to hemisphere
    ploti = ploti_sav(strcmp(EHemi,Hemi{hIdx}));
    
    % Skip hemisphere if no electrodes are found
    if (sum(ploti) ~= 0)
    

        % Condition electrode labels to hemisphere
        e_labels = e_labels_sav(strcmp(EHemi,Hemi{hIdx}),:);
        ENames = ENames_sav(strcmp(EHemi,Hemi{hIdx}),:);

        % Load surface
        [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,SURFACE_TYPE));

        % Load atlas
        C = ones(size(s_vert));
        if (readAtlas)
            %C = ones(size(s_vert));
            [a_vert,l,c] = read_annotation(sprintf('%s/%s/label/%sh.%s.annot',... 
                subjects_dir,sid,hemi,aname));
            c_norm = c.table(:,1:3) / max(max(c.table(:,1:3)));
            for i = 1:length(C)
                try
                    %fprintf('%i: %i - %i\n',tblIdx,l(i),c.table(tblIdx,5));
                    roi_color = c_norm(l(i) == c.table(:,5),:);
                    [nRtc,~] = size(roi_color);
                    if (nRtc > 1)
                        %roi_color = roi_color(1,:);
                        roi_color = mean(roi_color);
                    end

                    C(i,:) = roi_color;
                catch
                    C(i,:) = COLOR_UNKNOWN;
                end
            end
        end
        if strcmp(atlas,char([97,108,111,104,97]));C=jet(length(s_vert(:,1)));end

        % Saturaton
        Chsv = rgb2hsv(C);
        Chsv(:,2) = Chsv(:,2) * ATLAS_SATURATION;
        C = hsv2rgb(Chsv);
        C = C + (ATLAS_TINT - C) * ATLAS_TINT_FACTOR;

        % Plot surface
        p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
        h = p;
        %return
        set(p,'FaceVertexCdata',C);
        ax = gca;
        ax.Clipping = 'off';
        hold on;

        % Plot electrodes
        e_labels = e_labels(ploti,1);
        e_coords = s_vert(e_labels(:,1) + 1,:);
        %e_coords = Ecoords;
        for i = 1:length(e_labels(:,1))
            e_x = e_coords(i,1);
            e_y = e_coords(i,2);
            e_z = e_coords(i,3);

            %e_x = e_labels(i,2);
            %e_y = e_labels(i,3);
            %e_z = e_labels(i,4);
            %e_txt = ENames{e_labels(i,5)};
            e_txt = ENames{i};
            toff = 1.03;
            text(toff*e_x,toff*(e_y+ELEC_TXT_AOFFSET),toff*(e_z+ELEC_TXT_ZOFFSET),e_txt,'HorizontalAlignment',ELEC_HORIZ_ALIGN,'FontSize',ELEC_FONT_SIZE,'Color',ELEC_TXT_COLOR);

            if (readAtlas)
                % Find ROI for electrode
                [~,roiIdx] = min(sum((e_coords(i,1:3) - s_vert).^2,2));
                try
                    roiCIdx = c.table(:,5) == l(roiIdx);
                    roi_txt = c.struct_names{roiCIdx};
                    roi_txt = strsplit(roi_txt,'_');
                    if (strcmp(atlas,'mmp'))
                        roi_txt = roi_txt{2};
                    %elseif strcmp(atlas,'mk')
                    %    roi_txt = roi_txt{1};
                    %elseif strcmp(atlas,'fvea')
                    %    roi_txt = roi_txt{3};
                    elseif strcmp(atlas,'b')
                        roi_txt = strsplit(roi_txt{1},'.');
                        roi_txt = strjoin(roi_txt(~strcmp(roi_txt,'Brodmann')),'.');
                    else
                        roi_txt = strjoin(roi_txt,',');
                    end
                    roi_txt_color = c_norm(roiCIdx,:);
                catch
                    roi_txt = 'UNKNOWN';
                    roi_txt_color = c_norm(1,:);
                end
                if (SHOW_ROI_TXT)
                    [nRtc,~] = size(roi_txt_color);
                    if (nRtc > 1)
                        %roi_txt_color = roi_txt_color(1,:);
                        roi_txt_color = mean(roi_txt_color);
                    end
                    roi_txt_color = roi_txt_color*ROI_TXT_DARK_FACTOR;
                    text(toff*e_x,toff*(e_y+ELEC_TXT_AOFFSET),toff*(e_z-ELEC_TXT_ZOFFSET),roi_txt,'HorizontalAlignment',ELEC_HORIZ_ALIGN,'FontSize',ROI_FONT_SIZE,'Color',roi_txt_color);
                end
            end

            n = 21;
            s_color = zeros(n+1,n+1,3);
            if (readAtlas && COLOR_ELEC_WITH_ROI)
                for j = 1:3
                    s_color(:,:,j) = 0.35*roi_txt_color(j);
                end
            end
            [x,y,z] = sphere(n);
            radius = ELEC_RADIUS;
            s = surf(radius*x+e_x,radius*y+e_y,radius*z+e_z,s_color);
            hold on;
            s.FaceLighting = 'none';
            s.EdgeColor = 'none';
        end

        if (MOVE_CAM)
            % Calculate the geometric center
            geocenter = mean(s_vert);
            pt = e_coords(i,1:3);
            cam = pt - geocenter;
            [n_pt, ~] = (size(pt));
            if (n_pt ~= 0)
                %disp(size(geocenter))
                %cam = pt;
                cam_az = (atan2(cam(2),cam(1)) + pi/2)*(360/(2*pi));
                % the higher, the more likely to look from bottom
                view_down_thresh = 0;%pi/10;
                cam_el = (atan2(cam(3),sqrt(sum(cam.^2))) + view_down_thresh)*(360/(2*pi));
                
                cam_az_s = cam_az;
                cam_el_s = cam_el;
                % Round to nearest cardinal angle
                if (ROUND_CAM_ANGLE)
                    %cam_az = round((cam_az)/180)*180 + 90;
                    if ((cam_az >= 0) && (cam_az < 180))
                        cam_az = 90;
                    elseif ((cam_az >= 180) && (cam_az < 360))
                        cam_az = -90;
                    elseif ((cam_az < 0) && (cam_az >= -180))
                        cam_az = -90;
                    elseif ((cam_az < -180) && (cam_az >= -360))
                        cam_az = 90;
                    end
                    
                    if (cam_el < -75) %if (cam_el < -10)
                        cam_el = -90;
                    else
                        cam_el = 0;
                    end
                    %cam_el = round(cam_el/90)*90;
                end
            end
            %disp([cam_az,cam_el])
        end

        % Lighting
        axis tight;
        daspect([1 1 1]);
        view(0,0);
        if (SMOOTH_FACE)
            p.FaceLighting = 'gouraud';
        end
        % Intensity is reduced for 2 hemisphere plots
        p.AmbientStrength = 0.35 + 0.05*(n_hemi_for_light-1);
        if (cam_el < 0)
            p.DiffuseStrength = 0.33 - 0.25*(n_hemi_for_light-1);
        else
        p.DiffuseStrength = 0.4 - 0.25*(n_hemi_for_light-1);
        end
        p.SpecularStrength = 0;
        p.SpecularExponent = 1;
        p.BackFaceLighting = 'lit';
        
        cspread = 35;
%         if (cam_el >= 0)
        camlight(cam_az+180+cspread,cam_el * (-1) - 10);
        camlight(cam_az+180-cspread,cam_el * (-1) - 10);
%         else
%             
%         camlight(cam_az+180,cam_el * (-1));
%         end
%         cam_elev = 2;
%         ofc = 0;
%         camlight(-135,cam_elev+ofc);
%         camlight(45,cam_elev+ofc);
%         camlight(-225,cam_elev-ofc);
%         camlight(-45,cam_elev-ofc);
        %camlight(0,90);
        axis vis3d off;
        %lighting phong; % smooth lighting
        %material dull;
        if (BLUR_ROI)
            shading interp;
        else
            shading flat;
        end
    
    
    else
        n_hemi_for_light = 1;
    end
end

% Adjust final view
zoom(1)

if (MOVE_CAM)
    view([cam_az cam_el])
else
    if (length(Hemi) == 1)
        if strcmp(Hemi{1},'L')
            view(-90,0);
        end
    elseif (length(Hemi) == 2)
        view(90,90);
    else
        fprintf(2,'[!] Warning: I only understand 1 or 2 hemispheres.\n');
    end
end


end
