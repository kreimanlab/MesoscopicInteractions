close all;
clear;

dir_cor = '/media/klab/internal/data/coreg';

setenv('SUBJECTS_DIR',dir_cor);

% Patients
Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

for iSub = 1:length(Subjects)
    % Load localizations
    sid = Subjects{iSub};
    C = load(sprintf('%s/%s/label/all_parcellation.mat',dir_cor,sid));
    hemi = lower(C.EleHemi{1});
    
    % Load spherical
    fn_pial = sprintf('%s/%s/surf/%sh.sphere',dir_cor,sid,hemi);
    [vertex_coords_sph, faces_sph] = read_surf(fn_pial);
    TR_sph = triangulation(faces_sph+1,vertex_coords_sph);
    
    % Load pial
    fn_pial = sprintf('%s/%s/surf/%sh.pial',dir_cor,sid,hemi);
    [vertex_coords, faces] = read_surf(fn_pial);
    TR = triangulation(faces+1,vertex_coords);
    
    % Electrode pairs
    n_chan = length(C.EleLabels);
    for i1 = 1:(n_chan-1)
        for i2 = (i1+1):n_chan
            ei1 = find(C.EleCoords(:,5) == i1);
            ei2 = find(C.EleCoords(:,5) == i2);
            i1 = C.EleCoords(ei1,1);
            i2 = C.EleCoords(ei2,1);
            co_sph1 = vertex_coords_sph(i1,:);
            co_sph2 = vertex_coords_sph(i2,:);
            
            %co1 = vertex_coords(i1,:);
            %co2 = vertex_coords(i2,:);
            %coord1 = C.EleCoords(ei1,2:4);
            %coord2 = C.EleCoords(ei2,2:4);
            %diff1 = sqrt(sum((vertex_coords(i1,:) - coord1).^2));
            %diff2 = sqrt(sum((vertex_coords(i2,:) - coord2).^2));
            return
        end
    end
    
end