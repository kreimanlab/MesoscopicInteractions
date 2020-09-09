close all;
clear;

dir_cor = '/media/klab/internal/data/coreg';

setenv('SUBJECTS_DIR',dir_cor);

% Patients
Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

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