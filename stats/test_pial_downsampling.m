close all;
clear;

subjects_dir = '/media/jerry/internal/data/coreg';
sid_const = 'm00005'; %"fsaverage_sym";
hemi = 'r';
SURFACE_TYPE = 'pial';

[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,'pial'));
FV = triangulation(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3));
center = mean(s_vert,1);

[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,SURFACE_TYPE));
s_vert(:,1) = s_vert(:,1) - center(1);
s_vert(:,2) = s_vert(:,2) - center(2);
s_vert(:,3) = s_vert(:,3) - center(3);

nFV = reducepatch(faces+1,s_vert,0.1);

subplot(1,2,1)
p = trisurf(FV,'EdgeColor','none','FaceColor',[1 1 1]*0.8,'FaceAlpha',0.1);
daspect([1 1 1]);
view(90,0);
brainlight;
fprintf('[*] original faces:\t%i, verts: %i\n',length(faces(:,1)),length(s_vert(:,1)));

subplot(1,2,2)
p = trisurf(nFV.faces,nFV.vertices(:,1),nFV.vertices(:,2),nFV.vertices(:,3),'EdgeColor','none','FaceColor',[1 1 1]*0.8,'FaceAlpha',0.1);
daspect([1 1 1]);
view(90,0);
brainlight;
fprintf('[*] reduced faces:\t%i, verts: %i\n',length(nFV.faces(:,1)),length(nFV.vertices(:,1)));
