close all;
clear;

fn_paths = 'brainexport/m00124_5all.mat';
fprintf('[*] Loading %s ..\n',fn_paths)
load(fn_paths);

tube_radius = 0.1;
tube_n_edges = 3; % number of vertices of tube cross-section

% Construct mesh
for n = 1:n_counts
%for n = 13
    path = Paths{n};
    [n_path,~] = size(path);
    
%     Xc = [];
%     Yc = [];
%     Zc = [];
%     n_cylinders = n_path - 1;
%     for i = 1:n_cylinders
%         % Make cylinders
%         pt1 = path(i,:);
%         pt2 = path(i+1,:);
%         
%         % Scale by height
%         cyl_height = norm(pt2 - pt1);
% %         Xt = pt1(1) + X;
% %         Yt = pt1(2) + Y;
% %         Zt = pt1(3) + cyl_height*Z;
%         Xt = X;
%         Yt = Y;
%         Zt = cyl_height*Z;
%         
%         % Rotate
%         ptD = (pt2 - pt1)/norm(pt2 - pt1);
%         %rot = [atan2(ptD(2),ptD(1)), 0, acos(ptD(2)/1)];
%         rot = [atan2(ptD(3),ptD(2)), ...
%                atan2(ptD(1),ptD(3)), ...
%                atan2(ptD(2),ptD(1))];
%            
%         rotM = (rotationVectorToMatrix(rot));
%         %rotM = eye(3);
%         
%         pt10 = rotM  * [Xt(1,:); Yt(1,:); Zt(1,:);];
%         pt20 = rotM  * [Xt(2,:); Yt(2,:); Zt(2,:);];
%         
%         Xt = [pt10(1,:); pt20(1,:)];
%         Yt = [pt10(2,:); pt20(2,:)];
%         Zt = [pt10(3,:); pt20(3,:)];
%         
% %         fvc = surf2patch(Xt,Yt,Zt);
% %         trisurf(fvc.faces,fvc.vertices(:,1),fvc.vertices(:,2),fvc.vertices(:,3)); hold all
% %         daspect([1 1 1]);
% %         return
% %         qp = ptD ; % rotM * (ptD');
% %         quiver3(pt1(1),pt1(2),pt1(3),qp(1),qp(2),qp(3),5,'red'); hold on;
% 
%         Xc = [Xc; Xt+pt1(1)];
%         Yc = [Yc; Yt+pt1(2)];
%         Zc = [Zc; Zt+pt1(3)];
%     end
%     fvc = surf2patch(Xc,Yc,Zc);
%     trisurf(fvc.faces,fvc.vertices(:,1),fvc.vertices(:,2),fvc.vertices(:,3)); hold all
%     daspect([1 1 1]);

    
    tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]); hold on;
    
end


% debug
subjects_dir = '/media/jerry/internal/data/coreg';
SURFACE_TYPE = 'pial';
hemi = 'r';
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,SURFACE_TYPE));
p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
    'EdgeColor','none','FaceColor',0.9*[1 1 1],'FaceAlpha',0.5); 
daspect([1 1 1]);
p.AmbientStrength = 0.3 ;
p.DiffuseStrength = 0.4 ;
p.SpecularStrength = 0;
p.SpecularExponent = 1;
p.BackFaceLighting = 'lit';
cam_elev = 0;
camlight(-135,cam_elev);
camlight(45,cam_elev);
camlight(-225,cam_elev);
camlight(-45,cam_elev);

% 
% for n = 1:n_counts
% plot3(Paths{n}(:,1),Paths{n}(:,2),Paths{n}(:,3),'-'); hold all;
% end
% daspect([1 1 1])

fprintf('[!] Done.\n')