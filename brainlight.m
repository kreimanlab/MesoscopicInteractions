%
% Lighting config for brain surfaces plotted with trisurf
%
% p - surf handle
%

daspect([1 1 1]);
p.AmbientStrength = 0.3 ;
p.DiffuseStrength = 0.4 ;
p.SpecularStrength = 0;
p.SpecularExponent = 1;
p.BackFaceLighting = 'lit';
p.FaceLighting = 'gouraud';
cam_elev = 0;
camlight(-135,cam_elev);
camlight(45,cam_elev);
camlight(-225,cam_elev);
camlight(-45,cam_elev);