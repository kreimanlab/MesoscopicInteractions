close all;
clear;

dir_cache = './cache';
sid_const = 'm00124';
scale_mm_per_obj = 30;

mkdir('brainexport');
mkdir(sprintf('brainexport/%s',sid_const));

% Subject specific
fn_paths = sprintf('brainexport/%s_5all.mat',sid_const);
fprintf('[*] Loading %s ..\n',fn_paths)
load(fn_paths);
Paths = Paths(Paths_ind);
fn_cache = sprintf('%s/xsub_out_%s_1.mat',dir_cache,sid);
Ca = load(fn_cache);

tube_radius = 0.15;
tube_n_edges = 3; % number of vertices of tube cross-section
rotm = rotationVectorToMatrix([pi/2, 0, 0]);

% Construct mesh
n_col = 2^8;
cmap = corrcmap(n_col);
%cmap = jet(n_col);
pass_dist = (Ca.Dmats > Ca.dist_thresh);
pass_ct = (Ca.ct > Ca.ct_thresh);
mag_all = Ca.mag(pass_dist & pass_ct);
mag_max = max(mag_all);
mag_min = min(mag_all);
mag2col = @(x) cmap(round(((x - mag_min)/(mag_max - mag_min))*(n_col-1) + 1),:);


col_tube_spec = [1 1 1];
for n = 1:n_counts
    
    if (pass_dist(n) && pass_ct(n))
        path = Paths{n};
        [n_path,~] = size(path);

        
        col_tube = mag2col(Ca.mag(n));
        [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
        surf(X,Y,Z,'FaceColor',col_tube,'EdgeColor','none'); hold on;
        pa = surf2patch(X,Y,Z,'triangles');
        
        bchan1 = Count(Count(:,3) == n,1);
        bchan2 = Count(Count(:,3) == n,2);
        
        clear obj;
        s_vert_obj = pa.vertices;
        s_vert_obj = (rotm * (s_vert_obj'))';
        %s_vert_obj(:,1) = pa.vertices(:,1);
        %s_vert_obj(:,2) = pa.vertices(:,3);
        %s_vert_obj(:,3) = pa.vertices(:,2);
        mtl_name = sprintf('pair_%i_%i_%i_mag-%i',n,bchan1,bchan2,round(Ca.mag(n)*1000));
        obj.vertices = s_vert_obj / scale_mm_per_obj;
        obj.objects(1).type='usemtl'; % use material
        obj.objects(1).data=mtl_name;
        obj.objects(2).type='s'; % smooth shading
        obj.objects(2).data='off'; % 1 for on, "off" for off
        obj.objects(3).type='f';
        obj.objects(3).data.vertices = pa.faces;
        obj.material(1).type='newmtl';
        obj.material(1).data=mtl_name;
        obj.material(2).type='Ka'; % ambient color
        obj.material(2).data=col_tube;
        obj.material(3).type='Kd'; % diffuse color
        obj.material(3).data=col_tube;
        obj.material(4).type='Ks'; % specular color
        obj.material(4).data=col_tube; %col_tube_spec;
        obj.material(5).type='illum';
        obj.material(5).data=2;
        obj.material(6).type='Ns'; % 0 to 1000, specular exponent
        obj.material(6).data=0;
        currdir = pwd;
        fn_obj = sprintf('%s/brainexport/%s/%s_%s.obj',currdir,sid,sid,mtl_name);
        fprintf('[*] Saving pair to: %s\n',fn_obj);
        write_wobj(obj,fn_obj);
        
%         if (n > n_counts)
%             break
%         end
        
    end
end

% hemi
b1c1 = ecog.bip(1,1);
hemi_11 = elec_hemis(b1c1);
hemi = lower(hemi_11{1});

% Plot electrodes
n_elec_pts = 8;
r_elec = 2; % (mm)
col_elec = 0*[1 1 1];
col_elec_spec = [1 1 1];
for i = 1:(ecog.n_bchan)
    [x,y,z] = sphere(n_elec_pts);
    coord = ecog.bip(i,4:6);
    x = x + coord(1);
    y = y + coord(2);
    z = z + coord(3);
    mesh(x,y,z,'EdgeColor','none','FaceColor',col_elec); hold on;
    pa = surf2patch(x,y,z,'triangles');
    
    % Export to obj file
    clear obj;
    s_vert_obj = pa.vertices;
    s_vert_obj = (rotm * (s_vert_obj'))';
    %s_vert_obj(:,1) = pa.vertices(:,1);
    %s_vert_obj(:,2) = pa.vertices(:,3);
    %s_vert_obj(:,3) = pa.vertices(:,2);
    %mtl_name = sprintf('electrode%i',ecog.bip(i,1));
    mtl_name = sprintf('electrode%i',i);
    obj.vertices = s_vert_obj / scale_mm_per_obj;
    obj.objects(1).type='usemtl'; % use material
    obj.objects(1).data=mtl_name;
    obj.objects(2).type='s'; % smooth shading
    obj.objects(2).data='off'; % 1 for on, "off" for off
    obj.objects(3).type='f';
    obj.objects(3).data.vertices = pa.faces;
    obj.material(1).type='newmtl';
    obj.material(1).data=mtl_name;
    obj.material(2).type='Ka'; % ambient color
    obj.material(2).data=col_elec;
    obj.material(3).type='Kd'; % diffuse color
    obj.material(3).data=col_elec;
    obj.material(4).type='Ks'; % specular color
    obj.material(4).data=col_elec; %col_elec_spec;
    obj.material(5).type='illum';
    obj.material(5).data=2;
    obj.material(6).type='Ns'; % 0 to 1000, specular exponent
    obj.material(6).data=0;
    currdir = pwd;
    fn_obj = sprintf('%s/brainexport/%s/%s_%s.obj',currdir,sid,sid,mtl_name);
    fprintf('[*] Saving electrode to: %s\n',fn_obj);
    write_wobj(obj,fn_obj);

end

col_pial = 0.6*[1 1 1];
col_pial_spec = [1 1 1];
alpha_pial = 0.8;
subjects_dir = '/media/jerry/internal/data/coreg';
SURFACE_TYPE = 'pial';


% Save surface object
currdir = pwd;
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,SURFACE_TYPE));
fn_obj = sprintf('%s/brainexport/%s/%s_%s.obj',currdir,sid,sid,SURFACE_TYPE);
fprintf('[*] Saving brain to: %s\n',fn_obj);
clear obj;
s_vert_obj = s_vert;
s_vert_obj = (rotm * (s_vert_obj'))';
%s_vert_obj(:,1) = s_vert(:,1);
%s_vert_obj(:,2) = s_vert(:,3);
%s_vert_obj(:,3) = s_vert(:,2);
obj.vertices = s_vert_obj / scale_mm_per_obj;
obj.objects(1).type='usemtl'; % use material
obj.objects(1).data='pial';
obj.objects(2).type='s'; % smooth shading
obj.objects(2).data='1'; % 1 for on, "off" for off
obj.objects(3).type='f';
obj.objects(3).data.vertices = faces + 1;
obj.material(1).type='newmtl';
obj.material(1).data='pial';
obj.material(2).type='Ka'; % ambient color
obj.material(2).data=col_pial;
obj.material(3).type='Kd'; % diffuse color
obj.material(3).data=col_pial;
obj.material(4).type='Ks'; % specular color
obj.material(4).data=col_pial_spec;
obj.material(5).type='illum';
obj.material(5).data=4;
obj.material(6).type='Ns'; % 0 to 1000, specular exponent
obj.material(6).data=5;
obj.material(7).type='d';
obj.material(7).data=alpha_pial;
obj.material(7).type='Tr';
obj.material(7).data=(1-alpha_pial);
write_wobj(obj,fn_obj);



% Plot surface
p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
    'EdgeColor','none','FaceColor',col_pial,'FaceAlpha',alpha_pial); 
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



fprintf('[!] Done.\n')