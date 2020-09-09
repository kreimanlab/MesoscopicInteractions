close all;
clear;

dir_cache = './cache';

Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
   'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
   'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',... % 'sub23',
   'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',... % 'sub30',
   'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
   'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};


sid_const = 'fsaverage_sym';
path_ds_factor = 1;
scale_mm_per_obj = 30;
flag_export_obj = true;
flag_plot = false;

mkdir('brainexport');
mkdir(sprintf('brainexport/%s',sid_const));


% Calculate maximum and minimum colors
fn_colorca = 'cache/figure_brainexport_static_fsaverage.mat';

load('brainexport/red_6all_fsaverage.mat');
CaR = load(sprintf('%s/fig_cluster2_reduce.mat','cache'));
load('figure_brainexport_6all_fsaverage_red');


E(:,3) = E(:,3) - center(1);
E(:,4) = E(:,4) - center(2);
E(:,5) = E(:,5) - center(3);


% Remove electrodes in medial wall
A = CaR.A;
A(rem,:) = [];
A(:,rem) = [];

Ad = CaR.Ad;
Ad(rem,:) = [];
Ad(:,rem) = [];
n_A = length(A);

% colormap function
n_col = 2^8;
cmap = corrcmap(n_col);
mag_min = nanmin(A(:));
mag_max = nanmax(A(:));
mag2col = @(x) cmap(round(((x - mag_min)/(mag_max - mag_min))*(n_col-1) + 1),:);


tube_radius = 0.15;
tube_n_edges = 3; % number of vertices of tube cross-section
rotm = rotationVectorToMatrix([pi/2, 0, 0]);
n = 1;
Paths = Paths(Paths_ind);
NamesPair = {};
for i = 1:(n_A - 1)
    for j = (i+1):n_A
        mag = A(i,j);
        if (~isnan(mag))
            
            path = Paths{n};
            
            path = downsample(path,path_ds_factor);
            
            path(:,1) = path(:,1) - center(1);
            path(:,2) = path(:,2) - center(2);
            path(:,3) = path(:,3) - center(3);
            
            
            [n_path,~] = size(path);

            col_tube = mag2col(mag);
            [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
            
            if (flag_plot)
                surf(X,Y,Z,'FaceColor',col_tube,'EdgeColor','none'); hold on;
            end
            
            pa = surf2patch(X,Y,Z,'triangles');

            bchan1 = i; %Count(Count(:,3) == n,1);
            bchan2 = j; %Count(Count(:,3) == n,2);

            if (flag_export_obj)
                clear obj;
                s_vert_obj = pa.vertices;
                s_vert_obj = (rotm * (s_vert_obj'))';
                %s_vert_obj(:,1) = pa.vertices(:,1);
                %s_vert_obj(:,2) = pa.vertices(:,3);
                %s_vert_obj(:,3) = pa.vertices(:,2);
                mtl_name = sprintf('pair_%i_%i_%i_mag-%i',n,bchan1,bchan2,round(mag*1000));
                obj.vertices = s_vert_obj / scale_mm_per_obj;
                obj.objects(1).type='usemtl'; % use material
                obj.objects(1).data=mtl_name;
                obj.objects(2).type='s'; % smooth shading
                obj.objects(2).data='off'; % 1 for on, "off" for off
%                 obj.objects(3).type='o';
%                 obj.objects(3).data=mtl_name;
%                 obj.objects(4).type='g';
%                 obj.objects(4).data='pair';
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
                fn_obj = sprintf('%s/brainexport/%s/%s_%s.obj',currdir,sid_const,sid_const,mtl_name);
                fprintf('[*] Saving pair to: %s\n',fn_obj);
                write_wobj(obj,fn_obj);
                NamesPair = [NamesPair; {mtl_name}];
            end
            
        end
        
        n = n + 1;
    end
end


hemi = 'r';
%return

% Read projected electrodes
% label is zero-indexed, checked by freeview
% l = read_label(sprintf('%s/%s',subjects_dir,'fsaverage_sym'),sprintf('ielvis_%s',sid));
% if (isempty(l))
%     l = read_label(sprintf('%s','fsaverage_sym'),sprintf('ielvis_%s',sid));
% end
% [~,sIdx] = sort(l(:,end));
% l = l(sIdx,:);

[n_locations,~] = size(E);

% Plot electrodes
n_elec_pts = 4;
r_elec = 2; % (mm)
col_elec = 0.5*[1 1 1];
col_elec_spec = [1 1 1];
%for i = 1:(ecog.n_bchan)
%%
NamesElec = {};
for i = 1:n_locations
    [x,y,z] = sphere(n_elec_pts);
    %coord = ecog.bip(i,4:6);
    %coord = l(l(:,end) == ecog.bip(i,1),2:4);
    %coord = s_vert(l(l(:,end) == ecog.bip(i,1),1)+1,:);
    coord = E(i,3:5);
    x = x + coord(1);
    y = y + coord(2);
    z = z + coord(3);

    if (flag_plot)
        mesh(x,y,z,'EdgeColor','none','FaceColor',col_elec); hold on;
    end

    pa = surf2patch(x,y,z,'triangles');

    % List subject membership for location
    es = CaR.Es{i};
    usubs = sort(unique(es(:,1)));
    mstr = 'sub-';
    for iml = 1:length(usubs)
        if (iml ~= 1)
            mstr = [mstr,'-'];
        end
        mstr = [mstr,sprintf('%i',usubs(iml))];
    end

    % Export to obj file
    if (flag_export_obj)
        clear obj;
        s_vert_obj = pa.vertices;
        s_vert_obj = (rotm * (s_vert_obj'))';
        %s_vert_obj(:,1) = pa.vertices(:,1);
        %s_vert_obj(:,2) = pa.vertices(:,3);
        %s_vert_obj(:,3) = pa.vertices(:,2);
        %mtl_name = sprintf('electrode%i',ecog.bip(i,1));
        mtl_name = sprintf('elec-%i_%s',i,mstr);
        obj.vertices = s_vert_obj / scale_mm_per_obj;
        obj.objects(1).type='usemtl'; % use material
        obj.objects(1).data=mtl_name;
        obj.objects(2).type='s'; % smooth shading
        obj.objects(2).data='off'; % 1 for on, "off" for off
%         obj.objects(3).type='o';
%         obj.objects(3).data=mtl_name;
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
        fn_obj = sprintf('%s/brainexport/%s/%s_%s.obj',currdir,sid_const,sid_const,mtl_name);
        fprintf('[*] Saving electrode to: %s\n',fn_obj);
        write_wobj(obj,fn_obj);
        NamesElec = [NamesElec; {mtl_name}];
    end
end

% Export pair and elec names to json
%currdir = pwd;
%addpath('jsonlab-1.5');
%savejson('NamesElec',NamesElec,sprintf('%s/brainexport/%s/NamesElec.json',currdir,sid_const));
%savejson('NamesPair',NamesPair,sprintf('%s/brainexport/%s/NamesPair.json',currdir,sid_const));

col_pial = 0.6*[1 1 1];
col_pial_spec = [1 1 1];
alpha_pial = 0.8;
subjects_dir = '/media/jerry/internal/data/coreg';
SURFACE_TYPE = 'pial';

% Save surface object
currdir = pwd;
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,SURFACE_TYPE));
s_vert(:,1) = s_vert(:,1) - center(1);
s_vert(:,2) = s_vert(:,2) - center(2);
s_vert(:,3) = s_vert(:,3) - center(3);

if (flag_export_obj)
    fn_obj = sprintf('%s/brainexport/%s/%s_%s.obj',currdir,sid_const,sid_const,SURFACE_TYPE);
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
    obj.material(5).data=9;
    obj.material(6).type='Ns'; % 0 to 1000, specular exponent
    obj.material(6).data=5;
    obj.material(7).type='d';
    obj.material(7).data=alpha_pial;
    obj.material(7).type='Tr';
    obj.material(7).data=(1-alpha_pial);
    write_wobj(obj,fn_obj);
end



% Plot surface
if (flag_plot)
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
    
    view(90,0);
    axis off;
    set(gcf,'Position',[0 0 1920 1080]);
    colormap(cmap);
    cb = colorbar;
    set(cb,'Location','east');
    set(cb,'TickLength',0);
    set(cb,'Ticks',[mag_min mag_max]);
    caxis([mag_min mag_max]);
    print(gcf,'figures/brainexport_static_fsaverage_red','-depsc');
    print(gcf,'figures/brainexport_static_fsaverage_red','-dpng','-r300');
end


fprintf('[!] Done.\n')

