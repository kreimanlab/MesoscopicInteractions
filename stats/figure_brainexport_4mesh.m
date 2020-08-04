close all;
clear;

% Constants
flag_export_surf = false;
flag_debug = true;
SURFACE_TYPE = 'pial';
SPHERE_RADIUS = 100; % mm
SPHERE_CENTER = [0 0 0]; % mm
system('mkdir brainexport');
subjects_dir = '/media/jerry/internal/data/coreg';
dir_h5 = '/media/jerry/KLAB101/h5_notch20';
dir_res = '/media/jerry/KLAB101/results/coh_w10';
dir_cache = './cache';

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
   'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
   'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',... % 'm00045',
   'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',... % 'm00056',
   'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
   'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

% What to plot
%sid_i = randi([1 length(Subjects)]);
%sid = Subjects{sid_i};

%for sid_i = 1:length(Subjects)

%sid = 'm00005';
sid_i = 2;
sid = Subjects{sid_i};


% Read .h5
fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
ecog = H5eeg(fn_h5);
% Read ielvis electrodeNames
fn_enames = sprintf('%s/%s/elec_recon/%s.electrodeNames',subjects_dir,sid,sid);
fid = fopen(fn_enames,'r');
elec_names = textscan(fid,'%s','HeaderLines',2);
elec_names = elec_names{1};
elec_names = reshape(elec_names,3,[])';
elec_hemis = elec_names(:,3);
% elec_labels = elec_names(:,1);
fclose(fid);

% load cache
fn_be1 = sprintf('brainexport/%s.mat',sid);
if (~exist('Cb','var'))
    Cb = load(fn_be1);
end
fn_cache = sprintf('%s/xsub_out_%s_%i',dir_cache,sid,1);
Ca = load(fn_cache);

% load surface
hemi = Cb.hemi;
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,SURFACE_TYPE));

% dist thresh
pass_dist = (Ca.Dmats > Ca.dist_thresh);
pass_ctsig = (Ca.ct > Ca.ct_thresh);
mag_all = Ca.mag(pass_dist & pass_ctsig);

% normalize color
n_col = 100;
cmap = corrcmap(n_col);
mag_max = max(mag_all);
mag_min = min(mag_all);
mag2col = @(x) cmap(round(((x - mag_min)/(mag_max - mag_min))*(n_col-1) + 1),:);

count = 1;
for ii = 1:(ecog.n_bchan-1)
    for jj = (ii+1):ecog.n_bchan
        % Check significance
        %if (pass_dist(count) && pass_ctsig(count))
        if (true)
            Store = Cb.Stores{count};

            % get significance info
            mag = Ca.mag(count);
            coh_thresh = Ca.coh_thresh(count);
            ct = Ca.ct(count);
            ct_thresh = Ca.ct_thresh;

            %return
            if (flag_debug)
                path_coords = Store.path_coords;
                %col_line = mag2col(mag);
                col_line = 0.3*[1 1 1];
                plot3(path_coords(:,1),path_coords(:,2),path_coords(:,3),'-','Color',col_line); hold on;
            end
        
        end
        
        count = count + 1;
        
    end
end

if (flag_debug)
    % Plot electrodes
    
    for i = 1:ecog.n_bchan
        bcoord = ecog.bip(i,4:6);
        plot3(bcoord(1),bcoord(2),bcoord(3),'.','MarkerSize',40,'Color',[0 0 0]);
    end
    
    pt1 = Cb.pt1;
    pt2 = Cb.pt2;
    % debug
    p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
        'EdgeColor','none','FaceColor',0.9*[1 1 1],'FaceAlpha',0.5); 
    hold all;
    daspect([1 1 1]);
    %plot3(pt1(1),pt1(2),pt1(3),'black.'); hold on;
    %plot3(pt2(1),pt2(2),pt2(3),'black.'); hold on;
    %plot3(path_coords(:,1),path_coords(:,2),path_coords(:,3),'red-'); hold on;
    
    % plot center
    s_vert_c = mean(s_vert);
    plot3(s_vert_c(1),s_vert_c(2),s_vert_c(3),'red.','MarkerSize',80);
    
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
    if (strcmp(hemi,'r'))
        view(90,0);
    else
        view(-90,0);
    end
end

%return
% 
% Stores = {};
% %count = 1;
% counts = 1:(ecog.n_bchan-1);
% for ii = 1:(ecog.n_bchan-1)
%     count = counts(ii);
%     fprintf('%i of %i\n',count,length(counts));
%         
%     for jj = (ii+1):ecog.n_bchan
%         %fprintf('\t%i of %i\n',jj,length(counts));
%         
%         bchan1 = ii;
%         bchan2 = jj;
%         %bchan1 = randi([1 ecog.n_bchan]); %35;
%         %bchan2 = randi([1 ecog.n_bchan]); %84;
% 
%         b1c1 = ecog.bip(bchan1,1);
%         b1c2 = ecog.bip(bchan1,2);
%         b2c1 = ecog.bip(bchan2,1);
%         b2c2 = ecog.bip(bchan2,2);
% 
%         % get and check hemispheres
%         hemi_11 = elec_hemis(b1c1);
%         hemi_12 = elec_hemis(b1c2);
%         hemi_21 = elec_hemis(b2c1);
%         hemi_22 = elec_hemis(b2c2);
%         hemi = lower(hemi_11{1});
%         assert_hemis = all(strcmp({hemi_11{1},hemi_12{1},hemi_21{1}},hemi_22{1}));
%         if (~assert_hemis)
%             fprintf(2,'[W] Not all hemispheres are the same for subject: %s, %i, %i\n',sid,bchan1,bchan2)
%         end
% 
% 
% 
%         Store = struct;
%         Store.grad = {};
%         Store.grad_center = {};
%         Store.grad_umbrella = {};
% 
%         % Find vertex path between electrodes
%         % load surfs
%         [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,SURFACE_TYPE));
%         [s_vert_sphere, faces_sphere] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,'pial'));
%         center_sphere = mean(s_vert_sphere,1);
% 
% 
%         % get vertex numbers for electrodes
%         coord_b1c1 = ecog.bip(bchan1,4:6);
%         coord_b2c1 = ecog.bip(bchan2,4:6);
%         d_b1c1 = sqrt(sum((s_vert - coord_b1c1).^2,2));
%         d_b2c1 = sqrt(sum((s_vert - coord_b2c1).^2,2));
%         [~,v_index_b1c1] = min(d_b1c1);
%         [~,v_index_b2c1] = min(d_b2c1);
% 
%         % find path
%         %delta = 5; % step size (mm)
%         %delta_center = 4; % poke above the surface (mm)
%         break_mm_thresh = 4; % halt condition
%         collision_thresh = 2;
%         collision_probe = 6; % (mm)
%         collision_interval = 3; % (mm)
%         hitbox_ignore = 1; % this number times interval (mm) will be ignored
%         max_iter = 80;
%         max_iter_umbrella = 10;
%         stall_iter = 1; % stall threshold, number of iterations
%         stall_fac = 0.1;% stall delta multiplier
% 
%         umbrella_start = 0.2; % smaller values start shooting higher towards the "sky"
%         umbrella_stop = 0.8; % larger values allow each step to point to goal more
% 
%         pt1 = s_vert_sphere(v_index_b1c1,:); % pt1 is the start
%         pt2 = s_vert_sphere(v_index_b2c1,:); % pt2 is the goal
%         path = [v_index_b1c1];
%         path_coords = [pt1];
% 
%         delta = norm(pt2-pt1)/20; % step size (mm)
% 
%         v_path = v_index_b1c1;
%         pt1_sav = pt1;
%         max_iter_reached = false;
%         hitbox = linspace((1+hitbox_ignore),collision_probe,round((collision_probe-hitbox_ignore)/(collision_interval)));
%         TR = triangulation(faces_sphere+1,s_vert_sphere);
%         F = faceNormal(TR);
%         P = incenter(TR);
%         for i = 1:max_iter
%             if (flag_debug)
%                 fprintf('iter: %i / %i\n',i,max_iter)
%             end
% 
%             % check for stalls
%             if ((i > stall_iter) && (path(i) == path(i-stall_iter)))
%                 delta_sav = delta;
%                 delta = delta * (1 + stall_fac);
%                 collision_thresh_sav = collision_thresh;
%                 collision_thresh = collision_thresh * (1 - stall_fac);
%                 if (flag_debug)
%                     fprintf('\tdelta: %.2f -> %.2f\n',delta_sav,delta);
%                     fprintf('\tcollision_thresh: %.2f -> %.2f\n',collision_thresh_sav,collision_thresh);
%                 end
%             end
% 
%             theta_umbrella = linspace(umbrella_start,umbrella_stop,max_iter_umbrella);
% 
%             % tangent plane
%             grad = pt2 - pt1; % get initial gradient
%             grad = grad / (sqrt(sum(grad.^2))); % get unit gradient
%             %grad = grad * delta; % scale by step delta
% 
%             % get norms from all faces that have point pt1
%             grad_norm = mean(F((any((faces_sphere + 1) == v_path,2)),:));
%             grad_norm = grad_norm / (sqrt(sum(grad_norm.^2)));
% 
%             grad_center = pt1 - center_sphere;
%             grad_center = grad_center / (sqrt(sum(grad_center.^2))); % scale to unit vector
% 
%             % average norm and center vectors
%             grad_center = mean([grad_norm; grad_center]);
%             %grad_center = grad_norm;
%             grad_center = grad_center / (sqrt(sum(grad_center.^2))); 
%             %return
% 
% 
%             % debug plots
%             if (flag_debug)
%                 xp = grad*7;
%                 xx = pt1;
%                 quiver3(xx(1),xx(2),xx(3),xp(1),xp(2),xp(3),0.5,'blue'); hold on; % xx(1),xx(2),xx(3)
%                 xp = grad_center*5;
%                 xx = pt1;
%                 quiver3(xx(1),xx(2),xx(3),xp(1),xp(2),xp(3),0.5,'color',[1 0 1]); hold on; % xx(1),xx(2),xx(3)
%             end
%             Store.grad{i} = grad;
%             Store.grad_center{i} = grad_center;
% 
% 
%             for j = 1:max_iter_umbrella
%                 % start with pure normal/center vector
%                 grad_umbrella = (grad*theta_umbrella(j) + grad_center*(1-theta_umbrella(j)));
%                 grad_umbrella = grad_umbrella / (sqrt(sum(grad_umbrella.^2))); % scale to unit vector
%                 probe = (pt1') + (grad_umbrella') * hitbox;
%                 d_probe = zeros(1,length(hitbox));
%                 for k = 1:length(hitbox)
%                     d_probe(k) = min(sqrt(sum((s_vert_sphere - (probe(:,k)')).^2,2)));
%                 end
%                 cond_hit = any(d_probe < collision_thresh);
% 
%                 Store.grad_umbrella{i,j} = grad_umbrella;
%                 if (flag_debug)
%                     xp = grad_umbrella*3;
%                     xx = pt1;
%                     quiver3(xx(1),xx(2),xx(3),xp(1),xp(2),xp(3),0.5,'color',[0 1 0]); hold on; % xx(1),xx(2),xx(3)
%                     plot3(probe(1,:),probe(2,:),probe(3,:),'black.'); hold on;
%                 end
% 
%                 if (cond_hit)
%                     % On hit, recompute gradient, if necessary
%                     if (j > 1)
%                         grad_umbrella = (grad*theta_umbrella(j-1) + grad_center*(1-theta_umbrella(j-1)));
%                         grad_umbrella = grad_umbrella / (sqrt(sum(grad_umbrella.^2))); % scale to unit vector
%                     end
%                     break;
%                 end
% 
%             end
%             if (flag_debug)
%                 fprintf('umbrella: %i / %i\n',j,max_iter_umbrella)
%             end
%             pt_new = pt1 + grad_umbrella*delta;
%             d_grad = (sum((s_vert_sphere - pt_new).^2,2));
%             [~,v_path] = min(d_grad);
%             pt1 = s_vert_sphere(v_path,:);
% 
% 
%             % adaptively change umbrella params
%             if (j == max_iter_umbrella)
%                 umbrella_stop_sav = umbrella_stop;
%                 umbrella_stop = umbrella_stop * (1 + stall_fac);
%                 if (flag_debug)
%                     fprintf('\tumbrella_stop: %.2f -> %.2f\n',umbrella_stop_sav,umbrella_stop);
%                 end
%             end
%             if (j == 1)
%                 umbrella_start_sav = umbrella_start;
%                 umbrella_start = umbrella_start * (1 - stall_fac);
%                 if (flag_debug)
%                     fprintf('\tumbrella_start: %.2f -> %.2f\n',umbrella_start_sav,umbrella_start);
%                 end
%             end
% 
%             %return
%         %     grad = pt2 - pt1; % get initial gradient
%         %     grad = grad / (sqrt(sum(grad.^2))); % get unit gradient
%         %     grad = grad * delta; % scale by step delta
%         %     pt_new = pt1 + grad; % move towards pt2
%         %     % map point to nearest on spherical surface
%         %     d_grad = (sum((s_vert_sphere - pt_new).^2,2));
%         %     [~,v_path] = min(d_grad);
%         %     pt1 = s_vert_sphere(v_path,:); % update position
%         %     % project outside from the center
%         %     grad_center = pt1 - center_sphere;
%         %     grad_center = grad_center / (sqrt(sum(grad_center.^2))); % scale to unit vector
%         %     pt_new = pt1 + grad_center * delta_center; % poke outside
%         %     d_grad = (sum((s_vert_sphere - pt_new).^2,2));
%         %     [~,v_path] = min(d_grad);
%         %     pt1 = s_vert_sphere(v_path,:); % update position
% 
%             path = [path; v_path]; % save path
%             path_coords = [path_coords; pt1]; % save coordinates
% 
%             % check new position for break condition
%             break_cond = sqrt(sum((pt2 - pt1).^2,2)) < break_mm_thresh;
%             if (break_cond)
%                 max_iter_reached = true;
%                 break;
%             end
%         end
%         % if (max_iter_reached == false)
%         %     path = [v_index_b1c1; v_index_b2c1];
%         %     path_coords = [pt1; pt2];
%         % end
% 
%         % Close the loop on the path
%         path = [path; v_index_b2c1];
%         path_coords = [path_coords; pt2];
% 
%         pt1 = pt1_sav;
% 
% 
%         Store.path = path;
%         Store.path_coords = path_coords;
%         Store.pt1 = pt1;
%         Store.pt2 = pt2;
%         %Store.faces_sphere = faces_sphere;
%         %Store.s_vert_sphere = s_vert_sphere;
% 
%         if (flag_debug)
%             % debug
%             p = trisurf(faces_sphere + 1,s_vert_sphere(:,1),s_vert_sphere(:,2),s_vert_sphere(:,3),'EdgeColor','none','FaceColor',0.8*[1 1 1]); hold all;
%             plot3(pt1(1),pt1(2),pt1(3),'black.'); hold on;
%             plot3(pt2(1),pt2(2),pt2(3),'black.'); hold on;
%             plot3(path_coords(:,1),path_coords(:,2),path_coords(:,3),'red-'); hold on;
%             daspect([1 1 1]);
%             p.AmbientStrength = 0.3 ;
%             p.DiffuseStrength = 0.4 ;
%             p.SpecularStrength = 0;
%             p.SpecularExponent = 1;
%             p.BackFaceLighting = 'lit';
%             cam_elev = 0;
%             camlight(-135,cam_elev);
%             camlight(45,cam_elev);
%             camlight(-225,cam_elev);
%             camlight(-45,cam_elev);
%             if (strcmp(hemi,'r'))
%                 view(90,0);
%             else
%                 view(-90,0);
%             end
%         end
% 
%         Stores = [Stores, {Store}];
%         
%         %count = count + 1;
%     end
% end
% 
% save(['brainexport/',sid])



% Export surface at end, depends on hemi and surface files from above
% hemi = 'r'; % hemisphere override
if (flag_export_surf)
    %BRAIN_MESH_ALPHA = 1;
    %p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);

    % Save object
    fn_obj = sprintf('brainexport/%s_%s.obj',sid,SURFACE_TYPE);
    fprintf('[*] Saving brain to: %s\n',fn_obj);
    clear obj;
    obj.vertices = s_vert;
    obj.objects(1).type='f';
    obj.objects(1).data.vertices = faces + 1;
    write_wobj(obj,fn_obj);
end

%end

fprintf('[!] Done.\n')
