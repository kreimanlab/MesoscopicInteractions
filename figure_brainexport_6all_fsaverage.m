close all;
clear;

% Constants
flag_export_surf = false;
flag_debug = false;
SURFACE_TYPE = 'pial';
SPHERE_RADIUS = 100; % mm
SPHERE_CENTER = [0 0 0]; % mm
system('mkdir brainexport');
subjects_dir = '/media/jerry/internal/data/coreg';
dir_h5 = '/media/jerry/KLAB101/h5_notch20';

%subjects_dir = '/nas_share/cuenap_ssd/coregistration';
%dir_h5 = '/nas_share/cuenap/data/h5_notch20';
%subjects_dir = '/mnt/cuenap_ssd/coregistration';
%dir_h5 = '/mnt/cuenap/data/h5_notch20';
%Subjects = fliplr(Subjects);

Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
   'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
   'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',... % 'sub23',
   'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',... % 'sub30',
   'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
   'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};


[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,'fsaverage_sym','l','pial'));
center = mean(s_vert,1);
f = faces + 1;
P0 =  s_vert(f(:,1),:);
P1 =  s_vert(f(:,2),:);
P2 =  s_vert(f(:,3),:);

% What to plot
for sid_i = 4:length(Subjects)

    %sid_i = randi([1 length(Subjects)]);
    sid = Subjects{sid_i};
    %sid = 'sub3';
    fprintf('[%s] %i of %i\n',sid,sid_i,length(Subjects));
    

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
    
    % Read ielvis electrodes mapped onto average
    l = read_label(sprintf('%s/%s',subjects_dir,'fsaverage_sym'),sprintf('ielvis_%s',sid));
    if (isempty(l))
        l = read_label(sprintf('%s','fsaverage_sym'),sprintf('ielvis_%s',sid));
    end
    [~,sIdx] = sort(l(:,end));
    l = l(sIdx,:);
    
    n_counts = nchoosek(ecog.n_bchan,2);
    count = 1;
    Count = zeros(n_counts,3);
    for ii = 1:(ecog.n_bchan-1)
        for jj = (ii+1):ecog.n_bchan
            Count(count,1) = ii;
            Count(count,2) = jj;
            Count(count,3) = count;
            count = count + 1;
        end
    end
    
    %count = 1;
    %Paths = cell(1,n_counts);
    Paths = {};
    Paths_ind = [];
    
    % parfor
    parfor kk = 1:n_counts
%     for ii = 1:(ecog.n_bchan-1)
%         for jj = (ii+1):ecog.n_bchan

        ii = Count(kk,1);
        jj = Count(kk,2);
        count = Count(kk,3);
        
        
        fprintf('\t%i of %i\n',count,n_counts);
        tic;
        try
            bchan1 = ii;
            bchan2 = jj;
            %bchan1 = randi([1 ecog.n_bchan]); %35;
            %bchan2 = randi([1 ecog.n_bchan]); %84;

            b1c1 = ecog.bip(bchan1,1);
            b1c2 = ecog.bip(bchan1,2);
            b2c1 = ecog.bip(bchan2,1);
            b2c2 = ecog.bip(bchan2,2);

            if (b1c1 == b2c1)
                fprintf(2,'[W:] electrodes are the same: %i.\n',b1c1);
            end

            % get and check hemispheres
            hemi_11 = elec_hemis(b1c1);
            hemi_12 = elec_hemis(b1c2);
            hemi_21 = elec_hemis(b2c1);
            hemi_22 = elec_hemis(b2c2);
            hemi = lower(hemi_11{1});
            assert_hemis = all(strcmp({hemi_11{1},hemi_12{1},hemi_21{1}},hemi_22{1}));
            if (~assert_hemis)
                fprintf(2,'[W] Not all hemispheres are the same for subject: %s, %i, %i\n',sid,bchan1,bchan2)
            end

            % Store = struct;
            % Store.grad = {};
            % Store.grad_center = {};
            % Store.grad_umbrella = {};

            % Find vertex path between electrodes
            % load surfs
            
            %[s_vert_sphere, faces_sphere] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,'pial'));
            %center_sphere = mean(s_vert_sphere,1);
            

            delta = 2.5; % (mm)
            probe = [2, 5, 10, 15, 20]; % (mm)
            phis = linspace(1,170,32); % degrees
            %max_iter = 2*32;
            max_iter = 40;
            d_exit_thresh = 5; % (mm)

            % get vertex numbers for electrodes
            coord_b1c1 = s_vert(l(l(:,end) == b1c1,1)+1,:);
            coord_b2c1 = s_vert(l(l(:,end) == b2c1,1)+1,:);
            %coord_b1c1 = l(l(:,end) == b1c1,2:4);
            %coord_b2c1 = l(l(:,end) == b2c1,2:4);
            
            %coord_b1c1 = ecog.bip(bchan1,4:6);
            %coord_b2c1 = ecog.bip(bchan2,4:6);
            d_b1c1 = sqrt(sum((s_vert - coord_b1c1).^2,2));
            d_b2c1 = sqrt(sum((s_vert - coord_b2c1).^2,2));
            [~,v_index_b1c1] = min(d_b1c1);
            [~,v_index_b2c1] = min(d_b2c1);
            pt1 = s_vert(v_index_b1c1,:); % pt1 is the start
            pt2 = s_vert(v_index_b2c1,:); % pt2 is the goal
            pt1_sav = pt1;
            pt2_sav = pt2;
            
%             
%             path = [pt1;pt2];
% 
%             for j = 1:max_iter
%                 if (flag_debug)
%                     fprintf('%s: %i of %i\n',sid,j,max_iter);
%                 end
% 
%                 pta = pt1 - center;
%                 ptb = pt2 - center;
% 
%                 Q = cross(pta/norm(pta),ptb/norm(ptb));
%                 R = cross(Q/norm(Q),pta/norm(pta)); % orthogonal to pta
%                 R2 = cross(ptb/norm(ptb),Q/norm(Q)); % orthogonal to ptb
%                 probe_c_f = [];
%                 probe_c2_f = [];
%                 sat_1 = false;
%                 sat_2 = false;
%                 probe_c_before = NaN;
%                 probe_c2_before = NaN;
%                 c_before = NaN;
%                 c2_before = NaN;
%                 c_f = NaN;
%                 c2_f = NaN;
%                 for i = 1:length(phis)
%                     if (flag_debug)
%                         fprintf('\t%i of %i\n',i,length(phis));
%                     end
% 
%                     phi = (phis(i))*(pi/180);
%                     C = (pta/norm(pta)) * cos(phi) + R/(norm(R)) * sin(phi);
%                     C = C / norm(C);
%                     probe_c = (C') * probe + pt1';
%                     if (~sat_1)
%                         sat_1 = any(inpolyhedron(faces+1,s_vert,probe_c'));
%                     end
%                     if (sat_1 && (i > 1) && isempty(probe_c_f))
%                         probe_c_f = probe_c_before;
%                         c_f = c_before;
%                     end
% 
%                     phi2 = (phis(i))*(pi/180);
%                     C2 = (ptb/norm(ptb)) * cos(phi2) + R2/(norm(R2)) * sin(phi2);
%                     C2 = C2 / norm(C2);
%                     probe_c2 = (C2') * probe + pt2';
%                     if (~sat_2)
%                         sat_2 = any(inpolyhedron(faces+1,s_vert,probe_c2'));
%                     end
%                     if (sat_2 && (i > 1) && isempty(probe_c2_f))
%                         probe_c2_f = probe_c2_before;
%                         c2_f = c2_before;
%                     end
% 
% 
%                     %quiver3(pt1(1),pt1(2),pt1(3),C(1),C(2),C(3),10,'blue');
%                     %quiver3(pt2(1),pt2(2),pt2(3),C2(1),C2(2),C2(3),10,'green');
% 
%                     if (sat_1 && sat_2)
%                         break;
%                     end
% 
%                     probe_c_before = probe_c;
%                     probe_c2_before = probe_c2;
%                     c_before = C;
%                     c2_before = C2;
% 
%                 end
% 
%                 if (flag_debug && false)
%                     hold all;
%                     plot3(probe_c_f(1,:),probe_c_f(2,:),probe_c_f(3,:),'black.')
%                     plot3(probe_c2_f(1,:),probe_c2_f(2,:),probe_c2_f(3,:),'black.')
%                 end
% 
%                 % update
%                 pt1 = pt1 + c_f * delta;
%                 pt2 = pt2 + c2_f * delta;
% 
% 
% 
%                 % check for break
%                 d_between = norm(pt2 - pt1);
%                 cond_exit = d_between < d_exit_thresh;
%                 if (cond_exit)
%                     break
%                 end
% 
%                 % squeeze in path
%                 [n_path,~] = size(path);
%                 path = [path((1:j),:); pt1; pt2; path(((end-(j-1)):end),:)];
%             end




            path = [pt1];
            alphas = linspace(0,1,max_iter+2);
            alphas = alphas(2:(end-1));
            for j = 1:max_iter

                alpha = alphas(j);
                %phi = (phis(j))*(pi/180);

                pta = pt1 - center;
                ptb = pt2 - center;

                D = ptb * alpha + pta * (1 - alpha); % ray
                Dn = (-1) * D / norm(D);
                or = center + 120*Dn; % position ray 120 mm outside of center
                gD = gpuArray(Dn);

                [dist, flag] = arrayfun(@rayTriGPU, P0(:,1)', P0(:,2)', P0(:,3)', ...
                                    P1(:,1)', P1(:,2)', P1(:,3)', ...
                                    P2(:,1)', P2(:,2)', P2(:,3)', ...
                                    or(:,1), or(:,2), or(:,3), ...
                                    gD(:,1),gD(:,2),gD(:,3));

                distances = gather(dist);
                [~,mIdx] = min(distances);
                intpt = mean(s_vert(f(mIdx,:),:),1);

                %quiver3(or(1),or(2),or(3),Dn(1),Dn(2),Dn(3),5,'blue'); hold on;
                %plot3(intpt(:,1),intpt(:,2),intpt(:,3),'green.'); hold on;

                path = [path; intpt];
            end
            path = [path; pt2];


            %Paths{count} = path;
            Paths = [Paths, {path}];
            Paths_ind = [Paths_ind, count];

            if (flag_debug)
                % debug

                %quiver3(center(1),center(2),center(3),pta(1),pta(2),pta(3),1,'black');
                %quiver3(center(1),center(2),center(3),ptb(1),ptb(2),ptb(3),1,'color',0.33*[1 1 1]);

                %quiver3(center(1),center(2),center(3),R(1),R(2),R(3),norm(pta),'red');
                %quiver3(center(1),center(2),center(3),R2(1),R2(2),R2(3),norm(ptb),'color',[1 0 1]);

                p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
                    'EdgeColor','none','FaceColor',0.9*[1 1 1],'FaceAlpha',0.5); 
                hold all;
                pt1 = pt1_sav;
                pt2 = pt2_sav;
                plot3(pt1(1),pt1(2),pt1(3),'black.'); hold on;
                plot3(pt2(1),pt2(2),pt2(3),'black.'); hold on;
                plot3(path(:,1),path(:,2),path(:,3),'red-'); hold on;
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
                if (strcmp(hemi,'r'))
                    view(90,0);
                else
                    view(-90,0);
                end
            end
        
        catch e
            rethrow(e)
            fprintf(2,'E: skipped pair.\n');
        end

        t = toc;
        fprintf('\t\tTotal time: %.2f hrs.\n',(n_counts)*t/3600)
        %count = count + 1;
%         end
%     end
    end
    
    save_name = sprintf('brainexport/%s_6all_fsaverage',sid);
    %save(save_name,'-v7.3','-nocompression');
    save(save_name,'-v7.3');

end

% 
% % find path
% %delta = 5; % step size (mm)
% %delta_center = 4; % poke above the surface (mm)
% break_mm_thresh = 4; % halt condition
% collision_thresh = 2;
% collision_probe = 6; % (mm)
% collision_interval = 3; % (mm)
% hitbox_ignore = 1; % this number times interval (mm) will be ignored
% max_iter = 80;
% max_iter_umbrella = 10;
% stall_iter = 1; % stall threshold, number of iterations
% stall_fac = 0.1;% stall delta multiplier
% 
% umbrella_start = 0.2; % smaller values start shooting higher towards the "sky"
% umbrella_stop = 0.8; % larger values allow each step to point to goal more
% 
% pt1 = s_vert_sphere(v_index_b1c1,:); % pt1 is the start
% pt2 = s_vert_sphere(v_index_b2c1,:); % pt2 is the goal
% path = [v_index_b1c1];
% path_coords = [pt1];
% 
% delta = norm(pt2-pt1)/20; % step size (mm)
% 
% v_path = v_index_b1c1;
% pt1_sav = pt1;
% max_iter_reached = false;
% hitbox = linspace((1+hitbox_ignore),collision_probe,round((collision_probe-hitbox_ignore)/(collision_interval)));
% TR = triangulation(faces_sphere+1,s_vert_sphere);
% F = faceNormal(TR);
% P = incenter(TR);
% for i = 1:max_iter
%     if (flag_debug)
%         fprintf('iter: %i / %i\n',i,max_iter)
%     end
%     
%     % check for stalls
%     if ((i > stall_iter) && (path(i) == path(i-stall_iter)))
%         delta_sav = delta;
%         delta = delta * (1 + stall_fac);
%         collision_thresh_sav = collision_thresh;
%         collision_thresh = collision_thresh * (1 - stall_fac);
%         if (flag_debug)
%             fprintf('\tdelta: %.2f -> %.2f\n',delta_sav,delta);
%             fprintf('\tcollision_thresh: %.2f -> %.2f\n',collision_thresh_sav,collision_thresh);
%         end
%     end
%     
%     theta_umbrella = linspace(umbrella_start,umbrella_stop,max_iter_umbrella);
%     
%     % tangent plane
%     grad = pt2 - pt1; % get initial gradient
%     grad = grad / (sqrt(sum(grad.^2))); % get unit gradient
%     %grad = grad * delta; % scale by step delta
%     
%     % get norms from all faces that have point pt1
%     grad_norm = mean(F((any((faces_sphere + 1) == v_path,2)),:));
%     grad_norm = grad_norm / (sqrt(sum(grad_norm.^2)));
%     
%     grad_center = pt1 - center_sphere;
%     grad_center = grad_center / (sqrt(sum(grad_center.^2))); % scale to unit vector
%     
%     % average norm and center vectors
%     grad_center = mean([grad_norm; grad_center]);
%     %grad_center = grad_norm;
%     grad_center = grad_center / (sqrt(sum(grad_center.^2))); 
%     %return
%     
%     
%     % debug plots
%     if (flag_debug)
%         xp = grad*7;
%         xx = pt1;
%         quiver3(xx(1),xx(2),xx(3),xp(1),xp(2),xp(3),0.5,'blue'); hold on; % xx(1),xx(2),xx(3)
%         xp = grad_center*5;
%         xx = pt1;
%         quiver3(xx(1),xx(2),xx(3),xp(1),xp(2),xp(3),0.5,'color',[1 0 1]); hold on; % xx(1),xx(2),xx(3)
%     end
%     Store.grad{i} = grad;
%     Store.grad_center{i} = grad_center;
%     
%     
%     for j = 1:max_iter_umbrella
%         % start with pure normal/center vector
%         grad_umbrella = (grad*theta_umbrella(j) + grad_center*(1-theta_umbrella(j)));
%         grad_umbrella = grad_umbrella / (sqrt(sum(grad_umbrella.^2))); % scale to unit vector
%         probe = (pt1') + (grad_umbrella') * hitbox;
%         d_probe = zeros(1,length(hitbox));
%         for k = 1:length(hitbox)
%             d_probe(k) = min(sqrt(sum((s_vert_sphere - (probe(:,k)')).^2,2)));
%         end
%         cond_hit = any(d_probe < collision_thresh);
%         
%         Store.grad_umbrella{i,j} = grad_umbrella;
%         if (flag_debug)
%             xp = grad_umbrella*3;
%             xx = pt1;
%             quiver3(xx(1),xx(2),xx(3),xp(1),xp(2),xp(3),0.5,'color',[0 1 0]); hold on; % xx(1),xx(2),xx(3)
%             plot3(probe(1,:),probe(2,:),probe(3,:),'black.'); hold on;
%         end
%         
%         if (cond_hit)
%             % On hit, recompute gradient, if necessary
%             if (j > 1)
%                 grad_umbrella = (grad*theta_umbrella(j-1) + grad_center*(1-theta_umbrella(j-1)));
%                 grad_umbrella = grad_umbrella / (sqrt(sum(grad_umbrella.^2))); % scale to unit vector
%             end
%             break;
%         end
%         
%     end
%     if (flag_debug)
%         fprintf('umbrella: %i / %i\n',j,max_iter_umbrella)
%     end
%     pt_new = pt1 + grad_umbrella*delta;
%     d_grad = (sum((s_vert_sphere - pt_new).^2,2));
%     [~,v_path] = min(d_grad);
%     pt1 = s_vert_sphere(v_path,:);
%     
%     
%     % adaptively change umbrella params
%     if (j == max_iter_umbrella)
%         umbrella_stop_sav = umbrella_stop;
%         umbrella_stop = umbrella_stop * (1 + stall_fac);
%         if (flag_debug)
%             fprintf('\tumbrella_stop: %.2f -> %.2f\n',umbrella_stop_sav,umbrella_stop);
%         end
%     end
%     if (j == 1)
%         umbrella_start_sav = umbrella_start;
%         umbrella_start = umbrella_start * (1 - stall_fac);
%         if (flag_debug)
%             fprintf('\tumbrella_start: %.2f -> %.2f\n',umbrella_start_sav,umbrella_start);
%         end
%     end
%     
%     %return
% %     grad = pt2 - pt1; % get initial gradient
% %     grad = grad / (sqrt(sum(grad.^2))); % get unit gradient
% %     grad = grad * delta; % scale by step delta
% %     pt_new = pt1 + grad; % move towards pt2
% %     % map point to nearest on spherical surface
% %     d_grad = (sum((s_vert_sphere - pt_new).^2,2));
% %     [~,v_path] = min(d_grad);
% %     pt1 = s_vert_sphere(v_path,:); % update position
% %     % project outside from the center
% %     grad_center = pt1 - center_sphere;
% %     grad_center = grad_center / (sqrt(sum(grad_center.^2))); % scale to unit vector
% %     pt_new = pt1 + grad_center * delta_center; % poke outside
% %     d_grad = (sum((s_vert_sphere - pt_new).^2,2));
% %     [~,v_path] = min(d_grad);
% %     pt1 = s_vert_sphere(v_path,:); % update position
%     
%     path = [path; v_path]; % save path
%     path_coords = [path_coords; pt1]; % save coordinates
%     
%     % check new position for break condition
%     break_cond = sqrt(sum((pt2 - pt1).^2,2)) < break_mm_thresh;
%     if (break_cond)
%         max_iter_reached = true;
%         break;
%     end
% end
% % if (max_iter_reached == false)
% %     path = [v_index_b1c1; v_index_b2c1];
% %     path_coords = [pt1; pt2];
% % end
% 
% % Close the loop on the path
% path = [path; v_index_b2c1];
% path_coords = [path_coords; pt2];
% 
% pt1 = pt1_sav;
% 
% 
% Store.path = path;
% Store.path_coords = path_coords;
% Store.pt1 = pt1;
% Store.pt2 = pt2;
% Store.faces_sphere = faces_sphere;
% Store.s_vert_sphere = s_vert_sphere;
% 
% if (flag_debug)
%     % debug
%     p = trisurf(faces_sphere + 1,s_vert_sphere(:,1),s_vert_sphere(:,2),s_vert_sphere(:,3),'EdgeColor','none','FaceColor',0.8*[1 1 1]); hold all;
%     plot3(pt1(1),pt1(2),pt1(3),'black.'); hold on;
%     plot3(pt2(1),pt2(2),pt2(3),'black.'); hold on;
%     plot3(path_coords(:,1),path_coords(:,2),path_coords(:,3),'red-'); hold on;
%     daspect([1 1 1]);
%     p.AmbientStrength = 0.3 ;
%     p.DiffuseStrength = 0.4 ;
%     p.SpecularStrength = 0;
%     p.SpecularExponent = 1;
%     p.BackFaceLighting = 'lit';
%     cam_elev = 0;
%     camlight(-135,cam_elev);
%     camlight(45,cam_elev);
%     camlight(-225,cam_elev);
%     camlight(-45,cam_elev);
%     if (strcmp(hemi,'r'))
%         view(90,0);
%     else
%         view(-90,0);
%     end
% end
% 
% 
% 
% 
% 
% 
% % Export surface at end, depends on hemi and surface files from above
% % hemi = 'r'; % hemisphere override
% if (flag_export_surf)
%     %BRAIN_MESH_ALPHA = 1;
%     %p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),'EdgeColor','none','facealpha',BRAIN_MESH_ALPHA);
% 
%     % Save object
%     fn_obj = sprintf('brainexport/%s_%s.obj',sid,SURFACE_TYPE);
%     fprintf('[*] Saving brain to: %s\n',fn_obj);
%     clear obj;
%     obj.vertices = s_vert;
%     obj.objects(1).type='f';
%     obj.objects(1).data.vertices = faces + 1;
%     write_wobj(obj,fn_obj);
% end
% 

fprintf('[!] Done.\n')
