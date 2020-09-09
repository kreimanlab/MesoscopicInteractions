close all;
clear;

addpath('ndcrdspinterp');
load('brainexport/red_6all_fsaverage.mat');

ang_thresh = 90 + 45; % degrees
dist_thresh = 1; % millimeters
Tension = 0;
n = 1;
cmap = corrcmap(100);
alpha_pial = 0.3;
col_pial = 0.7*[1 1 1];
SURFACE_TYPE = 'pial';
hemi = 'r';
sid_const = 'fsaverage_sym';
subjects_dir = '/media/jerry/internal/data/coreg';
CaR = load('cache/fig_cluster2_reduce.mat');

[v_sph, fa_sph] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,'sphere'));
[v_pia, fa_pia] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,'pial'));
vf_sph = triangulation(fa_sph+1,v_sph(:,1),v_sph(:,2),v_sph(:,3));
vf_pia = triangulation(fa_pia+1,v_pia(:,1),v_pia(:,2),v_pia(:,3));


plen = [];
for i = (13*150 + (1:5*150)) %1:length(Paths) % randi([1 length(Paths)])
    path = Paths{i};
    [n_path,~] = size(path);
    path_mark = false(n_path,1);
    Px = path(:,1);
    Py = path(:,2);
    Pz = path(:,3);
    %K = convhull(Px,Py,Pz);
    %F = scatteredInterpolant(Px,Py,Pz);
    %Px = [path(1,1); path(:,1); path(end,1)]';
    %Py = [path(1,2); path(:,2); path(end,2)]';
    %Pz = [path(1,3); path(:,3); path(end,3)]';
    for j = 1:(length(Px)-1)
        %P0 = path(j,:);
        %P1 = path(j+1,:);
        %P2 = path(j+2,:);
        %n1 = (P0 - P1) / norm(P0 - P1);
        %n2 = (P2 - P1) / norm(P2 - P1);
        %ang = atan2(norm(cross(n1, n2)), dot(n1, n2)) * (180/pi);
        %fprintf('angle: %.2f\n',ang)
        %path_mark(j+1) = (ang > ang_thresh);
        
        % P1 & P2 are endpoints of spline.
        % P0 & P3 are used to calculate the slope of P1 & P2
        %[MatNbyNPlusOne]=crdatnplusoneval(P0,P1,P2,P3,T,n);
        
        %[XiYiZi]=crdatnplusoneval([Px(k),Py(k),Pz(k)],[Px(k+1),Py(k+1),Pz(k+1)],[Px(k+2),Py(k+2),Pz(k+2)],[Px(k+3),Py(k+3),Pz(k+3)],Tension,n);
        %plot3(XiYiZi(1,:),XiYiZi(2,:),XiYiZi(3,:),'blue-','linewidth',2); hold on;
        
        d = norm(path(j+1,:) - path(j,:));
        path_mark(j) = (d > dist_thresh);
        
        
    end
    
    pt1 = path(1,:);
    pt2 = path(end,:);
    
    id_0 = nearestNeighbor(vf_pia,pt1(1),pt1(2),pt1(3));
    id_f = nearestNeighbor(vf_pia,pt2(1),pt2(2),pt2(3));
    
    xsph_0 = vf_sph.Points(id_0,:);
    xsph_f = vf_sph.Points(id_f,:);
    xsph_0 = 100 * xsph_0/norm(xsph_0);
    xsph_f = 100 * xsph_f/norm(xsph_f);
    %n_res = round((atan2(norm(cross(pt1,pt2)),dot(pt1,pt2)) * 100) / 3);
    %path2 = zeros(n_res,3);
    %path2(1,:) = pt1;
    %path2(end,:) = pt2;
    
    path2 = pt1;
    %alphas = linspace(0,1,n_res);
    res_mm = 4;
    d_rate = 0.00120;
    alpha = 0;
    dalpha = 1/round((atan2(norm(cross(pt1,pt2)),dot(pt1,pt2)) * 100) / res_mm);
    dalpha_sav = dalpha;
    %n_alpha_adjust = 5;
    while (true)
        % Catch out of bound alpha deltas
        if ((dalpha <= 0) || (dalpha > dalpha_sav))
            dalpha = dalpha_sav;
        end
        
        % update alpha
        alphan = alpha + dalpha;
        
        % exit condition
        if(alphan > 1)
            break;
        end
        
        % -----------------------------------------------------------------
        xsph = xsph_0 * (1 - alphan) + xsph_f * alphan;
        xsph = 100*xsph/norm(xsph);
        id = nearestNeighbor(vf_sph,xsph(1),xsph(2),xsph(3));
        pt = vf_pia.Points(id,:);
        %------------------------------------------------------------------
        
        
%         % adjust alpha delta
%         for iadj = 1:n_alpha_adjust
%             error = norm(pt - path2(end,:)) - res_mm;
%             fprintf('%i error: %.14f\n',iadj,error)
%             dalpha = dalpha - d_rate * error;
%             alphan = alpha + dalpha;
%             
%             % -----------------------------------------------------------------
%             xsph = xsph_0 * (1 - alphan) + xsph_f * alphan;
%             xsph = 100*xsph/norm(xsph);
%             id = nearestNeighbor(vf_sph,xsph(1),xsph(2),xsph(3));
%             pt = vf_pia.Points(id,:);
%             %------------------------------------------------------------------
%         end
        
        error = norm(pt - path2(end,:)) - res_mm;
        %fprintf('\tabs error: %.14f\n',abs(error))
        %fprintf('alpha: %.4f\n',alpha)
        dalpha = dalpha - d_rate * error;
        
        % store point
        path2 = [path2; pt];
        
        % update alpha
        alpha = alphan;
    end
    path2 = [path2; pt2];
    plen = [plen; length(path2)];
    
    
    %ptCloud = pointCloud(path);
    %ptCloudA = pcdownsample(ptCloud,'nonuniformGridSample',12);
    
    %plot3(path(:,1),path(:,2),path(:,3),'black-'); hold all;
    %path2 = ptCloudA.Location;
    plot3(path2(:,1),path2(:,2),path2(:,3),'black-'); hold all;
    plot3(path2(end,1),path2(end,2),path2(end,3),'red.'); hold all;
    %plot3(path(path_mark,1),path(path_mark,2),path(path_mark,3),'red-'); hold on;
    
    
    
    
    
    
    
    

    

%     if (i == 149)
%         break
%     end
    
end
fprintf('mean path len: %.4f\n',mean(plen))

daspect([1 1 1]);
return

[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,SURFACE_TYPE));
p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
    'EdgeColor','none','FaceColor',col_pial,'FaceAlpha',alpha_pial);
hold all;
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