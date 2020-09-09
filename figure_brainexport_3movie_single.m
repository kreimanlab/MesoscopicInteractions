close all;
clear;


b1 = 35;
b2 = 84;


% Constants
flag_export_surf = false;
flag_debug = true;
SURFACE_TYPE = 'pial';
SPHERE_RADIUS = 100; % mm
SPHERE_CENTER = [0 0 0]; % mm
system('mkdir brainexport');
subjects_dir = '/media/jerry/internal/data/coreg';
dir_h5 = '/media/jerry/KLAB101/h5_notch20';
dir_art = '/media/jerry/KLAB101/h5_notch20/art_nosz2';
dir_res = '/media/jerry/KLAB101/results/coh_w10';
dir_cache = './cache';


Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
   'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
   'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',... % 'sub23',
   'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',... % 'sub30',
   'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
   'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

% What to plot
%sid_i = randi([1 length(Subjects)]);
%sid = Subjects{sid_i};

%for sid_i = 1:length(Subjects)

%sid = 'sub3';
sid_i = 3;
sid = Subjects{sid_i};
dir_mov = sprintf('brainexport/%s_movie_single',sid);
system(sprintf('mkdir %s',dir_mov));

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
fprintf('[*] Loading brainexport 1 cache..\n')
fn_be1 = sprintf('brainexport/%s_6all.mat',sid);
if (~exist('Cb','var'))
    Cb = load(fn_be1);
    Paths = Cb.Paths;
    Paths_ind = Cb.Paths_ind;
    Paths = Paths(Paths_ind);
end
fn_cache = sprintf('%s/xsub_out_%s_%i',dir_cache,sid,1);
fprintf('[*] Loading xsub_out_%s cache..\n',sid)
Ca = load(fn_cache);

% load graph
fprintf('[*] Loading coherence timeseries..\n')
fn_graph = sprintf('%s/%s_graph-pcBroadband.h5',dir_res,sid);
R = h5read(fn_graph,'/R');
[~,n_graph] = size(R);

% load artifacts
fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
art_idx = h5read(fn_art,'/art_idx') > 0;

% load surface
hemi = Cb.hemi;
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,SURFACE_TYPE));

%%

% dist thresh
pass_dist = (Ca.Dmats > Ca.dist_thresh);
pass_ctsig = (Ca.ct > Ca.ct_thresh);
mag_all = Ca.mag(pass_dist & pass_ctsig);

% normalize color
n_col = 200;
cmap = corrcmap(n_col);
%cmap = jet(n_col);
R_sig = R;
R_sig(art_idx) = NaN;
R_sig2 = R_sig(:,~all(art_idx));
R_sig = R_sig(pass_dist & pass_ctsig,~all(art_idx));
[~,n_graph_sig] = size(R_sig);

%mag_max = double(max(R_sig(:))); 
mag_max = double(max(R_sig(:))); %max(mag_all);
mag_min = double(min(R_sig(:))); %min(mag_all);
mag2col = @(x) cmap(round(((x - mag_min)/(mag_max - mag_min))*(n_col-1) + 1),:);

tube_n_edges = 9;
tube_radius = 1.5;
len_series = 8;
for k = 1:n_graph_sig

    h = figure('Units','Normalized','Position',[0 0 1 1],'visible','off');
    %[ha, pos] = tight_subplot(1,1,[.01 .01],[.01 .01],[.01 .01]);
    %axes(ha(1));
    
    count = 1;
    for ii = 1:(ecog.n_bchan-1)
        for jj = (ii+1):ecog.n_bchan
            
            cond_fwd = ( (ii == b1) && (jj == b2) );
            cond_rev = ( (ii == b2) && (jj == b1) );
            % Check significance
            %if (pass_dist(count) && pass_ctsig(count))
            if ((cond_fwd || cond_rev) && (pass_dist(count) && pass_ctsig(count)))
                %Store = Cb.Stores{count};
                path = Paths{count};

                % get significance info
                %mag = Ca.mag(count);
                mag = R_sig2(count,k);
                
                coh_thresh = Ca.coh_thresh(count);
                ct = Ca.ct(count);
                ct_thresh = Ca.ct_thresh;

                %return
                if (flag_debug && (~isnan(mag)))
                    %path_coords = Store.path_coords;
                    %path_coords = path;
                    if (mag > mag_max)
                        mag = mag_max;
                    end
                    
                    [n_path,~] = size(path);
                    [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
                    fv = surf2patch(X,Y,Z);
                    col_path = mag2col(Ca.mag(count));
                    p2 = trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),...
                        'EdgeColor','none','FaceColor',col_path); hold on;
                    p2.SpecularStrength = 0;
                    p2.SpecularExponent = 1;
                    
                    
                    %plot3(path_coords(:,1),path_coords(:,2),path_coords(:,3),'-','Color',mag2col(mag)); hold on;
                end
%                 if (mag > 0.3)
%                 fprintf('mag: %.3f\n',mag);
%                 end

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

        %pt1 = Cb.pt1;
        %pt2 = Cb.pt2;
        % debug
        p = trisurf(faces + 1,s_vert(:,1),s_vert(:,2),s_vert(:,3),'EdgeColor','none','FaceColor',0.8*[1 1 1]); hold all;
        daspect([1 1 1]);
        %plot3(pt1(1),pt1(2),pt1(3),'black.'); hold on;
        %plot3(pt2(1),pt2(2),pt2(3),'black.'); hold on;
        %plot3(path_coords(:,1),path_coords(:,2),path_coords(:,3),'red-'); hold on;

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
        if (strcmp(hemi,'r'))
            view(90,0);
        else
            view(-90,0);
        end
        % Just for sub3
        view(130,20);
        axis off;
        colormap(cmap);
        cb = colorbar;
        set(cb,'TickLength',0);
        set(cb,'Location','east');
        col_ticks = [mag_min mag_max];
        col_tickl = {sprintf('%.2f',mag_min),sprintf('%.2f',mag_max)};
        set(cb,'Ticks',col_ticks)
        set(cb,'TickLabels',col_tickl)
        caxis([mag_min mag_max])
    end
    
    series = sprintf('%i',k);
    while (length(series) < len_series)
        series = ['0',series];
    end
    
    print(h,sprintf('%s/%s',dir_mov,series),'-djpeg')
    close(h);
end


fprintf('[!] Done.\n')
