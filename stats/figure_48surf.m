close all;
clear;

subjects_dir = '/media/jerry/internal/data/coreg';
system('mkdir figures/48surf');

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
   'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
   'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',... % 'm00045',
   'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',... % 'm00056',
   'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
   'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

const_col_surface = 0.90*[1 1 1];
const_alpha_surface = 1;
const_elec_surface = 0*[1 1 1];
const_elec_low = 0.4*[1 1 1];

h = figure('visible','off','Units','pixels','Position',[0 0 1920 1080]);
[ha, pos] = tight_subplot(6,8,[.01 .01],[.01 .01],[.01 .01]);

hold all;

n_chan_all = [];

for i = 1:length(Subjects)

    sid = Subjects{i};
    axes(ha(i));
    
    
    % Get hemi
    fn_en = sprintf('%s/%s/elec_recon/%s.electrodeNames',subjects_dir,sid,sid);
    if (~exist(fn_en,'file'))
        fprintf(2,'E: Did not find file: %s\n',fn_en);
    end
    f_en = fopen(fn_en,'r');
    En = textscan(f_en,'%s','Delimiter','\n','HeaderLines',2);
    En = En{1};
    hemi = lower(En{1}(end));
    fclose(f_en);
    
    % Read surf
    surface_type = 'pial';
    [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,surface_type));
    
    % Get electrodes
    l = read_label(sprintf('%s/%s',subjects_dir,sid),sprintf('all_surf_ielvis'));
    if (isempty(l))
        l = read_label(sprintf('%s',sid),sprintf('all_surf_ielvis'));
    end
    [~,sIdx] = sort(l(:,end));
    l = l(sIdx,:); % zero-indexed
    [n_chan,~] = size(l);
    [Xe,Ye,Ze] = sphere(20);
    
    % Plot electrode
    const_elec_radius = 2;
    for i2 = 1:n_chan
        X = Xe*const_elec_radius + s_vert(l(i2,1)+1,1);
        Y = Ye*const_elec_radius + s_vert(l(i2,1)+1,2);
        Z = Ze*const_elec_radius + s_vert(l(i2,1)+1,3);
        elec_pa_i = surf2patch(X,Y,Z);
        vert = elec_pa_i.vertices;
        q = trisurf(elec_pa_i.faces,vert(:,1),vert(:,2),vert(:,3),...
            'EdgeColor','none','FaceColor',const_elec_surface);
        hold on;
        q.SpecularStrength = 0;
        q.SpecularExponent = 1;
    end
    
    

    p = trisurf(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3),...
        'EdgeColor','none','FaceColor',const_col_surface);
    brainlight;
%     daspect([1 1 1]);
%     p.AmbientStrength = 0.5 ;
%     p.DiffuseStrength = 0.5 ;
%     p.SpecularStrength = 0.2;
%     p.SpecularExponent = 1;
%     p.BackFaceLighting = 'lit';
%     p.FaceLighting = 'gouraud';
%     cam_elev = -15;
%     camlight(-135,cam_elev);
%     camlight(45,cam_elev);
%     %camlight(-225,cam_elev);
%     %camlight(-45,cam_elev);

    ax = gca;
    if (strcmp(hemi,'r'))
        x_t = ax.XLim(1);
        y_t = ax.YLim(2);
        view(90,0);
    else
        x_t = ax.XLim(1);
        y_t = (-1) * ax.YLim(2);
        view(-90,0);
    end
    axis off
    
    x_t = 0;%ax.Position(1);
    y_t = 0;%ax.Position(4);
    
%     text(x_t,y_t,sprintf('%i',str2double(sid(2:end))),'Units','Normalized',...
%         'HorizontalAlignment','left');
    
    if (i == 8)
        %break
    end
    
    n_chan_all = [n_chan_all, n_chan];

end


p_name = sprintf('figures/48surf/all_nchan-%i_mean-%i_std-%i_anon',...
    sum(n_chan_all),round(mean(n_chan_all)),round(std(n_chan_all)));
print(h,p_name,'-depsc');
print(h,p_name,'-djpeg','-r300');
print(h,p_name,'-dpng','-r600');
