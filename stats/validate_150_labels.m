close all;
clear;

toggle_plt = false;

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

AChans = cell(1,length(Subjects));
for iSub = 1:length(Subjects)

    % Load 150 labels
    fprintf('[%i of %i]\n',iSub,length(Subjects));
    sid = Subjects{iSub};
    C = load(sprintf('/media/jerry/internal/data/coreg/%s/label/all_parcellation_150.mat',sid));
    hemi = lower(C.EleHemi{1});

    % Load surface
    dir_corLp = '/media/klab/internal/data/coreg';
    [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',dir_corLp,sid,hemi,'pial'));

    if (toggle_plt)
        % Plot surface
        h = figure();
        trisurf(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3),'EdgeColor','None','FaceColor',[1 1 1]*0.8);
        hold all;
    end

    % Channel labels
    chan_labels = C.AtlLabels{1};
    c_coord = C.EleCoords(:,2:4);
    %chan_reg = C.AtlROIs{1}.RH.struct_names;
    Cai = load('cache/fig_cluster3_cluster_i.mat');
    chan_reg = chan_labels; %C.AtlLabels{1};
    %c_coord = s_vert(C.EleCoords(:,1)+1,:);

    iM = 1;
    CaT14 = load(sprintf('./cache/figure_t14_%i_150',iM));

    % Construct cluster_i
    cluster_i = load('cache/fig_cluster3_cluster_i.mat');
    cluster_i = cluster_i.cluster_i;
    cluster_i2 = zeros(size(cluster_i));
    for i = 1:length(cluster_i)
        idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(i))));
        cluster_i2(i) = idx;
    end
    cluster_i = cluster_i2;
    Ca = load(sprintf('%s/xsub_out_%s_%i.mat','./cache',sid,iM));
    CaRloc = load(sprintf('%s/fig_cluster2_reduce_%i_new.mat','cache',iM));

    atl_chans = zeros(1,Ca.ecog.n_bchan);
    for i = 1:(Ca.ecog.n_bchan)
        ic = Ca.ecog.bip(i,1) == C.EleCoords(:,end);

        bchan = C.EleCoords(i,end);

        % Build list of locations for subject
        sid_const_int_tmp = str2double(sid(2:end));
        Eb2avg = nan(Ca.ecog.n_bchan,1);
        ico = 0;
        for iet = 1:length(CaRloc.Es)
            est = CaRloc.Es{iet};
            est2 = CaRloc.Es2{iet};
            sidints = est(:,1);
            for ieu = 1:length(sidints)
                if ( (est(ieu,1)==sid_const_int_tmp) ) % (est2(ieu,6) > 0) &&
                    Eb2avg(est(ieu,2)) = iet;
                    ico = ico + 1;
                end
            end
        end

        b1 = Eb2avg(i);
        text_val = sprintf('%i',cluster_i(b1));
        
        % save
        atl_chans(i) = cluster_i(b1);
        % Add 150 custom parcellation text
        %desc = sprintf('%s<br>Area Custom: %i',desc,cluster_i(b1));

        %idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(i))));
        %text_val = CaT14.rois_plt{idx};
        %idx = find(strcmp(CaT14.rois_plt_all_parcellation,chan_reg{i}));
        %text_val = CaT14.rois_plt{idx};

        %text(c_coord(i,1),c_coord(i,2),c_coord(i,3),sprintf('%i',CaT14.cluster_i(Cai.cluster_i(idx))),'HorizontalAlignment','center');
        %text(c_coord(i,1)+10,c_coord(i,2),c_coord(i,3),sprintf('%s',chan_reg{i}),'HorizontalAlignment','center','FontSize',12);

        if (toggle_plt)
            if (strcmp(hemi,'r'))
                R_offset = 10;
            else
                R_offset = -10;
            end
            R_offset = 0;
            plot3(c_coord(ic,1),c_coord(ic,2),c_coord(ic,3),'black.','MarkerSize',12);
            text(c_coord(ic,1)+R_offset,c_coord(ic,2),c_coord(ic,3),sprintf('%s',text_val),'HorizontalAlignment','center','FontSize',12);
        end

    end

    if (toggle_plt)
        % Lighting
        brainlight;
        axis off;
        if (strcmp(hemi,'r'))
            view(90,0);
        else 
            view(-90,0);
        end
    end

    AChans{iSub} = atl_chans;
end

save('./cache/validate_150_labels','AChans');