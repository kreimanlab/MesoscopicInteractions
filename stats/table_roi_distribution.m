close all;
clear;

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

subjects_dir = '/media/klab/internal/data/coreg';
atlas_i = 1;
atlas_unknown_i = 1;

ofi = fopen('./table_roi_distribution_out.csv','w');

for atlas_i = 1:20

    % get full list of ROIs from first patient
    A1 = load(sprintf('%s/%s/label/all_parcellation.mat',subjects_dir,Subjects{1}));
    rois = A1.AtlROIs{atlas_i}.LH;
    roi_names = rois.struct_names;
    roi_names_count = zeros(length(roi_names),1);
    roi_names_count_l = zeros(length(roi_names),1);
    roi_names_count_r = zeros(length(roi_names),1);


    for i = 1:length(Subjects)
        sid = Subjects{i};
        A = load(sprintf('%s/%s/label/all_parcellation.mat',subjects_dir,sid));
        %fprintf('%s\t%s\n',sid,A.AtlNames{atlas_i});
        elec_roi_labels = A.AtlLabels{atlas_i};
        rois = A.AtlROIs{atlas_i}.LH;
        
        % Load bip
        Ca = load(sprintf('./cache/xsub_out_%s_%i_atl%i',sid,1,atlas_i));
        
        % Load hemisphere
        fn_hemi = sprintf('%s/%s/elec_recon/%s.electrodeNames',subjects_dir,sid,sid);
        He = textscan(fopen(fn_hemi,'r'),'%s','HeaderLines',2);
        He = reshape(He{1},3,[]);
        hemis = He(3,:)';
        
        % Filter only bip electrodes
        elec_roi_labels = elec_roi_labels(Ca.ecog.bip(:,1));
        hemis = hemis(Ca.ecog.bip(:,1));

        if (length(hemis) ~= length(elec_roi_labels))
            fprintf(2,'[E] Number of channels in hemisphere file: %i does not match number in .mat file: %i\n',length(hemis),length(elec_roi_labels));
        end

        % Increment roi count for each electrode
        for j = 1:length(elec_roi_labels)
            elec_roi = elec_roi_labels{j};
            roi_names_idx = find(strcmp(roi_names,elec_roi));
            if (isempty(roi_names_idx)) %unknown roi
                roi_names_count(atlas_unknown_i,:) = roi_names_count(atlas_unknown_i,:) + 1;
                if (strcmp(hemis{j},'L'))
                    roi_names_count_l(atlas_unknown_i,:) = roi_names_count_l(atlas_unknown_i,:) + 1;
                elseif (strcmp(hemis{j},'R'))
                    roi_names_count_r(atlas_unknown_i,:) = roi_names_count_r(atlas_unknown_i,:) + 1;
                else
                    fprintf(2,'[E] Hemisphere not recognized: %s\n',hemis{j});
                end
            else
                roi_names_count(roi_names_idx,:) = roi_names_count(roi_names_idx,:) + 1;
                if (strcmp(hemis{j},'L'))
                    roi_names_count_l(roi_names_idx,:) = roi_names_count_l(roi_names_idx,:) + 1;
                elseif (strcmp(hemis{j},'R'))
                    roi_names_count_r(roi_names_idx,:) = roi_names_count_r(roi_names_idx,:) + 1;
                else
                    fprintf(2,'[E] Hemisphere not recognized: %s\n',hemis{j});
                end
            end
        end
        %return
    end

    n_total = sum(roi_names_count);
    [roi_names_count_s,sidx] = sort(roi_names_count,'descend');
    roi_names = roi_names(sidx);
    roi_names_count_l = roi_names_count_l(sidx);
    roi_names_count_r = roi_names_count_r(sidx);
    
    for i = 1:length(roi_names)
        fprintf(ofi,'%s,%i,%s,%i,%.9f,%i,%.9f,%i,%.9f\n',...
            A.AtlNames{atlas_i},sidx(i)-1,roi_names{i},...
            roi_names_count_s(i),1*roi_names_count_s(i)/n_total,...
            roi_names_count_l(i),1*roi_names_count_l(i)/n_total,...
            roi_names_count_r(i),1*roi_names_count_r(i)/n_total);
    end
    

end

fclose(ofi);