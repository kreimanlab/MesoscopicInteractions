close all;
clear;

Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

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