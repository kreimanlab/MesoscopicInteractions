close all;
clear;

% Must run: figure_t14_allatl.m first
% Run figure_brainexport_static_fsaverage_red.m next

odir = '/home/jerry/Downloads/blend4web_ce/projects/brainview/assets/data';
cond_dryrun = false;

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

% --- copy paste ----------------------------------------------------------------
% deprecated
cluster_i = load('cache/fig_cluster3_cluster_i.mat');
cluster_i = cluster_i.cluster_i;
%cluster_i = 1:length(cluster_i); %does not work
% end deprecated

CaT14 = load(sprintf('./cache/figure_t14_%i_150',1));
cluster_i2 = zeros(size(cluster_i));
for i = 1:length(cluster_i)
    idx = find(strcmp(CaT14.rois_plt_all_parcellation,sprintf('%i',cluster_i(i))));
    cluster_i2(i) = idx;
end
cluster_i = cluster_i2; 
% --- end copy paste ----------------------------------------------------------------



for iM = 1:5 % 2:5
    for iSub = 0:48
        
        % Handle average brain
        if (iSub == 0)
            sid = 'fsaverage_sym';
            atlmax = 20;
        else
            sid = Subjects{iSub};
            atlmax = 21;
        end
        
        for atl = [0 2 3 7 atlmax] %0:atlmax
            cond_bypass = false;
            
            ifn = sprintf('./brainexport/controls_data_sub-%i_freq-%i_atl-%i',iSub,iM,atl);
            ofn = sprintf('%s/controls_data_sub-%i_freq-%i_atl-%i',odir,iSub,iM,atl);
            ofn_cb = sprintf('%s/d_controls_data_sub-%i_freq-%i_atl-%i',odir,iSub,iM,atl);
            
            fprintf('[In] %s.mat\n',ifn);
            fprintf('[Out] %s.json\n',ofn);
            toggle_cont = true;
            try
                C = load(ifn);
            catch
                fprintf(2,'[*] Skipped, no .mat found\n');
                toggle_cont = false;
            end
            
            if (toggle_cont)
            
                if ((iSub == 0) && (atl == 0))
                    % For fsaverage_sym
                    reg_dict = cell(length(C.S.atl_chans),1);
                    for ird = 1:length(reg_dict)
                        reg_dict{ird} = sprintf('%i',ird);
                    end
                else
                    % Match each point with index on adjacency matrix

                    % Build reg_dict
                    if (iSub == 0)

                        if (atl == 9)
                            % Exception for area STPa, which appeared twice in
                            % the correlation matrix due to overlap mapping
                            % non-uniqueness between 150-area custom regions
                            % and atlas regions. This only appeared for atl=9,
                            % FVEall

                            i_stpa = find(strcmp(C.S.labels,'STPa'));
                            i_pass = true(1,length(C.S.labels));
                            i_pass(i_stpa(end)) = false;

                            C.S.Img = C.S.Img(i_pass,i_pass,:);
                            C.S.A = C.S.A(i_pass,i_pass);
                            C.S.dendro_Z = C.S.dendro_Z(i_pass(1:(end-1)),:);
                            C.S.dendro_reorder = C.S.dendro_reorder(i_pass);
                            C.S.labels = C.S.labels(i_pass);
                        end

                        Ca = load(sprintf('./cache/xsub_out_%s_%i_atl%i','m00001',iM,atl));
                        reg_dict = Ca.C.AtlROIs{atl}.RH.struct_names;

                        % clean out empty region names
                        reg_dict_clean = {};
                        i_c = 1;
                        for ic = 1:length(reg_dict)
                            %if ((~isempty(reg_dict{ic})) && (~strcmp(lower(reg_dict{ic}),'unknown')))
                            if (~isempty(reg_dict{ic}))
                                reg_dict_clean{i_c} = reg_dict{ic};
                                i_c = i_c + 1;
                            end
                        end
                        reg_dict = reg_dict_clean;

                        % Check region dict size
                        cond_size_mismatch = (length(reg_dict) ~= length(C.S.atl_chansR));
                        if (cond_size_mismatch)
                            if (atl == 1)
                                reg_dict = reg_dict(2:end);
                                fprintf('[!] region mismatch fixed.\n');
                            else
                                fprintf(2,'[!] region dictionary size mismatch.\n');
                                return
                            end
                        end

                        % matrix 2 node definition
                        % - this is special because it is based on region
                        % overlap with 150 custom regions, as determined in
                        % custom_parcellation_atlas_overlap.m
                        idx_atl_chansR = false(size(C.S.atl_chansR));
                        reg_dict_clean = nettoyer_rois(reg_dict,atl);
                        for i_f = 1:length(idx_atl_chansR)
                            % was ich gesucht hab ich gefunden
                            idx_atl_chansR(i_f) = sum(strcmp(C.S.labels,reg_dict_clean{i_f})) > 0;
                        end
                        atl_chansR_bypass = C.S.atl_chansR(idx_atl_chansR);
                        %cond_bypass = true;
                        
                        % allatl labels
                        atl_labels = reg_dict(C.S.atl_chans);
                        
                        if (atl == 2)
                            % Format Desikan labels
                            for ia = 1:length(atl_labels)
                                atl_labels{ia} = convertRoiDK(atl_labels{ia});
                            end
                        elseif (atl == 7)
                            % Format Markov labels
                            atl_labels_mk = nettoyer_rois((atl_labels));
                            atl_labels = atl_labels_mk;
                        end
                        reg_ind = atl_labels;
                        %return;

                    else
                        % For single subjects
                        if (atl == 0)
                            % regions are simply bipolar electrodes
                            CaS = load(sprintf('./cache/xsub_out_%s_%i_atl%i',sid,iM,1));
                            reg_dict = cell(1,CaS.ecog.n_bchan);
                            for ird = 1:length(reg_dict)
                                reg_dict{ird} = sprintf('%i',ird);
                            end
                        elseif (atl == 21)
                            Ct = load(sprintf('/media/klab/internal/data/coreg/%s/label/all_parcellation_150.mat',Subjects{iSub}));
                            reg_dict = Ct.AtlROIs{1}.RH.struct_names;
                            %reg_dict = reg_dict(cluster_i);
                            %reg_ind = reg_dict(C.S.atl_chans);
                            atl_labels_bip = cell(length(C.S.atl_chans),1);
                            for ia = 1:length(atl_labels_bip)
                                atl_labels_bip{ia} = sprintf('%i',C.S.atl_chans(ia));
                            end
                            reg_ind = atl_labels_bip;
                        else
                            CaS = load(sprintf('./cache/xsub_out_%s_%i_atl%i',sid,iM,atl));
                            hemi = CaS.C.EleHemi{1};
                            if (strcmp(hemi,'R'))
                                reg_dict = CaS.C.AtlROIs{atl}.RH.struct_names;
                            else
                                reg_dict = CaS.C.AtlROIs{atl}.LH.struct_names;
                            end
                            atl_labels = CaS.C.AtlLabels{atl};
                            atl_labels_bip = atl_labels(CaS.ecog.bip(:,1));
                            reg_ind = atl_labels_bip;
                        end
                    end
                end
% 
%                 reg_index = C.S.atl_chans;
%                 
%                 % unknown indices
%                 reg_index(isnan(reg_index)) = 1;
                
                %
                % Set electrode-wise area labels
                % 
                if (atl ~= 0)
                    % For areas
%                     reg_index_str = reg_dict(reg_index);
                    
                    reg_index_str = reg_ind;
                    
                    % clean labels according to figure_t14_allatl
                    %reg_index_str = nettoyer_rois(reg_index_str,atl);
                else
                    % For bipolar electrodes
                    reg_index_str = reg_dict;
                end

                if (atl == 21)
                    % Custom parcellation
                    
                    % map nodes onto matrix
                    AC = load('./cache/validate_150_labels');
                    reg_index_A = nan(size(AC.AChans{iSub}));
                    for i_ris = 1:length(reg_index_A)
                        idx = find(strcmp(C.S.labels,sprintf('%i',AC.AChans{iSub}(i_ris))));
                        if (~isempty(idx))
                            reg_index_A(i_ris) = idx;
                        end
                    end
                    
                    % map matrix onto nodes
                    strList = cell(1,length(AC.AChans{iSub}));
                    for isl = 1:length(strList)
                        strList{isl} = sprintf('%i',AC.AChans{iSub}(isl));
                    end
                    reg_index_A_rev = cell(1,length(C.S.labels));
                    for i_ris = 1:length(reg_index_A_rev)
                        idx = find(strcmp(C.S.labels{i_ris},strList));
                        reg_index_A_rev{i_ris} = idx;
                    end
                    atl_chansR_bypass = reg_index_A_rev;
                    cond_bypass = true;
                else
                    % Other allatl parcellations
                    
                    reg_index_A = nan(size(reg_index_str));
                    % for every electrode, region label
                    for i_ris = 1:length(reg_index_str)
                        % area to search
                        suchen = reg_index_str{i_ris};
                        % get adj matrix index
                        a = find(strcmp(suchen,C.S.labels));
                        if (~isempty(a))
                            try
                                reg_index_A(i_ris) = a;
                            catch
                                fprintf(2,'Flag!\n');
                                %disp(a)
                                disp( C.S.labels(a) );
                                return
                            end
                        else
                            reg_index_A(i_ris) = NaN;
                        end
                    end
                end
                % reg_index_A contains Adj index for eah electrode

                % Match matrix node with point
                reg_index_A_rev = cell(size(C.S.labels));
                for i_ris = 1:length(reg_index_A_rev)
                    reg_index_A_rev{i_ris} = find(strcmp(C.S.labels{i_ris},reg_index_str));
                end

                % Clean labels
                rois_plt = C.S.labels;
                if ((atl == 2) && (iSub ~= 0))
                    for i = 1:length(rois_plt)
                        rois_plt{i} = convertRoiDK(rois_plt{i});
                    end
                end
                
                O = struct();
                O.A_bg_color_r = round(C.S.Img(:,:,1) * 255);
                O.A_bg_color_g = round(C.S.Img(:,:,2) * 255);
                O.A_bg_color_b = round(C.S.Img(:,:,3) * 255);
                O.A_coherence = C.S.A;
                if ((atl ~= 0) || (iSub == 0))
                    try
                        O.A_coherence_std = C.S.Astd;
                    catch
                        O.A_coherence_std = nan(size(C.S.A));
                    end
                end
                O.A_coherence_min = nanmin(C.S.A(C.S.A ~= 0));
                O.A_coherence_max = nanmax(C.S.A(C.S.A ~= 0));
%                 if ((iSub == 0) || ((iSub ~= 0) && (atl ~= 0)))
%                     O.A_npairs = C.S.A_npairs;
%                 end
                if (iSub == 0)
                    O.A_npairs = C.S.A_npairs;
                    O.A_npairs_sig = C.S.A_npairs_sig;
                    O.A_nusubs = C.S.A_nusubs;
                    O.A_nusubs_sig = C.S.A_nusubs_sig;
                else
                    if (atl ~= 0)
                        try
                            O.A_npairs = C.S.A_npairs;
                            O.A_npairs_sig = C.S.A_npairs_sig;
                        catch
                            O.A_npairs = nan(size(O.A_coherence));
                            O.A_npairs_sig = nan(size(O.A_coherence));
                        end
                    end
                end
                O.A_name = C.S.atl_name;
                O.A_labels = nettoyer_rois(rois_plt,atl);
                O.A_node2mat = reg_index_A; % value of row/col of A per node
                if (cond_bypass)
                    O.A_mat2node = atl_chansR_bypass;
                else
                    O.A_mat2node = reg_index_A_rev;
                end
                outStr = jsonencode(O);

                %return
                % Preview output
                disp(O.A_name);

    %             return
                if (~cond_dryrun)
                    
                    % Dendrogram
                    try
                        dendro_scale = 1 + length(C.S.labels)*(1/5);
                        h = figure('Position',round(dendro_scale*[0 0 55 300]),'visible','off');
                        set(gcf,'Color','None');
                        hD = dendrogram(C.S.dendro_Z,0,'Reorder',C.S.dendro_reorder,'Orientation','right','ColorThreshold',0);
                        for ihD = 1:length(hD)
                            hD(ihD).Color = 0.75*[1 1 1];
                            hD(ihD).LineWidth = 3;
                        end
                        %return
                        axis tight;
                        axis off;
                        im = frame2im(getframe(h));
                        im_alpha = all(im<=0,3) ~= 1;

                        % Crop left and right
                        idx_crop = ~all(im_alpha == 0);
                        im = im(:,idx_crop,:);
                        im_alpha = im_alpha(:,idx_crop);

                        % Crop top
                        cond_pass = all(im_alpha(1,:)==0);
                        while (cond_pass)
                            im = im(2:end,:,:);
                            im_alpha = im_alpha(2:end,:);
                            cond_pass = all(im_alpha(1,:)==0);
                        end

                        % Crop bottom
                        cond_pass = all(im_alpha(end,:)==0);
                        while (cond_pass)
                            im = im(1:(end-1),:,:);
                            im_alpha = im_alpha(1:(end-1),:);
                            cond_pass = all(im_alpha(end,:)==0);
                        end

                        %return
                        imwrite(im,[ofn_cb,'.png'],'png','Alpha',double(im_alpha));
                        close(h);
                    catch 
                        %rethrow(e);
                        fprintf('[!] dendrogram failed.\n');
                        h = figure('Position',[0,0,1,length(O.A_labels)],'visible','off');
                        imagesc(zeros(length(O.A_labels),1)); colormap gray;
                        im = frame2im(getframe(h));
                        im_alpha = false(size(im(:,:,1)));
                        imwrite(im,[ofn_cb,'.png'],'png','Alpha',double(im_alpha));
                        close(h);
                    end
                    
                    % Save json
                    sf = fopen([ofn,'.json'],'w');
                    fprintf(sf,outStr);
                    fclose(sf);
                end
            end
        end
    end
end



function [rois_out] = nettoyer_rois(rois_plt,atl)
%     for i = 1:length(rois_plt)
%         rois_plt{i} = replace(rois_plt{i},'/','_');
%     end

    for i = 1:length(rois_plt)
        % Single atlas exceptions
        if (startsWith(rois_plt{i},'L_'))
            rpt = strsplit(rois_plt{i},'L_');
            rpt = strsplit(rpt{end},'_ROI');
            rois_plt{i} = rpt{1};
        end
        if (endsWith(rois_plt{i},'_M132'))
            rpt = strsplit(rois_plt{i},'_M132');
            rois_plt{i} = rpt{1};
        end
        if (endsWith(rois_plt{i},'/M132'))
            rpt = strsplit(rois_plt{i},'/M132');
            rois_plt{i} = rpt{1};
        end
        if (startsWith(rois_plt{i},'FVE_all_'))
            rpt = strsplit(rois_plt{i},'FVE_all_');
            rois_plt{i} = rpt{end};
        elseif (startsWith(rois_plt{i},'FVE_'))
            rpt = strsplit(rois_plt{i},'FVE_');
            rois_plt{i} = rpt{end};
        end
        if (startsWith(rois_plt{i},'LVE00_'))
            rpt = strsplit(rois_plt{i},'LVE00_');
            rois_plt{i} = rpt{end};
        end
        if (startsWith(rois_plt{i},'PHT00_'))
            rpt = strsplit(rois_plt{i},'PHT00_');
            rois_plt{i} = rpt{end};
        end
        if (startsWith(rois_plt{i},'BoninBailey_'))
            rpt = strsplit(rois_plt{i},'BoninBailey_');
            rois_plt{i} = rpt{end};
        end
        if (startsWith(rois_plt{i},'FerryEtAl_00'))
            rpt = strsplit(rois_plt{i},'FerryEtAl_00');
            rois_plt{i} = rpt{end};
        end
        if (startsWith(rois_plt{i},'UD86_'))
            rpt = strsplit(rois_plt{i},'UD86_');
            rois_plt{i} = rpt{end};
        end
        if (startsWith(rois_plt{i},'PGR91_'))
            rpt = strsplit(rois_plt{i},'PGR91_');
            rois_plt{i} = rpt{end};
        end
        if (startsWith(rois_plt{i},'LyonKaas02_'))
            rpt = strsplit(rois_plt{i},'LyonKaas02_');
            rois_plt{i} = rpt{end};
        end
        if (startsWith(rois_plt{i},'BRL87_'))
            rpt = strsplit(rois_plt{i},'BRL87_');
            rois_plt{i} = rpt{end};
        end

        % Remove underscores
        rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '_',' ');
        rois_plt{i} = replace(rois_plt{i}(1:(end-0)), '.',' ');

%         if (atl == 2)
%             rois_plt{i} = convertRoiDK(rois_plt{i});
%         end
    end
    
    rois_out = rois_plt;
end
