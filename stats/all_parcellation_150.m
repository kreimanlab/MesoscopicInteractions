close all;
clear;

%
% a = load('/media/jerry/internal/data/coreg/m00005/label/all_parcellation.mat')
%
% a = 
% 
%   struct with fields:
% 
%     AtlLabels: {1×20 cell}
%      AtlNames: {1×20 cell}
%       AtlROIs: {1×20 cell}
%     EleCoords: [104×5 double]
%     EleLabels: {104×1 cell}
%       EleHemi: {104×1 cell}
%

% Patients
Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

Cai = load('cache/fig_cluster3_cluster_i.mat');
CaT14 = load(sprintf('./cache/figure_t14_%i_150',1));
cluster_i = Cai.cluster_i; %Cai.cluster_i;
iM = 1;
Ca5 = load(sprintf('/home/jerry/data/stats/cache/fig_cluster2_reduce_%i_new.mat',iM));
fn_ca1 = sprintf('/home/jerry/data/stats/cache/fig_cluster2_A_%i.mat',iM);
Ca1 = load(fn_ca1);
[n_rois,~] = size(Ca5.Es);
CaAnnot = load('cache/fig_cluster2_reduce_annot.mat');

% Construct AtlROIs
LH = struct();
LH.numEntries = n_rois;
LH.struct_names = cell(n_rois,1);
for i = 1:n_rois
    LH.struct_names{i} = sprintf('%i',cluster_i(i)); %cluster_i(i));
end
LH.table = nan(n_rois,5);
LH.table(:,1:3) = CaAnnot.cc * 255;
LH.orig_tab = '';
AtlROI = struct();
AtlROI.LH = LH;
AtlROI.RH = LH;
AtlROIs = {};
AtlROIs{1} = AtlROI;

% Construct AtlNames
atl_name = '150-Area_Parcellation';
AtlNames = {};
AtlNames{1} = atl_name;

for i = 1:length(Subjects)
    sid = Subjects{i};
    sid_int_i = str2double(sid(2:end));
    Ca = load(sprintf('./cache/xsub_out_%s_%i.mat',sid,iM));
    atl_label = cell(Ca.ecog.n_bchan,1);
    atl_label_raw = cell(Ca.ecog.n_chan,1);
    CaSub = load(sprintf('/media/jerry/internal/data/coreg/%s/label/all_parcellation.mat',sid));
    
    
    for primary = 1:2
        count = 1;
        % For each of 150 parcellations
        for j = 1:n_rois
            % List of sid_ints and bchans for each region
            E = Ca5.Es{j};
            [n_E2,~] = size(E);

            % Assign region to electrode
            for k = 1:n_E2
                sid_int = E(k,1);
                sid_bchan = E(k,2);

                if (sid_int == sid_int_i)
                    b1c1 = Ca.ecog.bip(sid_bchan,1);
                    b1c2 = Ca.ecog.bip(sid_bchan,2);

                    % Bipolar labels
                    atl_label{sid_bchan} = sprintf('%i',cluster_i(j));

                    % Raw electrode labels
                    if (primary == 1)
                        atl_label_raw{b1c2} = sprintf('%i',cluster_i(j));
                    elseif (primary == 2)
                        atl_label_raw{b1c1} = sprintf('%i',cluster_i(j));
                    end

                    fprintf('\t(%i) bchan %i (elec %i - elec %i) - area %i\n',count,sid_bchan,b1c1,b1c2,cluster_i(k));
                    count = count + 1;
                end
            end
        end
        %return
    end
    
    %return
%     % Map the rest of the raw electrodes
%     for j = 1:length(atl_label_raw)
%         if isempty(atl_label_raw{j})
%             % remap
%             b1 = find(Ca.ecog.bip(:,1) == j);
%             fprintf('[*] Remap electrode %i\n',j)
%             return
%         end
%     end
    
    AtlLabels = {};
    AtlLabels{1} = atl_label_raw;
    
    
    % Fix coordinates
    if (~strcmp(sid,'m00049'))
        [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s','/media/klab/internal/data/coreg',sid,lower(CaSub.EleHemi{1}),'pial'));
        c_coord = s_vert(CaSub.EleCoords(:,1)+1,:);
        EleCoords = CaSub.EleCoords;
        EleCoords(:,2:4) = c_coord;
    else
        EleCoords = CaSub.EleCoords;
    end
    EleLabels = CaSub.EleLabels;
    EleHemi = CaSub.EleHemi;
    
    of_name = sprintf('/media/jerry/internal/data/coreg/%s/label/all_parcellation_150.mat',sid);
    save(of_name,'AtlLabels','AtlNames','AtlROIs','EleCoords','EleLabels','EleHemi');
    
    %     AtlLabels: {1×20 cell}
    %      AtlNames: {1×20 cell}
    %       AtlROIs: {1×20 cell}
    %     EleCoords: [104×5 double]
    %     EleLabels: {104×1 cell}
    %       EleHemi: {104×1 cell}
    
    fprintf('[*] Wrote to: %s\n',of_name);
    
end

fprintf('[!] Done.\n')


