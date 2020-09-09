close all
clear
rng shuffle;

D_anat = load('../coreg/AdjMKFV.mat');


% Read
C = readcell('AdjMKFV_labels.xlsx');
labels_abc = cell(1,length(D_anat.MK_labels));
labels_pub = cell(1,length(D_anat.MK_labels));
for i = 1:length(D_anat.MK_labels)
    fprintf('%i\t%s\n',i,D_anat.MK_labels{i});
    %[~,~,labels] = xlsread('AdjMKFV_labels','B:B');
    
    l_abc = C{i,2};
    if (isnumeric(l_abc))
        labels_abc{i} = sprintf('%i',l_abc);
    else
        labels_abc{i} = l_abc;
    end
    
    l_pub = C{i,3};
    if (isnumeric(l_pub))
        labels_pub{i} = sprintf('%i',l_pub);
    else
        labels_pub{i} = l_pub;
    end
    %labels_pub{i} = l_pub;
end

% Convert index
pub_i = nan(length(labels_abc),1);
for i = 1:length(labels_abc)
    pub_i(i) = find(strcmp(labels_pub{i},labels_abc));
end

% Load matrix
A = D_anat.AdjMKflne(pub_i,pub_i);
rois_plt = labels_abc(pub_i);
imagesc(log10(A(:,~all(A==0))));
colormap hot;
yticks(1:length(rois_plt));
yticklabels(rois_plt);
xticks(1:length(rois_plt));
xticklabels(rois_plt);
xtickangle(90);
set(gca,'TickDir','out');