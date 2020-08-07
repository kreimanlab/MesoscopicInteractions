close all;

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
   'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
   'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',... % 'm00045',
   'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',... % 'm00056',
   'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
   'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

dir_coreg = '/media/jerry/internal/data/coreg';

for atl_i = [2, 7]

    Unknowns = zeros(1,length(Subjects));
    for i = 1:length(Subjects)
        sid = Subjects{i};
        AP = load(sprintf('%s/%s/label/all_parcellation.mat',dir_coreg,sid));
        label = AP.AtlLabels{atl_i};
        is_unknown = strcmp(label,'unknown');
        Unknowns(i) = sum(is_unknown);
    end
    fprintf('[%s]\n',AP.AtlNames{atl_i})
    fprintf('[*] total unknowns: %i\n',sum(Unknowns));
    fprintf('[*] subject average unknowns: %.4f\n',mean(Unknowns));
    fprintf('[*] subject std unknowns: %.4f\n',std(Unknowns));
    fprintf('Subjects with unknowns:\n')
    disp(find(Unknowns))


end