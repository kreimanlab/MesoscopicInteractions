close all;

Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
   'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
   'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',... % 'sub23',
   'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',... % 'sub30',
   'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
   'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

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