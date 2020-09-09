clear;
close all;

dir_h5 = '/nas_share/cuenap/data/h5_notch20';


Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

n250 = 0;
n256 = 0;
nUNK = 0;
nAll = 0;
for i = 1:length(Subjects)
    sid = Subjects{i};
    fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
    rate = round(h5readatt(fn_h5,'/h5eeg/eeg','rate'));
    if (rate == 250)
        fprintf('[*] %s rate: %i\n',sid,rate)
        n250 = n250 + 1;
    elseif (rate == 256)
        fprintf('[*] %s rate: %i\n',sid,rate)
        n256 = n256 + 1;
    else
        fprintf(2,'[!] Unknown rate: %i\n',rate)
        nUNK = nUNK + 1;
    end
    nAll = nAll + 1;
end

fprintf('[*] 250 Hz: %i of %i (%.5f)\n',n250,nAll,n250/nAll);
fprintf('[*] 256 Hz: %i of %i (%.5f)\n',n256,nAll,n256/nAll);
fprintf('[*] Unknown: %i of %i (%.5f)\n',nUNK,nAll,n256/nAll);
