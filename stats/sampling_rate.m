clear;
close all;

dir_h5 = '/nas_share/cuenap/data/h5_notch20';


Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

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
