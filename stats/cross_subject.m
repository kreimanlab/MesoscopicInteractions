close all;
clear;

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

iM = 1;
Nbchans = zeros(1,length(Subjects));
Nchans = zeros(1,length(Subjects));
for i = 1:length(Subjects)
    sid = Subjects{i};
    Ca = load(sprintf('./cache/xsub_out_%s_%i.mat',sid,iM));
    Nbchans(i) = Ca.ecog.n_bchan;
    Nchans(i) = Ca.ecog.n_chan;
end
fprintf('[*] n_chans: %.2f +- %.3f(%i - %i)\n',mean(Nchans),std(Nchans),min(Nchans),max(Nchans));

fprintf('[*] n_bchans: %.2f +- %.3f(%i - %i)\n',mean(Nbchans),std(Nbchans),min(Nbchans),max(Nbchans));
