close all;
clear;

Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

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
