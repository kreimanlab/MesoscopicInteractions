close all;
clear;

Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

dir_res = '/home/jerry/data/results/coh_w10';
atl = 2;
for iM = 1
    metric = metrics{iM};
    Ca14 = load(sprintf('./cache/figure_t14_%i_atl%i.mat',iM,atl));
    CaA = load(sprintf('./cache/xsub_out_all_%i_atl%i.mat',iM,atl));
    
    for iSub = 1:length(Subjects)
        sid = Subjects{iSub};
        Ca = load(sprintf('./cache/xsub_out_%s_%i_atl%i.mat',sid,iM,atl));
        fn_res = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
        R = double(h5read(fn_res,'/R'));
        return
    end
end
