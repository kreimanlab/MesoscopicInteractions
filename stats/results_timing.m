close all;
clear;

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

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
