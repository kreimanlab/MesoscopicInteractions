close all;
clear;

SubjectsX = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

iMX = 5;
ct_t = [];
c_t = [];
for i = 1:length(SubjectsX)
    sidX = SubjectsX{i};
    load(sprintf('./cache/xsub_out_%s_%i',sidX,iMX));
    coh_thresh = coh_thresh((Dmats>dist_thresh) & (ct > ct_thresh));
    bonf_denom = (ecog.n_comb/length(SubjectsX));
    ct_thresh_stats = binoinv(1-(p_val_ct/bonf_denom),n_graph,p_val/(bonf_seconds/w))/n_graph;
    fprintf('[%s] ct_t: %.4f coh_t: %.4f (+- %.4f)\n',sidX,ct_thresh_stats,mean(coh_thresh),std(coh_thresh));
    
    ct_t = [ct_t; ct_thresh_stats];
    c_t = [c_t; coh_thresh];
    %return
end

fprintf('[ALL] ct_t: %.6f (+- %.6f, %.6f-%.6f) coh_t: %.6f (+- %.6f)\n',...
    mean(ct_t),std(ct_t),min(ct_t),max(ct_t),mean(c_t),std(c_t));