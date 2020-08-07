close all;
clear;

SubjectsX = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

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