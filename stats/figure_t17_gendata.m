close all;

metricsp = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

for i = 1:5
%     metric = metricsp{i};
    CaT14 = load(sprintf('./cache/figure_t14_%i_150',i));
    %return
    data = CaT14.Adj_plt2_cl(CaT14.cluster_i,CaT14.cluster_i);
    %data(isnan(data)) = 0;
    [n_rows,~] = size(data);
    fn_out = sprintf('./NCestimation_V2/t14_%i',i);
    of = fopen(fn_out,'w');
    for j = 1:n_rows
        fprintf(of,'%.9f\t',data(j,:));
        fprintf(of,'\n');
    end
    fclose(of);
end