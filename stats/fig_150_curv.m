close all;

sid = 'fsaverage_sym';
hemi = 'r';
subjects_dir = '/media/klab/internal/data/coreg';
[s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,'pial'));
faces = faces + 1;
[curv, fnum] = read_curv(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid,hemi,'curv'));

% Correlate boundary and curvature
Caa2 = load('cache/fig_cluster2_reduce_new_annot_vertex_color2.mat');
vcolor = Caa2.vertex_color2;

% Check faces to find boundary
is_bound_v = false(length(curv),1);
for i = 1:length(curv)
    fac = vcolor(i,:);
    if (length(unique(fac)) == 1)
        % Face is boundary face
        is_bound_v(i) = true;
    end
end



Caa = load('cache/fig_cluster2_reduce_annot.mat');
cc_cluster = Caa.cc;

curv = sign(curv);

[rS,pS] = corr(is_bound_v,curv,'Type','Spearman');
fprintf('Correlation:\n')
fprintf('\tSpearman r: %.6f\n',rS);
fprintf('\tSpearman p-val: %.6d\n',pS);
% plot(is_bound_v, curv,'black.');
% xlabel('Is area boundary')
% ylabel('Curvature (radians)')
% set(gca,'TickDir','Out');

[P,H,STATS] = ranksum(curv(is_bound_v),curv(~is_bound_v));
%[Ht,Pt] = ttest(curv(is_bound_v),curv(~is_bound_v));


fprintf('[%s - %sh]\n',sid,hemi);
fprintf('Ranksum of curvature boundary (n=%i) vs no boundary (n=%i)\n',sum(is_bound_v),sum(~is_bound_v));
fprintf('\tp=%.9d\n',P)
fprintf('\tzval: %.9d\n',STATS.zval);
fprintf('\tranksum: %.9d\n',STATS.ranksum);


