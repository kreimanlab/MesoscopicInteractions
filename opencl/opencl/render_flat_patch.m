close all;
clear;


[pa] = read_patch_rev('coreg/lh.full.flat.patch.3d');
k = convhull(pa.x,pa.y);
[v_pial,tri_pial] = read_surf('coreg/lh.pial');
tri_pial = tri_pial - 1;
tri_flat = [];
tri_flat_i = 1;
for i = 1:length(tri_pial)
    tri_1 = find(pa.vno == tri_pial(i,1),1);
    tri_2 = find(pa.vno == tri_pial(i,2),1);
    tri_3 = find(pa.vno == tri_pial(i,3),1);
    if ((~isempty(tri_1)) && (~isempty(tri_2)) && (~isempty(tri_3)))
        tri_flat(tri_flat_i,:) = [tri_1, tri_2, tri_3];
        tri_flat_i = tri_flat_i + 1;
    end
    fprintf('%i of %i\n',i,length(tri_pial));
end

%plot(pa.x,pa.y,'black.');
trisurf(tri_flat,pa.x,pa.y,pa.z);
save('coreg/fsaverage_sym_flat');