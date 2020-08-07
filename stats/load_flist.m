Res = carveF(res_dir,'.h5');
Dist = carveF(res_dir,'.mat');
H5 = carveF(h5_dir,'.h5');
Art = carveF(art_dir,'.mat');

%   Initialize patient list
S = cell(size(H5));
S_isok = false(size(S));
S_art = cell(size(S));
S_graph = cell(size(S));
S_dist = cell(size(S));
for i = 1:length(H5)
    s = strsplit(H5{i},'/');
    s = strsplit(s{end},'.h5');
    S{i} = s{1};
    %   Check for required files for each patient
    cond_art = (sum(contains(Art,S{i})) ~= 0);
    graphStr = sprintf('%s_graph-%s',S{i},metric);
    distStr = sprintf('%s_dists-%s',S{i},metric);
    cond_graph = (sum(contains(Res,graphStr)) ~= 0);
    cond_dist = (sum(contains(Dist,distStr)) ~= 0);
    if (cond_art && cond_graph && cond_dist)
        S_isok(i) = true;
        S_art{i} = Art{contains(Art,S{i})};
        S_graph{i} = Res{contains(Res,graphStr)};
        S_dist{i} = Dist{contains(Dist,distStr)};
    end
end
S = S(S_isok);
S_art = S_art(S_isok);
S_graph = S_graph(S_isok);
S_dist = S_dist(S_isok);
S_h5 = H5(S_isok);
for i = 1:length(S)
    fprintf('S: %s\n',S{i});
    fprintf('\tS_art: %s\n',S_art{i});
    fprintf('\tS_graph: %s\n',S_graph{i});
    fprintf('\tS_dist: %s\n',S_dist{i});
    fprintf('\tS_h5: %s\n',S_h5{i});
end

clear Art cond_art cond_dist cond_graph Dist distStr graphStr H5 i Res s;