function [ O ] = getOutputFilenames( resultsDir, isFitted )
%getOutputFilenames returns a list of output files that are ready
%   isFitted - whether to return list of fitted or unfitted files

%fprintf('[*] Looking for results in: %s\n',resultsDir);
%fprintf('[*] Looking for h5eeg in: %s\n',h5Dir);
D = dir(resultsDir);
Disdir = [D.isdir];
Dname = {D.name};
Dname = Dname(~Disdir);
SubjectsAllPerm = Dname(startsWith(Dname,'m') & endsWith(Dname,'.h5') & contains(Dname,'perm'));
SubjectsAllGraph = Dname(startsWith(Dname,'m') & endsWith(Dname,'.h5') & contains(Dname,'graph'));
SubjectsAllDists = Dname(startsWith(Dname,'m') & endsWith(Dname,'.mat') & contains(Dname,'dists'));
n_f = length(SubjectsAllPerm);

% Trim graph files with no perm partner
SubjectsAllGraphN = {};
count = 1;
for i = 1:length(SubjectsAllGraph)
    sag = SubjectsAllGraph{i};
    sag = strsplit(sag,'.h5');
    %sad = replace(sag{1},'graph','dists');
    sag = replace(sag{1},'graph','perm');
    if ((sum(startsWith(SubjectsAllPerm,sag)) > 0))
        SubjectsAllGraphN{count} = SubjectsAllGraph{i};
        count = count + 1;
    end
end
SubjectsAllGraph = SubjectsAllGraphN;

% Trim perm files with no graph partner
SubjectsAllPermN = {};
count = 1;
for i = 1:length(SubjectsAllPerm)
    sap = SubjectsAllPerm{i};
    sap = strsplit(sap,'.h5');
    %sad = replace(sap{1},'perm','dists');
    sap = replace(sap{1},'perm','graph');
    sap = strsplit(sap,'-');
    sap = strjoin(sap(1:(end-1)),'-');
    if ((sum(startsWith(SubjectsAllGraph,sap)) > 0))
        SubjectsAllPermN{count} = SubjectsAllPerm{i};
        count = count + 1;
    end
end
SubjectsAllPerm = SubjectsAllPermN;

% Trim out dists
SubjectsAllPermN = {};
SubjectsAllGraphN = {};
count = 1;
for i = 1:length(SubjectsAllPerm)
    sap = SubjectsAllPerm{i};
    sap = strsplit(sap,'.h5');
    sad = replace(sap{1},'perm','dists');
    if (isFitted)
        dCond = (sum(startsWith(SubjectsAllDists,sad)) > 0);
    else
        dCond = (sum(startsWith(SubjectsAllDists,sad)) == 0);
    end
    
    if (dCond)
        SubjectsAllPermN{count} = SubjectsAllPerm{i};
        SubjectsAllGraphN{count} = SubjectsAllGraph{i};
        count = count + 1;
    end
end
SubjectsAllPerm = SubjectsAllPermN;
SubjectsAllGraph = SubjectsAllGraphN;

n_f = length(SubjectsAllPerm);
% for i = 1:n_f
%     if (isFitted)
%         fprintf('%s : %s : %s\n',SubjectsAllPerm{i},SubjectsAllGraph{i},SubjectsAllDists{i});
%     else
%         fprintf('%s : %s\n',SubjectsAllPerm{i},SubjectsAllGraph{i});
%     end
% end

O.SubjectsAllPerm = SubjectsAllPerm;
O.SubjectsAllGraph = SubjectsAllGraph;
O.SubjectsAllDists = SubjectsAllDists;
O.n_f = n_f;

end

