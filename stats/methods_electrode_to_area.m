close all
clear

Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

issame = 0;
isdrop = 0;
total = 0;
n_chan_tot = 0;
Nchans = [];
for i = 1:length(Subjects)
    atl = 2;
    ifn = sprintf('cache/xsub_out_%s_1_atl1.mat',Subjects{i});
    %fprintf('Read: %s\n',ifn);
    L = load(ifn);
    ele_rois = L.C.AtlLabels{atl};
    issame_bip = false(L.ecog.n_bchan,1);
    isdrop_bip = false(L.ecog.n_bchan,1);
    ekey = reshape(L.ecog.bip(:,1:2),[],1);
    for j = 1:L.ecog.n_bchan
        ee = L.ecog.bip(j,1:2);
        issame_bip(j) = strcmp(ele_rois{ee});
        isdrop_bip(j) = (sum(ekey == ee(1)) == 1) | (sum(ekey == ee(1)) == 1);
    end
    fprintf('%s\tsame: %i, total: %i, frac: %.2f %%\n',L.sid,...
        sum(issame_bip),numel(issame_bip),100*sum(issame_bip)/numel(issame_bip));
    dr = sum(isdrop_bip & (~issame_bip));
    fprintf('%s\tdrop: %i, total: %i, frac: %.2f %%\n',L.sid,...
        dr,numel(issame_bip),100*dr/numel(issame_bip));
    issame = issame + sum(issame_bip);
    total = total + numel(issame_bip);
    isdrop = isdrop + sum(isdrop_bip & (~issame_bip));
    n_chan_tot = n_chan_tot + L.ecog.n_chan;
    Nchans(i) = L.ecog.n_chan;
end
fprintf('Nchan total: %i\n',n_chan_tot);
fprintf('TOTAL\tsame: %i, total: %i, frac: %.2f %%\n',issame,total,100*issame/total);
fprintf('Number of dropped electrode area pairs: %i\n',isdrop);