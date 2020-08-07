close all
clear

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

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