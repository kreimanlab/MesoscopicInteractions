close all;
clear;

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
   'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
   'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',... % 'm00045',
   'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',... % 'm00056',
   'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
   'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

% Exclude subjects with HD grid
% Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
%            'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
%            'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',... % 'm00045',
%            'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',... % 'm00056',
%            'm00060','m00061','m00068','m00071','m00073'         ,'m00079','m00083',...
%            'm00084','m00095','m00096','m00097','m00100','m00107','m00122'       };
metrics = {'pcBroadband','pcTheta','pcAlpha','pcBeta','pcGamma'};

for iM = 1:5
    Ds = [];
    coh_thresh_all = [];
    for i = 1:length(Subjects)
        sid = Subjects{i};
        %fn_paths = sprintf('brainexport/%s_6all_%i.mat',sid,iM);
        fn_path2 = sprintf('cache/xsub_out_%s_%i_atl2.mat',sid,iM);
        %fprintf('[%i/%i] Load: %s\n',i,length(Subjects),fn_path2);
        Ca2 = load(fn_path2);
        %CaP = load(fn_paths);
        %interd = CaP.Ca.ecog.bip(:,3);
        interd = Ca2.ecog.bip(:,3);

        Ds = [Ds; interd];
        n_bip_pairs(i) = nchoosek(Ca2.ecog.n_bchan,2);

        d_idx = Ca2.Dmats > Ca2.dist_thresh;
        coh_thresh_all = [coh_thresh_all; Ca2.coh_thresh(d_idx)];
        %coh_thresh_all = [coh_thresh_all; Ca2.coh_thresh];
    end
    
    %coh_thresh_all = coh_thresh_all(coh_thresh_all ~= 1);

    fprintf('[*] %s Coh thresh\n',metrics{iM});
    fprintf('\tn: %i\n',length(coh_thresh_all));
    fprintf('\tmean: %.4f\n',mean(coh_thresh_all));
    fprintf('\tmedian: %.4f\n',median(coh_thresh_all));
    fprintf('\tstd: %.4f\n',std(coh_thresh_all));
    fprintf('\tmin: %.4f\n',min(coh_thresh_all));
    fprintf('\tmax: %.4f\n',max(coh_thresh_all));
end

% methods_bip_distance
% [*] pcBroadband Coh thresh
% 	n: 124865
% 	mean: 0.3473
% 	median: 0.2399
% 	std: 0.2671
% 	min: 0.1171
% 	max: 1.0000
% [*] pcDelta Coh thresh
% 	n: 124865
% 	mean: 0.8732
% 	median: 1.0000
% 	std: 0.2394
% 	min: 0.2367
% 	max: 1.0000
% [*] pcTheta Coh thresh
% 	n: 124865
% 	mean: 0.6526
% 	median: 0.8176
% 	std: 0.3530
% 	min: 0.0756
% 	max: 1.0000
% [*] pcAlpha Coh thresh
% 	n: 124865
% 	mean: 0.5323
% 	median: 0.4567
% 	std: 0.3291
% 	min: 0.1597
% 	max: 1.0000
% [*] pcBeta Coh thresh
% 	n: 124865
% 	mean: 0.2860
% 	median: 0.1598
% 	std: 0.2489
% 	min: 0.0903
% 	max: 1.0000