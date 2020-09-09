close all;
clear;

Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
   'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
   'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',... % 'sub23',
   'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',... % 'sub30',
   'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
   'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

% Exclude subjects with HD grid
% Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
%            'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
%            'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',... % 'sub23',
%            'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',... % 'sub30',
%            'sub33','sub34','sub35','sub36','sub37'         ,'sub39','sub40',...
%            'sub41','sub42','sub43','sub44','sub45','sub46','sub47'       };
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