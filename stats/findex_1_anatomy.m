% 2863 total left after ROI restriction

close all;
clear;
rng shuffle;

metric = 'pcBroadband';
metric_i = 1;
atl_i = 2;

% Dirs
dir_h5 = '/media/jerry/KLAB101/h5_notch20';
dir_cache = './cache';

% Go through all subjects
Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};


roi1 = {'parstriangularis','parsopercularis','parsorbitalis'};
roi2 = {'superiortemporal','middletemporal','inferiorparietal'};
d_thresh_mm = 30;

npairs = 0;
for i = (npairs+1):length(Subjects) % 1
    sid = Subjects{i};
    ecog = H5eeg(sprintf('%s/%s.h5',dir_h5,sid));
    fn_cache = sprintf('%s/xsub_out_%s_%i.mat',dir_cache,sid,metric_i);
    Ca = load(fn_cache);
    elec_rois = Ca.C.AtlLabels{atl_i};
    for j = 1:(ecog.n_bchan - 1)
        for k = (j+1):(ecog.n_bchan)
            b1c1 = ecog.bip(j,1);
            b1c2 = ecog.bip(j,2);
            b2c1 = ecog.bip(k,1);
            b2c2 = ecog.bip(k,2);
            roi_b1c1 = elec_rois{b1c1};
            roi_b2c2 = elec_rois{b1c2};
            roi_b2c1 = elec_rois{b2c1};
            roi_b2c2 = elec_rois{b2c2};
            b1hemi = Ca.C.EleHemi{b1c1};
            b2hemi = Ca.C.EleHemi{b2c1};
            d_b1b2 = sqrt(sum((ecog.bip(j,4:6) - ecog.bip(k,4:6)).^2));

            isin = any(strcmp(roi_b1c1,roi1)) && any(strcmp(roi_b2c1,roi2));
            isin_b = any(strcmp(roi_b1c1,roi2)) && any(strcmp(roi_b2c1,roi1));
            isfar = d_b1b2 >= d_thresh_mm;
            
            if ((isin || isin_b) && isfar)
                npairs = npairs + 1;
                pname = sprintf('%i_%s_%i%s_%i%s_%imm_%s_%s',npairs,sid,j,b1hemi,k,b2hemi,round(d_b1b2),roi_b1c1,roi_b2c1);
                fprintf('%s\n',pname);
                h = figure('visible','off');
                h2 = brainplot_one(sid,[b1c1 b2c1],'dk');
                copyobj(h2,gca);
                print(h,['figures/findex_1_anatomy/',pname],'-djpeg');
                close(h)
                %return
            end
            
        end
    end
    
    
end

