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
Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};


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

