% 70 left after cherry picking electrodes and enforcing fraction of arts
% 6 subjects

close all;
clear;
rng shuffle;

% Thresholds
d_thresh_mm = 30;
fracart_thresh = 0.01;

% Constants
metric = 'pcBroadband';
metric_i = 1;
atl_i = 2;

% Dirs
dir_h5 = '/media/jerry/KLAB101/h5_notch20';
dir_cache = './cache';
dir_art = '/media/klab/KLAB101/h5_notch20/art_nosz';

% Go through all subjects
Subjects = {'sub1','sub2','sub3','sub5','sub6','sub7',...
    'sub10',         'sub13','sub14',         'sub17',         'sub20',...
    'sub21','sub22','sub23','sub24','sub26','sub28',         'sub34',...
    'sub35','sub36','sub39','sub41','sub42','sub44','sub45',...
    };

Roi2 = { {'LT22','LT29','LT21','LT28','PT4','PT3','PT11','PT12'} ;... % sub1
    {'LT23','LT30'} ;... % sub2
    {'PT5','PT6','PT7'} ;... % sub3
    {'LLT23','LLT31','LPT6','LPT14'} ;... % sub5
    {'LT15','LT23','LT32','LT40'} ;... % sub6
    {'LT30','LT22','LT14','LT6','LT5','LT13'} ;... % sub7
    {'MT6','MT5','PT16','PT8','PT7'} ;... % sub10
    {'AT16','MT7','MT8','MT15','MT16','PT16'} ;... % sub13
    {'RAT8','RMT8','RMT7','RPT8','ROT8','ROT7'} ;... % sub14
    {'LT31','PT8','PT7','PT15','PT14','LT23'} ;... % sub17
    {'PA9','PT16','PA1','PA2','PA3','PA4','PA5'} ;... % sub20
    {'LF1','LF9','LA9'} ;... % sub21
    {'AT22','AT23','AT31','PT16','PT7'} ;... % sub22
    {'LT5','AF8','LT4','AF7'} ;... % sub23
    {'LT14','LT13','LT6','LT5','LP26','LP25','LP18','LP17','LP19','LP10','LP9','LP11','LP12'} ;... % sub24
    {'RI14','RI7','RI15','RI8','RP8','RI16'} ;... % sub26
    {'AT22','AT31','PT8','PT16','PT15','AT30'} ;... % sub28
    {'LLG35','LLG43','LLG44','LLG45','LLG38','LLF46','LLG39','LLG40','LLG32'} ;... % sub34
    {'AT23','AT31','AT30'} ;... % sub35
    {'TM7','TM6','TM15','TM14','TP16','TP15'} ;... % sub36
    {'BT15','BT23','BT22','BT30','BT29'} ;... % sub39
    {'LLT44','LLT45','LLT46','LLT53','LLT54','LLT47','LLT55','LLT39'} ;... % sub41
    {'LT30','LT38','LT37','LT45','LT53','LT61','LT52'} ;... % sub42
    {'AT32','AT31','PT8','PT7','PT15','PT14'} ;... % sub44
    {'LF41','LF33','LF25','LF17','LF18','LF10','LF11'} ;... % sub45
    };


roi1 = {'parsorbitalis','parstriangularis','parsopercularis'};
roi2 = {'superiortemporal','middletemporal','inferiorparietal'};


npairs = 0;
for i = (npairs+1):length(Subjects) % 1
    sid = Subjects{i};
    
    % Init H5eeg
    ecog = H5eeg(sprintf('%s/%s.h5',dir_h5,sid));
    
    % Load cache
    fn_cache = sprintf('%s/xsub_out_%s_%i.mat',dir_cache,sid,metric_i);
    Ca = load(fn_cache);
    elec_rois = Ca.C.AtlLabels{atl_i};
    elec_name = Ca.C.EleLabels;
    Roi2_t = Roi2{i};
    
    % Load artifacts
    fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
    art_idx = h5read(fn_art,'/art_idx');
    
    bchan_count = 1;
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
            %name_b1c1 = elec_name{b1c1};
            
            b1hemi = Ca.C.EleHemi{b1c1};
            b2hemi = Ca.C.EleHemi{b2c1};
            d_b1b2 = sqrt(sum((ecog.bip(j,4:6) - ecog.bip(k,4:6)).^2));

            
            
            % Passing conditions
            isin = any(strcmp(roi_b1c1,roi1)) && any(strcmp(roi_b2c1,roi2));
            isin_b = any(strcmp(roi_b1c1,roi2)) && any(strcmp(roi_b2c1,roi1));
            isfar = d_b1b2 >= d_thresh_mm;
            is_roi2 = any(strcmp(Roi2_t,elec_name{b1c1})) || any(strcmp(Roi2_t,elec_name{b2c1}));
            
            
            %return
            if ( (((isin || isin_b) && isfar) && is_roi2) )
                
                % Calculate fraction of artifacts
                art_idx_l = art_idx(bchan_count,:);
                fracart = sum(art_idx_l~=0)/length(art_idx_l);
                is_fracart = fracart <= fracart_thresh;
                
                % Show
                if is_fracart
                    npairs = npairs + 1;
                    pname = sprintf('%i_%s_%i%s_%i%s_%imm_%s_%s_fracart-%i',...
                        npairs,sid,j,b1hemi,k,b2hemi,round(d_b1b2),roi_b1c1,...
                        roi_b2c1,round(1000*fracart));
                    fprintf('%s\n',pname);
                    h = figure('visible','off');
                    h2 = brainplot_one(sid,[b1c1 b2c1],'dk');
                    copyobj(h2,gca);
                    print(h,['figures/findex_2_fracart/',pname],'-djpeg');
                    close(h)
                end

            end
            
            bchan_count = bchan_count + 1;
        end
    end
    
    
end

