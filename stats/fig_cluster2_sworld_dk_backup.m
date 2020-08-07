if strcmp(rois{i},'lateralorbitofrontal')
            rois2_pts(i,2) = rois2_pts(i,2) + 25; % left-right
            rois2_pts(i,3) = rois2_pts(i,3) - 6;
        elseif strcmp(rois{i},'parsorbitalis')
            rois2_pts(i,2) = rois2_pts(i,2) + 4; % left-right
            rois2_pts(i,3) = rois2_pts(i,3) + 2;
        elseif strcmp(rois{i},'superiorparietal')
            rois2_pts(i,2) = rois2_pts(i,2) + 5;
            rois2_pts(i,3) = rois2_pts(i,3) + 10;
        elseif strcmp(rois{i},'inferiorparietal')
            rois2_pts(i,2) = rois2_pts(i,2) - 8;
            rois2_pts(i,3) = rois2_pts(i,3) + 5;
        elseif strcmp(rois{i},'superiorfrontal')
            rois2_pts(i,2) = rois2_pts(i,2) + 30;
            rois2_pts(i,3) = rois2_pts(i,3) - 15;
        elseif strcmp(rois{i},'isthmuscingulate')
            rois2_pts(i,2) = rois2_pts(i,2) - 1;
            rois2_pts(i,3) = rois2_pts(i,3) + 7;
        elseif strcmp(rois{i},'fusiform')
            rois2_pts(i,2) = rois2_pts(i,2) - 12;
            rois2_pts(i,3) = rois2_pts(i,3) - 1;
        elseif strcmp(rois{i},'middletemporal')
            rois2_pts(i,2) = rois2_pts(i,2) + 28;
            rois2_pts(i,3) = rois2_pts(i,3) - 18;
        elseif strcmp(rois{i},'lingual')
            rois2_pts(i,2) = rois2_pts(i,2) - 4;
            rois2_pts(i,3) = rois2_pts(i,3) + 3;
        elseif strcmp(rois{i},'precuneus')
            rois2_pts(i,2) = rois2_pts(i,2) - 15;
            rois2_pts(i,3) = rois2_pts(i,3) + 8;
        elseif strcmp(rois{i},'lateraloccipital')
            rois2_pts(i,2) = rois2_pts(i,2) - 10;
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'rostralanteriorcingulate')
            rois2_pts(i,2) = rois2_pts(i,2) + 4;
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'caudalanteriorcingulate')
            rois2_pts(i,2) = rois2_pts(i,2) + 7;
            rois2_pts(i,3) = rois2_pts(i,3) - 4;
        elseif strcmp(rois{i},'frontalpole')
            rois2_pts(i,2) = rois2_pts(i,2) + 2;
            rois2_pts(i,3) = rois2_pts(i,3) + 5;
        elseif strcmp(rois{i},'parahippocampal')
            rois2_pts(i,2) = rois2_pts(i,2) - 6;
            rois2_pts(i,3) = rois2_pts(i,3) + 6;
        elseif strcmp(rois{i},'temporalpole')
            rois2_pts(i,2) = rois2_pts(i,2) + 0;
            rois2_pts(i,3) = rois2_pts(i,3) - 2;
        elseif strcmp(rois{i},'precentral')
            rois2_pts(i,2) = rois2_pts(i,2) - 10;
            rois2_pts(i,3) = rois2_pts(i,3) + 22;
        elseif strcmp(rois{i},'postcentral')
            rois2_pts(i,2) = rois2_pts(i,2) - 15;
            rois2_pts(i,3) = rois2_pts(i,3) + 30;
        elseif strcmp(rois{i},'supramarginal')
            rois2_pts(i,2) = rois2_pts(i,2) + 5;
            rois2_pts(i,3) = rois2_pts(i,3) - 2;
        end