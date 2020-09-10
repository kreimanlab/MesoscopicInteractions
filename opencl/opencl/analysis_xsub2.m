close all;
clear;

Metrics = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma','s','sP'};
Atlases = 1:20;

for me_i = 1:length(Metrics)
    for at_i = 1:length(Atlases)
        metric = Metrics{me_i};
        ATL_I = Atlases(at_i);

        try
            %metric = 'pcBroadband';
            n_sub = 51; % 51
            FWER = 0.01;

            % 1 - Destrieux
            % 2 - Desikan
            % 3 - PALS Brodmann
            % 4 - 7 networks
            % 5 - 17 networks
            % 6 - HCP-MMP1
            % 7 - M132
            % 8 - FVE91
            % 9 - FVEall
            %ATL_I = 6;

            % constants

            %min_n_sub = 2;
            %min_n_epair_is_interact = 2; % within each patient

            % read fsaverage_sym lh flattened
            % [pa] = read_patch_rev('coreg/lh.full.flat.patch.3d');
            % k = convhull(pa.x,pa.y);
            % [v_pial,tri_pial] = read_surf('coreg/lh.pial');
            % 
            % plot(pa.x,pa.y,'black.');

            load('coreg/fsaverage_sym_flat');

            %tri = delaunay(pa.x,pa.y);
            %trisurf(tri,pa.x,pa.y,pa.z);
            %plot3(pa.x,pa.y,pa.z,'black.','markersize',0.1)
            %view(0,90);
            %n_sub = 20;

            for i = 1:n_sub
                xfname = sprintf('./xsub/fitted/xsub-%s-%i',metric,i);
                fprintf('[%i/%i] Loading %s.\n',i,n_sub,xfname)
                load(xfname);
                fprintf('\tBuilding xsub matrices.\n')

                if (i == 1)
                    % constants
                    n_rois = length(AT.AdjCT{ATL_I});
                    %n_rois = 180;
                    p_value_cs = (FWER/nchoosek(n_rois,2))/20;

                    % init
                    adjct_xsub = cell(size(AT.AdjCT{ATL_I}));
                    adjct_xsubM = cell(size(AT.AdjCT{ATL_I}));
                    adjct_sub = cell(size(AT.AdjCT{ATL_I}));
                    adjct_bchan = cell(size(AT.AdjCT{ATL_I}));
                    adjct_dist = cell(size(AT.AdjCT{ATL_I}));
                    %adjct_ = cell(size(AT.AdjCT{ATL_I}));
                end

                % For loop through roi pairs
                for j = 1:length(adjct_xsub)
                    for k = 1:length(adjct_xsub)
                        adjcts = AT.AdjCT{ATL_I}{j,k};

                        % if there is coverage for this roi pair
                        if (~isempty(adjcts))
                            % note AdjCT_thresh is dependent on the patient
                            adjcts(adjcts < AdjCT_thresh) = nan;
                            adjct_xsub{j,k} = [adjct_xsub{j,k}, adjcts];
                            adjct_xsubM{j,k} = [adjct_xsubM{j,k}, AT.AdjMag{ATL_I}{j,k}];
                            adjct_sub{j,k} = [adjct_sub{j,k}, AT.Sub{ATL_I}{j,k}];
                            adjct_bchan{j,k} = [adjct_bchan{j,k}, AT.Bchan{ATL_I}{j,k}];
                            adjct_dist{j,k} = [adjct_dist{j,k}, AT.Dist{ATL_I}{j,k}];
                        end

                    end
                end

                fprintf('\tdone.\n')

            end

            AdjCS = zeros(size(adjct_xsub));
            AdjCS_mag = zeros(size(adjct_xsub));
            n_covered = 0;
            n_sig = 0;
            for j = 1:length(adjct_xsub)
                for k = 1:length(adjct_xsub)
                    adjcts = adjct_xsub{j,k};
                    adjctsM = adjct_xsubM{j,k};
                    if (isempty(adjcts))
                        roi_is_interact = -1;
                    else

                        % If at least n electrode pairs is significant between rois
                        theta = binoinv(1-p_value_cs,length(adjcts),A.const.P_VALUE_CT);
                        % if number of significant is above threshold value
                        if (sum(~isnan(adjcts)) >= theta)
                            roi_is_interact = 1;
                            n_sig = n_sig + 1;
                        else
                            roi_is_interact = 0;
                        end

                        n_covered = n_covered + 1;
                        AdjCS_mag(j,k) = nanmean(adjctsM);
                    end
                    AdjCS(j,k) = roi_is_interact;
                end
            end
            pct_sig = 100*(n_sig)/(n_covered);
            fprintf('Percent significant: %.2f\n',pct_sig);
            %save(sprintf('xsub/xsub-%s_%s',metric,AT.P.AtlNames{ATL_I}),'-v7.3');
            save(sprintf('xsub/xsub2-%s_atl-%i',metric,ATL_I),'-v7.3');

            % plotting
            %imagesc(AdjCS);
            %colormap(corrcmap(100))
            %colorbar;
        catch
            fprintf('> skip: %s, %i\n',metric,ATL_I);
        end

    end
end
