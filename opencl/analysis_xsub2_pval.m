close all;
clear;

if ismac
    resultsDir_pval = '/Volumes/RawData/data/results';
    h5Dir_pval = '/Volumes/RawData/scripts/synth/out';
elseif isunix
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        resultsDir_pval = '/media/klab/44/data/results';
        h5Dir_pval = '/media/klab/44/h5';
    else
        resultsDir_pval = '/mnt/cuenap2/data/results';
        h5Dir_pval = '/mnt/cuenap2/scripts/synth/out';
    end
end

A_sav = Analysis(resultsDir_pval,h5Dir_pval);

% whether to print adjacency matrices of all individual subjects
PLOT_ALL_IND_SUB = true;
% whether to print adjacency matrix for one randomly chosen subject
PLOT_LUCKY_SUB = true;

%Metrics = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma','s','sP'};
%Atlases = 1:20;
Metrics_pval = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};
Atlases_pval = 2;
for me_i = 1:length(Metrics_pval)
    for at_i = 1:length(Atlases_pval)
        metric_pval = Metrics_pval{me_i};
        ATL_I_CONST = Atlases_pval(at_i);

        try
            %metric = 'pcBroadband';
            n_sub = 51; % 51
            

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
            lucky_sub = randi([1 n_sub]);
            for i_sub_pval = 1:n_sub
                xfname_xsub = sprintf('./xsub/fitted/xsub-%s-%i',metric_pval,i_sub_pval);
                fprintf('[%i/%i] Loading %s.\n',i_sub_pval,n_sub,xfname_xsub)
                load(xfname_xsub,'-regexp','^(?!CT)\w');
                A = A_sav;
                % --- change p-value parameters ---------------------------
                A.const.P_VALUE = 0.05; % fixd based on looking at correlations over time
                A.const.P_VALUE_CT = 0.05; % 176 max electrodes/51 patients  
                %A.const.P_VALUE_CS = (0.05/nchoosek(31,2)); % 180 rois /20 atlases
                FWER = 0.05; % determines p_value_cs
                % ---------------------------------------------------------
                h5inf = h5info(sprintf('%s/%s_graph-%s.h5',resultsDir_pval,sid,metric_pval),'/R');
                n_w = h5inf.Dataspace.Size(2);
                AdjCT_thresh = binoinv((1-A.const.P_VALUE_CT),n_w,A.const.P_VALUE)/n_w;
                
                %fprintf('\tBuilding xsub matrices.\n')

                if (i_sub_pval == 1)
                    % constants
                    n_rois = length(AT.AdjCT{ATL_I_CONST});
                    %n_rois = 180;
                    
                    % ---------------------------------------------------------
                    %p_value_cs = (FWER/nchoosek(n_rois,2))/20;
                    p_value_cs = (FWER/nchoosek(n_rois,2));

                    % init
                    adjct_xsub = cell(size(AT.AdjCT{ATL_I_CONST}));
                    adjct_xsubM = cell(size(AT.AdjCT{ATL_I_CONST}));
                    adjct_sub = cell(size(AT.AdjCT{ATL_I_CONST}));
                    adjct_bchan = cell(size(AT.AdjCT{ATL_I_CONST}));
                    adjct_dist = cell(size(AT.AdjCT{ATL_I_CONST}));
                    %adjct_ = cell(size(AT.AdjCT{ATL_I}));
                end

                % For loop through roi pairs
                for j = 1:length(adjct_xsub)
                    for k = 1:length(adjct_xsub)
                        adjcts = AT.AdjCT{ATL_I_CONST}{j,k};

                        % if there is coverage for this roi pair
                        if (~isempty(adjcts))
                            % note AdjCT_thresh is dependent on the patient
                            adjcts(adjcts < AdjCT_thresh) = NaN;
                            adjct_xsub{j,k} = [adjct_xsub{j,k}, adjcts];
                            adjct_xsubM{j,k} = [adjct_xsubM{j,k}, AT.AdjMag{ATL_I_CONST}{j,k}];
                            adjct_sub{j,k} = [adjct_sub{j,k}, AT.Sub{ATL_I_CONST}{j,k}];
                            adjct_bchan{j,k} = [adjct_bchan{j,k}, AT.Bchan{ATL_I_CONST}{j,k}];
                            adjct_dist{j,k} = [adjct_dist{j,k}, AT.Dist{ATL_I_CONST}{j,k}];
                        end

                    end
                end

                %fprintf('\tdone.\n')

                % --- plot individual adjacency matrix --------------------
                if (((i_sub_pval == lucky_sub) && (PLOT_LUCKY_SUB))|| (PLOT_ALL_IND_SUB))
                    h_t = figure;
                    fig_w = 8.5;
                    fig_h = 11.0;
                    set(h_t,'Position',[0 0 fig_w*100 fig_h*100])
                    set(h_t, 'PaperUnits', 'Inches')
                    set(h_t, 'PaperPosition', [0 0 fig_w fig_h])    % can be bigger than screen 
                    set(h_t, 'PaperSize', [fig_w fig_h])    % Same, but for PDF output
                    ax1 = axes('Parent',h_t,'Units','Normalized','Position',[0.18 0.35 0.55 0.4]);
                    CT = A.getCT(i_sub_pval,metric_pval,false);
                    pct_sig_sub = sum(sum(~isnan(CT.AdjCT) & (CT.AdjCT ~= 0)))/numel(CT.AdjCT);
                    
                    %imagesc(CT.AdjCT);
                    %ax1 = axes('Parent',h_t,'Units','Normalized','Position',[0.18 0.35 0.55 0.4]);
                    ax1 = gca;
                    title_fig = sprintf('sub p=%.6f pCT=%.6f pCS=%.6f (%.2f %% sig)',...
                        A.const.P_VALUE,A.const.P_VALUE_CT,p_value_cs,100*pct_sig_sub);
                    % Correct for multiple hypothesis testing at individual level
                    AdjCT_thresh_ind = binoinv((1-(A.const.P_VALUE_CT/nchoosek(A.h5eeg.n_bchan{i_sub_pval},2))),n_w,(A.const.P_VALUE))/n_w;
                    fprintf('AdjCT_thresh: %.6f, AdjCT_thresh_ind: %.6f\n',AdjCT_thresh,AdjCT_thresh_ind)
                    ax2 = A.plotAdjElectrodeDistThresh(ax1, i_sub_pval, CT.AdjCT, ...
                        title_fig, dmat, CT.AdjCT, AdjCT_thresh_ind);
                    
                    %colormap(corrcmap(100));
                    c = colorbar(ax1,'Location','manual');
                    set(c,'Position',[0.86 0.35 0.02 0.4]);
                    set(c,'FontSize',8);
                    %colorbar;
                    
                    A.saveFig(h_t,sprintf('./figures/xsub_coh/test_sub-%s_%s-%i_%i_%i_%i',...
                        sid,metric_pval,ATL_I_CONST,1e6*A.const.P_VALUE,1e6*A.const.P_VALUE_CT,1e6*p_value_cs));
                    close(h_t);
                end
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
            %clear A_sav;
            clear CT;
            save(sprintf('xsub/xsub2-%s_atl-%i_AdjCTthresh2-%i',metric_pval,ATL_I_CONST,round(AdjCT_thresh*1e6)),'-v7.3');

            % plotting
            h_t = figure;
            imagesc(AdjCS);
            colormap(corrcmap(100))
            colorbar;
            title(sprintf('xsub p=%.6f pCT=%.6f pCS=%.6f (%.2f %% sig)',...
                A.const.P_VALUE,A.const.P_VALUE_CT,p_value_cs,pct_sig))
            A.saveFig(h_t,sprintf('./figures/xsub_coh/test_xsub_%s-%i_%i_%i_%i',...
                metric_pval,ATL_I_CONST,1e6*A.const.P_VALUE,1e6*A.const.P_VALUE_CT,1e6*p_value_cs));
            close(h_t);
        catch e
            fprintf('> skip: %s, %i\n',metric_pval,ATL_I_CONST);
            rethrow(e);
%             fprintf(2,'W> %s\n',e.identifier);
%             fprintf(2,'W> %s\n',e.message);
        end

    end
end

