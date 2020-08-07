close all;
%clear;

if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        %resultsDir = '/media/klab/44/data/results';
        resultsDir = '/media/klab/D0BA8713BA86F56E/data/results';
        h5Dir = '/media/klab/44/h5';
    elseif strcmp(strip(hname),'ubuntu_1604')
        fprintf('Found host: ubuntu_1604\n')
        resultsDir = '/nas_share/RawData/data/results';
        h5Dir = '/nas_share/RawData/scripts/synth/out';
    else
        resultsDir = '/mnt/cuenap2/data/results';
        h5Dir = '/mnt/cuenap2/scripts/synth/out';
    end
end

if (~exist('A','var'))
    A = Analysis(resultsDir,h5Dir);
end

%Metrics = {'p','s','pP','sP','sd','st','sa','sb','sg',...
%    'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma',...
%    'scBroadband','scDelta','scTheta','scAlpha','scBeta','scGamma'};
%Metrics = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma','s','sP'};
Metrics = {'pcDelta','pcTheta','pcGamma'};
%Metrics = {'s','sP'};

%Metrics = {'s','sP','scBroadband','scGamma'};
%Metrics = {'scDelta','scTheta','scAlpha','scBeta'};
%Metrics = {'sd','st','sa','sb','sg'};

%Metrics = {'s','scBroadband'};
%Metrics = {'s'};
%Metrics = {'sP','scGamma'};
NotEnvMetrics = {'p','s','pP','sP'};



for metric_i = 1:length(Metrics)
    
    
    metric = Metrics{metric_i};
    %CT_all = cell(1,A.h5eeg.n_sub);
    AdjCT_thresh_all = cell(1,A.h5eeg.n_sub);
    %AT_all = cell(1,A.h5eeg.n_sub);

    for i = 1:A.h5eeg.n_sub
    %for i = 3
        tic;
        %close all;
        try
            sid = A.h5eeg.subjects{i};
            dmat = A.getDistanceMatrix(i);

            % Load CT
%             try
%                 load(sprintf('xsub/fitted/xsub-%s-%i',metric,i),'CT');
%             catch
            CT = A.getCT(i,metric,false);
%             end
            %CT_all{i} = CT;

            % Map to atlas
%             try
%                 load(sprintf('xsub/fitted/xsub-%s-%i',metric,i),'AT')
%             catch
            AT = A.getAT(i,CT,dmat);
%             end
            %AT_all{i} = AT;

            % CT statistical threshold
            AdjCT_thresh = binoinv((1-A.const.P_VALUE_CT),CT.n_w,A.const.P_VALUE)/CT.n_w;
            %AdjCT_thresh_all{i} = AdjCT_thresh;
            
            % save intermediate
            save(sprintf('xsub/fitted/xsub-%s-%i',metric,i),'-v6');
        catch
            fprintf('Skip metric %i, sub %i.\n',metric_i,i)
        end
        
        t_sing = toc;
        fprintf('[%i/%i] (%i of %i) ETA for %s: %.2f min.\n',metric_i,length(Metrics),i,A.h5eeg.n_sub,metric,(A.h5eeg.n_sub-i)*t_sing/60)
        %save(sprintf('xsub-%s-%i',metric,i),'-v7.3');
    end

   %save(sprintf('xsub-%s',metric),'-v7.3')
    
end

