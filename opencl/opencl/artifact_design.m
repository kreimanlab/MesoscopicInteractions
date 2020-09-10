close all;

if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix    
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        resultsDir = '/media/klab/44/data/results';
        h5Dir = '/media/klab/44/h5';
    elseif strcmp(strip(hname),'ubuntu_1604')
        resultsDir = '/nas_share/RawData/data/results';
        h5Dir = '/nas_share/RawData/scripts/synth/out';
    else
        resultsDir = '/mnt/cuenap2/data/results';
        h5Dir = '/mnt/cuenap2/scripts/synth/out';
    end
end

clear A;
if (~exist('A','var'))
    A = Analysis(resultsDir,h5Dir);
end

w = 1;
n_rand = 1;
for i = 1:n_rand
    
    if (~ exist('Vchan','var'))
        sub_i = randi([1 A.h5eeg.n_sub]);
        chan_i = randi([1 A.h5eeg.n_chan{sub_i}]);
        Vchan = h5read(A.h5eeg.filenames{sub_i},'/h5eeg/eeg',[chan_i 1],[1 A.h5eeg.n_samples{sub_i}]);
        
        % normalize
        Vchan = (Vchan - mean(Vchan))/std(Vchan);
    end
    %samp = randi([1 A.h5eeg.n_samples{sub_i}-round(w*A.h5eeg.fs{sub_i})]);
    Starts = 1:round(w*A.h5eeg.fs{sub_i}):A.h5eeg.n_samples{sub_i};
    Std = zeros(1,length(Starts));
    Max = zeros(1,length(Starts));
    c = 1;
    for j = Starts
        start_i = j;
        end_i = start_i + round(w*A.h5eeg.fs{sub_i}) - 1;
        if (end_i > A.h5eeg.n_samples{sub_i})
            end_i = A.h5eeg.n_samples{sub_i};
        end
        Std(c) = std(Vchan(start_i:end_i));
        Max(c) = max(Vchan(start_i:end_i));
        c = c + 1;
    end
    T = linspace(0,length(Vchan)/(A.h5eeg.fs{sub_i} * 3600 * 24),length(Vchan)); %days
    plot(T,Vchan,'black.','MarkerSize',1);
    axis tight;
    xlabel('Days')
    ylabel('Z-score')
%     fprintf('%.4f\n',mean(Std));
%     [freq,x] = hist(Max,10000);
%     plot(x,freq/trapz(x,freq),'black-'); hold on;
%     [freqD,xD] = hist(Vchan,10000);
%     plot(xD,freqD/trapz(xD,freqD),'green-'); hold on;
%     plot(x,normpdf(x,mean(Vchan),std(Vchan)),'red'); hold on;
end

