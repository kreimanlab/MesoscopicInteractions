function [ X ] = carveFromAnnot( ecog, startidx, endidx, h5File, achan )
%carveFromAnnot takes in indices for annotations and returns the eeg data

    % Unit conversion
    UNIT = 1000; % 1000 for uV, 1 for mV

    % Read EEG
    eeg = ecog.readEEG({startidx,endidx});
    
    % CAR
    eegcar = double(eeg.data - mean(eeg.data,2));
    eegcar = eegcar/UNIT;

    % Adjust for sub21
    %if (str2double(h5File(5:6)) == 43)
    %    eegcar = eegcar / 1000;
    %end
    
    % Return data and features
    if (achan > 0)
        % Single case
        %eegArt = eegcar(:,achan);
        
        % Extract features
        Fs = ecog.fs;
        inddata = eegcar(:,achan);
        X = [(max(inddata)-min(inddata)) / median(max(eegcar) - min(eegcar)), ...
            mean(abs(diff(inddata)/(1/Fs))) / mean(median(abs(diff(eegcar)/(1/Fs))))];
    else
        % Global case
        %eegArt = eegcar;
        
        % Extract features
        n_chan = ecog.n_chan;
        Adj = ones(n_chan,n_chan);
        for j = 2:n_chan
            for k = 1:(n_chan-1)
                Adj(j,k) = fastcorr2(eegcar(:,j),eegcar(:,k));
                Adj(k,j) = Adj(j,k);
            end
        end
        X = [sum(sum(Adj.^2))/(n_chan^2) max(max(eegcar) - min(eegcar))];
    end
    
    %figure; plot(eegcar(:,uchans(chan)));

end

