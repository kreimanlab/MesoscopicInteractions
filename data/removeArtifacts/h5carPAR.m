function [ ] = h5carPAR( infilename, step )

    %STEP = 1/60; % minutes of data to process at once
    %step = round(60*STEP*ecog.fs);
    ecog = H5eeg(infilename);
    if (step > ecog.n_samples)
        fprintf('[!] Warning: step exceeds number of samples, setting step to n_samples.\n');
        %return
        step = ecog.n_samples;
    end
    
    try
        h5create(infilename,'/h5eeg/cavg',[1 Inf],'ChunkSize',[1 ecog.fs])
    catch
        fprintf('[!] Error: file already has a common average array, exiting h5car.\n');
        return
    end
    
    %c = 1;
    startIdxs = 1:step:ecog.n_samples;
    Nsi = length(startIdxs);
    parfor i = 1:Nsi
    %for start_idx = 1:step:ecog.n_samples
        start_idx = startIdxs(i);
        tic;
        end_idx = start_idx + step - 1;
        if (end_idx > ecog.n_samples)
            end_idx = ecog.n_samples;
        end
        %fprintf('(*) start: %i\tend: %i (of %i)\n',start_idx,end_idx,ecog.n_samples)
        
        % --- Main CAR code ---
        % Write CAR
        eeg = ecog.readEEG({start_idx end_idx});
        cavg = mean(eeg.data,2);
        h5write(infilename,'/h5eeg/eeg',(eeg.data-cavg)',[1 start_idx],[ecog.n_chan (end_idx-start_idx+1)]);
        
        % Save common average as new dataset
        h5write(infilename,'/h5eeg/cavg',cavg',[1 start_idx],[1 (end_idx-start_idx+1)])
        % --- Main CAR code ---
        
%         if (i == 1)
%             tprog = toc;
%         else
%             tprog = (tprog + toc)/2;
%         end
        fprintf('\t(*) %.2f%%, %.2f hrs left\n',100*i/ceil(ecog.n_samples/step),toc*(ceil(ecog.n_samples/step)-i)/3600);
        %c = c + 1;
    end
    
end