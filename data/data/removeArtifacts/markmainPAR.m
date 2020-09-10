% Last Edited - Feb 9, 2017
% Jerry Wang :: jwang04@g.harvard.edu

function markmain(infilearg)
% markmain() makrs artifacts and optionally performs common average
% referencing

    % Perform Common Average Referencing (CAR) ?
    DO_CAR = true;

    % Print splash
    fprintf('#------------------------------------------------------------#\n')
    fprintf('|   Artifact Removal Tool                                    |\n')
    fprintf('#------------------------------------------------------------#\n')

    load('QDAartifacts.mat')
    WINDOW = 1; %seconds
 
    % CAR filter
    fprintf('[$] CAR filtering...\n')
    if (DO_CAR)
        h5carPAR(infilearg, 100000);
        %h5car(infilearg, 100000);
        fprintf('[!] Done.\n');
    else
        fprintf('[!] Manually disabled. Edit markmain.m to turn on.\n');
    end

    % inputs
    infile = infilearg;
    outfile = strsplit(infile,'.hdf5');
    outfile = outfile{1};
    outfile = ['qda_',outfile,'.annot'];

    % load QDA model
    Ks = MdlQuadraticSing.Coeffs(2,1).Const+offsetSing;
    Ls = MdlQuadraticSing.Coeffs(2,1).Linear;
    Qs = MdlQuadraticSing.Coeffs(2,1).Quadratic;
    fsing = @(x1,x2) Ks + Ls(1)*x1 + Ls(2)*x2 + Qs(1,1)*x1.^2 + ...
        (Qs(1,2)+Qs(2,1))*x1.*x2 + Qs(2,2)*x2.^2;

    K = MdlQuadraticGlob.Coeffs(2,1).Const+offsetGlob;
    L = MdlQuadraticGlob.Coeffs(2,1).Linear;
    Q = MdlQuadraticGlob.Coeffs(2,1).Quadratic;
    fglob = @(x1,x2) K + L(1)*x1 + L(2)*x2 + Q(1,1)*x1.^2 + ...
        (Q(1,2)+Q(2,1))*x1.*x2 + Q(2,2)*x2.^2;

    % Read file
    fprintf('[-] Reading file %s...\n',infile)
    ecog = H5eeg(infile);

    % Precalculations
    width = round(WINDOW*ecog.fs);
    markedN = 0;
    Starts = 1:width:ecog.n_samples;
    cprog = round(length(Starts)/1000);
    N_SAMPLES = ecog.n_samples;
    N_CHAN = ecog.n_chan;
    counter = struct();
    markedSamps = 0;
    %progress = [];
    total = length(Starts);
    %for startidx = Starts
    
    % Create artifact dataset
    h5create(infile,'/h5eeg/artifacts',[1,total]);
    h5writeatt(infile,'/h5eeg/artifacts','width',width);
    %j = 0;
    parfor si = 1:total
        tic;
        
        % Adjust for the last chunk
        startidx = Starts(si);
        endidx = startidx + width - 1;
        if (endidx > N_SAMPLES)
            endidx = N_SAMPLES;
        end

        % Handle global events
        globExist = false;
        [x_glob] = carveFromAnnot(ecog, startidx, endidx, infile, 0);
        if (fglob(x_glob(1),x_glob(2)) > 0)
            % Write global artifact marks
            startGlob = [1 startidx]; % write position
            countGlob = [N_CHAN width]; % write block count
            artGlob = nan(N_CHAN,width); % block of NaNs to write
            h5write(infile,'/h5eeg/eeg',artGlob,startGlob,countGlob);
            h5write(infile,'/h5eeg/artifacts',ecog.n_chan,[1 si],[1 1]);
            markedSamps = markedSamps + width;
            
            % --- PRINT ---
            %fprintf('[%.2f%%\t] found global artifact %i\n',100*markedSamps/endidx,startidx);
            
            %eeg = ecog.readEEG({startidx,endidx});
            %for i = startidx:endidx
            %    fprintf(of,'%f\n',i);
            %end
            globExist = true;
        end

        singExist = 0;
        % Handle single events
        for eSing = 1:N_CHAN
            [x_sing] = carveFromAnnot(ecog, startidx, endidx, infile, eSing);
            if (fsing(x_sing(1),x_sing(2)) > 0)
                % Write single artifact marks
                startSing = [eSing startidx]; % write position
                countSing = [1 width]; % write block count
                artSing = nan(1,width); % block of NaNs to write
                h5write(infile,'/h5eeg/eeg',artSing,startSing,countSing);
                if (~globExist)
                    markedSamps = markedSamps + width/N_CHAN;
                end
                
                % --- PRINT ---
                %fprintf('[%.2f%%\t] found single artifact %i [%i]\n',100*markedSamps/endidx,startidx,eSing);
                
                %for i = startidx:endidx
                %    fprintf(of,'%i\n%f\n',(-1)*eSing,i);
                %end
                singExist = singExist + 1;
            end
        end
%         if ((singExist > 0) && (~globExist))
%             h5write(infile,'/h5eeg/artifacts',singExist,[1 si],[1 1]);
%         elseif ((singExist == 0) && (~globExist))
%             % Write blank to artifacts data if no events are found
%             h5write(infile,'/h5eeg/artifacts',singExist,[1 si],[1 1]);
%         end
        if (~globExist)
            h5write(infile,'/h5eeg/artifacts',singExist,[1 si],[1 1]);
        end
        
%         if (si == 1)
%             tprog = toc;
%         else
%             tprog = (tprog + toc)/2;
%         end
        
        % Print progress
%         if (mod(si-1,cprog) == 0)
%             fprintf('[%.2f%% marked] file: %s\tfinished %i of %i (%.2f%% finished), %.2f hrs left\n',...
%                 100*markedSamps/endidx,ecog.filename,startidx,N_SAMPLES,100*startidx/N_SAMPLES,tprog*(total-si)/3600);
%         end
        if (mod(si-1,cprog) == 0)
            fprintf('[%s]\tfinished %i of %i (%.2f%% finished), %.2f hrs left\n',...
                ecog.filename,startidx,N_SAMPLES,100*startidx/N_SAMPLES,toc*(total-si)/3600);
        end


        %j = j + 1;
    end

    markedFrac = 100*markedSamps/N_SAMPLES;
    fprintf('[!] Total data marked: %.2f%%\n', markedFrac);
    outtxt = strsplit(ecog.filename,'.hdf5');
    outtxt = [outtxt{1},'.txt'];
    fprintf('(*) Saving summary to: %s\n',outtxt);
    % write receipt
    ofID = fopen(outtxt,'w');
    fprintf(ofID,'[!] Total data marked: %.2f%%\n', markedFrac);
    fclose(ofID);
    
    %fclose(of);
    fprintf('[!] All Done.\n')

    exit()
end
