classdef H5eeg
    %[H5eeg]
    %   
    %   Implements a hdf5 data format for the efficient storage and
    %   retrieval of electrophysiological recordings.
    %
    %   Last Edited: Jerry Wang - jwang04@g.harvard.edu - 2017 Mar 10
    %
    %[Example Usage]
    %
    %   patient1 = H5eeg('sub1.hdf5')
    %   
    %   patient1.readAttributes() returns a struct of descriptive
    %   information about the data.
    %
    %   patient1.readEEG(<double>) will read at elapsed time
    %
    %   patient1.readEEG({<int>, <int>}) will read from first int index to
    %   second int index
    %
    %   patient1.readEEG({<int>, <int>, <int>, <int>},{<int>, <int>}) will 
    %   read from channels indicated by the second cell at beginning
    %   samples specified by the first and third ints of the first cell,
    %   with lengths specified by the second and fourth ints of the first
    %   cell
    %
    %   patient1.readEEG(<double>, <double>) will read the entire length of
    %   data for 2 channel indices.
    %
    
    properties
        filename; % file name
        n_chan; % number of channels
        n_samples; % number of samples
        start_time; % starting serial sample time
        end_time; % ending serial sample time
        fs; % sampling rate
        bip; % bipolar montage
        n_bchan; % number of bipolar channels
        n_comb; % n_bchan choose 2
    end
    
    %
    % H5eeg.findIndexWithElapsedTime
    %
    % @param etime
    %
    % @return index
    %
    % Finds the index within the data that contains the given elapsed time
    %
    methods (Hidden)
        function index = findIndexWithElapsedTime (self, etime)
            aux = h5read(self.filename,'/h5eeg/aux');
            labels = h5readatt(self.filename,'/h5eeg/aux','labels');
            isElapsedTime = strcmp(labels,'Elapsed time');
            if (any(isElapsedTime))
                preIndex = find(aux(isElapsedTime,:) == etime);
                if (~isempty(preIndex))
                    index = preIndex;
                else
                    fprintf(2, '[!] elapsed time was not found\n');
                end
            else
                fprintf(2,'[!] "Elapsed time" label not found in /h5eeg/aux\n');
            end
        end
    end
    
    methods
        %
        % H5eeg Constructor
        %
        % @param fn
        %
        % New object containing filename of the h5eeg file to manipulate
        %
        function self = H5eeg (fn)
            self.filename = fn;
            self.n_chan = NaN;
            self.n_samples = NaN;
            self.fs = NaN;
            self.start_time = NaN;
            self.end_time = NaN;
            self.bip = NaN;
            self.n_bchan = NaN;
            self.n_comb = NaN;
            if (exist(fn,'file'))
                %info = h5info(self.filename,'/h5eeg/eeg');
                %info = h5info(self.filename,'/h5eeg/aux');
                %self.n_chan = length(h5readatt(self.filename,'/h5eeg/eeg','gains'));
                %self.n_samples = info.Dataspace.Size(2);%info.Dataspace.Size(2);
                self.n_chan = h5readatt(self.filename,'/h5eeg/eeg','n_chan');
                self.n_samples = h5readatt(self.filename,'/h5eeg/eeg','n_samples');
                self.fs = h5readatt(self.filename,'/h5eeg/eeg','rate');
                tstamp = h5readatt(self.filename,'/h5eeg','timestamp');
                self.start_time = (tstamp + datenum(1970,1,1)*86400)*self.fs;
                self.end_time = self.start_time + self.n_samples - 1;
                self.bip = h5readatt(self.filename,'/h5eeg/eeg','bip');
                [n_bchan,~] = size(self.bip);
                self.n_bchan = n_bchan;
                self.n_comb = nchoosek(n_bchan,2);
            end
        end
        
        %
        % H5eeg.readAttributes
        %
        % @return att
        %
        % Reads attributes of the main group '/h5eeg' and saves it in a
        % single struct. The output attributes are: subject, experiment,
        % datetime, and timestamp.
        %
        function att = readAttributes(self)
            att.subject = h5readatt(self.filename,'/h5eeg','subject');
            att.experiment = h5readatt(self.filename,'/h5eeg','experiment');
            att.datetime = h5readatt(self.filename,'/h5eeg','datetime');
            att.timestamp = h5readatt(self.filename,'/h5eeg','timestamp');
            att.events = h5read(self.filename,'/h5eeg/events');
        end
        
        %
        % H5eeg.readEEG
        %
        % @return eeg
        %
        % Reads in electrophysiological recordings in /h5eeg/eeg
        %
        function eeg = readEEG(self, varargin)
            %eeg.labels = h5readatt(self.filename,'/h5eeg/eeg','labels');
            %eeg.offsets = h5readatt(self.filename,'/h5eeg/eeg','offsets');
            %eeg.gains = h5readatt(self.filename,'/h5eeg/eeg','gains');
            argc = length(varargin);
            switch (argc)
                case 0
                    eeg.data = h5read(self.filename,'/h5eeg/eeg');
                case 1
                    arg1in = varargin{1};
                    switch class(arg1in)
                        case 'double'
                            idx = self.findIndexWithElapsedTime(arg1in);
                            eeg.data = h5read(self.filename,'/h5eeg/eeg',[1 idx],[self.n_chan 1]');
                        case 'cell'
                            eeg.data = h5read(self.filename,'/h5eeg/eeg',[1 arg1in{1}],[self.n_chan (arg1in{2}-arg1in{1}+1)])';
                        otherwise
                            fprintf(2,'[!] unrecognized type: %s\n',...
                                class(arg1in));
                    end
                case 2
                    if (isa(varargin{1},'double') && isa(varargin{2},'double'))
                        % Two numbers used for two channels
                        eeg.data(:,1) = h5read(self.filename,'/h5eeg/eeg',[varargin{1} 1],[1 self.n_samples])';
                        eeg.data(:,2) = h5read(self.filename,'/h5eeg/eeg',[varargin{2} 1],[1 self.n_samples])';
                    elseif (isa(varargin{1},'cell') && isa(varargin{2},'cell'))
                        eeg.data(:,1) = h5read(self.filename,'/h5eeg/eeg',[varargin{2}{1} varargin{1}{1}],[1 varargin{1}{2}])';
                        eeg.data(:,2) = h5read(self.filename,'/h5eeg/eeg',[varargin{2}{2} varargin{1}{3}],[1 varargin{1}{4}])';
                    end
                otherwise
                    fprintf(2,'[!] unrecognized number of inputs: %i\n'...
                        ,argc);
            end
        end
        
        % Get distance matrix
        function Dmat = getDist(self)
            Dmat = zeros(self.n_bchan,self.n_bchan);
            for i2 = 1:(length(Dmat)-1)
                for j2 = (i2+1):length(Dmat)
                    a = reshape(self.bip(i2,4:end),2,[]);
                    b = reshape(self.bip(j2,4:end),2,[]);
                    dist2 = zeros(2,2);
                    for a1 = 1:2
                        for b1 = 1:2
                            dist2(a1,b1) = sqrt(sum((a(a1,:) - b(b1,:)).^2));
                        end
                    end
                    d = min(dist2(:));
                    Dmat(i2,j2) = d;
                    Dmat(j2,i2) = d;
                end
            end
        end
        
        % Returns the theoretical distances between electrodes
        function Dmat = getDist2(self)
            Dmat = self.getDist();
            try
                Dmat2 = h5read(self.filename,'/h5eeg/elec_recon/Dmat');
            catch
                fprintf('[!] Could not read /h5eeg/elec_recon/Dmat. This is an older version of the data file.\n')
            end
            Dmat(Dmat2 < 0) = -1;
        end
        
        % Return absolute time string from sample number
        function tiout = getTimeFromSample(self, stampsDir, sampleNum)
            LINES_HEADER = 15;
            files = h5readatt(self.filename,'/h5eeg','files');
            n_samples_l = h5readatt(self.filename,'/h5eeg','n_samples');
            linenum = 0;
            
            % find file that corresponds to sampleNum
            fi_i = find(sampleNum <= cumsum(n_samples_l),1);
            subfile = files{fi_i};
            sid = self.filename((regexp(self.filename,'(?<=/m0)(.*)(?=h5)')-2):(end-3));
            
            % convert stamp index to line number in text file
            stamp = double(h5read(self.filename,'/h5eeg/aux',[1 sampleNum],[1 1]));
            stamps_fn = sprintf('%s/%s_%s_stamps.txt',stampsDir,sid,subfile);
            times_fn = sprintf('%s/%s_%s_times.txt',stampsDir,sid,subfile);
            fprintf('%s\n',stamps_fn)
            
            % init starting stamp number
            f = fopen(stamps_fn,'r');
            for i = 1:LINES_HEADER
                tline = fgetl(f);
                linenum = linenum + 1;
            end
            fprintf('--- H5eeg.getTimeFromSample BEGIN ---\n')
            tline = fgetl(f);
            linenum = linenum + 1;
            start_stamp = str2double(tline);
            
            % beginning guess
            stamp_ln = stamp - start_stamp + 1;
            lines_t = stamp_ln + LINES_HEADER;
            linenum = lines_t;
            [~,tline2] = system(sprintf('sed -n %i,%ip %s',lines_t,lines_t,stamps_fn));
            est_stamp = str2double(tline2);
            if ((est_stamp < 0) || (round(est_stamp) ~= est_stamp))
                fprintf(2,'ERROR: invalid stamp %s in %s\n',tline2,stamps_fn);
                return;
            end
            
            % check guess
            stamp_guess = str2double(tline2);
            lines_times = 0;
%             fprintf('DEBUG\nStart stamp: %i\tGoal stamp: %i\tCurrent: %s\tLine num:%i\n',start_stamp,stamp,tline2,linenum);
            if (stamp_guess == stamp)
                lines_times = linenum;
%             elseif (stamp_guess < stamp)
%                 fprintf(2,'W: finished ascending search\nStart stamp: %i\tGoal stamp: %i\tCurrent: %s\tLine num:%i\n',start_stamp,stamp,tline2,linenum);
%                 lines_times = linenum;
%             elseif (stamp_guess > stamp)
%                 fprintf(2,'W: finished descending search\nStart stamp: %i\tGoal stamp: %i\tCurrent: %s\tLine num:%i\n',start_stamp,stamp,tline2,linenum);
%                 lines_times = linenum;
            else
%                 fprintf('DEBUG stage 2 start\nStart stamp: %i\tGoal stamp: %i\tCurrent: %s\tLine num:%i\n',start_stamp,stamp,tline2,linenum);
%                 fprintf('linenum: %i\n',linenum);
%                 fprintf('stamp_guess - stamp: %i\n',stamp_guess - stamp)
                stamp_guess_prev = stamp_guess;
                linenum = linenum - (stamp_guess - stamp);
%                 fprintf('linenum: %i\n',linenum);
                [~,tline2] = system(sprintf('sed -n %i,%ip %s',linenum,linenum,stamps_fn));
                stamp_guess = str2double(tline2);
                lines_times = linenum;
                fprintf(2,'[!] Corrected stamp %i -> %i, error: %i\n',stamp_guess_prev,stamp_guess,stamp_guess - stamp);
%                 fprintf('DEBUG stage 2 end\nStart stamp: %i\tGoal stamp: %i\tCurrent: %s\tLine num:%i\n',start_stamp,stamp,tline2,linenum);
            end
            
            % debug
            %fprintf('Start stamp: %i\tGoal stamp: %i\tCurrent: %s\tLine num:%i\n',start_stamp,stamp,tline2,linenum);
            
            fclose(f);

            [~,tline3] = system(sprintf('sed -n %i,%ip %s',lines_times,lines_times,times_fn));

            tiout = tline3;
            
            fprintf('--- H5eeg.getTimeFromSample EXIT_SUCCESS ---\n')
        end
        
        % Return video file, if exists, given time string
        function vidout = getVideo(self, videosDir, sampleNum, timein)
            
            % get subject and subfile strings
            files = h5readatt(self.filename,'/h5eeg','files');
            n_samples_l = h5readatt(self.filename,'/h5eeg','n_samples');
            fi_i = find(sampleNum <= cumsum(n_samples_l),1);
            subfile = files{fi_i};
            sid = self.filename((regexp(self.filename,'(?<=/m0)(.*)(?=h5)')-2):(end-3));
            
            % try to find video file
            indir = videosDir;
            d = dir(indir);
            disdir = [d.isdir];
            dname = {d.name};
            L = dname(~ disdir);
            outL = {};
            j = 1;
            for i = 1:length(L)
                if (contains(L{i},sid) && contains(L{i},subfile))
                    outL{j} = sprintf('%s/%s',indir,L{i});
                    j = j + 1;
                end
            end
            
            % default output
            vidout = '';
            if (length(outL) == 1)
                % calculate datenum of query time
                %datefmt = 'mmmm dd, yyyy HH:MM:SS.FFF AM';
                qtime = datenum(timein);
                
                % find query time in video spec file
                fid = fopen(outL{1},'r');
                D = textscan(fid,'%s','Delimiter','\n');
                D = D{1};
                fclose(fid);
                offset_secs = 0;
                for j = 1:length(D)
                    L = strsplit(D{j},'\t');
                    vfile = L{1};
                    start_t = datenum(L{2});
                    end_t = datenum(L{3});
                    if ( (qtime >= start_t) && (qtime <= end_t) )
                        offset_secs = (qtime - start_t)*(24*60*60);
                        break;
                    end
                end
                vidout = [sid,'_',subfile,'_',vfile((end-7):end),',',sprintf('%.0f',offset_secs)];
            elseif (length(outL) > 1)
                fprintf(2,'[E] more than 1 video file found for: %s, %s\n',sid,subfile);
                return;
            else
                fprintf(2,'[E] video file not found for: %s, %s\n',sid,subfile);
                return;
            end
        end
        
        % Return sample number from absolute time string
        function sout = getSampleFromTime(self, stampsDir, timeStr)
            LINES_HEADER = 15;
            files = h5readatt(self.filename,'/h5eeg','files');
            n_samples_l = h5readatt(self.filename,'/h5eeg','n_samples');
            ds_factors = h5readatt(self.filename,'/h5eeg','ds_factors');
            timestamp_l = double(h5readatt(self.filename,'/h5eeg','timestamp'));
            
            % convert to days
            timestamp_l = timestamp_l / (24*60*60) + datenum(1970,1,1);
            
            % convert to date num
            time_dn = datenum(timeStr);
            
            % find file that corresponds to timeStr
            fi_i = sum(time_dn >= timestamp_l);
            %fi_i = find(time_dn <= cumsum(timestamp_l),1);
            %fi_i = find(sampleNum <= cumsum(n_samples_l),1);
            subfile = files{fi_i};
            sid = self.filename((regexp(self.filename,'(?<=/m0)(.*)(?=h5)')-2):(end-3));
            
            % convert time string index to line number in text file
            stamps_fn = sprintf('%s/%s_%s_stamps.txt',stampsDir,sid,subfile);
            times_fn = sprintf('%s/%s_%s_times.txt',stampsDir,sid,subfile);
            fprintf('%s\n',times_fn)
            
            % init starting time
            f = fopen(times_fn,'r');
            linenum = 0;
            for i = 1:LINES_HEADER
                tline = fgetl(f);
                linenum = linenum + 1;
            end
            fprintf('--- H5eeg.getSampleFromTime BEGIN ---\n')
            tline = fgetl(f);
            linenum = linenum + 1;
            start_time_l = datenum(tline);
            
            
            % beginning guess
            dsf = ds_factors(fi_i);
            stamp_ln = round( (time_dn - start_time_l)*(24*60*60*self.fs*dsf) ) + 1;
            lines_t = stamp_ln + LINES_HEADER;
            linenum = lines_t;
            [~,tline2] = system(sprintf('sed -n %i,%ip %s',lines_t,lines_t,times_fn));
%             fprintf('--- tline2 ---\n %s\n --- end ---\n',tline2)
%             disp(times_fn)
%             est_stamp = str2double(tline2);
%             if ((est_stamp < 0) || (round(est_stamp) ~= est_stamp))
%                 fprintf(2,'ERROR: invalid stamp %s in %s\n',tline2,times_fn);
%                 return;
%             end

            if (isempty(tline2))
                fprintf(2,'[!] Could not map "%s" to %s\n',timeStr,sid);
                sout = '';
                fprintf('--- H5eeg.getSampleFromTime EXIT WITH ERROR ---\n')
                return
            end
            
            % check guess
            time_guess = datenum(tline2); %str2double(tline2);
            lines_times = 0;
            if (time_guess == time_dn)
                lines_times = linenum;
            else
                time_guess_prev = time_guess;
                linenum = linenum - round( (time_guess - time_dn)*(24*60*60*self.fs*dsf) );
                %disp(linenum)
                [~,tline2] = system(sprintf('sed -n %i,%ip %s',linenum,linenum,times_fn));
                time_guess = datenum(tline2); %str2double(tline2);
                lines_times = linenum;
                fprintf(2,'[!] Corrected time %s -> %s, error: %i days\n',datestr(time_guess_prev),datestr(time_guess),time_guess_prev-time_guess );
            end
            
            fclose(f);

            [~,tline3] = system(sprintf('sed -n %i,%ip %s',lines_times,lines_times,stamps_fn));

            % stamp
            stamp = str2double(tline3);
            
            % convert stamp to index number
            nsl = [0;cumsum(n_samples_l)] + 1; % starting samples
            stamp_start = double(h5read(self.filename,'/h5eeg/aux',[1 nsl(fi_i)],[1 1]));
            sampleNum_guess = nsl(fi_i) + floor((stamp - stamp_start)/dsf);
            % clip read bounds
            if (sampleNum_guess > self.n_samples)
                sampleNum_guess = self.n_samples;
            elseif (sampleNum_guess < 1)
                sampleNum_guess = 1;
            end
            stamp_guess = double(h5read(self.filename,'/h5eeg/aux',[1 sampleNum_guess],[1 1]));
            
            % check stamp guess
            if (abs(stamp_guess - stamp) <= dsf)
                sout = sampleNum_guess;
            else
                stamp_guess_prev = stamp_guess;
                sampleNum_guess = sampleNum_guess - floor((stamp_guess - stamp)/dsf);
                % clip read bounds
                if (sampleNum_guess > self.n_samples)
                    sampleNum_guess = self.n_samples;
                elseif (sampleNum_guess < 1)
                    sampleNum_guess = 1;
                end
                stamp_guess = double(h5read(self.filename,'/h5eeg/aux',[1 sampleNum_guess],[1 1]));
                fprintf(2,'[!] Corrected stamp %i -> %i, actual: %i\n',stamp_guess_prev,stamp_guess,stamp);
                sout = sampleNum_guess;
            end
            
            fprintf('--- H5eeg.getSampleFromTime EXIT_SUCCESS ---\n')
            
        end
        
    end
    
end

