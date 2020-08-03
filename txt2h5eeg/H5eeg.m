classdef H5eeg
    %[H5eeg]
    %   
    %   Implements a hdf5 data format for the efficient storage and
    %   retrieval of electrophysiological recordings.
    %
    %   Last Edited: Jerry Wang - jwang04@g.harvard.edu - 2016 Aug 28
    %
    %[Example Usage]
    %
    %   patient1 = H5eeg('m00001.hdf5')
    %   
    %   patient1.readAttributes() returns a struct of descriptive
    %   information about the data.
    %
    %   patient1.readEEG(<double>) will read at elapsed time
    %
    %   patient1.readEEG({<int>, <int>}) will read from first int index to
    %   second int index
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
            if (exist(fn,'file'))
                %info = h5info(self.filename,'/h5eeg/eeg');
                info = h5info(self.filename,'/h5eeg/aux');
                self.n_chan = length(h5readatt(self.filename,'/h5eeg/eeg','gains'));
                self.n_samples = info.Dataspace.Size(2);%info.Dataspace.Size(2);
                self.fs = h5readatt(self.filename,'/h5eeg/eeg','rate');
                tstamp = h5readatt(self.filename,'/h5eeg','timestamp');
                self.start_time = (tstamp + datenum(1970,1,1)*86400)*self.fs;
                self.end_time = self.start_time + self.n_samples - 1;
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
            eeg.labels = h5readatt(self.filename,'/h5eeg/eeg','labels');
            eeg.offsets = h5readatt(self.filename,'/h5eeg/eeg','offsets');
            eeg.gains = h5readatt(self.filename,'/h5eeg/eeg','gains');
            argc = length(varargin);
            switch (argc)
                case 0
                    eeg.data = h5read(self.filename,'/h5eeg/eeg');
                case 1
                    arg1in = varargin{1};
                    switch class(arg1in)
                        case 'double'
                            idx = self.findIndexWithElapsedTime(arg1in);
                            eeg.data = h5read(self.filename,'/h5eeg/eeg',[1 idx],[self.n_chan 1]);
                        case 'cell'
                            eeg.data = h5read(self.filename,'/h5eeg/eeg',[1 arg1in{1}],[self.n_chan (arg1in{2}-arg1in{1}+1)])';
                        otherwise
                            fprintf(2,'[!] unrecognized type: %s\n',...
                                class(arg1in));
                    end
                case 2
                    % Two numbers used for two channels
                    eeg.data(:,1) = h5read(self.filename,'/h5eeg/eeg',[varargin{1} 1],[1 self.n_samples])';
                    eeg.data(:,2) = h5read(self.filename,'/h5eeg/eeg',[varargin{2} 1],[1 self.n_samples])';
                otherwise
                    fprintf(2,'[!] unrecognized number of inputs: %i\n'...
                        ,argc);
            end
        end
    end
    
end

