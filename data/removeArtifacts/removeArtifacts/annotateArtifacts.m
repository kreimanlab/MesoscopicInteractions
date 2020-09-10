function annotateArtifacts

clear all; close all;% clc;

% Prompt for user information
fprintf('#------------------------------------------------------------#\n')
fprintf('|   Artifact Annotation                                      |\n')
fprintf('#------------------------------------------------------------#\n')
listing = dir();
dirnames = {listing.name};

% Prompt for dataset selection
fprintf('\n *  Available datasets:\n');
datasets = {};
l = 1;
for k = 1:length(listing)
    if ((~isempty(strfind(dirnames{k},'.hdf5'))) || (~isempty(strfind(dirnames{k},'.h5'))))
        fprintf('\t%s\n',dirnames{k});
        datasets{l} = dirnames{k};
        l = l + 1;
    end
end
%infilename = input(' $  Data file name (include file extension)\n: ','s');
%infilename = 'sub21_day3.hdf5';
%infilename = 'sub40_cut1.hdf5';
%infilename = 'sub33_3rd9h.hdf5';
%infilename = 'sub33_1st9h.hdf5';
infilename = 'sub21_day2.hdf5';
%infilename = 'm00125_1.hdf5';

% Extract dataset file name to include into annotation name
dset = strsplit(infilename,'.');
dset = dset{1};
fprintf('\n *  Existing annotations:\n');
for k = 1:length(listing)
    if (~isempty(strfind(dirnames{k},'.annot')))
        fprintf('\t%s\n',dirnames{k});
    end
end
%username = input(' $  Username (no extensions, will be appended by ''_'' followed by ''Data file name'')\n: ','s');
%username = 'jer';
%username = 'test';
username = 'qda';

% Parameters
INPUT_FILE_NAME = infilename;
UI_TIME_LENGTH = 7; % seconds
CHANNEL_OFFSET = 8; % standard deviations
CHANNEL_LABEL_FONT_SIZE = 8; % font size
CHANNEL_LABEL_COLOR = [0.2 0 0.4]; % color
EVENT_COLOR = [0.5 0 0]; % color
EVENTS_FONT_SIZE = 12; % event font size
EVENT_BELOW = 2; % how many lines below last channels to dock events
MARKER_SIZE = 6; % size of mark marker
MARKER_LINE_WIDTH = 1; % width of mark line
start_time = 1;
chan_offset = NaN;
chan_means = [];
Fs = NaN;
TIMESCALE_LIST = {'7s','15s','30s','1min','2min','4min'};
TIMESCALE_SECONDS = [7,15,30,60,120,240];
OFFSET_LIST = {'XXXFine','XXFine','XFine','Fine','Normal','Rough','XRough','XXRough','XXXRough','4XRough','5XRough'};
OFFSET_VALUES = [0.5,1,2,4,8,16,32,64,128,256,512];

% plotEEG properties
username = [username,'_',dset,'.annot'];
properties.outfilename = username;
properties.CAR = true; % Whether or not to subtract common average
properties.colored = 2; % color mode
properties.center = true; % whether or not to center voltages
properties.ymax = NaN;
properties.ymin = NaN;
properties.xmarked = [];
properties.data = [];
properties.showBrain = false;
properties.event_offset = 3; % how many channels between each event marker
colors.coral = ([248, 131, 121]-40)/255;
colors.singleMark = [0 1 1];
colors.seafoamgreen = ([114, 214, 155]-40)/255;
colors.background = [0.9 0.9 0.9];
colors.plotbackground = [1 1 1];
PLOT_LINE_BRIGHTNESS = 0.7;

% plotEEG persistence variables
%persist.eeg = [];
persist.start_time = 0;
persist.end_time = 0;

% read in previous annotations start times
try
    %Elist_t = textscan(fopen(username,'r'),'%f');
    %Elist_t = unique(Elist_t{1});
    [Elist_t_one, Elist_t] = readExcludelist(username);
catch
    Elist_t = [];
    Elist_t_one = [];
end
Elist_t = sort(unique(Elist_t))';
Elist_t_one = unique(sort(Elist_t_one)','rows');
%Elist_t_one = unique(sort([Elist_t_one(2,:);Elist_t_one(1,:)])','rows');

Elist = [];
ElistStr = {};
ElistOne = [];
ElistStrOne = {};
iii = 1;
iiiOne = 1;
previousEidx = 0;
previousEidxOne = 0;
for ii = 1:length(Elist_t)
    if ((Elist_t(ii) - previousEidx) > 1)
        Elist = [Elist, Elist_t(ii)];
        ElistStr{iii} = ['[',num2str(iii),'] ',num2str(Elist_t(ii))];
        iii = iii + 1;
    end
    previousEidx = Elist_t(ii);
end
for ii = 1:length(Elist_t_one)
    if ((Elist_t_one(ii,2) - previousEidxOne) > 1)
        ElistOne = [ElistOne, Elist_t_one(ii,2)];
        ElistStrOne{iiiOne} = ['[',num2str(iiiOne),'] ',...
            num2str(Elist_t_one(ii,2)),' (',num2str(Elist_t_one(ii,1)),')'];
        iiiOne = iiiOne + 1;
    end
    previousEidxOne = Elist_t_one(ii,2);
end
%[ElistOne,eloidx] = sort(ElistOne);
%ElistStrOne = ElistStrOne(eloidx);

% Read in data with H5eeg
ecog = H5eeg(INPUT_FILE_NAME);
att = ecog.readAttributes();
HILIGHT_LIST = h5readatt(ecog.filename,'/h5eeg/eeg','labels');
HILIGHT_VALUES = 1:length(HILIGHT_LIST);
HILIGHT = false(1,length(HILIGHT_LIST));
HILIGHT_FACTOR = 0.1; % fraction to multiply to RGB 1x3 color vector 

% make out file if it doesn't already exist
if (fopen(username,'r') <= 0)
    fprintf(fopen(username,'w'),'');
end

% Create window
f = figure('units','normalized','outerposition',[0 0 1 1],'Color',colors.background,'Name','ECoG Viewer :: Kreiman Lab 2016','NumberTitle','off');

% axes
ax = axes('units','normalized','outerposition',[-0.15,-0.005,1.245,1.06]);
set(ax,'Color',colors.plotbackground);
global ax2

% Buttons
btnNext = uicontrol('Style', 'pushbutton', 'String', '>', 'Units',...
    'Normalized','Position', [0.5 0.002 0.04 0.04],...
    'Callback', @onButtonNextClick);
btnPrev = uicontrol('Style', 'pushbutton', 'String', '<', 'Units',...
    'Normalized','Position', [0.46 0.002 0.04 0.04],...
    'Callback', @onButtonPrevClick);
btnNext = uicontrol('Style', 'pushbutton', 'String', '>>', 'Units',...
    'Normalized','Position', [0.54 0.002 0.04 0.04],...
    'Callback', @onButtonNext2Click);
btnPrev = uicontrol('Style', 'pushbutton', 'String', '<<', 'Units',...
    'Normalized','Position', [0.42 0.002 0.04 0.04],...
    'Callback', @onButtonPrev2Click);
btnCAROn = uicontrol('Style', 'pushbutton', 'String', 'CAR', 'Units',...
    'Normalized','Position', [0.01 0.002 0.04 0.04],...
    'Callback', @onButtonCAROnClick);
btnCAROff = uicontrol('Style', 'pushbutton', 'String', 'BRAIN', 'Units',...
    'Normalized','Position', [0.05 0.002 0.04 0.04],...
    'Callback', @onButtonCAROffClick);
btnColor = uicontrol('Style', 'pushbutton', 'String', 'COLOR', 'Units',...
    'Normalized','Position', [0.09 0.002 0.04 0.04],...
    'Callback', @onButtonColorClick);
btnColor = uicontrol('Style', 'pushbutton', 'String', 'EVENT', 'Units',...
    'Normalized','Position', [0.13 0.002 0.04 0.04],...
    'Callback', @onButtonEventClick);
btnErase = uicontrol('Style', 'pushbutton', 'String', 'ERASE', 'Units',...
    'Normalized','Position', [0.30 0.002 0.04 0.04],...
    'Callback', @onButtonEraseClick);
btnCenter = uicontrol('Style', 'pushbutton', 'String', 'CENTER', 'Units',...
    'Normalized','Position', [0.34 0.002 0.08 0.04],...
    'Callback', @onButtonCenterClick);
btnMark = uicontrol('Style', 'pushbutton', 'String', 'MARK ALL', 'Units',...
    'Normalized','Position', [0.58 0.002 0.04 0.04],...
    'Callback', @onButtonMarkClick);
btnMarkSingle = uicontrol('Style', 'pushbutton', 'String', 'MARK ONE', 'Units',...
    'Normalized','Position', [0.62 0.002 0.04 0.04],...
    'Callback', @onButtonMarkSingleClick);
popupTimescale = uicontrol('Style', 'popup', 'String', ...
    TIMESCALE_LIST, 'Units','Normalized','Position', ...
    [0.66 0.002 0.06 0.03],'Callback', @setmap);
popupOffset = uicontrol('Style', 'popup', 'String', ...
    OFFSET_LIST, 'Units','Normalized','Position', ...
    [0.72 0.002 0.06 0.03],'Callback', @setmapOffset);
popupHilight = uicontrol('Style', 'popup', 'String', ...
    HILIGHT_LIST, 'Units','Normalized','Position', ...
    [0.78 0.002 0.06 0.03],'Callback', @setmapHilight);
popupPrevAnnot = uicontrol('Style', 'popup', 'String', ...
    ElistStr, 'Units','Normalized','Position', ...
    [0.84 0.002 0.06 0.03],'Callback', @setmapPrevAnnot);
popupPrevAnnotOne = uicontrol('Style', 'popup', 'String', ...
    ElistStrOne, 'Units','Normalized','Position', ...
    [0.90 0.002 0.06 0.03],'Callback', @setmapPrevAnnotOne);
%txtLabel = uicontrol('Style','text','Units','Normalized', ...
%        'Position',[0.005 0.985 0.08 0.01],'String','Kreiman Lab :: 2016');
    
% Plot
plotEEG(start_time)

function [excludelistOne, excludelistAll] = readExcludelist(username)
    excludelist = textscan(fopen(username,'r'),'%f');
    excludelist = excludelist{1};
    % filter out single marks
    excludelistOne = zeros(2,sum(excludelist < 0));
    excludelistAll = zeros(1,length(excludelist) - 2*sum(excludelist < 0));
    skip = false;
    i2 = 1;
    i3 = 1;
    for i = 1:length(excludelist)
        if (~skip)
            if (excludelist(i) < 0)
                %excludelistOne = [excludelistOne, [excludelist(i+1); -1*excludelist(i)]];
                excludelistOne(:,i3) = [excludelist(i+1); -1*excludelist(i)];
                i3 = i3 + 1;
                skip = true;
            else
                %excludelistAll = [excludelistAll, excludelist(i)];
                excludelistAll(i2) = excludelist(i);
                i2 = i2 + 1;
            end
        else
            skip = false;
        end
    end
end

function plotEEG(start_time)
    
    % calculate end time
    end_time = start_time + ecog.fs*UI_TIME_LENGTH - 1;
    
    % Only read new data if end and start times are different
    repeat = (persist.start_time == start_time) && (persist.end_time == end_time);
    if (repeat)
        eeg = persist.eeg;
        %1+1;
    else
        tic;
        eeg = ecog.readEEG({start_time, end_time});
        readTime = toc;
        readSpeed = ((end_time-start_time+1)*ecog.n_chan*64/(8*1e6))/readTime;
        text(start_time,CHANNEL_OFFSET,sprintf('%i MB/s\n',round(readSpeed)));
    end
    eegsave = eeg;
    Fs = ecog.fs;
    
    % read exclude list
    properties.data = []; % reset data for plotting marked
    [excludelistOne, excludelistAll] = readExcludelist(username);
    excludelist = sort(unique(excludelistAll))';
    excludelistOne = unique(sort(excludelistOne)','rows');
    exclude_idx = excludelist((excludelist >= start_time) & (excludelist <= end_time));
    excludeOne_idx = excludelistOne((excludelistOne(:,2) >= start_time) & (excludelistOne(:,2) <= end_time),:);
    
    % convert data from integer to double, if necessary
    if (isa(eeg.data,'int16') )%&& (start_time == 1))
        eeg.data = double(eeg.data)./100000;
        for i = 1:ecog.n_chan
            eeg.data(:,i) = (eeg.data(:,i) - double(eeg.offsets(i)))/double(eeg.gains(i));
        end
    end
    
    % time
    time = (start_time:end_time);

    % Clear anatomy axes before plotting
    if (properties.showBrain)
        set(ax2,'visible','off');
        cla(ax2)
    end
    
    % only compute the offset during initialization
    if (start_time == 1)
        chan_offset = min(std(eeg.data));
    end
    
    % plot
    hold all;
    if (ecog.n_chan <= 100)
        n_color_repeat = ceil(ecog.n_chan/10)-1;
    else
        n_color_repeat = 10;
    end
    
    % Color control
    switch (properties.colored)
        case 2
            % Color by name of electrode prefix
            n_color_repeat = ecog.n_chan;
            
            colorLabels = cell(ecog.n_chan,1);
            for i = 1:ecog.n_chan
                entry = eeg.labels{i};
                colorLabels{i} = entry(isletter(entry));
            end
            [colorGroups,~,colorIdx] = unique(colorLabels);
            palette = jet(length(colorGroups));
            % Interweave colors
            palette = palette([1:2:length(colorGroups),2:2:length(colorGroups)]',:);
            cc = palette(colorIdx,:);
            for i = 1:ecog.n_chan
                if (HILIGHT(i))
                    cc(i,:) = HILIGHT_FACTOR * cc(i,:);
                end
            end
        case 1
            cc = jet(ceil(ecog.n_chan/n_color_repeat));
        case 0
            cc = 0.1*ones(ceil(ecog.n_chan/n_color_repeat),3);
    end
    
    if (properties.CAR)
        eeg.data = eeg.data - mean(eeg.data,2);
    end
    

    for i = 1:ecog.n_chan
        yoffset = -1*(i-1)*CHANNEL_OFFSET;
        xoffset = round(-0.015*ecog.fs*UI_TIME_LENGTH)+start_time-1;
        
        % only compute the voltage centers if set
        if (properties.center)
            chan_means(i) = median(eeg.data(:,i));
        end
        
        % Handle properties
        if (properties.CAR)
            if (properties.colored > 0)
                CHANNEL_LABEL_COLOR = 0.9;%[0 0.4 0.2];
            else
                CHANNEL_LABEL_COLOR = 5;
            end
%             common_avg = mean(eeg.data,2);
%             for j = 1:ecog.n_chan
%                 eeg.data(:,j) = eeg.data(:,j) - common_avg;
%             end
        else
            CHANNEL_LABEL_COLOR = 0.5;%[0.2 0 0.4];
        end
        
        % Plot
        data_i = (eeg.data(:,i) - chan_means(i))/chan_offset + yoffset;
        data_b = data_i;
        
        % Filter out marked samples
        % global marks
        eindex = exclude_idx - start_time + 1;
        data_i(eindex) = NaN;
        % single marks
        eOneSampIndex = excludeOne_idx(:,2) - start_time + 1; % samples
        eOneChanIndex = excludeOne_idx(:,1); % channels
        eOneIndex = find(eOneChanIndex == i);
        if (~isempty(eOneIndex))
            data_i(eOneSampIndex(eOneIndex)) = NaN;
        end
        properties.data(:,i) = data_i;
        % Downsample for high sampling rates
        %if (Fs >= 1000)
        %    plot(ax,downsample(time,2),downsample(data_i,2),'color',PLOT_LINE_BRIGHTNESS*cc(mod(i-1,n_color_repeat)+1,:));
        %elseif (Fs >= 2000)
        %    plot(ax,downsample(time,4),downsample(data_i,4),'color',PLOT_LINE_BRIGHTNESS*cc(mod(i-1,n_color_repeat)+1,:));
        %else
        %    plot(ax,time,data_i,'color',PLOT_LINE_BRIGHTNESS*cc(mod(i-1,n_color_repeat)+1,:));
        %end
        plot(ax,time,data_i,'color',PLOT_LINE_BRIGHTNESS*cc(mod(i-1,n_color_repeat)+1,:));
        
        t = text(ax,xoffset,yoffset,eeg.labels{i},'color',...
            cc(mod(i-1,n_color_repeat)+1,:)*CHANNEL_LABEL_COLOR,'FontWeight','bold');
        t.FontSize = CHANNEL_LABEL_FONT_SIZE;
        
        % plot marked segments
        plot(ax,time(eindex),data_b(eindex),'.','MarkerSize',1,'color',colors.coral);
        plot(ax,time(eOneSampIndex(eOneIndex)),data_b(eOneSampIndex(eOneIndex)),'.','MarkerSize',1,'color',colors.singleMark);
        
        % Save maximum and minimum y values for scale
        if (i == 1)
            properties.ymax = max(data_i);
        elseif (i == ecog.n_chan)
            properties.ymin = min(data_i);
        end
    end
    
    % reset center after computing
    properties.center = false;
    
    % Plot events
    event_idx = ismember(att.events.start_idx,start_time:end_time);
    switch class(att.events.name)
        case 'cell'
            event_names = att.events.name(event_idx);
        case 'char'
            event_names_char = att.events.name(:,event_idx)';
            [chars, ~] = size(event_names_char);
            event_names = cell(chars,1);
            for i = 1:chars
                event_names{i} = event_names_char(i,:);
            end
    end
    event_index = att.events.start_idx(event_idx);
    event_duration = att.events.duration(event_idx);
    % place below last channel
    yoffset = yoffset - EVENT_BELOW*CHANNEL_OFFSET;
    for i = 1:length(event_names)
        % Introduce drift to event labels
        yoffset= yoffset + properties.event_offset*CHANNEL_OFFSET;
        
        if (event_duration == 0)
            % instantaneous comments displayed as flags
            t = text(ax,event_index(i),yoffset,event_names{i},'color',EVENT_COLOR);
        else
            % length of event is plotted
            class(['[',event_names{i}])
            class(yoffset)
            class(event_index(i))
            t = text(ax,double(event_index(i)),yoffset,['[',event_names{i}],'color',EVENT_COLOR);
            t2 = text(ax,double(event_index(i))+double(event_duration(i)),yoffset,']','color',EVENT_COLOR);
            t2.FontSize = EVENTS_FONT_SIZE;
        end
        t.FontSize = EVENTS_FONT_SIZE;
    end
    
    % Set axis properties
    axis tight
    xlim(ax,[xoffset end_time])
    ax.XTick = (0:UI_TIME_LENGTH)*ecog.fs + (start_time-1);
    ax.XTickLabelRotation = 90;
    ax.YTick = [];
    ax.YTickLabelRotation = 90;
    set(ax,'TickDir','out')
    Xlabels = cell(1,UI_TIME_LENGTH+1);
    for i = 1:(UI_TIME_LENGTH+1)
        Xlabels{i} = datestr((ecog.start_time+(start_time-1)+(i-1)*...
            ecog.fs)/(ecog.fs*86400),'HH:MM:SS');
    end
    set(ax,'XTickLabel',Xlabels);
    set(ax,'XGrid','on','XMinorGrid','on');
    set(ax,'YTickLabel',{});
    if (start_time == 1)
        title(ax,strcat('Subject- ',strtrim(att.subject),', Experiment- ',...
            strtrim(att.experiment),', Start Time- ',datestr(...
            ecog.start_time/(ecog.fs*86400),'yyyy/mm/dd HH:MM:SS.FFF'),...
            ' (',num2str(ecog.fs),' Hz)',', End Time- ',...
            datestr(ecog.end_time/(ecog.fs*86400),...
            'yyyy/mm/dd HH:MM:SS.FFF')));
    end
    
    
    % Plot within plot
    
    if (properties.showBrain)
        ax2 = axes('Parent',f,'OuterPosition',[.3 .1 .8 .8]); %set(ax2,'visible','off');
        axes(ax2);
        paintColors = [0.8 0 0];
        % assume first 6 letters in filename is patient ID
        patientID = ecog.filename(1:6);
        
        % Default brain plot parameters
        hemi = 'both';
        view = 'lfront';
        if (properties.colored == 2)
            paintColors = cc;
        elseif (properties.colored == 1)
            paintColors = zeros(ecog.n_chan,3);
            for i = 1:ecog.n_chan
                paintColors(i,:) = cc(mod(i-1,n_color_repeat)+1,:);
            end
            cc = paintColors;
        end
        
        %------------------ PATIENT SPECIFIC SPECIFICATION ----------------
        switch (patientID)
            case 'sub40'
                if (properties.colored > 0)
                    %paintColors = cc([1:128,187:188,131:186]',:);
                    % Depth electrodes DE1-10 not plotted
                    paintColors = paintColors(1:176,:);
                end
                hemi = 'rh';
                view = 'rback';
            case 'sub21'
                hemi = 'rh';
                view = 'rfront';
            case 'sub33'
                hemi = 'rh';
                view = 'rfront';
        end
        %------------------------------------------------------------------
        
        plotBrainH(patientID,paintColors,hemi,view,eeg.labels);
    end
    % Persistence
    persist.eeg = eegsave;
    persist.start_time = start_time;
    persist.end_time = end_time;
end

% Next Button
function onButtonNextClick(~,~)
    start_time = start_time+UI_TIME_LENGTH*Fs;
    cla(ax)
    plotEEG(start_time);
end
function onButtonNext2Click(~,~)
    start_time = start_time+UI_TIME_LENGTH*Fs*10;
    cla(ax)
    plotEEG(start_time);
end

% Prev Button
function onButtonPrevClick(~,~)
    start_time = start_time-UI_TIME_LENGTH*Fs;
    cla(ax)
    plotEEG(start_time);
end
function onButtonPrev2Click(~,~)
    start_time = start_time-UI_TIME_LENGTH*Fs*10;
    cla(ax)
    plotEEG(start_time);
end

% CAR ON Button
function onButtonCAROnClick(~,~)
    properties.CAR = ~ properties.CAR;
    cla(ax)
    plotEEG(start_time);
end

% ANATOMY Button
function onButtonCAROffClick(~,~)
    properties.showBrain = ~ properties.showBrain;
    cla(ax)
    if (~properties.showBrain)
        set(ax2,'visible','off');
        cla(ax2)
    end
    plotEEG(start_time);
end

% Color Button
function onButtonColorClick(~,~)
    properties.colored = properties.colored + 1;
    if (properties.colored == 3)
        properties.colored = 0;
    end
    cla(ax)
    plotEEG(start_time);
end

% Event Button
function onButtonEventClick(~,~)
    if (properties.event_offset == 3)
        properties.event_offset = 0;
    elseif (properties.event_offset == 0)
        properties.event_offset = 1;
    elseif (properties.event_offset == 1)
        properties.event_offset = 3;
    end
    cla(ax)
    plotEEG(start_time);
end

% Center Button
function onButtonCenterClick(~,~)
    properties.center = true;
    cla(ax)
    plotEEG(start_time);
end

% Mark Button
function onButtonMarkClick(~,~)
    [x,y] = ginput(1);
    plot(ax,[x x],[properties.ymin properties.ymax],'color',[0 0.5 0],'linewidth',MARKER_LINE_WIDTH)
    plot(ax,x,y,'o','color',[0 0.5 0],'MarkerSize',MARKER_SIZE)
    [x2,y2] = ginput(1);
    plot(ax,[x2 x2],[properties.ymin properties.ymax],'color',[0.5 0 0],'linewidth',MARKER_LINE_WIDTH)
    plot(ax,x2,y2,'o','color',[0.5 0 0],'MarkerSize',MARKER_SIZE)
    of = fopen(properties.outfilename,'a');
    for i = x:x2
        fprintf(of,'%f\n',round(i));
        properties.xmarked = [properties.xmarked; round(i)];
    end
    fclose(of);
    
    for i = 1:ecog.n_chan
        plot(ax,properties.xmarked,properties.data(properties.xmarked-start_time+1,i),'color',colors.seafoamgreen);
    end
    properties.xmarked = [];
end

% Mark One Button
function onButtonMarkSingleClick(~,~)
    [x,y] = ginput(1);
    plot(ax,[x x],[properties.ymin properties.ymax],'color',[0 0.5 0],'linewidth',MARKER_LINE_WIDTH)
    plot(ax,x,y,'o','color',[0 0.5 0],'MarkerSize',MARKER_SIZE)
    [x2,y2] = ginput(1);
    plot(ax,[x2 x2],[properties.ymin properties.ymax],'color',[0.5 0 0],'linewidth',MARKER_LINE_WIDTH)
    plot(ax,x2,y2,'o','color',[0.5 0 0],'MarkerSize',MARKER_SIZE)
    
    %fprintf('x: %i x2: %i y: %i y2: %i\n', x, x2, y, y2)
    [~, mchan] = min(abs(properties.data(round(x)-start_time+1,:) - y));
    %fprintf('mchan: %i\n',mchan)
    
    of = fopen(properties.outfilename,'a');
    for i = x:x2
        fprintf(of,'-%i\n%f\n',mchan,round(i));
        properties.xmarked = [properties.xmarked; round(i)];
    end
    fclose(of);
    
    %for i = 1:ecog.n_chan
    plot(ax,properties.xmarked,properties.data(properties.xmarked-start_time+1,mchan),'color',colors.seafoamgreen);
    %end
    properties.xmarked = [];
end

% Erase Button
function onButtonEraseClick(~,~)
    [x,y] = ginput(1);
    plot(ax,[x x],[properties.ymin properties.ymax],'color',[0 0.5 0],'linewidth',MARKER_LINE_WIDTH)
    plot(ax,x,y,'o','color',[0 0.5 0],'MarkerSize',MARKER_SIZE)
    [x2,y2] = ginput(1);
    plot(ax,[x2 x2],[properties.ymin properties.ymax],'color',[0.5 0 0],'linewidth',MARKER_LINE_WIDTH)
    plot(ax,x2,y2,'o','color',[0.5 0 0],'MarkerSize',MARKER_SIZE)
    
%     for i = x:x2
%         fprintf(of,'%f\n',round(i));
%         properties.xmarked = [properties.xmarked; round(i)];
%     end

    % Marked region is x:x2;
    of = fopen(properties.outfilename,'r');
    eraseIdx = textscan(of,'%f'); eraseIdx = eraseIdx{1};

    %fprintf('x: %i, x2: %i\n',x,x2);
    % Filter out all marked indices
    eraseIdx = eraseIdx(~ismember(eraseIdx,round(x):round(x2)));
    
    % Filter out unpaired negative values
    for i = 1:(length(eraseIdx)-1)
        if (eraseIdx(i) < 0) && (eraseIdx(i+1) < 0)
            eraseIdx(i) = NaN;
        end
    end
    eraseIdx(isnan(eraseIdx)) = [];
    
    % final case
    if (eraseIdx(length(eraseIdx)) < 0)
        eraseIdx(length(eraseIdx)) = [];
    end
    
    fclose(of);
    
    % Save file
    of = fopen(properties.outfilename,'w');
    for i = 1:length(eraseIdx)
        fprintf(of,'%f\n',eraseIdx(i));
    end
    fclose(of);
    
    % Reset drawing area
    cla(ax)
    plotEEG(start_time);
    
%     for i = 1:ecog.n_chan
%         plot(ax,properties.xmarked,properties.data(properties.xmarked-start_time+1,i),'color',colors.seafoamgreen);
%     end
%     properties.xmarked = [];
end

% Popup
function setmap(source,~)
    val = source.Value;
    UI_TIME_LENGTH = TIMESCALE_SECONDS(val);
    cla(ax)
    plotEEG(start_time);
end
function setmapOffset(source,~)
    val = source.Value;
    CHANNEL_OFFSET = OFFSET_VALUES(val);
    cla(ax)
    plotEEG(start_time);
end
function setmapHilight(source,~)
    val = HILIGHT_VALUES(source.Value);
    HILIGHT(val) = ~ HILIGHT(val);
    cla(ax)
    plotEEG(start_time);
end
function setmapPrevAnnot(source,~)
    start_time = Elist(source.Value) - 4*ecog.fs;
    if (start_time < 1)
        start_time = 1;
    end
    cla(ax)
    plotEEG(start_time);
end
function setmapPrevAnnotOne(source,~)
    start_time = ElistOne(source.Value) - 4*ecog.fs;
    if (start_time < 1)
        start_time = 1;
    end
    cla(ax)
    plotEEG(start_time);
end

end
