close all; clear;

% -------------------------------------------------------------------------
% Run with 2 monitors: STDOUT on the left monitor
% -------------------------------------------------------------------------

% Location of h5 files to annotate from
h5dir = '/media/klab/44/h5';

% Output file prefix
of_pre = 'jw';
enable_write = false;

% Length of data to view
widthv = 10; % seconds

% Length of data for artifact annot
width = 1; % seconds

% Voltage offset for plotting (uV)
v_offset = 100; % uV
n_bchan_max = 154; % m00083
bkg_color = 0.33*[1 1 1];
mark_color = 0*[1 1 1];

% 20 Hz filtering
%[b,a] = butter(3,[0.2 0.6],'stop');

d = dir(h5dir);
dname = {d.name};
disdir = [d.isdir];
H5names = dname(~disdir);

% How many to annot
n_annot = 100;

for i = 1:n_annot
    % Choose subject at random
    h5_i = randi([1 length(H5names)],1,1);
    h5_i = 52;
    ecog = H5eeg(sprintf('%s/%s',h5dir,H5names{h5_i}));
    Fs = round(ecog.fs);
    
    % Power line frequencies (Hz)
    PLI_S = 56;
    PLI_E = 64;
    PL2_S = (ecog.fs-180)-2;
    PL2_E = (ecog.fs-180)+2;
    PL3_S = 117;
    PL3_E = 123;
    
    % Choose segment of data at random
    wv = round(ecog.fs * widthv);
    w = round(ecog.fs * width);
    start_i = randi([1 (ecog.n_samples - wv + 1)],1,1);
    end_i = start_i + wv - 1;
    
    % Read segment for view
    mb_read = (wv * ecog.n_chan * 32 * (1/8))*1e-6;
    tic;
    eeg = ecog.readEEG({start_i end_i});
    t_read = toc;
    fprintf('> Read %.2f MB/s\n',mb_read/t_read)
    
    % Define segment for marking
    start_j = start_i + floor(0.5 * wv);
    end_j = start_j + w - 1;
    
    % Plot
    h = figure;
    %set(h,'Units','Normalized','Position',[1 0 1 1])
    set(h,'Units','Normalized','Position',[0.35 0 0.65 0.65])
    
    % fill figure with axes
    AxesH = axes;
    drawnow;
    InSet = get(AxesH, 'TightInset');
    set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    
    
    bip = h5readatt(ecog.filename,'/h5eeg/eeg','bip');
    [n_bchan,~] = size(bip);
    T = linspace(0,widthv,wv);
    
    % Segment to mark
    start_p = (start_j-start_i+1);
    end_p = start_p + w - 1;
    
    % Plot bounds
    bound_color = 0.7*[1 0 0];
    plot([T(start_p) T(start_p)],[v_offset -v_offset*(n_bchan+1)],'Color',bound_color)
    hold on;
    plot([T(end_p) T(end_p)],[v_offset -v_offset*(n_bchan+1)],'Color',bound_color)
    hold on;
    
    % psd
    Pxx = [];
    Pxx2 = [];
    F = [];
    
    for j = 1:n_bchan
        
        % Plot eeg
        vbip = eeg.data(:,bip(j,1)) - eeg.data(:,bip(j,2));
        vbip = vbip - mean(vbip);
        vplot = vbip-v_offset*(j-0.5);
        plot(T,vplot,'Color',bkg_color);
        hold on;
        
        % Plot segment for mark
        
        plot(T(start_p:end_p),vplot(start_p:end_p),'Color',mark_color);
        hold on;
        
        % psd
        [pxx,f] = pwelch(vbip(start_p:end_p),Fs,[],[],Fs);
        [pxx2,~] = pwelch(vbip,Fs,[],[],Fs);
        Pxx = [Pxx,pxx];
        Pxx2 = [Pxx2,pxx2];
        F = [F, f];
    end

    axis([min(T) max(T) -v_offset*(n_bchan_max+1) v_offset]);
    set(AxesH,'TickLength',[0.005 0.005])
    xlabel('Seconds')
    ylabel('uV')
    
    % Plot psd of all visible
    h2 = figure;
    subplot(2,1,1)
    %set(h2,'Units','Normalized','Position',[0.5 0.5 0.5 0.5])
    set(h2,'Units','Normalized','Position',[0.5 0.75 0.25 0.25])
    plot(f,log10(max(Pxx2')),'-','Color',0.5*[1 1 1])
    hold on;
    plot(f,log10(mean(Pxx2,2)),'black');
    hold on;
    plot(f,log10(min(Pxx2')),'-','Color',0.5*[1 1 1])
    axis([min(f) max(f) -5 5]);
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    legend({'max','mean','min'})
    title('Whole region')
    
    % Plot psd of just marked region
    subplot(2,1,2)
    %set(h2,'Units','Normalized','Position',[0.5 0.5 0.5 0.5])
    plot(f,log10(max(Pxx')),'-','Color',0.5*[1 1 1])
    hold on;
    plot(f,log10(mean(Pxx,2)),'black');
    hold on;
    plot(f,log10(min(Pxx')),'-','Color',0.5*[1 1 1])
    axis([min(f) max(f) -5 5]);
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    legend({'max','mean','min'})
    title('Marked region')
    
    % Plot timeseries of max and min channel
    BRO_S = min(f);
    BRO_E = max(f);
    fmask = ((f >= BRO_S) & (f < PLI_S)) | ((f > PLI_E) & (f < PL2_S)) | ((f > PL2_E) & (f < PL3_S)) | ((f > PL3_E) & (f <= BRO_E));
    [~,i_min] = min(sum(Pxx(fmask,:)));%min(max(eeg.data) - min(eeg.data));
    [~,i_max] = max(sum(Pxx(fmask,:)));%max(max(eeg.data) - min(eeg.data));
    h3 = figure;
    %set(h3,'Units','Normalized','Position',[0 0.5 0.5 0.5])
    set(h3,'Units','Normalized','Position',[0.75 0.75 0.25 0.25])
    subplot(2,1,1)
    vbip = eeg.data(:,bip(i_max,1)) - eeg.data(:,bip(i_max,2));
    vbip = vbip - mean(vbip);
    plot(T(start_p:end_p),vbip(start_p:end_p),'Color',bkg_color);
    axis tight;
    ylabel('uV')
    xlabel('Seconds')
    title('Bipolar channel with minimum power')
    subplot(2,1,2)
    vbip = eeg.data(:,bip(i_min,1)) - eeg.data(:,bip(i_min,2));
    vbip = vbip - mean(vbip);
    plot(T(start_p:end_p),vbip(start_p:end_p),'Color',bkg_color);
    axis tight;
    ylabel('uV')
    xlabel('Seconds')
    title('Bipolar channel with maximum power')
    
    % compute bipolar of marked segment
%     Vbip = zeros(n_bchan,w);
%     for k = 1:n_bchan
%         Vbip(:,i) = eeg.data
%     end
%     [~,i_min] = min(max(eeg.data) - min(eeg.data));
%     [~,i_max] = max(max(eeg.data) - min(eeg.data));
    
    % Get input
    sid = strsplit(ecog.filename,'/');
    sid = sid{end};
    sid = strsplit(sid,'.');
    sid = sid{1};
    istr = input(sprintf('(%s:%i-%i)\tIs this an artifact? [y/N] ',sid,start_j,end_j),'s');
    isart = false;
    switch istr
        case ''
            isart = false;
        case 'y'
            isart = true;
    end
    close(h);
    close(h2);
    close(h3);
    
    % Append to file
    if (enable_write)
        of = fopen(sprintf('%s.txt',of_pre),'a');
        fprintf(of,'%s,%i,%i,%i\n',sid,start_j,end_j,isart);
        fclose(of);
    end
end

