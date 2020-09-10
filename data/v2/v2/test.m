close all;
clear;

load_annot;
n_features = 2;
Y_col = [];
Y_col(1,:) = [0 0 0];
%Y_col(2,:) = [1 0 0];
Y_col(2,:) = [0 1 1];
Y_col(3,:) = [158, 92, 0]/255;
Y_col(4,:) = [0 1 0];
Y_col(5,:) = [0 0 1];
Y_col(6,:) = [126, 0, 224]/255;
%no_show = true;
no_show = false;
out_fmt = '-dpng';

%for i = 1:length(Sub)
%102,869 - unplugged (,511,566,742)
%76,99 - 20 Hz
%172 - 20 Hz a lot unmarked
%95 - spike
%349 - missed artifact
%49 - missed green?
%854,176,950 - EKG
%1012 - missed power line (j = 41)

I_temp = randi([1 length(Sub)]);
%I_temp = 742;
%I_temp = 76;
system('mkdir test');
system('mkdir test_noshow');
fprintf('Annotation: %i\n',I_temp);
%for i = I_temp % 
%for i = 1:length(Sub)


for l = 1:(52*50)
    sid = USub{randi([1 length(USub)])};
    seconds_view = 10; % number of seconds to view
    h5fn = sprintf('%s/%s.h5',h5dir,sid);
    ecog = H5eeg(h5fn);
    Fs = round(ecog.fs);
    i = 1;
    Start = randi([1 (ecog.n_samples - (Fs * seconds_view))]);
    End = Start + (Fs * seconds_view) - 1;
    Label = 2;


    eeg = ecog.readEEG({Start(i) End(i)});
    bip = h5readatt(h5fn,'/h5eeg/eeg','bip');
    [n_bchan,~] = size(bip);

    h = figure;
    set(h,'Units','Normalized','Position',[0 0 1 1]);

    AxesH = axes;
    drawnow;
    InSet = get(AxesH, 'TightInset');
    InSet = InSet * 4;
    InSet(1) = InSet(1) * 4;
    set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    

    X = zeros(ecog.n_chan,n_features);
    T = linspace(0,seconds_view,(Fs * seconds_view));
    n_seg = 0;
    n_mark = 0;
    A = zeros(n_bchan,seconds_view);
    meanV = zeros(1,n_bchan);
    
    y_offset_uv = 150;
    for k = 1:seconds_view
        Pxx = [];
        F = [];
        Ys = zeros(1,n_bchan);
        for j = 1:n_bchan
            V = eeg.data(:,bip(j,1)) - eeg.data(:,bip(j,2));        
            k2 = (k-1)*Fs + 1;
            V = V(k2:(k2+Fs-1));
            if (k == 1)
                meanV(j) = mean(V);
            end
            Vsav = V;
            V = V - mean(V);
            Vp = Vsav - meanV(j);
            X(j,:) = features_v2(V, Fs);
            Y = classify_v2(X(j,:));

            Ys(j) = Y;
            
            yoffset = (j - 1) * y_offset_uv;
            col = Y_col(Y+1,:);
            if ~((Y ~= 0) && (no_show))
                plot(T(k2:(k2+Fs-1)),Vp + yoffset,'color',col,'LineWidth',0.01); hold on;
                %col = [1 1 1];
            end
            %plot(V + yoffset,'color',col); hold on;

            [pxx,f] = pwelch(V,round(Fs),[],[],Fs);
            Pxx = [Pxx, pxx];
            F = [F, f];
            
            A(j,k) = Y;
            if (Y ~= 0)
                n_mark = n_mark + 1;
            end
            n_seg = n_seg + 1;
        end

    end

    axis([min(T) max(T) -y_offset_uv 157*y_offset_uv]);
    As = sum(A ~= 0,2);
    title(sprintf('%s - %i BIP : %.0f %% marked (0 - %.0f %%, 1 - %.0f %%, 2 - %.0f %%)',...
        sid,n_bchan,100*(n_mark/n_seg),100*sum(As > 0)/n_bchan,100*sum(As > 1)/n_bchan,100*sum(As > 2)/n_bchan))
    xlabel('Seconds')
    %ylabel('Voltage (uV)')
    
    % Ylabel
    [n_bchan,~] = size(bip);
    chan_labels = h5readatt(h5fn,'/h5eeg/eeg','labels');
    try
        dk_labels = h5readatt(h5fn,'/h5eeg/eeg','labels_2');
    catch
        dk_labels = chan_labels;
    end
    bchan_labels = cell(1,n_bchan);
    for j = 1:n_bchan
        b1 = bip(j,1);
        b2 = bip(j,2);
        bchan_labels{j} = sprintf('%s-%s (%s,%s)',...
            chan_labels{b1},chan_labels{b2},dk_labels{b1},dk_labels{b2});
    end
    yticks(0:y_offset_uv:round(y_offset_uv*n_bchan))
    yticklabels(bchan_labels)
    set(gca, 'TickDir', 'out');
    set(gca, 'FontSize',6);
    
    return

    % 21 ekg
    %if (i == 35) % randi([1 length(Sub)])
    %     if (true)
    %         return
    %     end
    if (no_show)
        print(h,sprintf('test_noshow/%s_%i_%i_%i',sid,seconds_view,Start(i),End(i)),out_fmt)
    else
        print(h,sprintf('test/%s_%i_%i_%i',sid,seconds_view,Start(i),End(i)),out_fmt)
    end
    close(h);
    
    % plot psd
    h2 = figure;
    set(h2,'Units','Normalized','Position',[0.2 0 0.5 0.5]);
    plot(f,log10(max(Pxx')),'-','Color',0.5*[1 1 1])
    hold on;
    plot(f,log10(mean(Pxx,2)),'black');
    hold on;
    plot(f,log10(min(Pxx')),'-','Color',0.5*[1 1 1])
    for j2 = 1:n_bchan
        if (Ys(j2) == 3)
            plot(f,log10(Pxx(:,j2)),'Color',[0 0 1]);
            hold on;
        elseif (Ys(j2) == 4)
            plot(f,log10(Pxx(:,j2)),'Color',[1 0 1]);
            hold on;
        end
    end
    axis([min(f) max(f) -5 5]);
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    legend({'max','mean','min'})
    
    %return
    if (no_show)
        print(h2,sprintf('test_noshow/%s_%i_%i_%i_psd',sid,seconds_view,Start(i),End(i)),out_fmt)
    else
        print(h2,sprintf('test/%s_%i_%i_%i_psd',sid,seconds_view,Start(i),End(i)),out_fmt)
    end
    close(h2);
    
end
%end
