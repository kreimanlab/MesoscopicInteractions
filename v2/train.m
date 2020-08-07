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

%for i = 1:length(Sub)
%102,869,452,727 - unplugged
%76,99 - 20 Hz
%172 - 20 Hz a lot unmarked (fixed)
%95 - spike (fixed)
%349 - power line (fixed)
%49 - missed green? (fixed)
%854,176,950 - EKG
%1012 - missed power line (j = 41) (fixed)
%return
I_temp = randi([1 length(Sub)]);
I_temp = 727;
system('mkdir train');
system('mkdir train_noshow');
fprintf('Annotation: %i\n',I_temp);
for i = I_temp % 
%for i = 1:length(Sub)
    
    h5fn = sprintf('%s/%s.h5',h5dir,Sub{i});
    ecog = H5eeg(h5fn);
    eeg = ecog.readEEG({Start(i) End(i)});
    bip = h5readatt(h5fn,'/h5eeg/eeg','bip');
    [n_bchan,~] = size(bip);
    Fs = round(ecog.fs);
    
    h = figure;
    set(h,'Units','Normalized','Position',[0.8 0 0.2 1]);
    
    AxesH = axes;
    drawnow;
    InSet = get(AxesH, 'TightInset');
    %InSet = InSet * 4;
    set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    
    X = zeros(ecog.n_chan,n_features);
    Pxx = [];
    F = [];
    Ys = zeros(1,n_bchan);
    for j = 1:n_bchan
        V = eeg.data(:,bip(j,1)) - eeg.data(:,bip(j,2));
        V = V - mean(V);
        X(j,:) = features_v2(V, Fs);
        Y = classify_v2(X(j,:));
        Ys(j) = Y;
        
        yoffset = (j - 1) * 150;
        col = Y_col(Y+1,:);
        if ~((Y ~= 0) && (no_show))
            plot(V + yoffset,'color',col,'LineWidth',0.01); hold on;
            %col = [1 1 1];
        end
        %plot(V + yoffset,'color',col); hold on;
        
        [pxx,f] = pwelch(V,round(Fs),[],[],Fs);
       %pxx = log10(pxx);
        Pxx = [Pxx, pxx];
        F = [F, f];
    end
    axis([0 (End(i) - Start(i) + 1) -150 157*150]);
    title(sprintf('GTruth: %i',Label(i)))
    xlabel('Samples')
    ylabel('Voltage (uV)')
    
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
            plot(f,log10(Pxx(:,j2)),'Color',Y_col(4,:));
            hold on;
        elseif (Ys(j2) == 4)
            plot(f,log10(Pxx(:,j2)),'Color',Y_col(5,:));
            hold on;
        elseif (Ys(j2) == 2)
            plot(f,log10(Pxx(:,j2)),'Color',Y_col(3,:));
            hold on;
        elseif (Ys(j2) == 1)
            plot(f,log10(Pxx(:,j2)),'Color',Y_col(2,:));
            hold on;
        end
    end
    axis([min(f) max(f) -5 5]);
    xlabel('Frequency (Hz)')
    ylabel('Magnitude (dB)')
    legend({'max','mean','min'})
    
    return
    % 21 ekg
    %if (i == 35) % randi([1 length(Sub)])
%     if (true)
%         return
%     end
    if (no_show)
        print(h,sprintf('train_noshow/%s_%i_%i_%i',Sub{i},Start(i),End(i),Label(i)),'-depsc')
    else
        print(h,sprintf('train/%s_%i_%i_%i',Sub{i},Start(i),End(i),Label(i)),'-depsc')
    end
    close(h);
    close(h2);
end
