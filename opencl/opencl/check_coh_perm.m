close all;

if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        resultsDir = '/media/klab/D0BA8713BA86F56E/data/results';
        h5Dir = '/media/klab/44/h5';
    else
        resultsDir = '/mnt/cuenap2/data/results';
        h5Dir = '/mnt/cuenap2/scripts/synth/out';
    end
end

if (~exist('A','var'))
    A = Analysis(resultsDir,h5Dir);
end

system('mkdir figures/check_coh_perm');

% Pick subject
s_idx = 2;
Fs = round(A.h5eeg.fs{s_idx});

h5fname = A.h5eeg.filenames{s_idx};
w = round(60*A.h5eeg.fs{s_idx});
start_i1 = randi([1 (round(A.h5eeg.n_samples{s_idx}/2) - w)]);
start_i2 = start_i1 + round(A.h5eeg.n_samples{s_idx}/2) - 1;
V1 = h5read(h5fname,'/h5eeg/eeg',[1 start_i1],[A.h5eeg.n_chan{s_idx} w]);
V2 = h5read(h5fname,'/h5eeg/eeg',[1 start_i2],[A.h5eeg.n_chan{s_idx} w]);
bchan1 = randi([1 A.h5eeg.n_bchan{s_idx}]);
bchan2 = bchan1;
while (bchan1 == bchan2)
    bchan2 = randi([1 A.h5eeg.n_bchan{s_idx}]);
end
b1c1 = A.h5eeg.bip{s_idx}(bchan1,1);
b1c2 = A.h5eeg.bip{s_idx}(bchan1,2);
b2c1 = A.h5eeg.bip{s_idx}(bchan2,1);
b2c2 = A.h5eeg.bip{s_idx}(bchan2,2);
start_hour1 = start_i1 / (A.h5eeg.fs{s_idx} * 3600);
start_hour2 = start_i2 / (A.h5eeg.fs{s_idx} * 3600);
T = linspace(0,w/A.h5eeg.fs{s_idx},w);
fig_w = 1920;
fig_h = 1200;
h = figure;
set(h,'Position',[0 0 fig_w fig_h]);

subplot(6,2,1)
plot(T,V1(b1c1,:),'black');
xlabel('Time (sec)')
ylabel(sprintf('%s (uV)',A.h5eeg.chan_labels{s_idx}{b1c1}))
title(sprintf('%s - Hour %.4f',A.h5eeg.subjects{s_idx},start_hour1))

subplot(6,2,5)
v1x = xcorr(V1(b1c1,:),'coeff');
plot(((1:length(v1x)) - w - 1)/A.h5eeg.fs{s_idx},v1x,'black')
xlabel('Offset (sec)')
ylabel('Autocovariance')
title('Autocovariance')
% plot(T,V1(b1c2,:),'black');
% xlabel('Time (sec)')
% ylabel(sprintf('%s (uV)',A.h5eeg.chan_labels{s_idx}{b1c2}))
% title(sprintf('%s - Hour %.4f',A.h5eeg.subjects{s_idx},start_hour1))

subplot(6,2,3)
plot(T,V2(b2c1,:),'black');
xlabel('Time (sec)')
ylabel(sprintf('%s (uV)',A.h5eeg.chan_labels{s_idx}{b2c1}))
title(sprintf('%s - Hour %.4f',A.h5eeg.subjects{s_idx},start_hour2))

subplot(6,2,7)
v2x = xcorr(V2(b2c1,:),'coeff');
plot(((1:length(v2x)) - w - 1)/A.h5eeg.fs{s_idx},v2x,'black')
xlabel('Offset (sec)')
ylabel('Autocovariance')
title('Autocovariance')

subplot(6,2,9)
v12x = xcorr(V1(b1c1,:),V2(b2c1,:),'coeff');
plot(((1:length(v12x)) - w - 1)/A.h5eeg.fs{s_idx},v12x,'black')
xlabel('Offset (sec)')
ylabel('Cross Covariance')
title('Cross Covariance')

subplot(6,2,11)
[Cxx,fc] = mscohere(V1(b1c1,:),V2(b2c1,:),hamming(2*Fs),[],[],Fs);
plot(fc,sqrt(Cxx),'black')
xlabel('Frequency (Hz)')
ylabel('msCoherence')
axis tight
title('mscohere')

subplot(6,2,2)
Fs = round(A.h5eeg.fs{s_idx});
[pxx,f] = pwelch(V1(b1c1,:),2*Fs,[],[],Fs);
plot(f,10*log10(pxx),'black');
xlabel('Frequency (Hz)')
ylabel(sprintf('%s PSD (dB)',A.h5eeg.chan_labels{s_idx}{b1c1}))
axis([min(f) max(f) min(10*log10(pxx)) max(10*log10(pxx))]);
title('Power Spectral Density - Welch''s Method')

subplot(6,2,6)
[pd] = periodogram(xcorr(V1(b1c1,:),'coeff'));
fd = linspace(0,Fs/2,length(pd));
plot(fd,10*log10(pd),'color',0.5*[1 1 1]); hold on;
plot(fd,movmean(10*log10(pd),round(Fs/30)),'black');
xlabel('Frequency (Hz)')
ylabel(sprintf('%s PSD (dB)',A.h5eeg.chan_labels{s_idx}{b1c1}))
%axis([min(f) max(f) min(10*log10(pxx)) max(10*log10(pxx))]);
axis tight;
title('DFT of Autocovariance')
pd11 = pd;

subplot(6,2,4)
Fs = round(A.h5eeg.fs{s_idx});
[pxx,f] = pwelch(V2(b2c1,:),2*Fs,[],[],Fs);
plot(f,10*log10(pxx),'black');
xlabel('Frequency (Hz)')
ylabel(sprintf('%s PSD (dB)',A.h5eeg.chan_labels{s_idx}{b1c1}))
axis tight;
title('Power Spectral Density - Welch''s Method')

subplot(6,2,8)
[pd] = periodogram(xcorr(V2(b2c1,:),'coeff'));
fd = linspace(0,Fs/2,length(pd));
plot(fd,10*log10(pd),'color',0.5*[1 1 1]); hold on;
plot(fd,movmean(10*log10(pd),round(Fs/30)),'black');
xlabel('Frequency (Hz)')
ylabel(sprintf('%s PSD (dB)',A.h5eeg.chan_labels{s_idx}{b2c1}))
%axis([min(f) max(f) min(10*log10(pxx)) max(10*log10(pxx))]);
axis tight;
title('DFT of Autocovariance')
pd22 = pd;

subplot(6,2,10)
v12x = xcorr(V1(b1c1,:),V2(b2c1,:),'coeff');
[pd] = periodogram(v12x);
fd = linspace(0,Fs/2,length(pd));
plot(fd,10*log10(pd),'color',0.5*[1 1 1]); hold on;
plot(fd,movmean(10*log10(pd),round(Fs/30)),'black');
xlabel('Frequency (Hz)')
ylabel(sprintf('%s - %s CSD (dB)',...
    A.h5eeg.chan_labels{s_idx}{b1c1},...
    A.h5eeg.chan_labels{s_idx}{b2c1}))
%axis([min(f) max(f) min(10*log10(pxx)) max(10*log10(pxx))]);
axis tight
title('DFT of Cross Covariance')
pd12 = pd;

subplot(6,2,12)
% x12 = fft(v12x);
% x11 = fft(xcorr(V1(b1c1,:),'coeff'));
% x22 = fft(xcorr(V2(b2c1,:),'coeff'));
% c12 = (x12 .* conj(x12))./(sqrt(x11.^2) .* sqrt(x22.^2));
% c12 = c12(1:length(fd));
n_s = Fs*2;
for i = 1:(n_s/2):w
    i_end = i + n_s-1;
    if (i_end > w)
        break
    end
    win = hamming(n_s)';
    pd11 = fft(V1(b1c1,i:i_end).*win);
    pd22 = fft(V2(b2c1,i:i_end).*win);
    
    if (i == 1)
        x1 = pd11;
        x2 = pd22;
    else
        x1 = sum([x1; pd11]);
        x2 = sum([x2; pd22]);
    end
    
    %pd12 = fft(V1(b1c1,i:i_end).*win,V2(b2c1,i:i_end).*win);
    c12t = (pd11.* conj(pd22))./(sqrt((pd11 .* conj(pd11))) .* sqrt((pd22.* conj(pd22))));
    %c12t = (pd11.* conj(pd22))./(sqrt((pd11 .^ 2)) .* sqrt((pd22.^ 2)));
    %c12t = c12t.^2;
    if (i == 1)
        c12 = c12t(1:ceil(n_s/2));
    else
        c12 = sum( [c12 ; c12t(1:ceil(n_s/2))] );
    end
end

c12 = c12/length(1:(n_s/2):w);
x1 = x1/length(1:(n_s/2):w);
x2 = x2/length(1:(n_s/2):w);
c12t = (x1.* conj(x2))./(sqrt((x1 .* conj(x1))) .* sqrt((x2.* conj(x2))));
%c12t = (x1.* conj(x2))./(sqrt((x1 .^2)) .* sqrt((x2.^2)));


plot(linspace(0,Fs/2,length(c12)),abs(c12),'black') % fd,
xlabel('Frequency (Hz)')
ylabel('msCoherence')
axis tight

print(h,sprintf('figures/check_coh_perm/raw_voltage'),'-depsc');
close(h)


% ---
V1b = V1(b1c1,:) - V1(b1c2,:);
V2b = V2(b2c1,:) - V2(b2c2,:);

h = figure;
set(h,'Position',[0 0 fig_w fig_h]);

subplot(6,2,1)
plot(T,V1b,'black');
xlabel('Time (sec)')
ylabel(sprintf('%s (uV)',A.h5eeg.bip_labels{s_idx}{bchan1}))
title(sprintf('%s - Hour %.4f',A.h5eeg.subjects{s_idx},start_hour1))

subplot(6,2,5)
v1bx = xcorr(V1b,'coeff');
plot(((1:length(v1bx)) - w - 1)/A.h5eeg.fs{s_idx},v1bx,'black')
xlabel('Offset (sec)')
ylabel('Autocovariance')
title('Autocovariance')

subplot(6,2,3)
plot(T,V2b,'black');
xlabel('Time (sec)')
ylabel(sprintf('%s (uV)',A.h5eeg.bip_labels{s_idx}{bchan2}))
title(sprintf('%s - Hour %.4f',A.h5eeg.subjects{s_idx},start_hour2))

subplot(6,2,7)
v2bx = xcorr(V2b,'coeff');
plot(((1:length(v2bx)) - w - 1)/A.h5eeg.fs{s_idx},v2bx,'black')
xlabel('Offset (sec)')
ylabel('Autocovariance')
title('Autocovariance')

subplot(6,2,9)
v12bx = xcorr(V1b,V2b,'coeff');
plot(((1:length(v12bx)) - w - 1)/A.h5eeg.fs{s_idx},v12bx,'black')
xlabel('Offset (sec)')
ylabel('Cross Covariance')
title('Cross Covariance')

subplot(6,2,2)
Fs = round(A.h5eeg.fs{s_idx});
[pxx,f] = pwelch(V1b,2*Fs,[],[],Fs);
plot(f,10*log10(pxx),'black');
xlabel('Frequency (Hz)')
ylabel(sprintf('%s PSD (dB)',A.h5eeg.bip_labels{s_idx}{bchan1}))
axis tight;
title('Power Spectral Density - Welch''s Method')

subplot(6,2,6)
[pd] = periodogram(xcorr(V1b,'coeff'));
fd = linspace(0,Fs/2,length(pd));
plot(fd,10*log10(pd),'color',0.5*[1 1 1]); hold on;
plot(fd,movmean(10*log10(pd),round(Fs/30)),'black');
xlabel('Frequency (Hz)')
ylabel(sprintf('%s PSD (dB)',A.h5eeg.bip_labels{s_idx}{bchan1}))
%axis([min(f) max(f) min(10*log10(pxx)) max(10*log10(pxx))]);
axis tight;
title('DFT of Autocovariance')

subplot(6,2,4)
Fs = round(A.h5eeg.fs{s_idx});
[pxx,f] = pwelch(V2b,2*Fs,[],[],Fs);
plot(f,10*log10(pxx),'black');
xlabel('Frequency (Hz)')
ylabel(sprintf('%s PSD (dB)',A.h5eeg.bip_labels{s_idx}{bchan2}))
axis tight;
title('Power Spectral Density - Welch''s Method')

subplot(6,2,8)
[pd] = periodogram(xcorr(V2b,'coeff'));
fd = linspace(0,Fs/2,length(pd));
plot(fd,10*log10(pd),'color',0.5*[1 1 1]); hold on;
plot(fd,movmean(10*log10(pd),round(Fs/30)),'black');
xlabel('Frequency (Hz)')
ylabel(sprintf('%s PSD (dB)',A.h5eeg.bip_labels{s_idx}{bchan2}))
%axis([min(f) max(f) min(10*log10(pxx)) max(10*log10(pxx))]);
axis tight;
title('DFT of Autocovariance')


subplot(6,2,10)
v12bx = xcorr(V1b,V2b,'coeff');
[pd] = periodogram(v12bx);
fd = linspace(0,Fs/2,length(pd));
plot(fd,10*log10(pd),'color',0.5*[1 1 1]); hold on;
plot(fd,movmean(10*log10(pd),round(Fs/30)),'black');
xlabel('Frequency (Hz)')
ylabel(sprintf('%s - %s CSD (dB)',...
    A.h5eeg.chan_labels{s_idx}{b1c1},...
    A.h5eeg.chan_labels{s_idx}{b2c1}))
%axis([min(f) max(f) min(10*log10(pxx)) max(10*log10(pxx))]);
axis tight
title('DFT of Cross Covariance')

print(h,sprintf('figures/check_coh_perm/bip_voltage'),'-depsc');
close(h)