close all; clear;

fs = 250;
seconds = 5;
mu = 0;
sig = 100;
f1 = 120; %Hz
t = ((1:(seconds*fs)) - 1)/fs;
t = t';
%x = random('Normal',mu,sig,seconds*fs,1) + sig*sin((2*pi*f1)*t);

x = h5read('/mnt/cuenap2/scripts/synth/out/m00023.h5','/h5eeg/eeg',[17 500000],[1 fs*seconds]);

figure;
subplot(2,2,1)
plot(t,x,'black')
xlabel('Seconds')
ylabel('Microvolts')
axis([min(t) max(t) -6*sig 6*sig])

subplot(2,2,2)
WINDOW = fs;
NOVERLAP = round(0.1*fs);
NFFT = fs;
Fs = fs;
[Pxx,F] = pwelch(x,WINDOW,NOVERLAP,NFFT,Fs);
plot(F,10*log10(Pxx),'black')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
axis tight;

subplot(2,2,3)
x2 = xcov(x);
t2 = (-1) * flipud(t);
t2 = [t2(1:(end-1)); t];
plot(t2,x2,'black')
xlabel('Seconds')
ylabel('Autocovariance')
axis tight;

subplot(2,2,4)
WINDOW = fs;
NOVERLAP = round(0.1*fs);
NFFT = fs;
Fs = fs;
[Pxx,F] = pwelch(x2,WINDOW,NOVERLAP,NFFT,Fs);
plot(F,10*log10(Pxx),'black')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
axis tight;
