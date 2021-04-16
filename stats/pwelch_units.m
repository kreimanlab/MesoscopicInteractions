close all;
clear;

fs = 500; % Hz
width = 10; % seconds
freq_signal = 50; % Hz
t = linspace(0,10,fs*width);
V = 100 * ( sin(t*((freq_signal)*(2*pi))) + normrnd(0,1,[1 fs*width]) ); % uV

h = figure;
subplot(2,1,1)
plot(t,V,'black-');
xlabel('Time (seconds)');
ylabel('Electric Potential (uV)')

subplot(2,1,2)
[Pxx,F] = pwelch(V,fs,round(fs*0.9),[],fs);
plot(F,log(Pxx),'black-');
ylabel('PSD (dB/Hz)');
xlabel('Frequency (Hz)');