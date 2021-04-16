close all;
clear;

% n = 1000000;
% x = normrnd(10,1,1,n);
% y = normrnd(10,1,1,n);
% z = normrnd(10,1,1,n);
% X = normrnd(10,1,50,n);
% Y = prod(X);
% histogram(log10(Y),'Normalization','probability');

% h = figure;
% subplot(1,3,1);
% histogram(x,'Normalization','probability');
% subplot(1,3,2);
% histogram(y,'Normalization','probability');
% subplot(1,3,3);
% histogram((x.*y .*z),'Normalization','probability');

w = 10; % seconds
fs = 250; % Hz


% Pinknoise
n_chan = 2;
pinkNoise = dsp.ColoredNoise('Color','pink','SamplesPerFrame',fs*w,'NumChannels',n_chan);


% Linear model
alpha = 0.7;
coeff = 0.5805;

n_perm = 1000;
r = zeros(1,n_perm);
for i = 1:n_perm
    X = pinkNoise();
    x1 = X(:,1);
    x1 = (x1 - mean(x1))/std(x1);
    x2 = X(:,2);
    x2 = (x2 - mean(x2))/std(x2);
    offset = 0;
    error = x2;
    Y = (coeff * x1 + offset)*alpha + error*(1-alpha);
    r(i) = abs(corr(x1,Y));
end
fprintf('mu: %.6f, median: %.6f, std: %.6f\n',mean(r),median(r),std(r));



h = figure;
histogram(r);


% % Signal
% sig_freq = 40; % Hz
% xSig = sin(linspace(0,2*pi*w*sig_freq,w*fs));
% 
% % simulation
% n_perm = 6;
% snr = 0.5;
% Rbr = zeros(1,n_perm);
% Rbr_log = zeros(1,n_perm);
% Alpha = zeros(1,n_perm);
% for i = 1:n_perm
%     %tic;
%     X = pinkNoise();
%     x1 = xSig*(snr) + X(:,1);
%     alpha = rand();
%     %alpha = 1;
%     x2 = xSig*(snr) + X(:,2)*alpha;
%     [R] = coherence(x1,x2,fs);
%     Rbr(i) = R(6);
%     x2_log = xSig*(snr) + X(:,2)*(1/log10(alpha));
%     [R] = coherence(x1,x2_log,fs);
%     Rbr_log(i) = R(6);
%     Alpha(i) = alpha;
%     %t_sing = toc;
%     %fprintf('[*] Finished %i of %i, ETA: %.2f min\n',i,n_perm,(n_perm-i)*t_sing/60);
% end
% 
% % save
% fn_cache = './cache/lognormal_R';
% if (exist([fn_cache,'.mat'],'file'))
%     C = load(fn_cache);
%     Rbr = [C.Rbr, Rbr];
%     Rbr_log = [C.Rbr_log, Rbr_log];
%     Alpha = [C.Alpha, Alpha];
% end
% save(fn_cache,'Rbr','Rbr_log','Alpha');