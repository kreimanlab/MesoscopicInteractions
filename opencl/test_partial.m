close all;

%ifn = '/mnt/cuenap2/scripts/synth/out/m00043.h5';
ifn = '~/data/h5eeg/artifact_removed/test/opencl/h5/m00026.h5';
n_chan = h5readatt(ifn,'/h5eeg/eeg','n_chan');
fs = h5readatt(ifn,'/h5eeg/eeg','rate');
n_samples = h5readatt(ifn,'/h5eeg/eeg','n_samples');
bip = h5readatt(ifn,'/h5eeg/eeg','bip');
[n_bchan,~] = size(bip);

win = 60;


si = randi([1 n_samples-round(fs*win)]);
w = h5read(ifn,'/h5eeg/eeg',[1 si],[n_chan round(win*fs)])';

%si = randi([1 n_samples-round(fs*w)]);
si = 1;
v = h5read(ifn,'/h5eeg/eeg',[1 si],[n_chan round(win*fs)])';

vt = zeros(round(win*fs),n_bchan);
wt = zeros(round(win*fs),n_bchan);
for i = 1:n_bchan
    vt(:,i) = v(:,bip(i,1)) - v(:,bip(i,2));
    wt(:,i) = w(:,bip(i,1)) - w(:,bip(i,2));
end
v = vt;
w = wt;

r_bf = nan(n_bchan,n_bchan);
r_cov = nan(n_bchan,n_bchan);
r_bfP = nan(n_bchan,n_bchan);
r_covP = nan(n_bchan,n_bchan);

% Spearman
v = tiedrank(v);
w = tiedrank(w);

% covariance method
tic;
ci = inv(cov(v));
for i = 1:(n_bchan-1)
    for j = (i+1):n_bchan
        r_cov(i,j) = - ci(i,j)/sqrt(ci(i,i)*ci(j,j));
    end
end
t_cov = toc;
fprintf('covariance: %.4f seconds\n',t_cov)

% covariance method perm
tic;
ci = inv(cov([v,w]));
for i = 1:(n_bchan-1)
    for j = (i+1):n_bchan
        r_covP(i,j) = - ci(i,j+n_bchan)/sqrt(ci(i,i)*ci(j+n_bchan,j+n_bchan));
    end
end
t_cov = toc;
fprintf('covariance, perm: %.4f seconds\n',t_cov)

fprintf('%.8f\n',mean(r_covP(~isnan(r_covP))))
fprintf('%.8f\n',var(r_covP(~isnan(r_covP))))

return
% brute force method
tic;
for i = 1:(n_bchan-1)
    for j = (i+1):n_bchan
        cind = true(1,n_bchan);
        cind(i) = false;
        cind(j) = false;
        r_bf(i,j) = partialcorr(v(:,i),v(:,j),v(:,cind));
    end
end
t_bf = toc;
fprintf('brute force: %.4f seconds\n',t_bf)
fprintf('[*] Average corr coeff error: %.4f\n',sum(nansum(abs(r_bf - r_cov))) / nchoosek(n_bchan,2))


% brute force method
tic;
for i = 1:(n_bchan-1)
    for j = (i+1):n_bchan
        cind = true(1,n_bchan);
        cind2 = true(1,n_bchan);
        cind(i) = false;
        cind2(j) = false;
        r_bfP(i,j) = partialcorr(v(:,i),w(:,j),[v(:,cind),w(:,cind2)]);
    end
end
t_bf = toc;
fprintf('brute force: %.4f seconds\n',t_bf)
fprintf('[*] Average corr coeff error perm: %.4f\n',sum(nansum(abs(r_bfP - r_covP))) / nchoosek(n_bchan,2))