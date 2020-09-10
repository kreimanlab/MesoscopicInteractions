close all;

%ifn = '/mnt/cuenap2/scripts/synth/out/sub21.h5';
ifn = '~/data/h5eeg/artifact_removed/test/opencl/h5/sub11.h5';
n_chan = h5readatt(ifn,'/h5eeg/eeg','n_chan');
fs = h5readatt(ifn,'/h5eeg/eeg','rate');
n_samples = h5readatt(ifn,'/h5eeg/eeg','n_samples');
bip = h5readatt(ifn,'/h5eeg/eeg','bip');
[n_bchan,~] = size(bip);

win = 60;

si = randi([1 n_samples-round(fs*win)]);
w = h5read(ifn,'/h5eeg/eeg',[1 si],[n_chan round(win*fs)])';

si = randi([1 n_samples-round(fs*win)]);
%si = 1;
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

% sizes
[v_samp,~] = size(v);

% cross correlation method
tic;
r_xcov = zeros((v_samp - 1)*2 + 1);
c = cov(v);
v_std = sqrt(var(v));
iComb = 1;
n_comb = nchoosek(n_bchan,2);
r_coh = zeros();
for i = 1:(n_bchan-1)
    for j = (i+1):n_bchan
        r_cov(i,j) = c(i,j)/(v_std(i)*v_std(j));
        iComb = iComb + 1;
    end
end
t_cov = toc;
fprintf('graph: %.4f seconds\n',t_cov)

fprintf('\tmean: %.8f\n',mean(r_cov(~isnan(r_cov))))
fprintf('\tvariance: %.8f\n',var(r_cov(~isnan(r_cov))))

% Brute force
tic;
for i = 1:(n_bchan-1)
    for j = (i+1):n_bchan
        r_bf(i,j) = fastcorr2(v(:,i),v(:,j));
    end
end
t_bf = toc;
fprintf('brute force: %.4f seconds\n',t_bf)
fprintf('[*] Average corr coeff error: %.4f\n',sum(nansum(abs(r_bf - r_cov))) / nchoosek(n_bchan,2))


% covariance method perm
tic;
ci = inv(cov([v,w]));
for i = 1:(n_bchan-1)
    for j = (i+1):n_bchan
        r_covP(i,j) = - ci(i,j+n_bchan)/sqrt(ci(i,i)*ci(j+n_bchan,j+n_bchan));
    end
end
t_cov = toc;
fprintf('perm: %.4f seconds\n',t_cov)

fprintf('\tmean: %.8f\n',mean(r_covP(~isnan(r_covP))))
fprintf('\tvariance: %.8f\n',var(r_covP(~isnan(r_covP))))
