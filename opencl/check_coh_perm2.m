close all;

system('mkdir figures/check_coh_perm2');

h5fname = 'h5/m00030.h5';
ecog = H5eeg(h5fname);
bip = h5readatt(ecog.filename,'/h5eeg/eeg','bip');
w = round(ecog.fs * 60);
start1 = w * 60 + randi([1 120*w]);
V1 = ecog.readEEG({start1 (start1 + w - 1)});
V1 = V1.data;
start2 = start1 + w * 60;
V2 = ecog.readEEG({start2 (start2 + w - 1)});
V2 = V2.data;
fs = round(ecog.fs);

PLI_S = 56;
PLI_E = 64;
PL2_S = (fs - 180) - 2;
PL2_E = (fs - 180) + 2;
PL3_S = 117;
PL3_E = 123;
DEL_S = 0.5;
DEL_E = 3;
THE_S = 3;
THE_E = 8;
ALP_S = 8;
ALP_E = 12;
BET_S = 12;
BET_E = 25;
GAM_S = 25;
GAM_E = 100;
BRO_S = 0.5;
BRO_E = 125;

dist = 'Normal';
bchan1 = 10;
v1b = V1(:,bip(bchan1,1)) - V1(:,bip(bchan1,2));
bchan2 = 20;
intime = true;
if (intime)
    v2b = V1(:,bip(bchan2,1)) - V1(:,bip(bchan2,2));
else
    v2b = V2(:,bip(bchan2,1)) - V2(:,bip(bchan2,2));
end
pd1 = fitdist(v1b,dist);
pd2 = fitdist(v2b,dist);

% Resample random data
resample = false;
if (resample)
    v1b = random(pd1,size(v1b));
    v2b = random(pd2,size(v2b));
end

h = figure;
set(h,'Position',[0 0 1440 900]);

subplot(3,1,1)
T = linspace(start1/ecog.fs,(start1 + w - 1)/ecog.fs,w);
plot(T/3600,v1b,'black')
xlabel('Time (hours)')
ylabel('uV')
axis tight;

subplot(3,1,2)
T = linspace(start2/ecog.fs,(start2 + w - 1)/ecog.fs,w);
plot(T/3600,v2b,'black')
xlabel('Time (hours)')
ylabel('uV')
axis tight;

subplot(3,1,3)
nwindow = round(fs*0.1);
noverlap = round(nwindow*0.95);
win = hamming(nwindow);
%win = hann(nwindow);
%win = blackmanharris(nwindow);
%win = ones(nwindow,1);
%win = rectwin(nwindow);
[pxx,f] = mscohere(v1b,v2b,win,round(noverlap),[],fs);
pxx = sqrt(pxx);
rectangle('Position',[PLI_S 0 (PLI_E - PLI_S) 1],'FaceColor',0.8*[1 1 1],'LineStyle','none'); hold on;
rectangle('Position',[PL2_S 0 (PL2_E - PL2_S) 1],'FaceColor',0.8*[1 1 1],'LineStyle','none'); hold on;
rectangle('Position',[PL3_S 0 (PL3_E - PL3_S) 1],'FaceColor',0.8*[1 1 1],'LineStyle','none'); hold on;
plot(f,pxx,'black'); hold on;
xlabel('Frequency (Hz)')
ylabel('Coherence')
axis([min(f) max(f) 0 1]);
plot([BRO_S BRO_S],[0 1],'--black'); hold on;
plot([BRO_E BRO_E],[0 1],'--black'); hold on;
mask_bro = ((f > BRO_S) & (f < PLI_S)) | ((f > PLI_E) & (f < PL2_S)) | ((f > PL2_E) & (f < PL3_S)) | ((f > PL3_E) & (f < BRO_E));
r_bro = mean(pxx(mask_bro));
plot([BRO_S BRO_E],[r_bro r_bro],':red');

if (resample)
    print(h,sprintf('figures/check_coh_perm2/rand_nwindow-%i_noverlap-%i',nwindow,noverlap),'-depsc');
elseif (intime)
    print(h,sprintf('figures/check_coh_perm2/intime_nwindow-%i_noverlap-%i',nwindow,noverlap),'-depsc');
else
    print(h,sprintf('figures/check_coh_perm2/nwindow-%i_noverlap-%i',nwindow,noverlap),'-depsc');
end