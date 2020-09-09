close all;

%D = h5read('/media/jerry/internal/data/h5_notch20/sub3.h5','/h5eeg/eeg',[1 1],[1 500]);

ecog = H5eeg('/media/jerry/internal/data/h5_notch20/sub3.h5');

w = 10;
l = round(ecog.fs)*w;
R = ecog.readEEG({1,l});
X = R.data;
c1 = 1;
T = (0:(length(X(:,c1))-1))/(round(ecog.fs));

k = 0.7*w*round(ecog.fs);

% Show original
subplot(4,1,1)
V = X(:,c1);
plot(T,V,'black-');
xlabel('Time (seconds)');
ylabel('uV');
set(gca,'TickDir','out');
box off;
axis tight;

% Show filtered
subplot(4,1,2)
V = X(:,c1);
Vf = movmean(V,k);
plot(T,Vf,'black-');
xlabel('Time (seconds)');
ylabel('uV');
set(gca,'TickDir','out');
box off;
axis tight;


c2 = 13;
subplot(4,1,3)
V2 = X(:,c2);
plot(T,V2,'black-');
xlabel('Time (seconds)');
ylabel('uV');
set(gca,'TickDir','out');
box off;
axis tight;


subplot(4,1,4)
V2 = X(:,c2);
Vf2 = movmean(V2,k);
plot(T,Vf2,'black-');
xlabel('Time (seconds)');
ylabel('uV');
set(gca,'TickDir','out');
box off;
axis tight;

coherence(V,V2,round(ecog.fs))
coherence(Vf,Vf2,round(ecog.fs))