close all;
clear;

h5dir = './h5';
Sub1 = {'m00003','m00004','m00005','m00017','m00018','m00019','m00021',...
   'm00022','m00025','m00027','m00028','m00033'};

% Sub1 = {'m00043','m00044','m00045','m00049','m00052','m00053','m00055',...
%     'm00056','m00058','m00059','m00060','m00061'};
Sub2 = {'m00001','m00006','m00019','m00023','m00024','m00026','m00030',...
    'm00032','m00035','m00037','m00038','m00039'};

%n_samp = 10;
n_sub = length(Sub1);
w = 2*240;

Gain = zeros(1,n_sub);
Offs = zeros(1,n_sub);
VarR = zeros(1,n_sub);

for i = 1:n_sub
    s1 = sprintf('%s/%s.h5',h5dir,Sub1{i});
    s2 = sprintf('%s/%s.h5',h5dir,Sub2{i});
    n_samples1 = h5readatt(s1,'/h5eeg/eeg','n_samples');
    n_samples2 = h5readatt(s2,'/h5eeg/eeg','n_samples');
    n_chan1 = h5readatt(s1,'/h5eeg/eeg','n_chan');
    n_chan2 = h5readatt(s2,'/h5eeg/eeg','n_chan');
    w1 = round(w*h5readatt(s1,'/h5eeg/eeg','rate'));
    w2 = w1;%round(w*h5readatt(s2,'/h5eeg/eeg','rate'));
    
    %i1 = randi([1 n_samples1-w1]);
    %i2 = randi([1 n_samples2-w2]);
    
    i1 = randi([1 min([n_samples1-w1,n_samples2-w2])]);
    i2 = i1;
    n_chan = min([n_chan1,n_chan2]);
    
    v1 = h5read(s1,'/h5eeg/eeg',[1 i1],[n_chan w1]);
    v2 = h5read(s2,'/h5eeg/eeg',[1 i2],[n_chan w2]);
    % Subtract mean
    for j = 1:n_chan
        v1(j,:) = v1(j,:) - mean(v1(j,:));
        v2(j,:) = v2(j,:) - mean(v2(j,:));
    end
    
    % Linear regression
    V1 = [ones(length(v1(:)),1) v1(:)];
    b = V1\(v2(:));
    
    Gain(i) = b(2);
    Offs(i) = b(1);
    VarR(i) = range(v1(:))/range(v2(:));
    %return
end
fprintf('Average gain: %.6f uV\n',mean(Gain))
fprintf('Average offset: %.6f uV\n',mean(Offs))
fprintf('Average variance ratio: %.6f\n',mean(VarR))
fprintf('Done.\n')
