
[~,host] = system('hostname');
if contains(host,'ubuntu_1604')
    h5dir = '/nas_share/cuenap/data/h5_notch20';
else
    %h5dir = '/media/klab/44/h5';
    h5dir = '/media/klab/KLAB101/h5_notch20';
end

of_pre = 'jw';
widthv = 10; % seconds
width = 1; % seconds

% read annotations
infn = sprintf('%s.txt',of_pre);
of = fopen(infn,'r');
D = textscan(of,'%s');
D = D{1};
fclose(of);
n_annot = length(D);
Sub = cell(1,n_annot);
Start = zeros(1,n_annot);
End = zeros(1,n_annot);
Label = zeros(1,n_annot);
for i = 1:n_annot
    d = strsplit(D{i},',');
    Sub{i} = d{1};
    Start(i) = str2double(d{2});
    End(i) = str2double(d{3});
    Label(i) = str2double(d{4});
end

% Subject coverage
USub = unique(Sub);
USubCov = zeros(1,length(USub));
for i = 1:length(USub)
    USubCov(i) = sum(strcmp(USub{i},Sub));
    %fprintf('%s: %i of %i\n',USub{i},USubCov(i),length(Sub));
end
[cov_min,i_min] = min(USubCov);
[cov_max,i_max] = max(USubCov);

% Print info
fprintf('Number of labeled artifacts: %i\n',length(Label))
fprintf('Fraction labeled as artifact: %.3f\n',sum(Label)/length(Label))
fprintf('Number of subjects labeled: %i\n',length(USub))
fprintf('Worst coverage: %s with %i of %i annotations\n',USub{i_min},cov_min,length(Sub))
fprintf('Best coverage: %s with %i of %i annotations\n',USub{i_max},cov_max,length(Sub))