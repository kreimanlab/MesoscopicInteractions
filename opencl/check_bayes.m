close all;
clear;

if ismac
    resultsDir = '/Volumes/RawData/data/results';
    h5Dir = '/Volumes/RawData/scripts/synth/out';
elseif isunix
    [~,hname] = system('hostname');
    if strcmp(strip(hname),'hopper')
        resultsDir = '/media/klab/44/data/results';
        h5Dir = '/media/klab/44/h5';
    else
        resultsDir = '/mnt/cuenap2/data/results';
        h5Dir = '/mnt/cuenap2/scripts/synth/out';
    end
end

if (~exist('A','var'))
    A = Analysis(resultsDir,h5Dir);
end

%for i = 1:A.r.n_f
i = randi([1 A.r.n_f]);

% take sqrt for coherence
pr_D_H0 = h5read([A.r.results_dir,'/',A.r.perm_f{i}],'/R');
pr_D_H0_d = load([A.r.results_dir,'/',A.r.dist_f{i}],'d');
pr_D = h5read([A.r.results_dir,'/',A.r.graph_f{i}],'/R');

if (startsWith(A.r.metrics{i},'pc') || startsWith(A.r.metrics{i},'sc'))
    pr_D_H0 = sqrt(pr_D_H0);
    pr_D = sqrt(pr_D);
    xmin = 0;
    xmax = 1;
else
    xmin = -1;
    xmax = 1;
end

[n_comb,n_perm] = size(pr_D_H0);
[~,n_graph] = size(pr_D);
pair = randi([1 n_comb]);
[n,x] = hist(pr_D_H0(pair,:),round(n_perm/50));
[n2,x2] = hist(pr_D(pair,:),round(n_graph/50));
f = n/(trapz(x,n));
f2 = n2/(trapz(x2,n2));
figure;
plot(x,f,'-','color',0.5*[1 1 1]); hold on;
plot(x,pdf(pr_D_H0_d.d{pair},x),'-','color',0.3*[1 0.2 0.2]); hold on;
plot(x2,f2,'-','color',0*[1 1 1]);
axis([xmin xmax 0 max([f,f2])]);
xlabel(sprintf('Metric: %s',A.r.metrics{i}));
legend({'P[D|H = H0]','P[D|H = H0] fit','P[D]'})
ylabel('Probability Distribution Function')
title(sprintf('%s, bipolar channel pair %i',A.r.subjects{i},pair))
return
%end
