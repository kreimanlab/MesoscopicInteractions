close all;
clear;

mu = 0.2;
sigma = 0.1;

% Normal Sample
h = figure;
AdjV = normrnd(mu,sigma,[1 20000]);
[y,x] = hist(AdjV,130);
y_pdf = y / (trapz(x,y));
histogram(AdjV,'Normalization','pdf','facecolor',0.5*[1 1 1]); hold on;
%plot(x_cen,y_pdf,'black-');
x_cen = linspace(x(1),x(end),100);
y_norm = normpdf(x_cen,(nanmean((AdjV))),(nanstd((AdjV))));
plot(x_cen,y_norm,'black'); hold on;
y_norm_log = lognpdf(x_cen,(nanmean(log(AdjV))),(nanstd(log(AdjV))));
plot(x_cen,y_norm_log,'blue:');
xlabel('Random Variable')
ylabel('PDF')
legend({'Normal Sample','Normal Fit','Lognormal Fit'});

% LogNormal Sample
h = figure;
AdjV = normrnd(mu,sigma,[1 20000]);
AdjV = exp(10*AdjV)/1000;
[y,x] = hist(AdjV,130);
y_pdf = y / (trapz(x,y));
histogram(AdjV,'Normalization','pdf','facecolor',0.5*[1 1 1]); hold on;
%plot(x_cen,y_pdf,'black-');
x_cen = linspace(x(1),x(end),100);
y_norm = normpdf(x_cen,(nanmean((AdjV))),(nanstd((AdjV))));
plot(x_cen,y_norm,'blue:'); hold on;
y_norm_log = lognpdf(x_cen,(nanmean(log(AdjV))),(nanstd(log(AdjV))));
plot(x_cen,y_norm_log,'black-');
xlabel('Random Variable')
ylabel('PDF')
legend({'Lognormal Sample','Normal Fit','Lognormal Fit'});