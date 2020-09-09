close all;
clear;

n = 1000000;
x = normrnd(10,1,1,n);
y = normrnd(10,1,1,n);
z = normrnd(10,1,1,n);
X = normrnd(10,1,50,n);
Y = prod(X);
histogram(log10(Y),'Normalization','probability');


% h = figure;
% subplot(1,3,1);
% histogram(x,'Normalization','probability');
% subplot(1,3,2);
% histogram(y,'Normalization','probability');
% subplot(1,3,3);
% histogram((x.*y .*z),'Normalization','probability');