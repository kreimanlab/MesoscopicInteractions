USE_LAB_FILTER = 1;
Fs = 500;
Fpl = 60;

X = 1:100000;
Y = 500+10*sin(2*pi*Fpl*X/Fs+normrnd(0,0.1,1,length(X)));
[~,V] = filterMaster(Y',X',USE_LAB_FILTER,Fs);
close all;
plot(X,Y,'black');
hold on;
plot(X,V','red');
legend({'raw','filtered'})