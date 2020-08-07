close all;

h = figure;
set(h,'Position',[0 0 1920 1080])
%A = repmat(1:104,104,1)/104 + 0.2*(0.5*rand(104,104)-0.5);
A = repmat(1:104,104,1)/104 + normrnd(0,0.01,104,104);
imagesc(A);
colorbar;
caxis([0 1]);
%cc = norm_luminance(corrcmap(100));
%cc = norm_luminance(parula(100));
cc = corrcmap(100);
colormap(cc);