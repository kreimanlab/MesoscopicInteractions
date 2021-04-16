close all;
clear;

PN = load('pinknoise_tmp');

[n_bchan,n_w] = size(PN.Coh);
for i = 2:n_bchan
    PN.Coh(i,:) = PN.Coh(1,randperm(n_w));
end

h = figure('visible','off','position',[0 0 400 400]);
imagesc(PN.Coh);
xlabel('10 Seconds');
ylabel('Simulated Bipolar Pair');
cb = colorbar;
colormap(inferno(64));
cb.Label.String = 'Coherence';
print(h,'figures/pinkoise_smallworld_simulated_coh','-dpng');

c = 1;
A = nan(PN.n_nodes,PN.n_nodes);
for i = 1:(PN.n_nodes-1)
    for j = (i+1):PN.n_nodes
        cohval = mean(PN.Coh(c,:)) + normrnd(0,std(PN.Coh(1,:)));
        A(i,j) = cohval;
        A(j,i) = cohval;
        c = c + 1;
    end
end