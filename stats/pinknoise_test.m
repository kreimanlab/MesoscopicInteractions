close all;
clear;
rng default;

fs = 250; % Hz
w = 10; % sec
n_w = 1000;
n_nodes = 20;

x_starts = 1:(fs*w):(fs*w*n_w);
n_pairs = nchoosek(n_nodes,2);
Coh = zeros(n_pairs,n_w);

b = 1;
for c1 = 1:(n_nodes-1)
    for c2 = (c1+1):(n_nodes)
        
        % generate noise
        pinkNoise = dsp.ColoredNoise(1,fs*w*n_w,1);
        x1 = pinkNoise();
        x2 = pinkNoise();

        for i = 1:n_w
            % coherence
            idx_start = x_starts(i);
            idx_stop = idx_start + (fs*w) - 1;
            R = coherence( x1(idx_start:idx_stop), x2(idx_start:idx_stop), fs );
            Coh(b,i) = R(6);
        end
        
        if (b == 1)
            h = figure('Position',[0 0 400 150],'visible','off');
            t = (x_starts)/(fs * 3600);
            plot(t,Coh(1,:),'black.');
            hold all;
            plot([t(1) t(end)],[1*std(Coh(1,:))],'blue--');
            set(gca,'tickdir','out');
            box off;
            xlabel('Hours');
            ylabel('Coherence');
            print(h,'figures/pinknoise_coh_time','-r300','-dpng');
            close(h);
            
            h = figure('Position',[0 0 150 150],'visible','off');
            idx_sig = Coh(1,:) > (mean(Coh(1,:)) + 1*std(Coh(1,:)));
            histogram(Coh(1,:));
            hold on;
            histogram(Coh(1,idx_sig));
            xlabel('Coherence')
            ylabel('Count')
            box off
            set(gca,'tickdir','out');
            print(h,'figures/pinknoise_coh_time_hist','-r300','-dpng');
            close(h);
        end
        b = b + 1;
        %return
    end
end