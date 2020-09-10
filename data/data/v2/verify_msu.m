close all;

load_annot;

w = 10;

if (w > 10)
    fprintf(2,'W> Sampling rate drift affects artifact stride - graph stride length mismatch.\n')
end

% overwrite paths
[~,host] = system('hostname');
if contains(host,'ubuntu_1604')
    h5dir = '/nas_share/RawData/scripts/synth/out';
    %artdir = '/nas_share/RawData/scripts/synth/out_art';
    %artdir = '/nas_share/RawData/scripts/synth/h5_art_25uv';
    artdir = '/media/klab/KLAB101/h5_art_opti';
    resultsdir = '/nas_share/RawData/data/results/coh_w10';
elseif contains(host,'hopper')
    %h5dir = '/mnt/cuenap2/scripts/synth/out';
    h5dir = '/media/klab/KLAB101/h5_notch20';
    %artdir = '/mnt/cuenap2/scripts/synth/out_art_o2';
    %artdir = '/media/klab/D0BA8713BA86F56E/data/h5_art_50uv';
    artdir = '/mnt/cuenap2/scripts/synth/out_art';
    resultsdir = '/mnt/cuenap2/data/results/coh_w10';
    
elseif contains(host,'o2.rc.hms.harvard.edu')
    %h5dir = '/n/groups/kreiman/jerry/data/scripts/v2';
    h5dir = '/n/scratch2/jw324/data/h5_notch20';
    artdir = '/n/groups/kreiman/jerry/data/h5_art';
    resultsdir = '/n/groups/kreiman/jerry/data/results/coh_w10';
end

n_features = 2;
toggle_plot = true;

%USub = USub(randperm(length(USub)));
%USub = USub(1:(end-1));
USub = USub(end);

%parfor i = 1:length(USub)
for i = 1:1
%for i = 1:length(USub)
    
    Y_col = [];
    Y_col(1,:) = [0 0 0];
    %Y_col(2,:) = [1 0 0];
    Y_col(2,:) = [0 1 1];
    Y_col(3,:) = [158, 92, 0]/255;
    Y_col(4,:) = [0 1 0];
    Y_col(5,:) = [0 0 1];
    Y_col(6,:) = [126, 0, 224]/255;

    
    sid = USub{i};
    artfn = sprintf('%s/%s_art.h5',artdir,sid);
    h5fn = sprintf('%s/%s.h5',h5dir,sid);
    ecog = H5eeg(h5fn);
    graphBrfn = sprintf('%s/%s_graph_pcBroadband.h5',resultsdir,sid);
    %distBrfn = sprintf('%s/%s_dists-pcBroadband-10000.mat',resultsdir,sid);
    
    Art = h5read(artfn,'/artifacts');
    n_arts = h5readatt(artfn,'/artifacts','n_arts');
    art_width = h5readatt(artfn,'/artifacts','width');
    graph_width = round(ecog.fs*w);
    
    fprintf('(%i) %s] art_width: %.8f graph_width: %.8f\n',i,sid,art_width,graph_width)
    
    if (toggle_plot)
        system('mkdir verify');
        system('mkdir verify/art_idx');
        h = figure();
        set(h,'Units','Normalized','Position',[0 0 1 1]);
        set(h,'PaperUnits','inches');
        %set(h,'PaperSize',[20 16]);
        set(h,'PaperPosition',[0 0 30 10]);

        % Axis margins
        AxesH = axes;
        drawnow;
        InSet = get(AxesH, 'TightInset');
        InSet = InSet * 4;
        InSet(1) = InSet(1) * 4; %ytick label offset
        set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])

        % Plot
        imagesc(Art);
        colormap(Y_col(1:(max(Art(:))+1),:));

        % xticks
        T_hr = linspace(0,n_arts/3600,n_arts);
        T_tickval = 1:round(3600*24):n_arts;
        xticks(T_tickval);
        xlab = cell(1,length(T_tickval));
        for j = 1:length(T_tickval)
            xlab{j} = sprintf('%.0fh',T_hr(T_tickval(j)));
        end
        %xticklabels(T_hr(T_tickval));
        xticklabels(xlab);

        % yticks
        bip = h5readatt(h5fn,'/h5eeg/eeg','bip');
        [n_bchan,~] = size(bip);
        chan_labels = h5readatt(h5fn,'/h5eeg/eeg','labels');
        try
            dk_labels = h5readatt(h5fn,'/h5eeg/eeg','labels_2');
        catch
            dk_labels = chan_labels;
        end
        bchan_labels = cell(1,n_bchan);
        set(gca,'TickDir','out');

        for j = 1:n_bchan
            b1 = bip(j,1);
            b2 = bip(j,2);
            try
                bchan_labels{j} = sprintf('%s-%s (%s,%s)',...
                    chan_labels{b1},chan_labels{b2},dk_labels{b1},dk_labels{b2});
            catch
                bchan_labels{j} = sprintf('%s-%s',...
                    chan_labels{b1},chan_labels{b2});
            end
        end
        yticks(1:n_bchan)
        yticklabels(bchan_labels)
        title(sprintf('%s',sid))

        % Calculate usable electrodes
        n_comb = nchoosek(n_bchan,2);
        n_w = round(n_arts/w) - 1;
        art_idx = false(n_comb,n_w);
        for l = 1:n_w
            l_start = (l-1)*w + 1;
            l_end = l_start + w - 1;
            
            k3 = 1;
            for k1 = 1:(n_bchan-1)
                for k2 = (k1+1):n_bchan
                    cond1 = sum(Art(k1,l_start:l_end) ~= 0) > 0;
                    cond2 = sum(Art(k2,l_start:l_end) ~= 0) > 0;
                    if ( cond1 || cond2 )
                    %if ( (Art(k1,l) > 0) || (Art(k2,l) > 0) )
                        art_idx(k3,l) = true;
                    else
                        art_idx(k3,l) = false;
                    end
                    k3 = k3 + 1;
                end
            end
        end
        h2 = figure;
        set(h2,'Units','Normalized','Position',[0 0 1 1]);
        set(h2,'Units','Normalized','Position',[0 0 1 1]);
        set(h2,'PaperUnits','inches');
        %set(h2,'PaperSize',[20 16]);
        set(h2,'PaperPosition',4*[0 0 30 10]);
        % Axis margins
        AxesH = axes;
        drawnow;
        InSet = get(AxesH, 'TightInset');
        InSet = InSet * 4;
        InSet(1) = InSet(1) * 4; %ytick label offset
        set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
        
        imagesc(art_idx);
        colormap bone;
        xticks(T_tickval/w);
        xticklabels(xlab);

        comb_labels = cell(1,n_bchan);
        j = 1;
        for k1 = 1:(n_bchan-1)
            for k2 = (k1+1):n_bchan
                c1b1 = bip(k1,1);
                c1b2 = bip(k1,2);
                c2b1 = bip(k2,1);
                c2b2 = bip(k2,2);
                comb_labels{j} = sprintf('%s-%s %s-%s',...
                chan_labels{c1b1},chan_labels{c1b2},chan_labels{c2b1},chan_labels{c2b2});
                j = j + 1;
            end
        end
        yticks(1:n_comb)
        yticklabels(comb_labels)
        set(gca,'TickDir','out');
        frac_art = sum(sum(art_idx ~= 0))/numel(art_idx);
        title(sprintf('%s - %.2f%% artifacts, w = %i',sid,100*frac_art,w))
        
        % Save figs
        out_fmt = '-depsc';
        print(h,sprintf('verify/%s_artifacts',sid),out_fmt,'-r400')
        print(h2,sprintf('verify/%s_artifacts_comb',sid),out_fmt,'-r400')
        close(h);
        close(h2);
        save_verify(sid,Art,Y_col,art_idx,n_bchan,w,art_width,graph_width,bip,chan_labels,bchan_labels,comb_labels,T_tickval,xlab);
        %save(sprintf('verify/art_idx/%s',sid),'-v6','art_idx','n_bchan','w','art_width','graph_width','bip','chan_labels','bchan_labels','comb_labels','T_tickval','xlab')
        %return
    end
    
end

function [ ] = save_verify ( sid,Art,Y_col,art_idx,n_bchan,w,art_width,graph_width,bip,chan_labels,bchan_labels,comb_labels,T_tickval,xlab )
    save(sprintf('verify/art_idx/%s',sid),'-v6','Art','Y_col','art_idx','n_bchan','w',...
        'art_width','graph_width','bip','chan_labels','bchan_labels',...
        'comb_labels','T_tickval','xlab');
end

% 
% for l = 1:1
%     sid = USub{randi([1 length(USub)])};
%     seconds_view = 10; % number of seconds to view
%     h5fn = sprintf('%s/%s.h5',h5dir,sid);
%     ecog = H5eeg(h5fn);
%     Fs = round(ecog.fs);
%     i = 1;
%     Start = randi([1 (ecog.n_samples - (Fs * seconds_view))]);
%     End = Start + (Fs * seconds_view) - 1;
%     Label = 2;
% 
%     eeg = ecog.readEEG({Start(i) End(i)});
%     bip = h5readatt(h5fn,'/h5eeg/eeg','bip');
%     [n_bchan,~] = size(bip);
% 
%     h = figure;
%     set(h,'Units','Normalized','Position',[0 0 1 1]);
% 
%     AxesH = axes;
%     drawnow;
%     InSet = get(AxesH, 'TightInset');
%     InSet = InSet * 4;
%     set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
% 
%     X = zeros(ecog.n_chan,n_features);
%     T = linspace(0,seconds_view,(Fs * seconds_view));
%     n_seg = 0;
%     n_mark = 0;
%     A = zeros(n_bchan,seconds_view);
%     meanV = zeros(1,n_bchan);
%     for k = 1:seconds_view
%         for j = 1:n_bchan
%             V = eeg.data(:,bip(j,1)) - eeg.data(:,bip(j,2));        
%             k2 = (k-1)*Fs + 1;
%             V = V(k2:(k2+Fs-1));
%             if (k == 1)
%                 meanV(j) = mean(V);
%             end
%             Vsav = V;
%             V = V - mean(V);
%             Vp = Vsav - meanV(j);
%             X(j,:) = features_v2(V, Fs);
%             Y = classify_v2(X(j,:));
% 
%             yoffset = (j - 1) * 150;
%             col = Y_col(Y+1,:);
%             if ~((Y ~= 0) && (no_show))
%                 plot(T(k2:(k2+Fs-1)),Vp + yoffset,'color',col,'LineWidth',0.01); hold on;
%                 %col = [1 1 1];
%             end
%             %plot(V + yoffset,'color',col); hold on;
% 
%             A(j,k) = Y;
%             if (Y ~= 0)
%                 n_mark = n_mark + 1;
%             end
%             n_seg = n_seg + 1;
%         end
% 
%     end
% 
%     axis([min(T) max(T) -150 157*150]);
%     As = sum(A ~= 0,2);
%     title(sprintf('%s - %i BIP : %.0f %% marked (0 - %.0f %%, 1 - %.0f %%, 2 - %.0f %%)',...
%         sid,n_bchan,100*(n_mark/n_seg),100*sum(As > 0)/n_bchan,100*sum(As > 1)/n_bchan,100*sum(As > 2)/n_bchan))
%     xlabel('Seconds')
%     ylabel('Voltage (uV)')
%     %return
% 
%     % 21 ekg
%     %if (i == 35) % randi([1 length(Sub)])
%     %     if (true)
%     %         return
%     %     end
%     if (no_show)
%         print(h,sprintf('test_noshow/%s_%i_%i_%i',sid,seconds_view,Start(i),End(i)),out_fmt)
%     else
%         print(h,sprintf('test/%s_%i_%i_%i',sid,seconds_view,Start(i),End(i)),out_fmt)
%     end
%     close(h);
%     
% end
