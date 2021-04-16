% Generate training data for sleep

close all; clear;
rng('shuffle');

[~,hname] = system('hostname');

% Paths
if contains(hname,'kraken')
    dir_h5 = '/media/jerry/KLAB101/h5_notch20';
    dir_art = '/media/jerry/KLAB101/h5_notch20/art_nosz';
    dir_video = '/media/jerry/internal/data/videos';%'/media/klab/KLAB101/videos';
    dir_stamp = '/media/jerry/internal/data/stamps';%'/media/klab/KLAB101/stamps';
elseif contains(hname,'ubuntu_1604')
    dir_h5 = '/nas_share/cuenap/data/h5_notch20';
    dir_art = '/nas_share/cuenap/data/h5_notch20/art_nosz';
    dir_video = '/nas_share/cuenap/data/videos';
    dir_stamp = '/nas_share/cuenap/data/stamps';
end

% Output path
dir_out = './sleeptrain';
system(sprintf('mkdir %s',dir_out));

% Constants
n_per_sub = 10*5*12; % number of random samples per patient
spectral_res = 0.5; % Hz
overlap = 0.8; % overlap fraction
voff = 250; % uV - offset voltage for plotting
f_cutoff_min = 2; % hz
f_cutoff_max = 50; % Hz

% program constants
n_cmap = 100;

% Subjects with video list (n = 37)
Subjects = { ...
    'm00019';...
    'm00023';...
    'm00024';...
    'm00026';...
    'm00030';...
    'm00032';...
    'm00035';...
    'm00037';...
    'm00038';...
    'm00039';...
    'm00043';...
    'm00044';...
    'm00045';...
    'm00047';...
    'm00048';...
    'm00049';...
    'm00052';...
    'm00053';...
    'm00055';...
    'm00058';...
    'm00059';...
    'm00060';...
    'm00061';...
    'm00068';...
    'm00071';...
    'm00073';...
    'm00075';...
    'm00079';...
    'm00083';...
    'm00084';...
    'm00095';...
    'm00096';...
    'm00097';...
    'm00100';...
    'm00107';...
    'm00122';...
    'm00124';...
};

%Subjects = {'m00030','m00060'};
% 61
Subjects = {'m00019'};
%Subjects = {'m00023'};

% Shuffle
Subjects = Subjects(randperm(length(Subjects)));

% Subjects loop
for i = 1:length(Subjects)
    sid = Subjects{i};
    fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
    fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
    
    % H5eeg object
    ecog = H5eeg(fn_h5);
    fs = round(ecog.fs);
    bip = ecog.bip;
    [n_bchan,~] = size(bip);
    n_chan = ecog.n_chan;
    
    % read artifacts
    fprintf('[%s] Read artifacts file %s ..\n',sid,fn_art);
    art_idx = h5read(fn_art,'/art_idx');
    art_idx_any = any(art_idx,1); % collapse across pairs
    w = double(h5readatt(fn_art,'/art_idx','w'));
    w_samp = w * round(ecog.fs);
    
    % check data has enough no-artifact segments
    if (sum(~art_idx_any) < n_per_sub)
        fprintf(2,'[ERR: %s] too many artifacts, lower n_per_sub.\n',sid);
    end
    
    % get baseline power spectrum
%     interval_hr = 7; % hours
%     interval = (round(60/w) * interval_hr); 
%     psd_b = [];
%     f_b = [];
%     J = 1:interval:length(art_idx_any);
%     fprintf('[*] Computing baseline power spectrum (n=%i)..\n',length(J));
%     fac = 1/spectral_res;
%     parfor j2 = 1:length(J)
%         j = J(j2);
%         samp_i = (j - 1) * w_samp  + 1;
%         V = h5read(fn_h5,'/h5eeg/eeg',[1 samp_i],[n_chan,w_samp]);
%         % bip
%         Vbip = zeros(n_bchan,w_samp);
%         for k = 1:n_bchan
%             Vbip(k,:) = V(bip(k,1)) - V(bip(k,2));
%         end
%         [Pxx,f] = pwelch(Vbip',round(fs/fac),round(overlap*fs/fac),round(fs/fac),fs);
%         psd_b = [psd_b, mean(Pxx,2)];
%         f_b = [f_b, f];
%         %fprintf('[*] baseline psd %.2f%%\n',100*j/length(art_idx_any));
%     end
%     psd_b1 = mean(psd_b,2);
%     psd_b2 = std(psd_b')';
%     f_b = f_b(:,1);
    
    % trim
%     psd_b1 = psd_b1(2:(end-1));
%     psd_b2 = psd_b2(2:(end-1));
%     f_b = f_b(2:(end-1));

    % sample ignoring artifacts
    samp_bank = nan(1,n_per_sub);
    for j = 1:n_per_sub
        sat = true;
        while (sat)
            % try random sample
            samp = randi([1 length(art_idx_any)]);

            % has artifact and not repeated sample?
            sat = art_idx_any(samp) || any(samp_bank == samp);
        end
        % store sample to avoid repeat sample
        samp_bank(j) = samp;
        
        % convert w-second index to h5 index 
        samp_i = (samp - 1) * w_samp  + 1;
        
        % get time string
        try
            realTime = ecog.getTimeFromSample(dir_stamp,samp_i);
        catch
            fprintf(2,'[!] break: could not get time from sample.\n');
            exit();
        end
        
        % get video file name
        vidStr = 'VIDEO_NOT_FOUND';
        vidStr = ecog.getVideo(dir_video,samp_i,realTime);
        
        fprintf('\t* %s at %s\t(%s)\n',sid,realTime,vidStr);
        
        % plot EEG
        h = figure('visible','off');
        set(h,'Position',2*[0 0 1920 1080])
        V = h5read(fn_h5,'/h5eeg/eeg',[1 samp_i],[ecog.n_chan,w_samp]);
        t = linspace(0,w,w_samp);
        %voff = mean(std(V'));
        
        % gate frequency
%         f_idx = (f_b > f_cutoff_min) & (f_b < f_cutoff_max);
%         cax_min = min(f_b(f_idx));
%         cax_max = max(f_b(f_idx));
        cmap = jet(n_cmap);

        chan_labels = h5readatt(ecog.filename,'/h5eeg/eeg','labels');
        bchan_labels = cell(n_bchan,1);
        for k = 1:n_bchan
            bchan_labels{k} = [chan_labels{bip(k,1)},'-',chan_labels{bip(k,2)}];
            
            v = V(ecog.bip(k,1),:) - V(ecog.bip(k,2),:);
            v = v - mean(v);
            
            
%             % Color by frequency with largest positive deviation
%             [pxx,f] = pwelch(v,round(fs/fac),round(overlap*fs/fac),round(fs/fac),fs);
%             pxx = pxx(2:(end-1));
%             f = f(2:(end-1));
%             
%             % Get color for largest spectral difference
%             ldiff = log(pxx) - log(psd_b1);
%             ldiff = ldiff(f_idx);
%             f2 = f(f_idx);
%             [~,midx] = max(ldiff);
%             fmax = f2(midx);
%             col_idx = round((fmax-cax_min)/(cax_max-cax_min))*(n_cmap - 1) + 1;
            
            % Plot with color
            %plot(t,v - voff*(k-1),'-','Color',cmap(col_idx,:));
            plot(t,v - voff*(k-1),'-black','LineWidth',1)
            hold all;
            text(t(1),(-1)*voff*(k-1),bchan_labels{k},'fontsize',7);
            
        end
        
        
        
        % plot scale
        plot([0 1],[1.5*voff 1.5*voff],'red-','LineWidth',2); hold on;
        text(0.5,2*voff,'1s')
        plot([0 0],[1.5*voff (1.5*voff+4*voff)],'red-','LineWidth',2); hold on;
        text(0,1.5*voff+2*voff,sprintf('%.1fmV',(4*voff/1000)))
        
        % plot frequency scale
        t1 = linspace(2,2.5,w_samp);
        plot(t1,voff*sin(1.5*2*pi*t1) + 3 * voff,'black');
        text(2,4.8*voff,'Delta(1.5)')
        
        t2 = linspace(3,3.5,w_samp);
        plot(t2,voff*sin(5.5*2*pi*t2) + 3 * voff,'black');
        text(3,4.8*voff,'Theta(5.5)')
        
        t3 = linspace(4,4.5,w_samp);
        plot(t3,voff*sin(10*2*pi*t3) + 3 * voff,'black');
        text(4,4.8*voff,'Alpha(10)')
        
        t4 = linspace(5,5.5,w_samp);
        plot(t4,voff*sin(24*2*pi*t4) + 3 * voff,'black');
        text(5,4.8*voff,'Beta(24)')
        
        t5 = linspace(6,6.5,w_samp);
        plot(t5,voff*sin(35*2*pi*t5) + 3 * voff,'black');
        text(6,4.8*voff,'Gamma(35)')
        
        
        set(gca,'TickDir','Out')
        xlabel('Seconds')
%        colormap(cmap)
%        cb = colorbar;
        set(gca,'box','off');
%        set(cb,'TickLength',0);
        set(gca,'Xgrid','on');
        set(gca,'XTick',linspace(0,w,w*2+1));
        axis([min(t) max(t) (-1*(n_bchan+2)*voff) (6.5)*voff]);
        set(gca,'Ytick',[]); % remove y axis
%        caxis([cax_min cax_max]);

        % stretch axes
        InSet = get(gca, 'TightInset');
        set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
        
        
        % output name
        vstr = replace(vidStr,',','_');
        vstr = replace(vstr,'.','-');
        ofname = sprintf('%ssec_samp-%i',vstr,samp_i);
        
        % print picture
        system(sprintf('mkdir sleeptrain/%s',sid));
        print(h,sprintf('sleeptrain/%s/%s',sid,ofname),'-djpeg','-r250');
        fprintf('* wrote to: sleeptrain/%s/%s\n',sid,ofname)
        close(h);
        
        % write to log
        system(sprintf('echo "%s" >> sleeptrain/%s/log.txt',ofname,sid));
        
    end
    

end


