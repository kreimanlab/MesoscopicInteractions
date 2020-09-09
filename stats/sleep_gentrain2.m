% Generate training data for sleep
close all; clear;
rng('shuffle');

% EEG trace param
XAXIS_SEC = 30;
FLAG_COPY_VIDEO = true;
FLAG_EXIT_ON_ERR = false;


[~,hname] = system('hostname');

% Paths
if contains(hname,'kraken')
    %dir_h5 = '/media/jerry/KLAB101/h5_notch20';
    %dir_art = '/media/jerry/KLAB101/h5_notch20/art_nosz';
    dir_h5 = '/mnt/cuenap/data/h5_notch20';
    dir_art = '/mnt/cuenap/data/h5_notch20/art_nosz';
    dir_video = '/media/jerry/internal/data/videos';%'/media/klab/KLAB101/videos';
    dir_stamp = '/media/jerry/internal/data/stamps';%'/media/klab/KLAB101/stamps';
    dir_rawdata = '/mnt/symphony_1/RawData';
elseif contains(hname,'ubuntu_1604')
    dir_h5 = '/nas_share/cuenap/data/h5_notch20';
    dir_art = '/nas_share/cuenap/data/h5_notch20/art_nosz';
    dir_video = '/nas_share/cuenap/data/videos';
    dir_stamp = '/nas_share/cuenap/data/stamps';
    dir_rawdata = '/mnt/symphony_1/RawData';
end

% Output path
dir_out = './sleeptrain2';
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
    'sub5';...
    'sub8';...
    'sub9';...
    'sub11';...
    'sub14';...
    'sub15';...
    'sub17';...
    'sub18';...
    'sub19';...
    'sub20';...
    'sub21';...
    'sub22';...
    'sub23';...
    'sub24';...
    'sub25';...
    'sub26';...
    'sub27';...
    'sub28';...
    'sub29';...
    'sub31';...
    'sub32';...
    'sub33';...
    'sub34';...
    'sub35';...
    'sub36';...
    'sub37';...
    'sub38';...
    'sub39';...
    'sub40';...
    'sub41';...
    'sub42';...
    'sub43';...
    'sub44';...
    'sub45';...
    'sub46';...
    'sub47';...
    'sub48';...
};

%Subjects = {'sub21'};
%Subjects = {'sub33'};
Subjects = {'sub40'};

% Shuffle
%Subjects = Subjects(randperm(length(Subjects)));

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
    
    if isnan(XAXIS_SEC)
        w = double(h5readatt(fn_art,'/art_idx','w'));
    else
        w = round(XAXIS_SEC);
    end
    w_samp = w * round(ecog.fs);
    
    % check data has enough no-artifact segments
    if (sum(~art_idx_any) < n_per_sub)
        fprintf(2,'[ERR: %s] too many artifacts, lower n_per_sub.\n',sid);
    end

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
            gotString = true;
        catch e
            fprintf(2,'[!] break: could not get time from sample.\n');
            if (FLAG_EXIT_ON_ERR)
                exit();
            end
            gotString = false;
        end
        
        % Get video
        if (gotString)
            % get video file name
            vidStr = 'VIDEO_NOT_FOUND';
            vidStr = ecog.getVideo(dir_video,samp_i,realTime);
            
            % output name
            vstr = replace(vidStr,',','_');
            vstr = replace(vstr,'.','-');
            ofname = sprintf('%i_%ssec_samp-%i',j,vstr,samp_i);

            fprintf('\t* %s at %s\t(%s)\n',sid,realTime,vidStr);
            
            % Copy video
            if (FLAG_COPY_VIDEO && (~strcmp(vidStr,'VIDEO_NOT_FOUND')))
                %vidStr = 'sub21_4ce6814b_1173.avi,23'
                vidL = strsplit(vidStr,',');
                vidFn = vidL{1}; %sprintf('%s/%s/%s',dir_video,sid,vidL{1});
                vidS = str2double(vidL{2});
                
                % find video file location
                % find /mnt/symphony_1/RawData/sub21 | grep 4ce6814b | grep 1173 | grep "\.avi"
                vidL2 = strsplit(vidFn,'.avi');
                vidL2 = strsplit(vidL2{1},'_');
                
                fprintf('[*] Fetching video file..\n');
                btxt = sprintf('find %s/%s | grep %s | grep %s | grep avi',dir_rawdata,vidL2{1},vidL2{2},vidL2{3});
                [~,vidPath] = system(btxt);
                vidPath = strsplit(vidPath,'\n');
                vidPath = vidPath{1};
                fprintf('\t(-) found: %s\n',vidPath);
                
                % attach video
                vPtr = VideoReader(vidPath);
                vidLeft = vPtr.Duration - vidS;
                
                % Check if second video needs loading
                isVidOverflow = vidLeft <= XAXIS_SEC;
                if (isVidOverflow)
                    fprintf('\t[!] Video cut short: %.1f of %.1f sec remains\n',vidLeft,XAXIS_SEC);
                    %vidSfinish = vptr.Duration;
                else
                    %vidSfinish = vidS + XAXIS_SEC;
                end
                
                out_fmt = 'avi';
                fn_out_vid = sprintf('%s/%s/%s.%s',dir_out,sid,ofname,out_fmt);
                fn_out_aud = sprintf('%s/%s/%s.%s',dir_out,sid,ofname,'wav');

                % Read video
                vPtr.CurrentTime = vidS;
                %currAxes = axes;
                vidIn = {};
                duration = 0;
                while hasFrame(vPtr)
                    vidFrame = readFrame(vPtr);
                    vidIn = [vidIn; {vidFrame}];
                    %image(vidFrame, 'Parent', currAxes);
                    %currAxes.Visible = 'off';
                    %pause(1/v.FrameRate);
                    duration = duration + 1/vPtr.FrameRate;
                    if (duration > XAXIS_SEC)
                        break;
                    end
                end
%                 
%                 NFrames = round(Vptr.FrameRate*Vptr.Duration);
%                 vFrames = 1:NFrames;
%                 vSeconds = vFrames / vPtr.FrameRate;
%                 vIdx = (vSeconds >= vidS) & (vSeconds <= vidSfinish);

                % attach video to write
                fprintf('[!] Writing to: %s..\n',fn_out_vid);
                
                % -- tested working method -------------------------------
%                 Wptr = VideoWriter(fn_out_vid,'Motion JPEG AVI');
%                 Wptr.FrameRate = vPtr.FrameRate;
%                 open(Wptr);
%                 for iF = 1:length(vidIn)
%                     writeVideo(Wptr,vidIn{iF});
%                 end
%                 close(Wptr);

                % --- ffmpeg method ---
                ftxt = sprintf('ffmpeg -v 8 -i "%s" -vcodec copy -ss %.3f -t %.3f %s',vidPath,vidS,duration,fn_out_vid);
                fprintf(['\t',ftxt,'\n']);
                system(ftxt);
                
                ftxtA = sprintf('ffmpeg -v 8 -i "%s" -vn -acodec copy -ss %.3f -t %.3f %s',vidPath,vidS,duration,fn_out_aud);
                fprintf(['\t',ftxtA,'\n']);
                system(ftxtA);
                
                fprintf('\t(-) Success!\n');
            end

            % plot EEG
            h = figure('visible','off');
            set(h,'Position',2*[0 0 1920 1080])
            V = h5read(fn_h5,'/h5eeg/eeg',[1 samp_i],[ecog.n_chan,w_samp]);
            t = linspace(0,w,w_samp);
            cmap = jet(n_cmap);
            chan_labels = h5readatt(ecog.filename,'/h5eeg/eeg','labels');
            bchan_labels = cell(n_bchan,1);
            for k = 1:n_bchan
                bchan_labels{k} = [chan_labels{bip(k,1)},'-',chan_labels{bip(k,2)}];

                v = V(ecog.bip(k,1),:) - V(ecog.bip(k,2),:);
                v = v - mean(v);

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
            f_plot = 1.5;
            t1 = linspace(2,2.5,w_samp);
            plot(t1,voff*sin(f_plot*2*pi*t1) + 3 * voff,'black');
            text(2,4.8*voff,sprintf('Delta(%.1fHz)',f_plot));

            t2 = linspace(3,3.5,w_samp);
            plot(t2,voff*sin(5.5*2*pi*t2) + 3 * voff,'black');
            text(3,4.8*voff,'Theta(5.5Hz)')

            t3 = linspace(4,4.5,w_samp);
            plot(t3,voff*sin(10*2*pi*t3) + 3 * voff,'black');
            text(4,4.8*voff,'Alpha(10Hz)')

            t4 = linspace(5,5.5,w_samp);
            plot(t4,voff*sin(24*2*pi*t4) + 3 * voff,'black');
            text(5,4.8*voff,'Beta(24Hz)')

            t5 = linspace(6,6.5,w_samp);
            plot(t5,voff*sin(35*2*pi*t5) + 3 * voff,'black');
            text(6,4.8*voff,'Gamma(35Hz)')

            set(gca,'TickDir','Out')
            xlabel('Seconds')
            set(gca,'box','off');
            set(gca,'Xgrid','on');
            set(gca,'XTick',linspace(0,w,w*2+1));
            axis([min(t) max(t) (-1*(n_bchan+2)*voff) (6.5)*voff]);
            set(gca,'Ytick',[]); % remove y axis

            % stretch axes
            InSet = get(gca, 'TightInset');
            set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])

            

            % print picture
            system(sprintf('mkdir %s/%s',dir_out,sid));
            print(h,sprintf('%s/%s/%s',dir_out,sid,ofname),'-djpeg','-r250');
            fprintf('* wrote to: %s/%s/%s\n',dir_out,sid,ofname)
            close(h);

            % write to log
            system(sprintf('echo "%s" >> %s/%s/log.txt',ofname,dir_out,sid));
            
            %return
        end
        
        
    end
    
end


