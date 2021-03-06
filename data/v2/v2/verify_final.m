close all;
clear;

%load_annot;

% Constants
w = 10;

% Parameters
X1_t1 = [10,25,50,75,100];
X1_t2 = [1000,1250,1500,1750,2000,2500,3000];
X2_t = [200,250,300,350,400,500]/4;
n_opti = length(X1_t1)*length(X1_t2)*length(X2_t);

% Final thresholds
f_x1_t1 = 1; % uV
f_x1_t2 = 5; % uV
f_x2_t = 5; % uV/ms

% overwrite paths
%[~,host] = system('hostname');
%if contains(host,'ubuntu_1604')
%    h5dir = '/nas_share/RawData/scripts/synth/out';
%    %artdir = '/nas_share/RawData/scripts/synth/out_art';
%    %artdir = '/nas_share/RawData/scripts/synth/h5_art_25uv';
%    artdir = '/media/klab/KLAB101/h5_art_opti';
%    resultsdir = '/nas_share/RawData/data/results/coh_w10';
%elseif contains(host,'hopper')
%    %h5dir = '/mnt/cuenap2/scripts/synth/out';
%    h5dir = '/media/klab/KLAB101/h5_notch20';
%    %artdir = '/mnt/cuenap2/scripts/synth/out_art_o2';
%    %artdir = '/media/klab/D0BA8713BA86F56E/data/h5_art_50uv';
%    artdir = '/media/klab/KLAB101/h5_art_opti';
%    resultsdir = '/mnt/cuenap2/data/results/coh_w10';
%elseif contains(host,'o2.rc.hms.harvard.edu')
%    %h5dir = '/n/groups/kreiman/jerry/data/scripts/v2';
%    h5dir = '/n/scratch2/jw324/data/h5_notch20';
%    artdir = '/n/groups/kreiman/jerry/data/h5_art_opti';
%    resultsdir = '/n/groups/kreiman/jerry/data/results/coh_w10';
%end
h5dir = '../h5_notch20';
artdir = '../h5_notch20/art';
%resultsdir = '../results/coh_w10';

n_features = 4;
toggle_sav = true;
toggle_plot = true;
toggle_drawplt = true;
toggle_savplot = true;

%for i = 1:1
%for i = 1:length(USub)
%USub = USub(randperm(length(USub)));
%USub = USub(1:(end-1));

USub = {'example'};
%USub = {'sub48','sub30','sub45','mSu'};
%for i = length(USub) %mSu
for i = 1:length(USub) % par

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
    %graphBrfn = sprintf('%s/%s_graph_pcBroadband.h5',resultsdir,sid);
    %distBrfn = sprintf('%s/%s_dists-pcBroadband-10000.mat',resultsdir,sid);
    n_arts = h5readatt(artfn,'/artifacts','n_arts');
    art_width = h5readatt(artfn,'/artifacts','width');
    graph_width = round(ecog.fs*w);
    bip = h5readatt(h5fn,'/h5eeg/eeg','bip');
    [n_bchan,~] = size(bip);
    system('mkdir verify');
    system('mkdir verify/art_idx');
    
    fprintf('(%i) %s] art_width: %.8f graph_width: %.8f\n',i,sid,art_width,graph_width)
    
    c_i = 1;
    for i1 = 1:length(X1_t1)
        for i2 = 1:length(X1_t2)
            for i3 = 1:length(X2_t)
    
                % select final artifact threshold
                if ((i1 == f_x1_t1) && (i2 == f_x1_t2) && (i3 == f_x2_t))
                    tic;

                    x1_t1 = X1_t1(i1);
                    x1_t2 = X1_t2(i2);
                    x2_t = X2_t(i3);

                    %Art = h5read(artfn,'/artifacts',[(1+(c_i-1)*n_bchan) 1],[n_bchan n_arts]);
                    Art = h5read(artfn,'/artifacts',[1 1],[n_bchan n_arts]);

                    c_i = c_i + 1;

                    if (toggle_plot)


                        if (toggle_drawplt)
                            h = figure('visible','off');
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
                                    bchan_labels{j} = sprintf('%s-%s (%s,%s)',...
                                        chan_labels{b1},chan_labels{b2},chan_labels{b1},chan_labels{b2});
                                end
                            end
                            yticks(1:n_bchan)
                            yticklabels(bchan_labels)
                            title(sprintf('%s',sid))

                        end

                        % Calculate usable electrodes
                        n_comb = nchoosek(n_bchan,2);
                        n_w = round(n_arts/w) - 1;
                        art_idx = false(n_comb,n_w);
                        for l = 1:n_w
                            l_start = (l-1)*w + 1;
                            l_end = l_start + w - 1;

                            % Precompute conditions
                            condp = zeros(n_bchan,1);
                            for k4 = 1:n_bchan
                                condp(k4) = sum(Art(k4,l_start:l_end) ~= 0) > 0;
                            end

                            k3 = 1;
                            for k1 = 1:(n_bchan-1)
                                %ak1 = Art(k1,l_start:l_end);

                                %c1t = Art(k1,l_start:l_end);
                                %cond1 = sum( c1t ~= 0) > 0;

                                for k2 = (k1+1):n_bchan

                                    %cond2 = sum(Art(k2,l_start:l_end) ~= 0) > 0;
                                    %if ( cond1 || cond2 )
                                    %if ( (Art(k1,l) > 0) || (Art(k2,l) > 0) )
    %                                     art_idx(k3,l) = true;
    %                                 else
    %                                     art_idx(k3,l) = false;
    %                                 end

                                    art_idx(k3,l) = (condp(k1) || condp(k2));

                                    %art_idx(k3,l) = sum([ak1,Art(k2,l_start:l_end)]) > 0;
                                    k3 = k3 + 1;
                                end
                            end
                        end
                        frac_art = sum(sum(art_idx ~= 0))/numel(art_idx);
                        t_single = toc;
                        eta_hr = t_single*(n_opti - (c_i-1))/3600;

                        fprintf('%s\t%i\t%i\t%i\t%.1f\t%.2f\t%.2f hr\n',...
                            sid,w,x1_t1,x1_t2,x2_t,100*frac_art,eta_hr)



                        %Debug
    %                     hz = figure;
    %                     imagesc(art_idx);
    %                     colormap gray;

                        if (toggle_drawplt)
                            h2 = figure('visible','off');
                            set(h2,'Units','Normalized','Position',[0 0 1 1]);
                            set(h2,'PaperUnits','inches');
                            %set(h2,'PaperSize',[20 16]);
                            set(h2,'PaperPosition',3*[0 0 30 10]);
                            set(gca,'Fontsize',1);
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
                            title(sprintf('%s - %.2f%% artifacts, w = %i',sid,100*frac_art,w))
                        end

                        % Save figs
                        if (toggle_savplot)
                            out_fmt = '-depsc';
                            %print(h,sprintf('verify/%s_artifacts',sid),out_fmt,'-r400')
                            print(h,sprintf('verify/%s_artifacts',sid),out_fmt)
                            %print(h2,sprintf('verify/%s_artifacts_comb',sid),out_fmt,'-r400')
                            print(h2,sprintf('verify/%s_artifacts_comb',sid),out_fmt)
                        end

                        if (toggle_drawplt)
                            close(h);
                            close(h2);
                        end

                        if (toggle_sav)
                            save_verify(sid,Art,Y_col,art_idx,n_bchan,w,art_width,graph_width,bip,chan_labels,bchan_labels,comb_labels,T_tickval,xlab);
                        end
                    end
    
                else
                    %fprintf('c_i: %i\n',c_i);
                    c_i = c_i + 1;
                end
            end
        end
    end
            
end

function [ ] = save_verify ( sid,Art,Y_col,art_idx,n_bchan,w,art_width,graph_width,bip,chan_labels,bchan_labels,comb_labels,T_tickval,xlab )
    save(sprintf('verify/art_idx/%s',sid),'-v6','Art','Y_col','art_idx','n_bchan','w',...
        'art_width','graph_width','bip','chan_labels','bchan_labels',...
        'comb_labels','T_tickval','xlab');
end
