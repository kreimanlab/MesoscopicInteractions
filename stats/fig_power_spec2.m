close all;
clear;

system('mkdir figures/power_spec2');


Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};
Subjects = {'m00032'};

metric = 'pcBroadband';
metrici = 1;
dir_res = '/media/jerry/untitled/results/coh_w10';
dir_art = '/media/jerry/untitled/h5_notch20/art_nosz';
dir_cache = './cache';
dir_h5 = '/media/klab/untitled/h5_notch20';


P = [];
P0 = [];
F = [];
P_allsub = [];
for i = 1:length(Subjects)
    sid = Subjects{i};
    fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
    fprintf('[*] Reading: %s ..\n',fn_graph);
    fn_cache = sprintf('%s/xsub_out_%s_%i.mat',dir_cache,sid,metrici);
    R = h5read(fn_graph,'/R');
    R = R(:,1:(end-1));
    
    fn_perm = sprintf('%s/%s_perm-%s-10000.h5',dir_res,sid,metric);
    Rp = h5read(fn_perm,'/R');
    Rp = Rp(:,1:(end-1));
    
    Ca = load(fn_cache);
    fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
    fprintf('[*] Reading: %s ..\n',fn_art);
    art_idx = h5read(fn_art,'/art_idx');
    
    fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
    ecog = H5eeg(fn_h5);
    n_samples = h5readatt(fn_h5,'/h5eeg','n_samples');
    ds_factors = h5readatt(fn_h5,'/h5eeg','ds_factors');
    starts = cumsum([1; n_samples]);
    starts = starts(1:(end-1));
    
    % Remove artifacts
    R(art_idx>0) = NaN;
    %Rp(art_idx>0) = NaN;
    
    % Remove nonsignificant
    cond_all = (Ca.Dmats > Ca.dist_thresh) & (Ca.ct > Ca.ct_thresh);
    R = R(cond_all,:);
    Rp = Rp(cond_all,:);
    
    [n_R,n_R2] = size(R);
    %fs = 1/(Ca.w /(3600));
    fs = 1/(Ca.w /(60));
    window = hamming(fs);
    
    for k = 1:length(starts)
        start_idx = starts(k);
        end_idx = start_idx + n_samples(k) - 1;
        if (end_idx > ecog.n_samples)
            end_idx = ecog.n_samples;
        end
        % convert index to w=10 index
        start_idx = floor(start_idx/(Ca.w*round(ecog.fs))) + 1;
        end_idx = floor(end_idx/(Ca.w*round(ecog.fs))) + 1;
        if (end_idx > n_R2)
            end_idx = n_R2;
        end
        
        % Check that there are enough samples in subfile
        cond_enough = length(start_idx:end_idx) >= length(window);
        
        if (cond_enough)
            %return
            for j = 1:(n_R)
                % Interp nans
                x = double(R(j,start_idx:end_idx));
                nanx = isnan(x);
                t = 1:numel(x);
                x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
    %             if (isnan(x(end)))
    %                 x(end) = x(end-1);
    %             end
    %             if (isnan(x(1)))
    %                 x(1) = x(2);
    %             end
                if (~isnan(nanmean(x)))
                    x(isnan(x)) = nanmean(x);
                else
                    x(isnan(x)) = 0;
                    fprintf(2,'[!] Warn, all coherence are NaNs\n.');
                end

                [pxx,f] = pwelch(x,window,[],[],fs);

                P = [P, pxx];

                % perm
                x0 = double(Rp(j,:));
                nanx = isnan(x0);
                t = 1:numel(x0);
                x0(nanx) = interp1(t(~nanx), x0(~nanx), t(nanx));
    %             if (isnan(x0(end)))
    %                 x0(end) = x0(end-1);
    %             end
    %             if (isnan(x0(1)))
    %                 x0(1) = x0(2);
    %             end
                %x0(isnan(x0)) = nanmean(x0);
                if (~isnan(nanmean(x0)))
                    x0(isnan(x0)) = nanmean(x0);
                else
                    x0(isnan(x0)) = 0;
                    fprintf(2,'[!] Warn, all coherence are NaNs\n.');
                end
                [pxx0,f0] = pwelch(x0,window,[],[],fs);

                P0 = [P0, pxx0];

                fprintf('\t(%s) %i of %i , %i of %i\n',sid,k,length(starts),j,n_R)
            end
        
        end
    end
    
    try
        % plot per subject
        h = figure('visible','off','Position',[0,0,600,300]);
        Pp = P; %P - P0;
        pxx = nanmean(Pp,2);
        pxx_s = nanstd(Pp')';
        plot(f,10*log10(real(pxx)),'black-'); hold all;
        plot(f,10*log10(real(pxx + pxx_s)),'black:'); hold on;
        plot(f,10*log10(real(pxx - pxx_s)),'black:'); % )10*log10(
        axis tight;
        xlabel('Frequency (min^-^1)')
        ylabel('Power Spectral Density (dB)')
        box off;
        set(gca,'TickDir','out');
        print(h,sprintf('./figures/power_spec2/%i_%s',i,Subjects{i}),'-depsc');
        print(h,sprintf('./figures/power_spec2/%i_%s',i,Subjects{i}),'-dpng','-r900');
        close(h);

        P = [];
        P_allsub = [P_allsub, P];
    end
end


%%

Pp = P_allsub; %P - P0;

pxx = nanmean(Pp,2);
%pxx = P(:,5*13);
pxx_s = nanstd(Pp')';

% semilogy(f,pxx,'black-'); hold all;
% semilogy(f,pxx+pxx_s,'black:'); hold all;
% semilogy(f,pxx-pxx_s,'black:'); hold all;

h = figure('visible','on','Position',[0,0,600,300]);

plot(f,10*log10(real(pxx)),'black-'); hold all;
plot(f,10*log10(real(pxx + pxx_s)),'black:'); hold on;
plot(f,10*log10(real(pxx - pxx_s)),'black:'); % )10*log10(

axis tight;

xlabel('Frequency (min^-^1)')
ylabel('Power Spectral Density (dB)')
box off;
set(gca,'TickDir','out');

