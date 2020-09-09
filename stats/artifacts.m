% Table S3: Percentage of artifacts removed according to each criterion

close all;
clear;

%dir_art = '/media/klab/internal/data/h5_notch20/art_nosz';
dir_art = '/home/jerry/data/v2/art_szs';

Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

trig_plot = true;
ST = load('SubjectsTimes.mat');
plot_n_page = 8;

iM = 1;
Nbchans = zeros(1,length(Subjects));
Nchans = zeros(1,length(Subjects));
f1 = [];
n_sub = length(Subjects);
n_p = 1;
n_pall = ceil(n_sub/plot_n_page);

n_cmap = 100;
cmap = jet(n_cmap);

for i = 1:n_sub
    sid = Subjects{i};
    Ca = load(sprintf('./cache/xsub_out_%s_%i.mat',sid,iM));
    fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
    a = h5read(fn_art,'/artifacts');
    [~,n_sec] = size(a);
    %art_idx = h5read(fn_art,'/art_idx');
    n = [sum(sum(a == 1))/numel(a);
        sum(sum(a == 2))/numel(a);
        sum(sum(a == 3))/numel(a);
        sum(sum(a == 4))/numel(a);
        sum(sum(a == 5))/numel(a);
        sum(sum(a == 6))/numel(a)];
    f1 = [f1, n];
    
    if (trig_plot)
        n_i = mod(i,plot_n_page);
        if (n_i == 1)
            h = figure('visible','off','PaperUnits','inches','PaperPosition',[0 0 8.5 11]);
        end
        
        
%         starti = ST.SubjectsIndex{i};
%         endi = cumsum(starti-1);
%         endi = [endi(2:end); Inf];
%         startt = ST.SubjectsTimes{i};
%         
%         stride = 60;
%         night = zeros(1,n_sec);
%         for j = 1:length(startt)
%             if (~ isnan(startt{j}))
%                 dt = datetime(startt{j},'format','MM/dd/uuuu HH:mm:ss');
%                 art_range = ceil([starti(j),endi(j)]/Ca.ecog.fs);
%                 art_range(2) = art_range(2) - 1;
%                 if (art_range(2) > n_sec)
%                     art_range(2) = n_sec;
%                 end
%                 for k = art_range(1):stride:art_range(2)
%                     sec_elapsed = k - art_range(1);
%                     currtime = dt + seconds(sec_elapsed);
%                     
%                     nidx = (currtime.Hour / 24);
%                     %ncolor = cmap(round(nidx*n_cmap),:);
%                     
%                     night(k:(k+stride-1)) = nidx;
%                 end
%             end
%         end
%         return
        
        % show
%         subplot(plot_n_page,1,mod(i-1,plot_n_page)+1);
%         imagesc([a; night(1:n_sec)]);
%         colormap('hot');
%         ylabel('Chan');
%         xlabel('Days');
%         title(sid);
        
        % xlabel ticks
        xl = 1:(60*60*24):n_sec;
        xticks(xl);
        xticklabels(1:length(xl));
        
        if (n_i == 0)
            print(h,sprintf('figures/artifacts_%i_of_%i',n_p,n_pall),'-dpng','-r400');
            close(h);
            n_p = n_p + 1;
        end
        
    end
    
end

for i = 1:length(f1(:,1))
    f = f1(i,:);
    fprintf('[%i] mean: %.4f std: %.4f max: %.4f min: %.4d\n',i,mean(f),std(f),max(f),min(f));
end