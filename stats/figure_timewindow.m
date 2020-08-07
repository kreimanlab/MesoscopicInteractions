close all;
%clear;

system('mkdir timewindow');
metric = 'pcBroadband';
metrici = 1;
dir_w5 = '/media/jerry/KLAB101/results_res5hz_w5';
dir_w10 = '/media/jerry/KLAB101/results/coh_w10';
dir_w15 = '/media/jerry/KLAB101/results_res5hz_w15';
dir_art = '/media/jerry/KLAB101/h5_notch20/art_nosz2';
dir_cache = 'cache';

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
   'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
   'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',... % 'm00045',
   'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',... % 'm00056',
   'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
   'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

trig_D = true;
D_0 = [];
G = [];
Gs = [];
P = [];
Ps = [];
Diff = [];
count = 1;
Sids = {};
for i = 1:7
    sid = Subjects{i};
    
    try
        
        fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
        fprintf('[*] reading %s\n',fn_art);
        artifacts = h5read(fn_art,'/artifacts');
        fn_cache = sprintf('%s/xsub_out_%s_%i.mat',dir_cache,sid,metrici);
        fprintf('[*] reading %s\n',fn_cache);
        Ca = load(fn_cache);
    
        if (true || (~exist('D_w5','var')))
            if (trig_D)
                % Load fitted dists
                fn_dist_w5 = sprintf('%s/%s_dists-%s-10000.mat',dir_w5,sid,metric);
                fn_dist_w10 = sprintf('%s/%s_dists-%s-10000.mat',dir_w10,sid,metric);
                fn_dist_w15 = sprintf('%s/%s_dists-%s-10000.mat',dir_w15,sid,metric);

                fprintf('[*] reading %s ..\n',fn_dist_w5);
                D_w5 = load(fn_dist_w5);
                fprintf('[*] reading %s ..\n',fn_dist_w10);
                D_w10 = load(fn_dist_w10);
                fprintf('[*] reading %s ..\n',fn_dist_w15);
                D_w15 = load(fn_dist_w15);
            end

            % Load perm
            fn_p_w5 = sprintf('%s/%s_perm-%s-10000.h5',dir_w5,sid,metric);
            fn_p_w10 = sprintf('%s/%s_perm-%s-10000.h5',dir_w10,sid,metric);
            fn_p_w15 = sprintf('%s/%s_perm-%s-10000.h5',dir_w15,sid,metric);

            fprintf('[*] reading %s ..\n',fn_p_w5);
            P_w5 = h5read(fn_p_w5,'/R');
            fprintf('[*] reading %s ..\n',fn_p_w10);
            P_w10 = h5read(fn_p_w10,'/R');
            fprintf('[*] reading %s ..\n',fn_p_w15);
            P_w15 = h5read(fn_p_w15,'/R');

            % Load graphs
            fn_gr_w5 = sprintf('%s/%s_graph-%s.h5',dir_w5,sid,metric);
            fn_gr_w10 = sprintf('%s/%s_graph-%s.h5',dir_w10,sid,metric);
            fn_gr_w15 = sprintf('%s/%s_graph-%s.h5',dir_w15,sid,metric);

            fprintf('[*] reading %s ..\n',fn_gr_w5);
            G_w5 = h5read(fn_gr_w5,'/R');
            fprintf('[*] reading %s ..\n',fn_gr_w10);
            G_w10 = h5read(fn_gr_w10,'/R');
            fprintf('[*] reading %s ..\n',fn_gr_w15);
            G_w15 = h5read(fn_gr_w15,'/R');

            [~,n_gr_w5] = size(G_w5);
            [~,n_gr_w10] = size(G_w10);
            [~,n_gr_w15] = size(G_w15);
            
            % Remove artifacts
            G_w5(artifacts>0) = NaN;
            G_w10(artifacts>0) = NaN;
            G_w15(artifacts>0) = NaN;
            
            % Remove neighbors
            G_w5(Ca.Dmats<=Ca.dist_thresh,:) = NaN;
            G_w10(Ca.Dmats<=Ca.dist_thresh,:) = NaN;
            G_w15(Ca.Dmats<=Ca.dist_thresh,:) = NaN;
            P_w5(Ca.Dmats<=Ca.dist_thresh,:) = NaN;
            P_w10(Ca.Dmats<=Ca.dist_thresh,:) = NaN;
            P_w15(Ca.Dmats<=Ca.dist_thresh,:) = NaN;
            %return
        end

        n_comb = length(D_w5.d);
        %Mus = zeros(3,n_comb);
        for j = 1:n_comb
            if (trig_D)
                D_0(1,count) = D_w5.d{j}.mu;
                D_0(2,count) = D_w10.d{j}.mu;
                D_0(3,count) = D_w15.d{j}.mu;
            end

            P(1,count) = nanmean(P_w5(j,:));
            P(2,count) = nanmean(P_w10(j,:));
            P(3,count) = nanmean(P_w15(j,:));
            Ps(1,count) = nanstd(P_w5(j,:));
            Ps(2,count) = nanstd(P_w10(j,:));
            Ps(3,count) = nanstd(P_w15(j,:));

            G(1,count) = nanmean(G_w5(j,:));
            G(2,count) = nanmean(G_w10(j,:));
            G(3,count) = nanmean(G_w15(j,:));
            Gs(1,count) = nanstd(G_w5(j,:));
            Gs(2,count) = nanstd(G_w10(j,:));
            Gs(3,count) = nanstd(G_w15(j,:));

            count = count + 1;
        end
        
        
    catch e
        fprintf('\t* Skipped: %s\n',sid)
        %rethrow(e)
    end
    
    Sids = [Sids, {sid}];
end
%%

n_bins = 400;
n_xtix = 5;

if (trig_D)
    [fd5,xd5] = hist(D_0(1,:),n_bins); % w5
    [fd10,xd10] = hist(D_0(2,:),n_bins); % w10
    [fd15,xd15] = hist(D_0(3,:),n_bins); % w15
end


[fp5,xp5] = hist(P(1,:),n_bins); % w5
[fp10,xp10] = hist(P(2,:),n_bins); % w10
[fp15,xp15] = hist(P(3,:),n_bins); % w15

[fg5,xg5] = hist(G(1,:),n_bins); % w5
[fg10,xg10] = hist(G(2,:),n_bins); % w10
[fg15,xg15] = hist(G(3,:),n_bins); % w15

dp5 = (G(1,:)-P(1,:)) ./ sqrt(0.5*(Gs(1,:).^2 + Ps(1,:).^2));
dp10 = (G(2,:)-P(2,:)) ./ sqrt(0.5*(Gs(2,:).^2 + Ps(2,:).^2));
dp15 = (G(3,:)-P(3,:)) ./ sqrt(0.5*(Gs(3,:).^2 + Ps(3,:).^2));
[f5,x5] = hist(dp5,n_bins); % w5
[f10,x10] = hist(dp10,n_bins); % w10
[f15,x15] = hist(dp15,n_bins); % w15



%xmin = -0;
%xmax = 0.1;
%xtix = linspace(xmin,xmax,5);

h = figure('visible','on','Position',[0 0 0.5*1080 0.5*1080]);
%[ha, pos] = tight_subplot(3,1,[.1 .1],[.01 .01],[.01 .01]);

%axes(ha(1));
subplot(3,1,1)
% if (trig_D)
%     plot(xd5,fd5/trapz(xd5,fd5),'black:'); hold on;
% else
%     plot(xp5,fp5/trapz(xp5,fp5),'black:'); hold on;
% end
% plot(xg5,fg5/trapz(xg5,fg5),'black-'); hold on;
plot(x5,f5/trapz(x5,f5),'black-'); hold on;
ylabel({'Probability';'w = 5'})
%xticks(xtix);
set(gca,'TickDir','out');
box off;
ax = gca;
xticks(linspace(ax.XLim(1),ax.XLim(2),n_xtix));
%ax.XLim(1) = xmin;
%ax.XLim(2) = xmax;


subplot(3,1,2)
% if (trig_D)
%     plot(xd10,fd10/trapz(xd10,fd10),'black:'); hold on;
% else
%     plot(xp10,fp10/trapz(xp10,fp10),'black:'); hold on;
% end
% plot(xg10,fg10/trapz(xg10,fg10),'black-'); hold on;
plot(x10,f10/trapz(x10,f10),'black-'); hold on;
ylabel({'Probability';'w = 10'})
%xticks(xtix);
set(gca,'TickDir','out');
box off;
ax = gca;
xticks(linspace(ax.XLim(1),ax.XLim(2),n_xtix));
%ax.XLim(1) = xmin;
%ax.XLim(2) = xmax;


subplot(3,1,3)
% if (trig_D)
%     plot(xd15,fd15/trapz(xd15,fd15),'black:'); hold on;
% else
%     plot(xp15,fp15/trapz(xp15,fp15),'black:'); hold on;
% end
% plot(xg15,fg15/trapz(xg15,fg15),'black-'); hold on;
plot(x15,f15/trapz(x15,f15),'black-'); hold on;
xlabel(sprintf('d'''))
ylabel({'Probability';'w = 15'})
%xticks(xtix);
set(gca,'TickDir','out');
box off;
ax = gca;
xticks(linspace(ax.XLim(1),ax.XLim(2),n_xtix));
%ax.XLim(1) = xmin;
%ax.XLim(2) = xmax;

p_name = sprintf('figures/timewindow/nbpairs-%i_mean5-%i_mean10-%i_mean15-%i',...
    length(G(1,:)),round(1000*nanmean(dp5)),round(1000*nanmean(dp10)),round(1000*nanmean(dp15)) );

print(h,p_name,'-depsc');
print(h,p_name,'-dpng','-r300');




%%
% Plot statistics

hp = figure('visible','on','Position',[0 0 0.5*1080 0.5*1080]);

mean5 = nanmean(dp5);
mean10 = nanmean(dp10);
mean15 = nanmean(dp15);
std5 = nanstd(dp5);
std10 = nanstd(dp10);
std15 = nanstd(dp15);
col_bar = 0.5*ones(1,3);
size_cap = 12;

means = [mean5, mean10, mean15];
stds = [std5, std10, std15];
ba = bar(1:3,means); hold on;
ba(1).FaceColor = col_bar;
xticklabels({'w = 5','w = 10','w = 15'});

n_means = length(means);
n_comp = nchoosek(n_means,2);
Dat = [dp5; dp10; dp15];
Dat(:,any(isnan(Dat))) = [];
pval = 0.01;
for i = 1:(n_means-1)
    for j = (i+1):n_means
        %[~,p] = ttest(Dat(i,:),Dat(j,:));
        [~,p] = ranksum(Dat(i,:),Dat(j,:));
        h = (p < (pv, ranal/n_comp));
        fprintf('[*] t-test: %i, %i\t%i\t%.4d\n',i,j,h,p)
    end
end

er = errorbar(1:3,means,stds,stds,...
    'CapSize',size_cap); hold on;

er.LineStyle = 'none';
er.Color = [0 0 0];

set(gca,'TickDir','out');
box off;
ylabel(sprintf('d''')) % Normalized \\DeltaCoherence
ax = gca;
yticks(linspace(ax.YLim(1),ax.YLim(2),3))

p_name = sprintf('figures/timewindow/nbpairs-%i_mean5-%i_mean10-%i_mean15-%i_bar',...
    length(G(1,:)),round(1000*nanmean(dp5)),round(1000*nanmean(dp10)),round(1000*nanmean(dp15)) );

print(hp,p_name,'-depsc');
print(hp,p_name,'-dpng','-r300');
fprintf('[!] Done.\n')

fprintf('[!] Done.\n')