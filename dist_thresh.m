close all;
clear;
rng shuffle;
fig_fmt = '-depsc';

if (~ exist('dist_thresh.mat','file'))

metric = 'pcBroadband';
n_perm = 10000;

% Fast i/o definitions
dir_art = '../data/h5_notch20/art_nosz';
dir_res = '../opencl/results';
dir_cor = '../data/coregistration';
dir_sub = '../data/coregistration';

% Slow i/o definitions
dir_h5 = '/media/klab/KLAB101/h5_notch20';

metrics = {'pcBroadband','pcDelta','pcTheta','pcAlpha','pcBeta','pcGamma'};

% Patients
Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48',...
    'mSu'};

% Exclude monkey
Subjects = Subjects(1:(end-1));

% Shuffle
Subjects = Subjects(randperm(length(Subjects)));

system('mkdir figures');
system('mkdir figures/dist_thresh');

Elec = {};
Status = [];
Distance = [];
count = 1;

% Main loop
for i = 1:length(Subjects)
    sid = Subjects{i};
    fn_art = sprintf('%s/%s_art.h5',dir_art,sid);
    fn_dist = sprintf('%s/%s_dists-%s-%i.mat',dir_res,sid,metric,n_perm);
    fn_graph = sprintf('%s/%s_graph-%s.h5',dir_res,sid,metric);
    fn_h5 = sprintf('%s/%s.h5',dir_h5,sid);
    fn_coreg = sprintf('%s/%s/label/all_parcellation.mat',dir_cor,sid);
    
    % Check if files exist
    ckf = {fn_art,fn_dist,fn_graph,fn_h5,fn_coreg};
    for j = 1:length(ckf)
        if (~exist(ckf{j},'file'))
            fprintf(2,'E> File not found: %s\n',ckf{j});
            return
        end
    end
    
    % Init H5eeg
    ecog = H5eeg(fn_h5);
    %Dmat = ecog.getDist();
    
    % Get labels
    %chan_labels = h5readatt(fn_h5,'/h5eeg/eeg','labels');
    C1 = load(sprintf('%s/%s/label/all_parcellation.mat',dir_cor,sid));
    chan_labels = C1.EleLabels;
    bchan_labels = cell(ecog.n_bchan,1);
    for ii0 = 1:ecog.n_bchan
        bchan_labels{ii0} = sprintf('%s-%s',chan_labels{ecog.bip(ii0,1)},chan_labels{ecog.bip(ii0,2)});
    end
    
    % Build istance matrix
    Dmat = zeros(ecog.n_bchan,ecog.n_bchan);
    for i2 = 1:(length(Dmat)-1)
        for j2 = (i2+1):length(Dmat)
            a = reshape(ecog.bip(i2,4:end),2,[]);
            b = reshape(ecog.bip(j2,4:end),2,[]);
            dist2 = zeros(2,2);
            min_current = Inf;
            for a1 = 1:2
                for b1 = 1:2
                    dist2(a1,b1) = sqrt(sum((a(a1,:) - b(b1,:)).^2));
                    if (min(dist2(a1,b1)) < min_current)
                        min_current = dist2(a1,b1);
                        c1 = ecog.bip(i2,a1);
                        c2 = ecog.bip(j2,b1);
                    end
                end
            end
            d = min_current;
            %d = min(dist2(:));
            Dmat(i2,j2) = d;
            Dmat(j2,i2) = d;
            
            % Get grid name
            idx1 = regexp(chan_labels{c1},'\d');
            idx2 = regexp(chan_labels{c2},'\d');
            
            c1_grid = chan_labels{c1}(1:(idx1(1) - 1));
            c2_grid = chan_labels{c2}(1:(idx2(1) - 1));
            c1_num = str2num(chan_labels{c1}(idx1));
            c2_num = str2num(chan_labels{c2}(idx2));
            
            % Check if same grid
            if (strcmp(c1_grid, c2_grid))
                % compute whether they are neighbors, diagonal neighbors, or not
                
                % load grid dimensions
                ifn = sprintf('%s/%s/elec_recon/%s.mgrid',dir_sub,sid,sid);
                fh = fopen(ifn,'r');
                dat = textscan(fh,'%s');
                fclose(fh);
                dat = dat{1};
                grids = {};
                dims = [];
                for k = 1:length(dat)
                    if (strcmp(dat{k},'#Description'))
                        gd = strsplit(dat{k+1},'_');
                        gd = gd{end};
                        grids = [grids {gd}]; 
                    end
                    if (strcmp(dat{k},'#Dimensions'))
                        dim1 = str2num(dat{k+1});
                        dim2 = str2num(dat{k+2});
                        dims = [dims; [dim1 dim2]];
                    end
                end
                grids = grids(2:end)';
                
                % find grid
                idx_grid = find(strcmp(grids,c1_grid),1);
                grid = grids{idx_grid};
                dim = dims(idx_grid,:);
                
                % Are they immediate neighbors?
                % neighbor if: electrodes are consecutive, and the
                % preceding electrode [min] does is not the last electrode
                % in the long dimension [dim(2)]
                is_n1 = (abs(c1_num - c2_num) == 1) & (mod(min([c1_num,c2_num]),dim(2)) ~= 0);
                % or:
                is_n1 = is_n1 | ((max([c1_num,c2_num]) - dim(2)) == min([c1_num,c2_num]));
                
                % Are they diagonal neighbors?
                % check top-right, bottom-left
                is_n2 = ((c1_num + dim(2) + 1) == c2_num);
                is_n2 = is_n2 | ((c2_num + dim(2) + 1) == c1_num);
                % check top-left, bottom-right
                is_n2 = is_n2 | ((c1_num + dim(2) - 1) == c2_num);
                is_n2 = is_n2 | ((c2_num + dim(2) - 1) == c1_num);
                
                status = 0;
                if (is_n1)
                    status = 1;
                end
                if (is_n2)
                    status = 2;
                end
                if (is_n1 && is_n2)
                    status = 3;
                end
                
                %fprintf('%s - %s\t%i\t%.2f mm\n',chan_labels{c1},chan_labels{c2},status,d)
                
                Elec{count} = sprintf('%s_%s_%s',sid,chan_labels{c1},chan_labels{c2});
                Status(count) = status;
                Distance(count) = d;
                count = count + 1;
            end
            
            %return
            
        end
    end
    
    fprintf('(%i/%i) %s\n',i,length(Subjects),sid)
end

% Clear loop indices
clear i;
clear j;

% Print finish message
fprintf('Done.\n')

else
    load('dist_thresh.mat');
end

[c,x] = hist(Distance(Status == 0),300);
[c_n1,~] = hist(Distance(Status == 1),x);
[c_n2,~] = hist(Distance(Status == 2),x);
ymax = 600;
final_thresh = 17;
h = figure;
set(h,'Position',[0 0 0.5*1080 0.25*1080]);
% plot(x,c/trapz(x,c),'black'); hold on;
% plot(x,c_n1/trapz(x,c_n1),'red'); hold on;
% plot(x,c_n2/trapz(x,c_n2),'blue'); hold on;
rectangle('Position',[0 0 final_thresh ymax],'LineStyle','none','FaceColor',0.85*[1 1 1]);
hold all;
plot(x,c,'black'); hold on;
plot(x,c_n1,'black:'); hold on;
plot(x,c_n2,'black--'); hold on;
%plot(final_thresh,ymax,'blackx');
xlabel('Distance (mm)');
ylabel('Number of Pairs')
yticks(0:200:ymax)
axis([min(Distance) max(Distance) 0 ymax])
set(gca,'TickDir','out')
box off;
print(h,sprintf('figures/dist_thresh/distance_histogram_pairs-%i',length(Distance)),fig_fmt);
close(h);

n_rocs = 8*60;
dthreshs = linspace(0,30,n_rocs + 1);
FP = zeros(n_rocs,1);
TP = zeros(n_rocs,1);
for i_d = 1:length(dthreshs)
    %dthresh = 20;
    dthresh = dthreshs(i_d);
    fracTP = sum(Distance((Status==0)) > dthresh)/sum(Status == 0);
    fracFN = sum(Distance((Status==0)) <= dthresh)/sum(Status == 0);

    fracFP = sum(Distance((Status>0)) > dthresh)/sum(Status > 0);
    fracTN = sum(Distance((Status>0)) <= dthresh)/sum(Status > 0);
    
    FP(i_d) = fracFP;
    TP(i_d) = fracTP;
end

% find 20 mm
new_thresh = 17;
[~,i_20] =min(abs(dthreshs - new_thresh));

fprintf('dist_thresh = %.2f\n',new_thresh);
fprintf('TP: %.3f\n',TP(i_20));
fprintf('FP: %.3f\n',FP(i_20));
h = figure;
set(h,'Position',[0 0 0.25*1080 0.25*1080]); 
plot([FP; 0],[TP; 0],'black-'); hold on;
plot([0 1],[0 1],'--black'); hold on;
plot(FP(i_20),TP(i_20),'black.','MarkerSize',10)
xlabel('False positive rate');
ylabel('True positive rate');
box off;
axis([0 1 0 1])
set(gca,'TickDir','out');
print(h,sprintf('figures/dist_thresh/distance_ROI_thresh-%i',new_thresh),fig_fmt);
close(h);

% get final count of dist thresh removed
iM = 1;
n_total_pairs = 0;
n_dist_neg = 0;
n_dist_bthresh = 0;
n_remain = 0;
for i = 1:(length(Subjects))
    %sid = Subjects{i};
    sid = sprintf('sub%i',i);
    fn_cache = sprintf('%s/xsub_out_%s_%i_atl2.mat','./cache',sid,iM);
    Ca = load(fn_cache);
    
    
    dist_neg = Ca.Dmats < 0;
    dist_bthresh = (Ca.Dmats <= Ca.dist_thresh) & (~dist_neg);
    n_total_pairs = n_total_pairs + length(Ca.Dmats);
    n_dist_neg = n_dist_neg + sum(dist_neg);
    n_dist_bthresh = n_dist_bthresh + sum(dist_bthresh);
    
    n_remain = n_remain + sum(Ca.Dmats > Ca.dist_thresh);
    %return
end
fprintf('n_total_pairs: %i\n',n_total_pairs);
fprintf('n_dist_neg: %i\n',n_dist_neg);
fprintf('n_dist_bthresh: %i\n',n_dist_bthresh);
fprintf('n_dist_neg or n_dist_bthresh: %i\n',n_dist_neg + n_dist_bthresh);
fprintf('fraction n_dist_neg: %.6f\n',(n_dist_neg)/n_total_pairs);
fprintf('fraction n_dist_bthresh: %.6f\n',(n_dist_bthresh)/n_total_pairs);
fprintf('fraction both: %.6f\n',(n_dist_neg + n_dist_bthresh)/n_total_pairs);

% 
% dist_thresh
% dist_thresh = 17.00
% TP: 0.921
% FP: 0.018
% n_total_pairs: 148404
% n_dist_neg: 16451
% n_dist_bthresh: 7088
% n_dist_neg or n_dist_bthresh: 23539
% fraction n_dist_neg: 0.110853
% fraction n_dist_bthresh: 0.047762
% fraction both: 0.158614
