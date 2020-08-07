close all;
clear;

subjects_dir = '/media/jerry/internal/data/coreg';
sid_const = 'fsaverage_sym';
hemi = 'r';

n_col = 2^12;
distort_frac = 0.5; % From 0 to 1, 0 being no distortion, 1 being perfectly distorted

% =========
metrici = 1;
load(sprintf('brainexport/red_6all_fsaverage_%i.mat',metrici));
CaT14 = load(sprintf('./cache/figure_t14_%i_150',metrici));
A = CaT14.Adj_plt2;
A(A==0) = NaN;
            
% Fit to t dist
X = A(~isnan(A));
[f,x] = hist(X,64);
f_pdf = f/(trapz(x,f));
% DistNames = {'Beta','Binomial','BirnbaumSaunders','Burr','Exponential','ExtremeValue',...
%     'Gamma','GeneralizedExtremeValue','GeneralizedPareto','HalfNormal','InverseGaussian',...
%     'Kernel','Logistic','Loglogistic','Lognormal','Nakagami','NegativeBinomial',...
%     'Normal','Poisson','Rayleigh','Rician','Stable','tLocationScale','Weibull'};
% DistSSE = zeros(size(DistNames));
% DistPD = cell(size(DistNames));
% for i = 1:length(DistNames)
%     pd = fitdist(A(~isnan(A)),DistNames{i});
%     res = sum((f_pdf - pdf(pd,x)).^2);
%     DistSSE(i) = res;
%     DistPD{i} = pd;
%     fprintf('[%s] SSE: %.9f\n',DistNames{i},res);
% end


pd = fitdist(A(~isnan(A)),'Kernel');
%pd = DistPD{i};
%load('corrcmap_test-pd_np');
%pd = pd_np;

h = figure;
plot(x,f_pdf,'black-'); hold all;
plot(x,pdf(pd,x),'red-'); hold all;
xlabel('Coherence');
ylabel('Probability');


% colormap function
Eraw = E;
E(:,3) = E(:,3) - center(1);
E(:,4) = E(:,4) - center(2);
E(:,5) = E(:,5) - center(3);

n_A = length(A);


% experimental colormap distortion
distort_x = linspace(min(X),max(X),n_col);
distort_pdf = pdf(pd,distort_x);
%distort_pdf = distort_pdf/trapz(distort_x,distort_pdf); % sum(distort_x.*distort_pdf); %
%plot(distort_x,distort_pdf,'blue--');
distort_idx = ones(1,n_col);
for i = 2:n_col
    sc = trapz(distort_x(1:i),distort_pdf(1:i)) / trapz(distort_x(1:n_col),distort_pdf(1:n_col));
    fprintf('(%i) %.0f\n',i,sc*n_col);
    distort_idx(i) = round(sc*n_col*distort_frac + (1-distort_frac)*i);
end

plot(distort_x(distort_idx),distort_pdf,'blue-'); hold all;
set(gca,'TickDir','out');
box off;
axis tight;
print(h,'figures/corrcmap_test_colorshift','-depsc');
print(h,'figures/corrcmap_test_colorshift','-dpng','-r800');

cmap = corrcmap(n_col);
cmap = cmap(distort_idx,:);

save('./cache/corrcmap_test-cmap','cmap','n_col');

mag_min = nanmin(A(:));
mag_max = nanmax(A(:));
mag2col = @(x) cmap(round(((x - mag_min)/(mag_max - mag_min))*(n_col-1) + 1),:);

h2 = figure;
path_ds_factor = 1;
tube_radius = 0.15;
tube_n_edges = 3; % number of vertices of tube cross-section
rotm = rotationVectorToMatrix([pi/2, 0, 0]);
n = 1;
Paths = Paths(Paths_ind);
NamesPair = {};
for i = 1:(n_A - 1)
    for j = (i+1):n_A
        mag = A(i,j);
        if (~isnan(mag))

            path = Paths{n};

            path = downsample(path,path_ds_factor);

            path(:,1) = path(:,1) - center(1);
            path(:,2) = path(:,2) - center(2);
            path(:,3) = path(:,3) - center(3);


            [n_path,~] = size(path);

            try
                col_tube = mag2col(mag);
            catch
                col_tube = cmap(end,:);
            end
            [X,Y,Z] = tubeplot(path(:,1),path(:,2),path(:,3),tube_radius,ones(1,n_path),tube_n_edges,[0 0 1]);
            pa = surf2patch(X,Y,Z,'triangles');

            surf(X,Y,Z,'FaceColor',col_tube,'EdgeColor','none'); hold on;

        end

        n = n + 1;
    end
end
% =========


% 
% 
% col_pial = 0.6*[1 1 1];
% col_pial_spec = [1 1 1];
% alpha_pial = 0.2;
% 
% [s_vert, faces] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,'pial'));
%             FV = triangulation(faces+1,s_vert(:,1),s_vert(:,2),s_vert(:,3));
%             
% p = trisurf(faces + 1,s_vert(:,1)- center(1),s_vert(:,2)- center(2),s_vert(:,3)- center(3),...
%     'EdgeColor','none','FaceColor',col_pial,'FaceAlpha',alpha_pial); 
% p.AmbientStrength = 0.3 ;
% p.DiffuseStrength = 0.4 ;
% p.SpecularStrength = 0;
% p.SpecularExponent = 1;
% p.BackFaceLighting = 'lit';
% cam_elev = 0;
% camlight(-135,cam_elev);
% camlight(45,cam_elev);
% camlight(-225,cam_elev);
% camlight(-45,cam_elev);

view(90,0);
daspect([1 1 1]);
axis off;
set(gcf,'Position',[0 0 1920 1080]);
colormap(cmap);
cb = colorbar;
set(cb,'Location','east');
set(cb,'TickLength',0);
set(cb,'Ticks',linspace(mag_min,mag_max,5));
cb.Label.String = 'Coherence';
caxis([mag_min mag_max]);