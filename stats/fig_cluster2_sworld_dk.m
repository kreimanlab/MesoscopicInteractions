close all;
clear;

% The Watts-strogatz definition for small world network is:
%   L >= L_random but C >> C_random
%           ws    - Watts-Strogatz 1998
%           sigma - Humphries-Gurney,2008 (sigma > 1 for small world network)
%           omega - Telesford-Laurienti,2011 (0 to 1: 1 is most small world)
%           I     - Neal,2017 (0 to 1: 1 is most small world)


% Number of random initializations for small world index (default: 12)
n_MC = 12;

% Atlas index
atl = 2; % 2 - Desikan-Killiany

% Load human adjacency matrix
%fn_ca3 = './cache/fig_cluster2_fill_1.mat';
%Ca3 = load(fn_ca3);

fn_ca = './cache/xsub_out_all_1.mat';
Ca = load(fn_ca);
Ahs = nan(Ca.n_rois,Ca.n_rois);
An = Ahs;
Ausub = Ahs;
Ad = Ahs;
for i = 1:Ca.n_rois
    for j = 1:Ca.n_rois
        if (~isempty(Ca.AdjAtl{i,j}))
            lt = Ca.AdjAtl{i,j};
            li = isnan(lt) | (lt == 0);
            Ahs(i,j) = mean(lt(~li));
            An(i,j) = length(Ca.AdjAtlN{i,j});
            Ausub(i,j) = length(unique(Ca.AdjAtl_sid{i,j}));
            Ad(i,j) = mean(Ca.adjct_dist{i,j});
        end
    end
end

% Remove unknown node
rois = Ca.rois;
unk_i = strcmpi(rois,'UNKNOWN');
rois = rois(~unk_i);
Ahs = Ahs(~unk_i,~unk_i);
An = An(~unk_i,~unk_i);
Ausub = Ausub(~unk_i,~unk_i);
Ad = Ad(~unk_i,~unk_i);


% load cache
iM = 1;
fn_ca4 = sprintf('./cache/figure_t14_%i',iM);
Ca4 = load(fn_ca4);
rois2 = Ca4.rois_plt(Ca4.cluster_i);
Ahs = Ca4.Adj_plt2(Ca4.cluster_i,Ca4.cluster_i);

% 
% 
% % Label nodes
% rois2 = cell(size(rois));
% for i = 1:length(rois2)
%     rois2{i} = convertRoiDK(rois{i});
% end

NodeTable = table(rois2,'VariableNames',{'Name'});

% Build node colors
NodeColors = Ca.C.AtlROIs{atl}.LH.table(:,1:3) / (2^8 - 1);
NodeRois = Ca.C.AtlROIs{atl}.LH.struct_names;
NodeRoisDK = cell(size(NodeRois));
for i = 1:length(NodeRois)
    NodeRoisDK{i} = convertRoiDK(NodeRois{i});
end
rois2_col = zeros(length(rois2),3);
for i = 1:length(rois2)
    % Find roi
    sIdx = strcmpi(NodeRoisDK,rois2{i});
    rois2_col(i,:) = NodeColors(sIdx,:);
end

% make graph
Agr = Ahs;
Agr(isnan(Agr)) = 0;
G = graph(Agr,NodeTable,'upper','omitselfloops');


% % network plot
% h = figure;
% Lcolor = [0 0 0];
% hg = plot(G,'Layout','force','EdgeColor',Lcolor,'NodeColor',rois2_col,'NodeLabel',[]); %'Layout','force','EdgeLabel',G.Edges.Weight
% %layout(hg,'force','UseGravity',true);
% %layout(hg,'force','WeightEffect','direct');
% axis off;


% % circular plot
% h = figure;
% addpath('circularGraph');
% circularGraph(Agr,'Label',rois2,'Colormap',rois2_col);


% Get rois2 original names
rois2_raw = cell(size(rois2));
for i = 1:length(rois2)
    sIdx = strcmpi(NodeRoisDK,rois2{i});
    rois2_raw{i} = NodeRois{sIdx};
end

% Brain plot
dir_coreg = '/media/jerry/internal/data/coreg';
h = figure('Position',[0 0 1800 900]); %[0 0 1800 900]
set(h,'PaperUnits','inches');
set(h,'PaperPosition',[0 0 12 6]);

[v,f] = read_surf(sprintf('%s/%s/surf/rh.pial',dir_coreg,'fsaverage_sym'));
[v_a, l, ct] = read_annotation(sprintf('%s/%s/label/rh.aparc.annot',dir_coreg,'fsaverage_sym'));

% Build roi list
col_unknown = 0.7*[1 1 1];
col_nocov = 0.7*[1 1 1];
col_mix = 1*[1 1 1];
col_mix_coeff = 0.8;
rois2_pts = nan(length(rois),3);
VertexColor = repmat(col_unknown,length(v_a),1);
for i = 1:length(rois)
    t_roi = rois{i};
    
    sIdx = strcmpi(ct.struct_names,t_roi);
    a_id = ct.table(sIdx,5);
    
   
    if (sum(strcmp(rois2_raw,t_roi)) == 0)
        % no coverage areas
        a_col = col_nocov;
    else
        % mix color
        a_col = ct.table(sIdx,1:3)/(2^8 - 1);
        a_col = col_mix_coeff * col_mix + a_col * (1 - col_mix_coeff);
        
        % roi midpoint
        rois2_pts(i,:) = median(v(l == a_id,:),1);
        
        % Manually adjust roi center
        if strcmp(rois{i},'lateralorbitofrontal')
            rois2_pts(i,2) = rois2_pts(i,2) + 35; % left-right
            rois2_pts(i,3) = rois2_pts(i,3) - 7;
        elseif strcmp(rois{i},'parsorbitalis')
            rois2_pts(i,2) = rois2_pts(i,2) - 10.5; % left-right
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'parsopercularis')
            rois2_pts(i,2) = rois2_pts(i,2) + 2.5; % left-right
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'parstriangularis')
            rois2_pts(i,2) = rois2_pts(i,2) + 17; % left-right
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'caudalmiddlefrontal')
            rois2_pts(i,2) = rois2_pts(i,2) - 3; % left-right
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'superiorparietal')
            rois2_pts(i,2) = rois2_pts(i,2) - 8;
            rois2_pts(i,3) = rois2_pts(i,3) + 10;
        elseif strcmp(rois{i},'inferiorparietal')
            rois2_pts(i,2) = rois2_pts(i,2) - 19;
            rois2_pts(i,3) = rois2_pts(i,3) + 1;
        elseif strcmp(rois{i},'superiorfrontal')
            rois2_pts(i,2) = rois2_pts(i,2) + 20;
            rois2_pts(i,3) = rois2_pts(i,3) - 15;
        elseif strcmp(rois{i},'isthmuscingulate')
            rois2_pts(i,2) = rois2_pts(i,2) - 6;
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'fusiform')
            rois2_pts(i,2) = rois2_pts(i,2) - 18;
            rois2_pts(i,3) = rois2_pts(i,3) + 4;
        elseif strcmp(rois{i},'superiortemporal')
            rois2_pts(i,2) = rois2_pts(i,2) - 8;
            rois2_pts(i,3) = rois2_pts(i,3) + 6;
        elseif strcmp(rois{i},'middletemporal')
            rois2_pts(i,2) = rois2_pts(i,2) + 22;
            rois2_pts(i,3) = rois2_pts(i,3) - 13;
        elseif strcmp(rois{i},'inferiortemporal')
            rois2_pts(i,2) = rois2_pts(i,2) + 8;
            rois2_pts(i,3) = rois2_pts(i,3) - 15;
        elseif strcmp(rois{i},'lingual')
            rois2_pts(i,2) = rois2_pts(i,2) - 22;
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'precuneus')
            rois2_pts(i,2) = rois2_pts(i,2) - 5;
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'lateraloccipital')
            rois2_pts(i,2) = rois2_pts(i,2) - 10;
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'rostralanteriorcingulate')
            rois2_pts(i,2) = rois2_pts(i,2) + 0;
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'caudalanteriorcingulate')
            rois2_pts(i,2) = rois2_pts(i,2) + 4;
            rois2_pts(i,3) = rois2_pts(i,3) - 0;
        elseif strcmp(rois{i},'frontalpole')
            rois2_pts(i,2) = rois2_pts(i,2) + 0;
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'parahippocampal')
            rois2_pts(i,2) = rois2_pts(i,2) - 9;
            rois2_pts(i,3) = rois2_pts(i,3) + 6;
        elseif strcmp(rois{i},'entorhinal')
            rois2_pts(i,2) = rois2_pts(i,2) + 2.5;
            rois2_pts(i,3) = rois2_pts(i,3) - 1;
        elseif strcmp(rois{i},'temporalpole')
            rois2_pts(i,2) = rois2_pts(i,2) + 2;
            rois2_pts(i,3) = rois2_pts(i,3) - 0;
        elseif strcmp(rois{i},'precentral')
            rois2_pts(i,2) = rois2_pts(i,2) + 3;
            rois2_pts(i,3) = rois2_pts(i,3) + 12;
        elseif strcmp(rois{i},'postcentral')
            rois2_pts(i,2) = rois2_pts(i,2) - 12;
            rois2_pts(i,3) = rois2_pts(i,3) + 30;
        elseif strcmp(rois{i},'supramarginal')
            rois2_pts(i,2) = rois2_pts(i,2) - 16;
            rois2_pts(i,3) = rois2_pts(i,3) + 5;
        elseif strcmp(rois{i},'bankssts')
            rois2_pts(i,2) = rois2_pts(i,2) - 1;
            rois2_pts(i,3) = rois2_pts(i,3) + 0;
        elseif strcmp(rois{i},'transversetemporal')
            rois2_pts(i,2) = rois2_pts(i,2) - 5;
            rois2_pts(i,3) = rois2_pts(i,3) + 5;
        end
    end
    
    VertexColor(l == a_id,:) = repmat(a_col,sum(l == a_id),1);
end


% parameters
TR = triangulation(f+1,v(:,1),v(:,2),v(:,3));
mag_x = max(v(:,1)) - min(v(:,1));
thresh_x = min(v(:,1)) + mag_x * (0.3410); % larger numbers favor bottom
off_z = 0;%150; % 146
off_y = 145;%-30;
layer_x = 100;
col_surf = 0.8*[1 1 1];


% Bottom brain
TR2 = triangulation(f+1,v(:,1)*(-1) + mag_x,v(:,2)*(-1) + off_y,v(:,3) - off_z);
p = trisurf(TR,'EdgeColor','none','FaceVertexCData',VertexColor); % 'FaceColor',col_surf
hold on;

% local lighting
p.AmbientStrength = 0.7 ;
p.DiffuseStrength = 0.2 ;
p.SpecularStrength = 0;
p.SpecularExponent = 0.1;
p.BackFaceLighting = 'lit';
p.FaceLighting = 'gouraud';


% Top brain
hold on;
p = trisurf(TR2,'EdgeColor','none','FaceVertexCData',VertexColor);

% local lighting
p.AmbientStrength = 0.7 ;
p.DiffuseStrength = 0.2 ;
p.SpecularStrength = 0;
p.SpecularExponent = 0.1;
p.BackFaceLighting = 'lit';
p.FaceLighting = 'gouraud';
cam_elev = 0;
camlight(-135,cam_elev);
camlight(45,cam_elev);
camlight(-225,cam_elev);
camlight(-45,cam_elev);

% Plot ROI centers
rois2_front = [1,1,1,0,1, 0,0,1,0,1, 0,1,1,0,0, 1,0,1,0,0,... % paracentral
               0,0,0,1,0, 0,1,1,0,0, 1,1]; % precentral
rois2_degree = zeros(size(rois2));
rois2_pts2 = zeros(length(rois2),3);
AhsD = Ahs;
for i = 1:length(AhsD)
    AhsD(i,i) = NaN;
end
for i = 1:length(rois2)
    sIdx = strcmp(rois,rois2_raw{i});
    %if (rois2_pts(i,1) > thresh_x)
    if (rois2_front(i) == 1)
        hold on;
        rois2_pts2(i,1) = rois2_pts(sIdx,1)+layer_x;
        rois2_pts2(i,2) = rois2_pts(sIdx,2);
        rois2_pts2(i,3) = rois2_pts(sIdx,3);
%         plot3(rois2_pts(sIdx,1)+layer_x,rois2_pts(sIdx,2),rois2_pts(sIdx,3),...
%             '.','MarkerSize',2,'Color',rois2_col(i,:));
    else
        hold on;
        rois2_pts2(i,1) = rois2_pts(sIdx,1)+layer_x;
        rois2_pts2(i,2) = rois2_pts(sIdx,2)*(-1) + off_y;
        rois2_pts2(i,3) = rois2_pts(sIdx,3) - off_z;
%         plot3(rois2_pts(sIdx,1)+layer_x,rois2_pts(sIdx,2)*(-1) + off_y,rois2_pts(sIdx,3) - off_z,...
%             '.','MarkerSize',2,'Color',rois2_col(i,:));
    end
    
    % Calculate node degrees
    rois2_degree(i) = sum((~isnan(AhsD(i,:))) & (AhsD(i,:) ~= 0));
end

bus_offset = 1; % lateral width (not from center)
hbank = [];

% Plot region marker
reg_offset = -0.5; % left-right offset
reg_h = 2 / 2;
col_border_mix = [1 1 1];
col_border_mix_coeff = 0.3;
for i = 1:length(rois2)
    c = rois2_pts2(i,:);
    ca = c(2) + (reg_offset - (rois2_degree(i)/2))*bus_offset;
    cb = c(2) + (rois2_degree(i) - (rois2_degree(i)/2) - 0.5)*bus_offset;
    hold on;
    col_border = rois2_col(i,:) * (1 - col_border_mix_coeff) + col_border_mix * col_border_mix_coeff;
    %plot3([c(1) c(1)],[ca cb],[c(3) c(3)],'-','Color',rois2_col(i,:),'LineWidth',2);

%     plot3([c(1) c(1)],[ca cb],[c(3)-reg_h c(3)-reg_h],'-','Color',col_border,'LineWidth',1);
%     plot3([c(1) c(1)],[ca cb],[c(3)+reg_h c(3)+reg_h],'-','Color',col_border,'LineWidth',1);
%     plot3([c(1) c(1)],[ca ca],[c(3)-reg_h c(3)+reg_h],'-','Color',col_border,'LineWidth',1);
%     plot3([c(1) c(1)],[cb cb],[c(3)-reg_h c(3)+reg_h],'-','Color',col_border,'LineWidth',1);

    p1 = [c(1)+layer_x ca c(3)-reg_h];
    p2 = [c(1)+layer_x ca c(3)+reg_h];
    p3 = [c(1)+layer_x cb c(3)+reg_h];
    p4 = [c(1)+layer_x cb c(3)-reg_h]; 

    x = [p1(1) p2(1) p3(1) p4(1)];
    y = [p1(2) p2(2) p3(2) p4(2)];
    z = [p1(3) p2(3) p3(3) p4(3)];

    fill3(x, y, z, 1, 'EdgeColor','none','FaceColor',col_border);
    
    % Add to horizontal bank
    hbank = [hbank; c(3)-reg_h];
    hbank = [hbank; c(3)+reg_h];
    hbank = [hbank; c(3)];
end

% Plot lines
rois2_busc = zeros(size(rois2));

hstep = 0.8*bus_offset; % all horizontal lines must be this far apart
count = 1;

n_col = 100;
cmap = inferno(n_col);
mag_all = Ahs((~isnan(Ahs)) & (Ahs ~= 0));
mag_max = max(mag_all);
mag_min = min(mag_all);
mag2col = @(x) cmap(round(((x - mag_min)/(mag_max - mag_min))*(n_col-1) + 1),:);

% line width
lw = 1.5;

for i = 1:(length(rois2)-1)
    for j = (i+1):length(rois2)
        % For each line
        mag = Ahs(i,j);
        
        if ((~isnan(mag)) && (mag ~= 0))
            col = mag2col(mag);
%             sIdx1 = strcmp(rois,rois2_raw{i});
%             sIdx2 = strcmp(rois,rois2_raw{j});
            
            % Start and stop coordinates
            c1 = rois2_pts2(i,:);
            c2 = rois2_pts2(j,:);
            
            % Offset by bus count
            c1(2) = c1(2) + (rois2_busc(i) - rois2_degree(i)/2)*bus_offset;
            c2(2) = c2(2) + (rois2_busc(j) - rois2_degree(j)/2)*bus_offset;
            
            % vertical difference
            diff = c2(3) - c1(3);
            
            % horizontal break point;
            horiz = c1(3)+(diff*0.5);
            sat = min(abs(horiz - hbank)) > hstep;
%             if ((rois2_front(i) == 1) && (rois2_front(j) == 1))
%                 rsign = 1;
%             elseif ((rois2_front(i) == 0) && (rois2_front(j) == 0))
%                 rsign = -1;
%             else
%                 % cross-brain behavior
%                 rsign = -1;
%             end
            if (strcmp(rois2_raw{i},'precentral') || strcmp(rois2_raw{j},'precentral'))
                rsign = 1;
            elseif (strcmp(rois2_raw{i},'postcentral') || strcmp(rois2_raw{j},'postcentral'))
                rsign = 1;
            elseif (strcmp(rois2_raw{i},'precuneus') || strcmp(rois2_raw{j},'precuneus'))
                rsign = 1;
            elseif (strcmp(rois2_raw{i},'bankssts') || strcmp(rois2_raw{j},'bankssts'))
                rsign = -1;
            elseif (strcmp(rois2_raw{i},'parsorbitalis') || strcmp(rois2_raw{j},'parsorbitalis'))
                rsign = -1;
            elseif ((c1(3) > 0) && (c1(3) > 0))
                rsign = 1;
            elseif ((c1(3) <= 0) && (c1(3) <= 0))
                rsign = -1;
            else
                rsign = 1;
            end

            while (~sat)
                %rsign = 1 - 2*rand(); %1 - 2*randi([0 1]);
                %horiz = horiz + hstep * rsign;
                horiz = horiz + rsign * hstep * 0.01;
                sat = min(abs(horiz - hbank)) > hstep;
            end
            
            % save horizontal break points
            hbank = [hbank; horiz];
            
            hold on;
            % First vertical line
            v1x = [c1(1) c1(1)];
            v1y = [c1(2) c1(2)];
            v1z = [c1(3) horiz];
            plot3(v1x,v1y,v1z,'-','color',col,'LineWidth',lw);
            
            % Horizontal line
            hx = [c1(1) c2(1)];
            hy = [c1(2) c2(2)];
            hz = [horiz horiz];
            plot3(hx,hy,hz,'.-','color',col,'LineWidth',lw,'MarkerSize',lw*2)
            
            % Second vertical line
            v2x = [c2(1) c2(1)];
            v2y = [c2(2) c2(2)];
            v2z = [c2(3) horiz];
            plot3(v2x,v2y,v2z,'-','color',col,'LineWidth',lw);
            % plot straight line across
            %plot3([coord_1(1) coord_2(1)],[coord_1(2) coord_2(2)],[coord_1(3) coord_2(3)],'-black');
            
            
            % increment bus count
            rois2_busc(i) = rois2_busc(i) + 1;
            rois2_busc(j) = rois2_busc(j) + 1;
            
            if (count > Inf)
                daspect([1 1 1]);
                view(90,0);
                axis off;
                return
            end
        end
        count = count + 1;
    end
end


% global camera
daspect([1 1 1]);
view(90,0);
axis off;
cb = colorbar;
colormap(cmap);
caxis([mag_min mag_max]);
n_cticks = 4;
cticks = linspace(mag_min,mag_max,n_cticks);
cticks_lab = cell(1,n_cticks);
for i = 1:n_cticks
    cticks_lab{i} = sprintf('%.2f',cticks(i));
end
set(cb,'Location','west');
set(cb,'TickLength',0);
set(cb,'Ticks',cticks);
set(cb,'TickLabels',cticks_lab);


% Expand space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% save
print(h,'./figures/fig_cluster2_sworld_dk','-depsc','-r400');
print(h,'./figures/fig_cluster2_sworld_dk','-dpng','-r400');
close(h);

% Small world index calculations
trig_swi = false;

if (trig_swi)

    %Ahs = Ca3.A2;
    n_Ahs = length(Ahs);


    fprintf('[*] Computing small world index for humans..\n')
    [ws,sigma,omega,I] = swi(Ahs,n_MC);
    fprintf('\tn = %i\n',n_Ahs);
    fprintf('\t%s\n',ws)
    fprintf('\tsigma = %.6f\n',sigma);
    fprintf('\tomega = %.6f\n',omega);
    fprintf('\tSWI = %.6f\n',I);

    % Load individual adjacency matrices
    Ca = load('./cache/xsub_out_all_1.mat');
    iM = 1;
    for i = 1:length(Ca.Subjects)
        sid = Ca.Subjects{i};
        Cas = load(sprintf('./cache/xsub_out_%s_%i.mat',sid,iM));
    end

    % fprintf('[*] Calculating Watts-Strogatz metrics for H sapiens..\n')
    % [Lhs,Chs] = watts_strogatz(Ahs);
    % fprintf('\tn = %i\n',length(Ahs));
    % fprintf('\tL = %.6f\n\tC = %.6f\n',Lhs,Chs);
    % 
    % % Compute random graph metrics
    % fprintf('[*] Computing random graph metrics..\n')
    % n_MC = 12;
    % [Lhs_rand,Chs_rand] = watts_strogatz_random(Ahs, n_MC);
    % 
    % p_Lhs = 1 - sum((Lhs > sort(Lhs_rand)) | (Lhs < sort(Lhs_rand)))/length(Lhs_rand);
    % p_Chs = 1 - sum((Chs > sort(Chs_rand)) | (Chs < sort(Chs_rand)))/length(Chs_rand);
    % fprintf('\tL_rand = %.6f +- %.6f\n\tC_rand = %.6f +- %.6f\n',...
    %     mean(Lhs_rand),std(Lhs_rand),mean(Chs_rand),std(Chs_rand));
    % fprintf('[*] Permutation test with n: %i, 2-sided p_L = %.3d, p_C = %.3d\n',...
    %     n_MC,p_Lhs,p_Chs);
    % 
    % % small world metric
    % C = Chs;
    % Cr = mean(Chs_rand);
    % L = Lhs;
    % Lr = mean(Lhs_rand);
    % sigma = (C/Cr)/(L/Lr); % if sigma > 1 then small world
    % fprintf('\tsigma: %.6f\n',sigma);






    % Load c elegans adjacency matrix
    % L.R. Varshney, B.L. Chen, E. Paniagua, D.H. Hall and D.B. Chklovskii (2011)
    fn_ce = './connectome/NeuronConnect.csv';
    f_ce = fopen(fn_ce,'r');
    Dce = textscan(f_ce,'%s','Delimiter','\n');
    Dce = Dce{1};
    fclose(f_ce);
    Ece = {};
    j = 1;
    for i = 2:length(Dce)
        line = strsplit(Dce{i},',');
        if (~strcmp(line{3},'NMJ'))
            Ece{j,1} = line{1};
            Ece{j,2} = line{2};
            Ece{j,3} = line{3};
            Ece{j,4} = str2double(line{4});
            j = j + 1;
        end
    end
    Ace_labels = unique({Ece{:,1}; Ece{:,2}});
    n_Ace = length(Ace_labels);
    Ace = zeros(n_Ace,n_Ace);
    % build adjacency from edges
    [n_Ece,~] = size(Ece);
    for i = 1:n_Ece
        neur1 = Ece{i,1};
        neur2 = Ece{i,2};
        idx_neur1 = find(strcmp(Ace_labels,neur1),1);
        idx_neur2 = find(strcmp(Ace_labels,neur2),1);
        Ace(idx_neur1,idx_neur2) = Ece{i,4};
        Ace(idx_neur2,idx_neur1) = Ece{i,4};
    end

    fprintf('[*] Computing small world index for C elegans..\n')
    [ws,sigma,omega,I] = swi(Ace,n_MC);
    fprintf('\tn = %i\n',n_Ace);
    fprintf('\t%s\n',ws)
    fprintf('\tsigma = %.6f\n',sigma);
    fprintf('\tomega = %.6f\n',omega);
    fprintf('\tSWI = %.6f\n',I);


    % fprintf('\n[*] Calculating Watts-Strogatz metrics for C elegans..\n')
    % [Lce,Cce] = watts_strogatz(Ace);
    % fprintf('\tn = %i\n',length(Ace));
    % fprintf('\tL = %.6f\n\tC = %.6f\n',Lce,Cce);
    % 
    % % Compute random graph metrics
    % fprintf('[*] Computing random graph metrics..\n')
    % %n_MC = 12;
    % [Lce_rand,Cce_rand] = watts_strogatz_random(Ahs, n_MC);
    % 
    % p_Lce = 1 - sum((Lce > sort(Lce_rand)) | (Lce < sort(Lce_rand)))/length(Lce_rand);
    % p_Cce = 1 - sum((Cce > sort(Cce_rand)) | (Cce < sort(Cce_rand)))/length(Cce_rand);
    % fprintf('\tL_rand = %.6f +- %.6f\n\tC_rand = %.6f +- %.6f\n',...
    %     mean(Lce_rand),std(Lce_rand),mean(Cce_rand),std(Cce_rand));
    % fprintf('[*] Permutation test with n: %i, 2-sided p_L = %.3d, p_C = %.3d\n',...
    %     n_MC,p_Lce,p_Cce);
    % 
    % % small world metric
    % C = Cce;
    % Cr = mean(Cce_rand);
    % L = Lce;
    % Lr = mean(Lce_rand);
    % sigma = (C/Cr)/(L/Lr); % if sigma > 1 then small world
    % fprintf('\tsigma: %.6f\n',sigma);





    % Macaque
    D_anat = load('../coreg/AdjMKFV.mat');
    D = load('../coreg/AdjSu.mat');

    % convert to symmetric
    AdjMKneurons = D_anat.AdjMKneurons;
    AdjMKneurons = AdjMKneurons + AdjMKneurons';
    AdjMKneuronst = AdjMKneurons;
    AdjMKneuronst(AdjMKneuronst == 0) = NaN;
    AdjMK = log10(AdjMKneuronst);
    AdjMK_bin = ~isnan(AdjMK);
    AdjMK_bin = double(AdjMK_bin);

    Amaca = AdjMK;
    Amacf = D.AdjMag;
    isn = all(isnan(Amacf));
    Amacf(isn,:) = [];
    Amacf(:,isn) = [];
    Amacf(isnan(Amacf)) = 0;

    fprintf('[*] Computing small world index for Markov-Kennedy..\n')
    [ws,sigma,omega,I] = swi(Amaca,n_MC);
    fprintf('\tn = %i\n',length(Amaca));
    fprintf('\t%s\n',ws);
    fprintf('\tsigma = %.6f\n',sigma);
    fprintf('\tomega = %.6f\n',omega);
    fprintf('\tSWI = %.6f\n',I);

    fprintf('[*] Computing small world index for macaque functional interactions..\n')
    [ws,sigma,omega,I] = swi(Amacf,n_MC);
    fprintf('\tn = %i\n',length(Amacf));
    fprintf('\t%s\n',ws);
    fprintf('\tsigma = %.6f\n',sigma);
    fprintf('\tomega = %.6f\n',omega);
    fprintf('\tSWI = %.6f\n',I);

    % 
    % %Compute
    % fprintf('\n[*] Calculating Watts-Strogatz metrics for Markov-Kennedy..\n')
    % [Lmaca,Cmaca] = watts_strogatz(Amaca);
    % fprintf('\tn = %i\n',length(Amaca));
    % fprintf('\tL = %.6f\n\tC = %.6f\n',Lmaca,Cmaca);
    % 
    % % Compute random graph metrics
    % fprintf('[*] Computing random graph metrics..\n')
    % %n_MC = 12;
    % [Lmaca_rand,Cmaca_rand] = watts_strogatz_random(Amaca, n_MC);
    % 
    % p_Lmaca = 1 - sum((Lmaca > sort(Lmaca_rand)) | (Lmaca < sort(Lmaca_rand)))/length(Lmaca_rand);
    % p_Cmaca = 1 - sum((Cmaca > sort(Cmaca_rand)) | (Cmaca < sort(Cmaca_rand)))/length(Cmaca_rand);
    % fprintf('\tL_rand = %.6f +- %.6f\n\tC_rand = %.6f +- %.6f\n',...
    %     mean(Lmaca_rand),std(Lmaca_rand),mean(Cmaca_rand),std(Cmaca_rand));
    % fprintf('[*] Permutation test with n: %i, 2-sided p_L = %.3d, p_C = %.3d\n',...
    %     n_MC,p_Lmaca,p_Cmaca);
    % 
    % % small world metric
    % C = Cmaca;
    % Cr = mean(Cmaca_rand);
    % L = Lmaca;
    % Lr = mean(Lmaca_rand);
    % sigma = (C/Cr)/(L/Lr); % if sigma > 1 then small world
    % fprintf('\tsigma: %.6f\n',sigma);
    % 
    % 
    % 
    % %Compute
    % fprintf('\n[*] Calculating Watts-Strogatz metrics for mSu..\n')
    % [Lmacf,Cmacf] = watts_strogatz(Amacf);
    % fprintf('\tn = %i\n',length(Amacf));
    % fprintf('\tL = %.6f\n\tC = %.6f\n',Lmacf,Cmacf);
    % 
    % % Compute random graph metrics
    % fprintf('[*] Computing random graph metrics..\n')
    % %n_MC = 12;
    % [Lmacf_rand,Cmacf_rand] = watts_strogatz_random(Amacf, n_MC);
    % 
    % p_Lmacf = 1 - sum((Lmacf > sort(Lmacf_rand)) | (Lmacf < sort(Lmacf_rand)))/length(Lmacf_rand);
    % p_Cmacf = 1 - sum((Cmacf > sort(Cmacf_rand)) | (Cmacf < sort(Cmacf_rand)))/length(Cmacf_rand);
    % fprintf('\tL_rand = %.6f +- %.6f\n\tC_rand = %.6f +- %.6f\n',...
    %     mean(Lmacf_rand),std(Lmacf_rand),mean(Cmacf_rand),std(Cmacf_rand));
    % fprintf('[*] Permutation test with n: %i, 2-sided p_L = %.3d, p_C = %.3d\n',...
    %     n_MC,p_Lmacf,p_Cmacf);
    % 
    % % small world metric
    % C = Cmacf;
    % Cr = mean(Cmacf_rand);
    % L = Lmacf;
    % Lr = mean(Lmacf_rand);
    % sigma = (C/Cr)/(L/Lr); % if sigma > 1 then small world
    % fprintf('\tsigma: %.6f\n',sigma);
    % 

end