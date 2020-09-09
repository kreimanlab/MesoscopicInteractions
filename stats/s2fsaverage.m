clear;
close all;

addpath(genpath('/home/jerry/Documents/MATLAB/iELvis/iELVis_MAIN'));
addpath(genpath('/home/jerry/Documents/MATLAB/iELvis/iELVis_MATLAB_ADMIN'));

%subjects_dir = '/media/jerry/internal/data/coreg';
subjects_dir = '../data/coregistration';
setenv('SUBJECTS_DIR',subjects_dir);

groupAvgCoords = [];
groupLabels = [];
groupIsLeft = [];
cfg = [];
cfg.plotEm = 0;
%sid_const = 'sub6';
%subs = {'sub6'}; %{'sub1','sub3'};
subs = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
        'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
        'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',... % 'sub23',
        'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',... % 'sub30',
        'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
        'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

label_header = '#!ascii label  , from subject ';
label_header2 = '  vox2ras=TkReg';
Etable = {};
for a = 1:length(subs)
    sid_const = subs{a};
    fprintf('Working on Participant %s\n',subs{a});
    [avgCoords, elecNames, isLeft]=sub2AvgBrain(subs{a},cfg);
    groupAvgCoords = avgCoords; %[groupAvgCoords; avgCoords];
    groupLabels = [groupLabels; elecNames];
    groupIsLeft = [groupIsLeft; isLeft];
    %close all;
    
    if (isLeft(1))
        hemi = 'l';
    else
        hemi = 'r';
    end
    groupAvgCoords(:,1) = groupAvgCoords(:,1) .* (-1 * 2*(0.5-isLeft));
    avgCoords(:,1) = avgCoords(:,1) .* (-1 * 2*(0.5-isLeft));
    avgCoords(:,1) = avgCoords(:,1) * (-1);
    
    
    %groupAvgCoords(:,1)
    
    % Read electrode index
    [v_pia, fa_pia] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,'fsaverage_sym','r','pial'));
    vf_pia = triangulation(fa_pia+(1-min(fa_pia(:))),v_pia(:,1),v_pia(:,2),v_pia(:,3));
    %l = read_label(sprintf('%s/%s',subjects_dir,sid_const),'all_surf_ielvis');
    %if (isempty(l))
    l = read_label(sprintf('%s',subs{a}),'all_surf_ielvis');
    l(:,1) = l(:,1)+1;
    [~,sIdx] = sort(l(:,end));
    l = l(sIdx,:);
    %    fprintf('[!] Error fixed.\n')
    %end
    
    % ----- sort --------------------------------------------
    groupAvgCoords = groupAvgCoords(sIdx,:);
    avgCoords = avgCoords(sIdx,:);
    elecNames = elecNames(sIdx,:);
    isLeft = isLeft(sIdx,:);
    
    
    % Check number of channels match
    n_chan = length(l(:,end));
    n_chan_ielvis = length(avgCoords(:,1));
    if (n_chan ~= n_chan_ielvis)
        fprintf(2,'[E] number of channels mismatch.\n')
        return
    end
    
    % Output to fsaverage_sym label
    fn_out = sprintf('%s/fsaverage_sym/label/ielvis_%s.label',subjects_dir,subs{a});
    of = fopen(fn_out,'w');
    fprintf(of,'%s%s%s\n',label_header,subs{a},label_header2);
    fprintf('[*] Wrote label to: %s\n',fn_out);
    fprintf(of,'%i\n',n_chan);
    for i = 1:n_chan
        coord = avgCoords(i,:);
        idx = vf_pia.nearestNeighbor(coord);
        distance = norm(vf_pia.Points(idx,:) - coord);
        if (distance > 0)
            fprintf(2,'\t(*) remapped vertex onto mesh: %s, dist: %.2f mm\n',elecNames{i},distance);
        end
        fprintf(of,'%i %.3f %.3f %.3f %i\n',idx-1,coord(1),coord(2),coord(3),l(i,end));
    end
    fclose(of);
    
    
    
    
    
    % Plotting ----------------------------------------------------------
    col_pial = [1 1 1]*0.8;
    al_pial = 0.8; % 0.8 for other figures
    print_res = '-r1000';
    lcolors = rand(n_chan,3)*0.8 + [1 1 1]*0.2;
    
    h = figure('Position',[0 0 1080 1080],'visible','off','PaperUnits','inches','PaperPosition',[0 0 6 6]);
    ha = tight_subplot(2,2,[.01 .03],[.1 .01],[.01 .01]);

    emark_size = 24;
    emark_fontsize = 6;
    
    %subplot(2,2,1);
    for ii = 1:2
        axes(ha(ii));
        [v_pia, fa_pia] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,'fsaverage_sym','r','pial'));
        vf_pia = triangulation(fa_pia+(1-min(fa_pia(:))),v_pia(:,1),v_pia(:,2),v_pia(:,3));
        p = trisurf(fa_pia + (1-min(fa_pia(:))),v_pia(:,1),v_pia(:,2),v_pia(:,3),'EdgeColor','none','FaceColor',col_pial,'FaceAlpha',al_pial); hold on;
        hold all;

        for i = 1:length(groupAvgCoords(:,1))
            % Plot marker
            px = (-1)*groupAvgCoords(i,1);
            py = groupAvgCoords(i,2);
            pz = groupAvgCoords(i,3);
            plot3(px,py,pz,'black.','MarkerSize',emark_size,'Color',lcolors(i,:));
            
            % Plot text
            offset_const = 0;
            px = px + (3-2*ii)*offset_const;
            text(px,py,pz,sprintf('%i',i),'fontsize',emark_fontsize,'HorizontalAlignment','center')
        end
        brainlight;
        sign = (2*(1.5 - ii));
        view(90 * sign,0);
        axis tight;
        axis off;
        title('fsaverage\_sym','fontsize',10);
    end

    %copyobj(ha(1),ha(2));

    %%
    % Plot electrodes old style
    % h = figure;


    %subplot(2,2,3);
    for ii = 3:4
        axes(ha(ii));
        % Read patient
        [v_pia, fa_pia] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,'pial'));
        vf_pia = triangulation(fa_pia+(1-min(fa_pia(:))),v_pia(:,1),v_pia(:,2),v_pia(:,3));
        l = read_label(sprintf('%s/%s',subjects_dir,sid_const),'all_surf_ielvis');
        if (isempty(l))
            l = read_label(sprintf('%s',sid_const),'all_surf_ielvis');
            fprintf('[!] Error fixed.\n')
        end
        l(:,1) = l(:,1)+1;
        [~,sIdx] = sort(l(:,end));
        l = l(sIdx,:);
        hold all;
        p = trisurf(fa_pia + (1-min(fa_pia(:))),v_pia(:,1),v_pia(:,2),v_pia(:,3),'EdgeColor','none','FaceColor',col_pial,'FaceAlpha',al_pial); hold on;
        [n_l,~] = size(l);
        for i = 1:n_l
            pt = v_pia(l(i,1),:);
            plot3(pt(1),pt(2),pt(3),'.','Color',lcolors(i,:),'MarkerSize',emark_size); hold on;
            text(pt(1),pt(2),pt(3),sprintf('%i',i),'fontsize',emark_fontsize,'HorizontalAlignment','center')
        end
        brainlight;
        sign = (2*(1.5 - ii));
        if (strcmp(hemi,'r'))
            view(90 * sign,0);
        else
            view(-90 * sign,0);
        end
        axis tight;
        axis off;
        title(sprintf('Subject %i',a),'fontsize',10);
    end

%     fn_plot = sprintf('%s/fsaverage_sym/label/figures/%s_to_fsaverage_sym_rh',subjects_dir,sid_const);
%     print(h,fn_plot,'-depsc');
    %fn_plot = sprintf('%s/fsaverage_sym/label/s2fsaverage_sym_figures/subject_%i_fsaverage_sym_rh',subjects_dir,a);
    fn_plot = sprintf('%s/fsaverage_sym/label/s2fsaverage_sym_figures/Figure_W1-%i_%s_fsaverage_sym_rh',subjects_dir,a,sid_const);
    
    print(h,fn_plot,'-dpng',print_res);
    close(h);
    
    
    
    
    
    % ============= Plot bipolar labels ===================================
    Ca = load(sprintf('./cache/xsub_out_%s_1.mat',sid_const));
    h = figure('Position',[0 0 1080 1080],'visible','off','PaperUnits','inches','PaperPosition',[0 0 6 3]);
    ha = tight_subplot(1,2,[.01 .03],[.1 .01],[.01 .01]);
    lcolors = ones(n_chan,3)*0.7 + [1 1 1]*0.2;
    for ii = 1:2
        axes(ha(ii));
        % Read patient
        [v_pia, fa_pia] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,sid_const,hemi,'pial'));
        vf_pia = triangulation(fa_pia+(1-min(fa_pia(:))),v_pia(:,1),v_pia(:,2),v_pia(:,3));
        l = read_label(sprintf('%s/%s',subjects_dir,sid_const),'all_surf_ielvis');
        if (isempty(l))
            l = read_label(sprintf('%s',sid_const),'all_surf_ielvis');
            fprintf('[!] Error fixed.\n')
        end
        l(:,1) = l(:,1)+1;
        [~,sIdx] = sort(l(:,end));
        l = l(sIdx,:);
        hold all;
        p = trisurf(fa_pia + (1-min(fa_pia(:))),v_pia(:,1),v_pia(:,2),v_pia(:,3),'EdgeColor','none','FaceColor',col_pial,'FaceAlpha',al_pial); hold on;
        [n_l,~] = size(l);
        
        % Load bipolar info
        %Ca.ecog.bip
        for i = 1:Ca.ecog.n_bchan %n_l
            b1c1 = Ca.ecog.bip(i,1);
            l_idx = find(l(:,end) == b1c1);
            pt = v_pia(l(l_idx,1),:);
            plot3(pt(1),pt(2),pt(3),'.','Color',lcolors(i,:),'MarkerSize',emark_size); hold on;
            text(pt(1),pt(2),pt(3),sprintf('%i',i),'fontsize',emark_fontsize,'HorizontalAlignment','center')
        end
        brainlight;
        sign = (2*(1.5 - ii));
        if (strcmp(hemi,'r'))
            view(90 * sign,0);
        else
            view(-90 * sign,0);
        end
        axis tight;
        axis off;
        title(sprintf('Subject %i',a),'fontsize',10);
    end
    fn_plot = sprintf('%s/fsaverage_sym/label/s2fsaverage_sym_figures/Figure_W1-%iB_%s_fsaverage_sym_rh_n-%i',subjects_dir,a,sid_const,Ca.ecog.n_bchan);
    print(h,fn_plot,'-dpng',print_res);
    close(h);
    % ============= END Plot bipolar labels ===============================
    
    
    
    
    % Build electrode location table
    fn_ap = sprintf('%s/%s/label/all_parcellation.mat',subjects_dir,sid_const);
    ap = load(fn_ap);
    Ca = load(sprintf('./cache/xsub_out_%s_1.mat',sid_const));
    %Atls = [1 2];
    Atls = 1:20;
    length(Atls);
    et = {};
    for atI = Atls
        dict = ap.AtlROIs{atI}.LH.struct_names;
        lbl = ap.AtlLabels{atI};
        dict_n = zeros(size(dict));
        for i = 1:length(lbl)
            sIdx = strcmp(dict,lbl{i});
            dict_n(sIdx) = dict_n(sIdx) + 1;
        end
        et = [et, {dict, dict_n}];
    end
    
    Etable = [Etable; {et}];
%     if (a == 3)
%         return
%     end
end

%%
for ia = 1:length(Atls)
    fprintf('\n%s\n',ap.AtlNames{ia});
    dict_ns = []; %zeros(size(et{ia*2}));
    for i = 1:length(Etable)
        et = Etable{i};
        dict = et{ia*2-1};
        dict_n = et{ia*2};
        
        dict_ns = [dict_ns, dict_n];
    end
    
    dict_ns_all = sum(dict_ns,2);
    n_all = sum(dict_ns_all);
    d = dict_ns;
    d(d==0) = NaN;
    dict_ns_std = nanstd(d')';
    
    n_cov = 0;
    for j = 1:length(dict)
        if (dict_ns_all(j) > 0 )
            fprintf('%s,%i,%.1f,%.1f%%\n',dict{j},dict_ns_all(j),dict_ns_std(j),...
                100*dict_ns_all(j)/n_all);
            n_cov = n_cov + 1;
        end
    end
    fprintf('coverage %.2f%% total %i\n',100*n_cov/length(dict),length(dict));
    
end



% show current map





%% Plot electrodes ielvis style
% cfg=[];
% cfg.view='l';
% cfg.elecCoord=[groupAvgCoords groupIsLeft];
% cfg.elecNames=groupLabels;
% cfg.showLabels='n';
% cfg.title='PTs on Avg. Brain';
% cfgOut=plotPialSurf('fsaverage',cfg);



% h = figure;
% [v_pia, fa_pia] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,'fsaverage_sym','l','pial'));
% vf_pia = triangulation(fa_pia+(1-min(fa_pia(:))),v_pia(:,1),v_pia(:,2),v_pia(:,3));
% p = trisurf(fa_pia + (1-min(fa_pia(:))),v_pia(:,1),v_pia(:,2),v_pia(:,3),'EdgeColor','none','FaceColor',col_pial,'FaceAlpha',al_pial); hold on;
% hold all;
% 
% for i = 1:length(groupAvgCoords(:,1))
%     plot3(groupAvgCoords(i,1),groupAvgCoords(i,2),groupAvgCoords(i,3),'black.','MarkerSize',20,'Color',lcolors(i,:));
%     text(groupAvgCoords(i,1),groupAvgCoords(i,2),groupAvgCoords(i,3),sprintf('%i',i))
% end
% brainlight;
% view(-90,0);
% axis off;


return
%%

% Read average
[v_pia, fa_pia] = read_surf(sprintf('%s/%s/surf/%sh.%s',subjects_dir,'fsaverage_sym',hemi,'pial'));
vf_pia = triangulation(fa_pia+(1-min(fa_pia(:))),v_pia(:,1),v_pia(:,2),v_pia(:,3));
l = read_label(sprintf('%s/%s',subjects_dir,'fsaverage_sym'),sprintf('ielvis_%s',sid_const));
if (isempty(l))
    l = read_label(sprintf('%s','fsaverage_sym'),sprintf('ielvis_%s',sid_const));
    fprintf('[!] Error fixed.\n')
end
[~,sIdx] = sort(l(:,end));
l = l(sIdx,:);
subplot(1,2,1)
hold all;
p = trisurf(fa_pia + (1-min(fa_pia(:))),v_pia(:,1),v_pia(:,2),v_pia(:,3),'EdgeColor','none','FaceColor',col_pial,'FaceAlpha',al_pial); hold on;
brainlight;
view(90,0);
axis off;
%v_pia(l(1:5,1),:)

% cfg=[];
% cfg.view='l';
% cfg.elecCoord=[groupAvgCoords groupIsLeft];
% cfg.elecNames=groupLabels;
% cfg.showLabels='n';
% cfg.title='PTs on Avg. Brain';
% cfgOut=plotPialSurf('fsaverage',cfg);
