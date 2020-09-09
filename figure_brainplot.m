close all;
clear;

Subjects = {'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',...
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',...
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48',...
    };
fig_fmt = '-dpng';
system('mkdir figures');
system('mkdir figures/brainplot');
atl = '';
for i = 1:length(Subjects)
    sid = Subjects{i};
    hb = subaxis(6,8,i, 'SpacingHoriz', 0.01, 'SpacingVert', 0.01, 'Padding', 0, 'Margin', 0);
    %hb = subplot(6,8,i);
    h = brainplot_c1(sid,atl);
    copyobj(h,hb);
    %print(h,sprintf('figures/brainplot/%s-%s',sid,atl),fig_fmt);
    %close(h);
end
%tightfig;