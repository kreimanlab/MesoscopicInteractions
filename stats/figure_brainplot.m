close all;
clear;

Subjects = {'m00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',...
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',...
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124',...
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