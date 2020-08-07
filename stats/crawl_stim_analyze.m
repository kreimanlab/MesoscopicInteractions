close all;
clear;

Subjects = { ...
    'm00001','m00003','m00005','m00006','m00019','m00021','m00022','m00023',...
    'm00024','m00025','m00026','m00027','m00028','m00030','m00032','m00033',...
    'm00035','m00037','m00038','m00039','m00043','m00044','m00045','m00047',... % 'm00045',
    'm00048','m00049','m00052','m00053','m00055','m00056','m00058','m00059',... % 'm00056',
    'm00060','m00061','m00068','m00071','m00073','m00075','m00079','m00083',...
    'm00084','m00095','m00096','m00097','m00100','m00107','m00122','m00124'};

keyword_speech = {'speech','speak','answer','talk','word','language','say','said','voice','convers','dialogue','discuss','express','vocal','diction','utter'};

f = fopen('crawl_stim.log','r');
D = textscan(f,'%s','Delimiter','\n');
D = D{1};
fclose(f);

of = fopen('crawl_stim.speech','w');

cc = 1;
for i = 1:length(D)
    ifname = D{i};
    cond_in = contains(ifname,Subjects);
    if (cond_in)
        fprintf('[%03d] %s\n',cc,ifname);
        f = fopen(ifname,'r');
        D1 = textscan(f,'%s','Delimiter','\n');
        D1 = D1{1};
        for j = 1:length(D1)
            d1 = D1{j};
            if (contains(lower(d1),keyword_speech))
                fprintf(of,'%s\t%s\n',ifname,d1);
            end
        end
        fclose(f);
        cc = cc + 1;
    end
end

fclose(of);
