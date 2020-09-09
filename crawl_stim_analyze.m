close all;
clear;

Subjects = { ...
    'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8',...
    'sub9','sub10','sub11','sub12','sub13','sub14','sub15','sub16',...
    'sub17','sub18','sub19','sub20','sub21','sub22','sub23','sub24',... % 'sub23',
    'sub25','sub26','sub27','sub28','sub29','sub30','sub31','sub32',... % 'sub30',
    'sub33','sub34','sub35','sub36','sub37','sub38','sub39','sub40',...
    'sub41','sub42','sub43','sub44','sub45','sub46','sub47','sub48'};

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
