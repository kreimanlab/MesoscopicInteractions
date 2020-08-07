close all;
clear;

dir_h5 = '/media/jerry/KLAB101/h5_notch20';
dir_art = '/media/jerry/KLAB101/h5_notch20/art_nosz';
dir_video = '/media/jerry/internal/data/videos';%'/media/klab/KLAB101/videos';
dir_stamp = '/media/jerry/internal/data/stamps';%'/media/klab/KLAB101/stamps';


% Output path
dir_out = './sleeptrain';
d = dir(dir_out);
dname = {d.name};
disdir = [d.isdir];
dirs = dname(disdir);
dirs = dirs(3:end);
for i = 1:length(dirs)
    fnames = carveF(sprintf('%s/%s',dir_out,dirs{i}),'.jpg');
    Annot = {};
    aIdx = 1;
    outf = fopen(sprintf('./sleeptrain/%s.csv',dirs{i}),'w');
    for j = 1:length(fnames)
        fname = fnames{j};
        fl = strsplit(fname,'.jpg');
        fl = fl{1};
        fl = strsplit(fl,'/');
        fl = fl{end};
        fl = strsplit(fl,'-');
        if (length(fl) == 4)
            fl1 = strsplit(fl{1},'_');
            fl3 = strsplit(fl{3},'_');
            fl4 = strsplit(fl{4},'_');
            sid = fl1{1};
            fname = fl1{2};
            samp = str2double(fl3{1});
            label = str2double(fl4{1});
            if (~isnan(label))
                Annot{aIdx,1} = sid;
                Annot{aIdx,2} = fname;
                Annot{aIdx,3} = samp;
                Annot{aIdx,4} = label;
                
                fprintf(outf,'%s,%s,%i,%i\n',sid,fname,samp,label);
                
                aIdx = aIdx + 1;
            end
        end
    end
    fclose(outf);
end