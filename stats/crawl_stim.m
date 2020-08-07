clear;
close all;

crawl_path = '/mnt/symphony_1/RawData';

fprintf('[*] generating subpaths for: %s\n',crawl_path);
paths = genpath(crawl_path);
paths = strsplit(paths,':');

for i = 1:length(paths)
    path = paths{i};
    ifns = carveF(path,'.ent.txt');
    for j = 1:length(ifns)
        ifn = ifns{j};
        fprintf('%s\n',ifn);
    end
end
