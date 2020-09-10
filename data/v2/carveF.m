function [ outL ] = carveF( indir, suffix )
    d = dir(indir);
    disdir = [d.isdir];
    dname = {d.name};
    L = dname(~ disdir);
    outL = {};
    j = 1;
    for i = 1:length(L)
        if (endsWith(L{i},suffix))
            outL{j} = sprintf('%s/%s',indir,L{i});
            j = j + 1;
        end
    end
end