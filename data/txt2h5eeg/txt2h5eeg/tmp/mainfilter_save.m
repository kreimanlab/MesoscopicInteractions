function [ ] = mainfilter_save( path, data )
%mainfilter_save saves data into path
%   this method exists so that mainfilter can be parallelized at the file
%   level

%data = d(:,2:end);
%data = d;
save(path, '-v6', 'data');

end

