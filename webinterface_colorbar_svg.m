close all;
clear;

Generate_colorbar();

function Generate_colorbar(filename)
if nargin == 0
    filename = 'colormap.svg';
end

T = {'inferno'};
%   T={ ...
%     'jet'
%     'hsv'
%     'hot'
%     'cool'
%     'spring'
%     'summer'
%     'autumn'	
%     'winter'
%     'gray'
%     'bone'
%     'copper'
%     'pink'
%     'lines'
%     'colorcube'
%     'flag'
%     'prism'
% 	};    

fid = fopen(filename,'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n');
fprintf(fid,'<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n');
fprintf(fid,'<svg\n');
fprintf(fid,'\txmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n');
fprintf(fid,'\txmlns:svg="http://www.w3.org/2000/svg"\n');
fprintf(fid,'\txmlns="http://www.w3.org/2000/svg"\n');
fprintf(fid,'\twidth="512px"\n');
fprintf(fid,'\theight="%dpx">\n',12*numel(T)); % 32
fprintf(fid,'\tpreserveAspectRatio="none"\n');
fprintf(fid,'\t\n');
fprintf(fid,'\t\n');

for j=1:numel(T)
    cmap = colormap(T{j});
    [n_cmap,~] = size(cmap);
    if (n_cmap > 64)
        cmap = downsample(cmap,round(n_cmap/64));
        [n_cmap,~] = size(cmap);
    end
    fprintf(fid,'\t<g id="g_%s">\n',T{j});
    for i=1:size(cmap,1)
        colorname = sprintf('%02x%02x%02x',...
            fix(cmap(i,1)*255),fix(cmap(i,2)*255),fix(cmap(i,3)*255));
        fprintf(fid,'\t\t<rect\n');
        fprintf(fid,'\t\t\tid="rect%d"\n',1000+n_cmap*(j-1)+i);
        fprintf(fid,'\t\t\tx="%d"\n',8*(i-1));
        fprintf(fid,'\t\t\ty="%d"\n',36*(j-1));
        fprintf(fid,'\t\t\theight="32"\n');
        fprintf(fid,'\t\t\twidth="9"\n'); % 9
        fprintf(fid,'\t\t\tstyle="fill:#%s" />\n',colorname);   
    end
    fprintf(fid,'\t\t</g>\n\n');
end
fprintf(fid,'</svg>\n');
fclose(fid);
end