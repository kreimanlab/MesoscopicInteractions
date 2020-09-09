function [ cc ] = corrcmap( n )

%addpath('lchconversions')

% hsl
%c1 = [0.30 1 0.5];
%c2 = [0.04 1 0.5];

%c1 = [77 -110 110];
%c2 = [93 100 110];

%c1 = [80 78 203];
%c2 = [80 95 58];
% 
%c2 = [10 28 -90];
%c1 = [100 0 90];




% ----------------------------------------------------------------------------
%cc = flipud(plasma(n));
%cc = plasma(n);

%cc = viridis(n);
%cc = inferno(n);
cc = getPyPlot_cMap('nipy_spectral', n, false, 'python3'); %'nipy_spectral'
%cc = getPyPlot_cMap('gist_ncar', n, false, 'python3'); %'nipy_spectral'
%cc = getPyPlot_cMap('rainbow', n, false, 'python3'); %'nipy_spectral'
%cc = getPyPlot_cMap('gist_rainbow', n, false, 'python3'); %'nipy_spectral'
%cc = getPyPlot_cMap('gnuplot2', n, false, 'python3'); %'nipy_spectral'
%cc = getPyPlot_cMap('cubehelix', n, false, 'python3'); %'nipy_spectral'
%cc = cc * 0.9 + [1 1 1] * 0.1;
%cc = cc(4:(end-3),:);
% gist_ncar

% iseven = (mod(n,2) == 0);
% if (iseven)
%     n1 = round(n/2);
%     n2 = round(n/2);
% else
%     n1 = floor(n/2);
%     n2 = ceil(n/2);
% end
% cc = [magma(n1); ...
%     flipud(viridis(n2))];

% c2 = [10        50      -45     ];
% c1 = [98       40       (45)      ];
% cc = zeros(n,3);
% count = 1;
% for i = linspace(0,1,n)
%     cT = i*c2 + (1-i)*c1;
%     
%     im_lch = cT;
%     C_lch2lab= makecform('lch2lab');
%     im_lab=applycform(im_lch,C_lch2lab);
%     im_rgb = lab2rgb(im_lab,'OutputType','uint8');
%     %C_lab2srgb= makecform('lab2srgb');
%     %im_rgb=applycform(im_lab,C_lab2srgb);
%     
%     cc(count,:) = im_rgb; %lch2rgb(cT);
%     count = count + 1;
% end
% cc = cc/255;
% ----------------------------------------------------------------------------



% 
% c3 = [10 50 -45];
% c2 = [55 40 0];
% c1 = [100 0 90];
% 
% if (mod(n,2) == 1)
%     n = n + 1;
% end
% 
% cc = zeros(n,3);
% count = 1;
% for i = linspace(0,1,n/2)
%     cT = i*c2 + (1-i)*c1;
%     
%     im_lch = cT;
%     C_lch2lab= makecform('lch2lab');
%     im_lab=applycform(im_lch,C_lch2lab);
%     im_rgb = lab2rgb(im_lab,'OutputType','uint8');
%     %C_lab2srgb= makecform('lab2srgb');
%     %im_rgb=applycform(im_lab,C_lab2srgb);
%     
%     cc(count,:) = im_rgb; %lch2rgb(cT);
%     count = count + 1;
% end
% for i = linspace(0,1,n/2)
%     cT = i*c3 + (1-i)*c2;
%     
%     im_lch = cT;
%     C_lch2lab= makecform('lch2lab');
%     im_lab=applycform(im_lch,C_lch2lab);
%     im_rgb = lab2rgb(im_lab,'OutputType','uint8');
%     %C_lab2srgb= makecform('lab2srgb');
%     %im_rgb=applycform(im_lab,C_lab2srgb);
%     
%     cc(count,:) = im_rgb; %lch2rgb(cT);
%     count = count + 1;
% end
% 
% cc = cc/255;
%wL = [0.2126;0.7152;0.0722];


% 
% c5 = [10        50      -45     ];
% c4 = [32.5      35.5    -11.25  ];
% c3 = [55        26      22.5    ];
% c2 = [77.5      13.5    56.25   ];
% c1 = [100       0       90      ];
% 
% if (mod(n,2) == 1)
%     n = n + 1;
% end
% 
% cc = zeros(n,3);
% count = 1;
% for i = linspace(0,1,n/4)
%     cT = i*c2 + (1-i)*c1;
%     im_lch = cT;
%     C_lch2lab= makecform('lch2lab');
%     im_lab=applycform(im_lch,C_lch2lab);
%     im_rgb = lab2rgb(im_lab,'OutputType','uint8');
%     cc(count,:) = im_rgb; %lch2rgb(cT);
%     count = count + 1;
% end
% for i = linspace(0,1,n/4)
%     cT = i*c3 + (1-i)*c2;
%     im_lch = cT;
%     C_lch2lab= makecform('lch2lab');
%     im_lab=applycform(im_lch,C_lch2lab);
%     im_rgb = lab2rgb(im_lab,'OutputType','uint8');
%     cc(count,:) = im_rgb; %lch2rgb(cT);
%     count = count + 1;
% end
% for i = linspace(0,1,n/4)
%     cT = i*c4 + (1-i)*c3;
%     im_lch = cT;
%     C_lch2lab= makecform('lch2lab');
%     im_lab=applycform(im_lch,C_lch2lab);
%     im_rgb = lab2rgb(im_lab,'OutputType','uint8');
%     cc(count,:) = im_rgb; %lch2rgb(cT);
%     count = count + 1;
% end
% for i = linspace(0,1,n/4)
%     cT = i*c5 + (1-i)*c4;
%     im_lch = cT;
%     C_lch2lab= makecform('lch2lab');
%     im_lab=applycform(im_lch,C_lch2lab);
%     im_rgb = lab2rgb(im_lab,'OutputType','uint8');
%     cc(count,:) = im_rgb; %lch2rgb(cT);
%     count = count + 1;
% end
% cc = cc/255;

end

