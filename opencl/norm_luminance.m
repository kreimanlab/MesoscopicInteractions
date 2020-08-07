function [ cc2 ] = norm_luminance( cc )
%norm_luminance normalizes luminance of a colormap

[n,nColor] = size(cc);
if (nColor ~= 3)
    fprintf(2,'[!] Input must be colormap.\n')
end

% Luminance weights
wL = [0.2126;0.7152;0.0722];

% initialize
N_sim = 200000;
step = 0.01;
%x = ones(1,nColor);
x = ones(n,1);
%x = ones(n,nColor);
% objective function to min
val = -Inf;
for i = 1:N_sim
    xN = x + step/2 - step*rand(size(x));
    ccN = (xN + cc);
    ccN = ccN - min(ccN(:));
    %ccN = abs(ccN);
    ccN = ccN / max(ccN(:));
    valN = 0/(0.5 - mean(ccN*wL)) + 0.1/(var(ccN*wL)) + 1/sum(ccN*wL - ones(n,1)) + sum(abs(xN(:)))/(1*n);
    if (valN > val)
        x = xN;
        val = valN;
        fprintf('[*] objective fn val: %.8f\n',valN)
    end
end

cc2 = cc + x;
cc2 = cc2 - min(cc2(:));
%cc2 = abs(cc2);
cc2 = cc2 / max(cc2(:));
end

