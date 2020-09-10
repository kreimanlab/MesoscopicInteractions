function [ X ] = features_v2( V, Fs )
% V - 1 by n vector

X(1) = max(V) - min(V);
X(2) = max(abs(diff(V)));

% [pxx,f] = pwelch(V,Fs,[],[],Fs);
% Pxx = log10(pxx);
%X(3) = mean(Pxx((f > 19) & (f < 21))) / max(Pxx((f > 0) & (f < 15)));
%X(3) = mean(Pxx((f > 19) & (f < 21))) / mean(Pxx((f > 5) & (f < 15)));
%X(3) = max(Pxx((f > 18) & (f < 22)));
%X(3) = max(Pxx((f > 19) & (f < 21))) / min([Pxx((f > 18) & (f <= 19)),Pxx((f >= 21) & (f < 22))]);
%X(3) = max(Pxx((f > 19) & (f < 21))) / mean(Pxx((f >= 21) & (f <= 23)));
%X(3) = mean(Pxx((f >= 19) & (f <= 21))) / mean(Pxx((f > 0) & (f <= 15)));

%X(3) = mean(Pxx((f >= 19.5) & (f <= 20.5))) / mean(Pxx((f >= 0) & (f <= 5)));
%X(3) = mean(Pxx((f >= 19) & (f <= 21))) - mean(Pxx((f >= 25) & (f <= 35)));

% Final
% X(3) = mean(Pxx((f >= 18) & (f <= 22))) - mean(Pxx((f >= 5) & (f <= 14)));
% X(4) = mean(Pxx((f >= 57.5) & (f <= 62.5))) - mean(Pxx((f >= 0) & (f <= 50)));

%DEBUG
% fprintf('%.4f\n',X(3));
% if (X(3) > 1)
%     h = figure;
%     plot(f,Pxx,'black');
%     axis tight;
%     a = input('s');
% end

end

