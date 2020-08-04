function [ D2 ] = nanDist4clusteringThresh(XI, XJ)
%nanDist4clustering Custom "pdist" function for computing distance
%   Marginalized distance (ignoring nans) between two input vectors XI and
%   XJ, as opposed to imputed distances, the matlab default

    n = size(XI,2);
    sqdx = (XI-XJ).^2;
    nstar = sum(~isnan(sqdx),2); % Number of pairs that do not contain NaNs
    nstar(nstar == 0) = NaN; % To return NaN if all pairs include NaNs
    D2squared = nansum(sqdx,2).*n./nstar; % Correction for missing coordinates
    D2 = sqrt(D2squared);

    %D2 = D2 * (10 * n/max(nstar));
%     if (max(nstar) <= (n/2))
%         D2 = D2 * 5;
%     end
    
end

