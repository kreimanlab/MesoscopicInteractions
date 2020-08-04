function [ D2 ] = nanDist4clustering_neg_inf(XI, XJ)
%nanDist4clustering Custom "pdist" function for computing distance
%   Marginalized distance (ignoring nans) between two input vectors XI and
%   XJ, as opposed to imputed distances, the matlab default


    % Do not accept negative values. Used only for
    % mainClusterValidationNC_compute.m
    XI(XI < 0) = NaN;
    XJ(XJ < 0) = NaN;

    n = size(XI,2);
    sqdx = (XI-XJ).^2;
    nstar = sum(~isnan(sqdx),2); % Number of pairs that do not contain NaNs
    nstar(nstar == 0) = NaN; % To return NaN if all pairs include NaNs
    D2squared = nansum(sqdx,2).*n./nstar; % Correction for missing coordinates
    D2 = sqrt(D2squared);

    if (any(isnan(D2)))
        fprintf(2,'[!WARN!] nanDist2clustering_neg output NaN, setting dist to Inf.\n');
        D2(isnan(D2)) = Inf;
    end
end

