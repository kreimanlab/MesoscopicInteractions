function [l_WS,c_WS] = watts_strogatz(A2)
% Calculates the Watts-Strogatz metrics for small-world networks:
%   L - the typical path length
%   C - the clustering coefficient, local average
%
%   Usage:
%       [L,C] = watts_strogatz(AdjacencyMatrix)
%
    Adj_proc = (A2 > 0);

    % clustering coefficient calculation
    ci = clusteringcoef(Adj_proc);
    c_WS = nanmean(ci);
    %fprintf('[!] Watts and Strogatz C_actual: %.6f\n',c_WS);

    %[n_A2,~] = size(A2);
    G = graph(Adj_proc);
    
    dG = triu(distances(G));
    l_WS = nanmean(dG((dG ~=0) & (~isinf(dG)) & (~isnan(dG))));
%     
%     
%     A2_pathl = nan(size(A2));
%     %A2_paths = cell(size(A2));
%     pathl = nan(nchoosek(n_A2,2),1);
%     ka2 = 1;
%     %E = nan(nchoosek(n_A2,2),1);
%     for ia2 = 1:(n_A2-1)
%         for ja2 = (ia2+1):n_A2
%             path = shortestpath(G,ia2,ja2);
%             spl = length(path);
%             A2_pathl(ia2,ja2) = spl;
%             A2_pathl(ja2,ia2) = spl;
%             %A2_paths{ia2,ja2} = path;
%             %A2_paths{ja2,ia2} = path;
%             pathl(ka2) = spl;
%             %E(ka2) = A2(ia2,ja2);
%             ka2 = ka2 + 1;
%         end
%     end
%     l_WS = nanmean(pathl(:));
    %fprintf('[!] Watts and Strogatz L_actual: %.6f\n',l_WS);
end

