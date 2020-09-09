function [l_WS,c_WS] = watts_strogatz_lattice(A2)
% Calculates the Watts-Strogatz metrics for the theoretical lattice case
%   L - the typical path length
%   C - the clustering coefficient, local average
%
%   Usage:
%       [L,C] = watts_strogatz_lattice(AdjacencyMatrix)
%

    n_nodes = length(A2);
    n_edges = sum(sum(triu(A2)>0));
    %n_eposs = nchoosek(n_nodes,2);
    A = zeros(n_nodes,n_nodes);
    k = 1; % edge count
    n = 1; % node position count
    nn = 1; % nearest neighbor count
    while (true)
        % break when number of edges run out
        if (k > n_edges)
            break;
        end
        
        % set lattice edges
        n1 = n;
        n2 = n + nn;
        % wrap around
        if (n2 > n_nodes)
            n2 = n2 - n_nodes;
        end
        A(n1,n2) = 1;
        A(n2,n1) = 1;
        
        % increment nearest neighbor count on full circle
        if (n == n_nodes)
            n = 0;
            nn = nn + 1;
        end
        
        k = k + 1;
        n = n + 1;
    end
    
    [l_WS,c_WS] = watts_strogatz(A);
% 
%     % clustering coefficient calculation
%     ci = clusteringcoef(A2 > 0);
%     c_WS = nanmean(ci);
%     %fprintf('[!] Watts and Strogatz C_actual: %.6f\n',c_WS);
% 
%     [n_A2,~] = size(A2);
%     G = graph(A2 > 0);
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
% 
%     l_WS = nanmean(pathl(:));
%     %fprintf('[!] Watts and Strogatz L_actual: %.6f\n',l_WS);
end

