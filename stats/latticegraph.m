function [A] = latticegraph(A2)
    n_nodes = length(A2);
    n_edges = sum(sum(triu(A2)>0));
    
    % Vectorize nonzero values
    A2_vec = zeros(n_edges,1);
    nvec = 1;
    for i = 1:n_nodes
        for j = i:n_nodes
            a2v = A2(i,j);
            if (a2v > 0)
                A2_vec(nvec) = a2v;
                nvec = nvec + 1;
            end
        end
    end
    % Shuffle weights
    A2_vec = A2_vec(randperm(length(A2_vec)));
    
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
        A(n1,n2) = A2_vec(k); %1;
        A(n2,n1) = A2_vec(k); %1;
        
        % increment nearest neighbor count on full circle
        if (n == n_nodes)
            n = 0;
            nn = nn + 1;
        end
        
        k = k + 1;
        n = n + 1;
    end
end

