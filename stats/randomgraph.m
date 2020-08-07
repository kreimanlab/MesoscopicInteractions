function [Arand] = randomgraph(A)
    n_nodes = length(A);
    n_edges = sum(sum(triu(A)>0));
    n_eposs = nchoosek(n_nodes,2);
    
    % Vectorize nonzero values
    A2_vec = zeros(n_edges,1);
    nvec = 1;
    for i = 1:(n_nodes)
        for j = (i):n_nodes
            a2v = A(i,j);
            if (a2v > 0)
                A2_vec(nvec) = a2v;
                nvec = nvec + 1;
            end
        end
    end
    % Shuffle weights
    A2_vec = A2_vec(randperm(length(A2_vec)));
    
    Esamp = zeros(n_eposs,1);
    Esamp(1:n_edges) = A2_vec; %1;
    E = Esamp(randperm(n_eposs));
    Arand = zeros(n_nodes,n_nodes);
    kk = 1;
    for ii = 1:(n_nodes-1)
        for jj = (ii+1):n_nodes
            e = E(kk);
            Arand(ii,jj) = e;
            Arand(jj,ii) = e;
            kk = kk + 1;
        end
    end

end

