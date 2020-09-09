function [Lhs_rand,Chs_rand] = watts_strogatz_random(A, n_MC)
    %n_MC = 12;
    n_nodes = length(A);
    n_edges = sum(sum(triu(A)>0));
    n_eposs = nchoosek(n_nodes,2);
    Lhs_rand = [];
    Chs_rand = [];
    
    for i = 1:n_MC
        Esamp = zeros(n_eposs,1);
        Esamp(1:n_edges) = 1;
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
        [Ltmp,Ctmp] = watts_strogatz(Arand);
        Lhs_rand = [Lhs_rand; Ltmp];
        Chs_rand = [Chs_rand; Ctmp];
    end

end

