function [L,C,Ll,Cl,Lr,Cr] = swi_raw(A,n_MC)
%swi - Small world network index
%
%       [sigma,omega,SWI] = swi(A,n_MC)
%
%       inputs:
%           A     - Adjacency matrix
%           n_MC  - number of random network samples
%
%       outputs:
%           ws    - Watts-Strogatz 1998
%           sigma - Humphries-Gurney,2008 (sigma > 1 for small world network)
%           omega - Telesford-Laurienti,2011 (0 to 1: 1 is most small world), (originally -1 to 1: 0 is most small world, + is more random, - is lattice)
%           I     - Neal,2017 (0 to 1: 1 is most small world)

    % 49.9 % shortestpath - 97.4%
    [L,C] = watts_strogatz(A);
    % 25.6 %
    [Ll,Cl] = watts_strogatz_lattice(A);
    % 24.5 %
    [Lr,Cr] = watts_strogatz_random(A,n_MC);
    Lr = nanmean(Lr);
    Cr = nanmean(Cr);
    
%     fprintf('L: %.2f\n',L);
%     fprintf('C: %.2f\n',C);
%     fprintf('Ll: %.2f\n',Ll);
%     fprintf('Cl: %.2f\n',Cl);
%     fprintf('Lr: %.2f\n',Lr);
%     fprintf('Cr: %.2f\n',Cr);
    
    ws = sprintf('(%.6f >> %.6f and %.6f =~ %.6f)',C,Cr,L,Lr);
    sigma = (C/Cr)/(L/Lr);
    omega = (Lr/L) - (C/Cl);
    %omega = 1 - abs(omega);
    I = ((L - Ll)/(Lr - Ll)) * ((C - Cr)/(Cl - Cr));

end

