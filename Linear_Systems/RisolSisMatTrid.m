function [x] = RisolSisMatTrid(M, b)
% SOLVE TRIDIAGONAL MATRIX SYSTEM
% Solves systems with tridiagonal matrix using LU factorization method
% Input:
%        M    tridiagonal matrix associated with the system
%        b    right-hand side vector
% Output:
%        x    solution vector of the system

    n = length(M);
    % extract main, lower, and upper diagonals from matrix M
    ad = diag(M,0);     % main diagonal
    bd = [0; diag(M,-1)]; % lower diagonal with leading zero
    cd = [diag(M,1); 0];  % upper diagonal with trailing zero
    
    % by definition we derive alpha, beta and c (coincides with cd)
    alfa = zeros(n,1);
    beta = zeros(n,1);
    
    alfa(1) = ad(1);
    for i=2:n
        beta(i) = bd(i)/alfa(i-1);
        alfa(i) = ad(i) - beta(i)*cd(i-1);
    end
    
    % remove extra elements
    beta(1) = [];
    cd(end) = [];
    
    % L: lower bidiagonal matrix (specific for tridiagonal)
    L = eye(n) + diag(beta,-1);
    % U: upper bidiagonal matrix (specific for tridiagonal)  
    U = diag(alfa,0) + diag(cd,1);
    
    % solve system L*z=b
    z = RSL_SI(L, b);
    % solve system U*x=z
    x = RSL_SA(U, z);
end