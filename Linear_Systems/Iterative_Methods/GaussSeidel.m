function [x, iter] = GaussSeidel(A, b, x0, tol, max_iter)
% Gauss-Seidel iterative method for solving linear systems Ax = b
% Input:
%        A        coefficient matrix
%        b        right-hand side vector
%        x0       initial guess
%        tol      tolerance for convergence
%        max_iter maximum number of iterations
% Output:
%        x        solution vector
%        iter     number of iterations performed

    if nargin < 5
        max_iter = 100;
    end
    if nargin < 4
        tol = 1e-6;
    end
    
    n = size(A, 1);
    D = diag(diag(A));
    L = -tril(A, -1);
    U = -triu(A, 1);
    
    % Iteration matrix
    T = (D - L) \ U;
    c = (D - L) \ b;
    
    % Check convergence condition
    rho = max(abs(eig(T)));
    if rho >= 1
        warning('Gauss-Seidel method may not converge (spectral radius = %.4f)', rho);
    end
    
    x = x0;
    for iter = 1:max_iter
        x_new = T * x + c;
        
        if norm(x_new - x, inf) < tol
            x = x_new;
            return;
        end
        
        x = x_new;
    end
    
    warning('Maximum number of iterations reached without convergence');
end

