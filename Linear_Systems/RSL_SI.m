function x = RSL_SA(U,b)
% SOLVE LINEAR SYSTEM WITH UPPER TRIANGULAR MATRIX
% BACKWARD SUBSTITUTION
% Input:
%        U    upper triangular matrix
%        b    right-hand side vector
% Output:
%        x    solution vector

    n = length(b);
    x = zeros(n,1);
    % solve for the last unknown
    x(n) = b(n)/U(n,n);
    % solve for unknowns i=n-1,...,1
    for i=n-1:-1:1
        % solve for the i-th unknown
        x(i) = (b(i) - U(i,i+1:n)*x(i+1:n)) / U(i,i);
    end
end