function x = RSL_SI(L,b)
% SOLVE LINEAR SYSTEM WITH LOWER TRIANGULAR MATRIX
% FORWARD SUBSTITUTION
% Input:
%        L    lower triangular matrix
%        b    right-hand side vector
% Output:
%        x    solution vector

    n = length(b);
    x = zeros(n,1);
    % solve for the first unknown
    x(1) = b(1)/L(1,1);
    % solve for unknowns i=2,...,n
    for i=2:n
        % solve for the i-th unknown
        x(i) = (b(i) - L(i,1:i-1)*x(1:i-1)) / L(i,i);
    end
end