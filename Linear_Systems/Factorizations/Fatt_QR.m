function [Q,R] = Fatt_QR(A) 
% QR FACTORIZATION
% Input:
%        A    matrix to factorize
% Output:
%        Q    matrix product of all elementary reflectors
%        R    upper triangular matrix

% Apply QR factorization algorithm to matrix A. We iterate to obtain
% elementary reflectors of each submatrix of A and then multiply them
% to the previous A. At each iteration A will vary according to the
% elementary reflector it was multiplied by, and finally in the last
% step of the loop matrix A will be the upper triangular matrix R.

    [m, n]=size(A);
    Q = eye(m);
    
    for i=1:min(m-1,n)
        % vector containing first column of current submatrix
        x = A(i:m,i);
        
        % A is not uniquely determined, so it's possible to choose
        % the norm of x with positive or negative sign. The most
        % appropriate choice to calculate v is that x+norm_x has
        % the same sign so that the sum is not 0.
        if x(1) >= 0
            nor_x = norm(x,2);
        else
            nor_x = -norm(x,2);
        end
        
        e1 = zeros(length(x),1);
        e1(1) = 1;
        v = x + nor_x * e1;
        
        if norm(v) > 0
            % elementary reflector: Householder matrix (symmetric, orthogonal, idempotent)
            % build a matrix that transforms matrix A into upper triangular
            P_sub = eye(length(v)) - (2*v*v')/(v'*v);
            
            % Build the full P matrix of dimension m x m
            P = eye(m);
            P(i:m,i:m) = P_sub;
            
            % Modify matrix A for next iteration
            A = P*A;
            
            % Build matrix Q by multiplying previous elementary reflectors
            Q = Q*P;
        end
    end
    R = A;
end