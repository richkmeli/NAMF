function [L, U] = Fatt_LU_NOPIV(A)
% LU FACTORIZATION WITHOUT PIVOTING
% Input:
%        A    matrix to factorize
% Output:
%        L    lower triangular matrix
%        U    upper triangular matrix

    % matrix dimensions n rows, m columns
    [n,m] = size(A);
    % iterations over n rows
    L=zeros(n,m);
    U=A;
    
    for i = 1:n-1
        % SCALING
        for j=i+1:n
            if U(i,i) ~=0
                % Build matrix of Gauss coefficients values
                tmp = U(j,i)/U(i,i);
                % build Gauss reduction matrix
                U(j,:) = U(j,:) - tmp*U(i,:);
                L(j,i) = tmp;
            end
        end
    end
    L = L + eye(n);
end