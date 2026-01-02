function [y] = val_pol_carat_trid(A, lambda)
% CHARACTERISTIC POLYNOMIAL EVALUATION FOR TRIDIAGONAL MATRICES
% Calculates the values of the characteristic polynomial of matrix A at
% point lambda, p(lambda)=det(A-lambda*I).
% Input:
%        A       tridiagonal matrix
%        lambda  any point at which to evaluate the characteristic polynomial,
%                not an eigenvalue, because eigenvalues are obtained by setting
%                the characteristic polynomial equal to zero and seeing for
%                which values it vanishes.
% Output:
%        y       value of characteristic polynomial at point lambda

    n=length(A);
    D=zeros(n+1,1);
    % application of Thomas algorithm for characteristic polynomial calculation
    D(1)=1;
    D(2)=A(1,1)-lambda;
    for i=2:n
        D(i+1)=(A(i,i)-lambda)*D(i)-A(i-1,i)*A(i,i-1)*D(i-1);
    end
    y=D(end);
end