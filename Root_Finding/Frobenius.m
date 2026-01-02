function [Mat_Frob] = Frobenius(coeffp)
% FROBENIUS MATRIX
% Constructs the Frobenius matrix of the polynomial with coefficients 'coeffp'
% Input:
%        coeffp    polynomial coefficients
% Output:
%        Mat_Frob  Frobenius matrix

    n=length(coeffp)-1;
    Mat_Frob=diag(ones(n-1,1),-1);
    Mat_Frob(:,n)=-coeffp(n+1:-1:2)'/coeffp(1);
end