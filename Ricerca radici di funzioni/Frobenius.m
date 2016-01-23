function [ Mat_Frob ] = Frobenius( coeffp )
% MATRICE FROBENIUS
% costruisce la matrice di frobenius del polinomio con i coefficienti
% 'coeffp'
%-----------------INPUT
%   coeffp : coefficienti del polinomio
%-----------------OUTPUT
%   Mat_Frob : matrice di frobenius
%-----------------------
%   
    n=length(coeffp)-1;
    Mat_Frob=diag(ones(n-1,1),-1);
    Mat_Frob(:,n)=-coeffp(n+1:-1:2)'/coeffp(1);

end

