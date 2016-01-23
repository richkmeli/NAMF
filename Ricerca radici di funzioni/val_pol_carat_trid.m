function [ y ] = val_pol_carat_trid( A , lambda )
% VAUTAZIONE POLINOMIO CARATTERISTICO TRIDIAGONALI
%   Calcola le immagini del polinomio caratteristico della matrice A nei
%   punto lambda, p(lambda)=det(A-lambda*Id).
%------------------INPUT
%   A : matrice tridiagonale
%   lambda: punto qualsiasi in cui valutare il polinomio caratteristo, non è
%           autovalore, perche l'autovalore si ottiene ponendo uguale a zero il
%           polinomio caratteristico e vedendo per quali valori si annulla.
%------------------OUTPUT
%   y : immagine del polinomio caratteristico nel punto lambda
%-----------------------
    n=length(A);
    D=zeros(n+1,1);
    % applicazione dell'algoritmo di Thomas per il calcolo del polinomio
    % caratteristico.
    D(1)=1;
    D(2)=A(1,1)-lambda;
    for i=2:n
        D(i+1)=(A(i,i)-lambda)*D(i)-A(i-1,i)*A(i,i-1)*D(i-1);
    end
    y=D(end);
end

