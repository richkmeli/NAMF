function [ x ] = RisolSisMatTrid( M,b )
% Risoluzione di sistemi con matrice tridiagonale utilizzando il metodo di
% fattorizzazione LU
% input: M: matrice tridiagonale associata al sistema
%        b: vettore termnini noti
% x: vettore soluzione del sistema
    n=length(M);
    % estraiamo la diagonale principale, inferiore e superiore dalla
    % matrice M
    ad=diag(M,0);
    bd=[0 diag(M,-1)'];
    cd=[diag(M,1)' 0];
    % per definizione ricaviamo alfa, beta e c(coicide con cd)
    alfa(1)=ad(1);
    for i=2:n
        beta(i)=bd(i)/alfa(i-1);
        alfa(i)=ad(i)-beta(i)*cd(i-1);
    end
    % togliamo l'elemento in piu in testa al vettore beta
    beta(1)=[];
    cd(end)=[];
    % L: matrice bidiagonale inferiore (specifica per tridiagonale)
    L=zeros(n)+diag(ones(1,n),0)+diag(beta,-1);
    % U: matrice bidiagonale superiore (specifica per tridiagonale)
    U=zeros(n)+diag(alfa,0)+diag(cd,1);
    
    %risolviamo il sistema L*z=b
    z=L\b;
    %risolviamo il sistema U*x=z
    x=U\z;
    
end

