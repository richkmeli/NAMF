function x = RSL_SI(U,b)
%   RISOLUZIONE SISTEMA LINEARE DI UNA MATRICE TRIANGOLARE SUPERIORE
%   SOSTITUZIONE ALL'INDIETRO

        n = length(b);
        % ricaviamo l'incognita dell'ultima riga
        x(n) = b(n)/U(n,n);
        % ricaviamo le i incognite con i=n-1,...,1
        for i=n-1: -1 :1
        % ricaviamo l'elemento i dell'incognita
        x(i) = (b(i)-U(i,i+1:n)*x(i+1:n)') / U(i,i);
        end
        x=x';
    end