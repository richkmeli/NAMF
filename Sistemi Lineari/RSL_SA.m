function x = RSL_SA(L,b)
%   RISOLUZIONE SISTEMA LINEARE DI UNA MATRICE TRIANGOLARE INFERIORE
%   SOSTITUZIONE IN AVANTI

    n = length(b);
    % ricaviamo l'incognita della prima riga
    x(1) = b(1)/L(1,1);
    % ricaviamo le i incognite con i=2,...,n
    for i=2:n
    % ricaviamo l'elemento i dell'incognita
    x(i) = (b(i)-L(i,1:i-1)*x(1:i-1)') / L(i,i);
    end
    x=x';
end