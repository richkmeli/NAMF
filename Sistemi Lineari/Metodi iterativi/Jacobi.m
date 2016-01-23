function [err, xnew] = Jacobi(A, b, x0, iter_max, toll)
% valutiamo l'errore del metodo di Jacobi iterando fino a raggiungere una
% determinata tolleranza(toll).
% ------------INPUT----------------
% iter_max : numero massimo iterazioni
% toll : tolleranza
% ---------------------------------

    D = diag(diag(A));
    B = D-tril(A);
    C = D-triu(A);
    
    invD = diag(1./diag(A));
    % matrice di iterazione 
    J = invD*(B+C);
    % esaminiamo la convergenza del metodo valutando la norma infinita
    % degli autovalori di J
    rag_spet = norm(eig(J),inf);
    if ( rag_spet < 1 )
        disp('Convergente');
    else
        disp('Non convergente');
    end
    
    xold = x0;
    i = 1;
    test = 1;
    err = zeros(1, iter_max);
    while i < iter_max && test > toll
        % ricaviamo il vettore soluzione di quella spefica iterazione i
        xnew = J*xold+invD*b;
        err(i) = norm(xnew-xold);
        test = err(i);
        xold = xnew;
        i = i+1;
    end
    err = err(1:i-1);
    % numero di iterazioni effettuate
    if ( rag_spet < 1 )
        disp(' in '); i 
        disp(' iterazioni');
    end
end

