function [ rad , num_iter ] = rad_newton_vett( funz,der_funz,x0, toll, max_iter)
% RADICE NEWTON
% calcola la radice e in quante iterate si riesce ad ottenere una radice
% secondo una determinata tolleranza entro un numero massimo di iterazioni,
% iniziando dal punto x0
%-----------------INPUT
% funz: funzione su cui trovare la radice
% der_funz: derivata della funzione
% x0: punto iniziale
% toll: minima variazione accettabile tra l'approssimazione della radice e quella
%       della iterata precedente
% max_iter: numero massimo di iterazioni
%-----------------OUTPUT
% rad: approssimazione della radice
% num_iter: numero di iterate in cui viene trovata la radice
%-----------------------
    num_iter=0;
    rad = zeros(max_iter,1);
    rad(1) = x0;
    test = toll+1;
    while ( test > toll && num_iter < max_iter && der_funz(rad(num_iter+1))~=0)
        num_iter=num_iter+1;
        rad(num_iter+1)=rad(num_iter)-(funz(rad(num_iter))/der_funz(rad(num_iter)));
        % valuti la tolleranza come la differenza tra la radice precedente e
        % la radice dell'iteratazione precedente
        test = abs(rad(num_iter+1)-rad(num_iter));
    end
end

