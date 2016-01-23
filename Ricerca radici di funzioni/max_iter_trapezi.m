function [ max_iter ] = max_iter_trapezi( der2_funz , a , b , toll )
%	MASSIME ITERATE TRAPEZI
%   Calcolo delle massime iterazioni necessarie per approssimare con una
%   determinata tolleranza 'toll' l'integrale della funzione la cui derivata
%   seconda è 'der2_funz', si utilizza il metodo dei trapezi.
%   l'errore ammissibile deve essere maggiorato dalla tolleranza
%------------------INPUT
%   der2_funz : derivata seconda della funzione integranda
%   a , b : estremi di integrazione
%   toll : tolleranza richiesta per approssimare l'integrale
%------------------OUTPUT
%   max_iter : massime iterate necessarie
%------------------------

    % punti in cui valutare la derivata seconda, 101 per punto medio
    x= linspace(a,b,101);
    % passo, caso specifico dei trapezi n=1
    h=(b-a)/1;
    % dalla formula dell'errore di integrazione numerica possiamo ricavare:
    % maggiorazione del valore della derivata seconda
    maggioraz = norm(der2_funz(x),inf);
    H = sqrt(12*toll/(maggioraz * h));
    max_iter = h / H;
end

