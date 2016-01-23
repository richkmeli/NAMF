function [ rad_min, rad_max ] = rad_minmax( coeffp , toll, max_iter)
% RADICE MINIMA E MASSIMA
% Calcola la radice reale minima e massima con una tolleranza 'toll' entro 
% un determinato numero di iterazioni di un  polinomio i cui coefficienti 
% sono 'coeffp' con il metodo di newton.
%------------------INPUT
%   coeffp : coefficienti del polinoio
%   toll : tolleranza per il metodo di newton
%   max_iter : massimo numero di iterazioni
%------------------OUTPUT
%   rad_min : minima radice reale del polinomio
%   rad_max : massima radice reale del polinomio
%------------------------

    n=length(coeffp)-1;
    % METODO DI HIRSH
    % norma inf ,massimo della somma degli elementi sulle righe
    % vettore contenente la somma delle righe, primo elemento trattato
    % singolarmente perche non vi sono 1 sulla prma riga della mat di frob.
    sum_rig=[abs((coeffp(n+1)/coeffp(1))) 1+abs(coeffp(n:-1:2)/coeffp(1))];
    max_cerc_rig = max(sum_rig);
    % norma 1, massimo della somma degli elementi sulle colonne
    % vettore contenente la somma delle colonne
    sum_col=[ 1 sum(abs(coeffp(2:n+1)/coeffp(1)))];
    max_cerc_col = max(sum_col);
    V = min(max_cerc_rig, max_cerc_col);
    
    rad_min = rad_pol_newton(coeffp,-V, toll,max_iter);
    rad_max = rad_pol_newton(coeffp,V, toll,max_iter);
end

