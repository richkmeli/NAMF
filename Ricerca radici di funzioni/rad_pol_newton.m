function [ x_new ] = rad_pol_newton( coeffp, x0, toll, max_iter)
% RADICE POLINOMIO CON METODO DI NEWTON
% calcola la radice reale con una certa tolleranza 'toll' entro 'max_iter'
% iterate , partendo dal 'x0' con il metodo di Newton
%------------------INPUT
%   coeffp : coefficienti del polinoio
%   x0 : punto iniziale
%   toll : tolleranza, differenza tra l'approssimazione della radice tra un
%       iterata(x_old) e la successiva(x_nex)
%   max_iter : massimo numero di iterazioni
%------------------OUTPUT
%   x_new : radice reale del polinomio
%------------------------

    x_old=x0;
    test=toll+1;
    num_iter = 0;
    while (test > toll && num_iter < max_iter)
        num_iter=num_iter+1;
        [y_x dy_x] = Horner(coeffp,x_old);
        % controllo sulla derivata prima, se è un punto di massimo o minimo
        % , si ferma
        if dy_x == 0
            test = 0;
        else
            x_new = x_old -(y_x/dy_x);
            test = abs(x_old - x_new);
        end
        x_old = x_new;
    end
end

