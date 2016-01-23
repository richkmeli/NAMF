function [ y, dy ] = Horner( coeffp, x)
% HORNER
% algoritmo che permette di valutare un polinomio x^n+a_1*x^n-1+...+a_n
% svolgendo n addizioni e n moltiplicazioni. quindi riusciamo a valutare un
% polinomio con meno moltiplicazioni di un polinomio tradizionale (n(n+1)/2
% ). Algoritmo ideale quando si cercano radici reali con metodi iterativi.
%-----------------INPUT
%   coeffp : coefficienti del polinomio
%   x : punto in cui calcolare polinomio e derivata del polinomio
%-----------------OUTPUT
%   y : immagine del polinomio nel punto x
%   dy : immagine della derivata del polinomio nel punto x
%----------------------
%   Calcola i valori del polinomio e della sua derivata nel punto x,
%   sfruttando la definizione ricorsiva
    n=length(coeffp)-1;
    y= coeffp(1);
    for i=2:n+1
        y= y*x + coeffp(i);
    end
    % il metodo di horner permette di calcolare con meno operazioni anche
    % la derivata prima e la derivata seconda di p_n(x). la derivata prima
    % sarà della forma x^n-1+b_1*x^n-2+...+b_n-1
    % numero di argomenti richiesti in output, per evitare calcolo derivata
    % se non richiesto
    if nargout > 1
        % calcolo derivata del polinomio
        dy= coeffp(1)*n;
        for i= 2: n
            dy=dy*x+coeffp(i)*(n-i+1);
        end
    end
end

