function quadr = trap_comp_integr_dop(a,b,c,d,m,n,funct)
% Calcolo dell'integrale doppio tramite prodotto di due quadratura 
% con metodo dei trapezi composite.
% ------------INPUT----------------
% a,b: Estremi integrazione su x
% c,d: Estremi integrazione su y
% m: iterate della formula dei trapezi per l'integrale su x, cioè
%    sottintervalli in dividere l'intervallo di x
% n: iterate della formula dei trapezi per l'integrale su y, cioè
%    sottintervalli in dividere l'intervallo di y
% funct: Funzione da integrare su x,y
% ---------------------------------
    quadr = 0;
    % passo dei sottintervalli riferiti all'integrale di x
    H = (b-a)/m;
    % passo dei sottintervalli riferiti all'integrale di y
    K = (d-c)/n;
    % punti dei sottointervalli
    x = linspace(a,b,m+1);
    y = linspace(c,d,n+1);
    % sommatoria di tutti i valori delle quadrature calcolate nei
    % sottointervalli
    for i = 1:m
       for j = 1:n
          quadr = quadr + funct(x(i),y(j)) + funct(x(i),y(j+1)) + funct(x(i+1),y(j)) + funct(x(i+1),y(j+1)) ;
       end
    end
    % complietiamo la formula dei trapezi composita con h/2 e k/2
    quadr = (H/2)*(K/2)*quadr;
end