function Q=NewCot_chiu_semp(fun,a,b,n)
% NEWTON COTES CHIUSA SEMPLICE (su un intervallo)
%------------------INPUT
% fun : funzione integranda
% a,b : estremi di integrazione
% n+1 : numero di nodi
%-----------------OUTPUT
% Q : quadratura risultante
%-----------------------

    h=(b-a)/n;
    % intervallo su cui costruire la quadratura
    x=linspace(a,b,n+1)';
    
    switch n % grado
        case 1  % formula dei trapezi
            alpha=[1/2 1/2];
        case 2  % formula di Cavalieri-Simpson
            alpha=[1/3 4/3 1/3];
        case 3  % formula dei tre ottavi
            alpha=[3/8 9/8 9/8 3/8];
        case 4
            alpha=[14/45 64/45 24/45 64/45 14/45];
        case 5
            alpha=[95/288 375/288 250/288 250/288 375/288 95/288];

    end
    Q=h*alpha*fun(x);
end