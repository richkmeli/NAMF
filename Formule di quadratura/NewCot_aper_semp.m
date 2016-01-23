function Q=NewCot_aper_semp(fun,a,b,n)
% NEWTON COTES APERTE SEMPLICE (su un intervallo)
%------------------INPUT
% fun : funzione integranda
% a,b : estremi di integrazione
% n+1 : numero di nodi
%-----------------OUTPUT
% Q : quadratura risultante
%-----------------------

    h=(b-a)/(n+2);
    % intervallo su cui costruire la quadratura
    x=linspace(a+h,b-h,n+1)';

    switch n % grado
        case 0  % formula dei rettangoli o del punto medio
            alpha=2;
        case 1
            alpha=[3/2 3/2];
    end
    Q=h*alpha*fun(x);
end