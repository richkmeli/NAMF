function Q=NewCot_chiu_comp(fun,a,b,n,N)
% NEWTON COTES CHIUSE COMPOSITE
%------------------INPUT
% fun : funzione integranda
% a,b : estremi di integrazione
% n+1 : numero di nodi
% N numero di suddivisioni di [a,b]
%-----------------OUTPUT
% Q : quadratura risultante
%-----------------------
    
    H=(b-a)/N;
    % intervallo su cui costruire la quadratura
    X=linspace(a,b,N+1)';
   Q=0;
   for i=1:N
       Q=Q + NewCot_chiu_semp(fun,X(i),X(i+1),n);
   end
end