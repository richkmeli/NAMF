function [quadr, num_iter] = NCcc_Num_iter(funz,a,b,n,toll,max_iter)
% NEWTON COTES CHIUSE COMPOSITE E NUMERO ITERAZIONI
% numero iterate per approssimare l'integrale da a a b della funzione funz 
% su n+1 nodi,ad una tolleranza di toll entro un massimo di iterate max_iter
%------------------INPUT
%   funz : funzione integranda
%   a,b : estremi di integrazione
%   n : grado della quadratura, n+1 nodi
%   toll : tolleranza, differenza tra il valore esatto e la quadratura
%   max_iter : massimo numero di iterazioni
%------------------OUTPUT
%   quadr : valore della quadratura
%   num_iter : numero di iterate effettuate per raggiungere un valore di
%              tollerenza minore di toll
%------------------------
val_esa=integral(funz,a,b);
quadr = NewCot_chiu_comp(funz,a,b,n,1);
% errore commesso inizialmente
test=abs(val_esa-quadr);
num_iter=1;
while (test>toll && num_iter < max_iter)
    num_iter=num_iter+1;
    % calcolo la quadratura con (num_iter) sottointervalli
    quadr = NewCot_chiu_comp(funz,a,b,n,num_iter);
    test=abs(val_esa-quadr);
end