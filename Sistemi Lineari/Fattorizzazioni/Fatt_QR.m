function [Q,R] = Fatt_QR (A) 
% FATTORIZZAZIONE QR
%------------------INPUT
% A : matrice da fattorizzare
%-----------------OUTPUT
% Q : Matrice prodotto di tutti i riflettori elementari
% R : Matrice triangolare superiore
%-----------------------
% Applichiamo l'algoritmo di fattorizzazione QR alla matrice A. iteriamo
% per ottenere i riflettori elementari di ogni sottomatrice di A per poi
% moltiplicarli alla A precendente, infatti a ogni iterazione A varierà a
% seconda del riflettore elementare a cui è stata moltiplicata, ed infine
% nell'ultimo passo del ciclo la matrice A sarà la matrice triangolare
% superiore R.
    [m n]=size(A);
    % matrice rispetto alla base canonica, utilizzata per ricavare il
    % vettore v
    Me=zeros(m,n)+diag(diag(ones(m,n)));
    for i=1:m-1
          % vettore contenente prima riga di A
          x=A(i:m,i);
          % A non è univocamete determinata, quindi è possibile scegliere
          % la norma di x di segno positivo o negativo, la scelta piu
          % opportuna per calcolare v è che x+norm_x sia di segno concorde
          % in modo che la somma non sia 0.
          if x(1) >0
            nor_x=norm(x,2);
          else x(1) <0
            nor_x=-norm(x,2);
          end
          v=x + nor_x.*Me(i:m,i);
          % riflettore elementare: matrice di dimensione della i-esima
          % sottomatrice ( è matrice di hauseholder: simm,ortog,idempoten),
          % cioè costruiamo una matrice che trasformi la matrice A in una
          % triangolare superiore.
          P=eye(m+1-i,n+1-i) - ((2*v*v')/(norm(v,2))^2);
          % Costruiamo la matrice P di dimensione sempre m,n
          tmpMe=Me;
          tmpMe(i:m,i:n)=zeros(m+1-i,n+1-i) + P;
          P=tmpMe;
          % Modifichiamo la matrice A per l'iterata successiva
          A = P*A;
          % Costruisce la matrice Q facendo il prodotto tra matrici del
          % riflettore elementare precedente(Q)(mantiene ortogononalità 
          % e il riflettore elementare secondo la i sottomatrice. utilizziamo
          % il costrutto if per %inizializzare nella prima iterata Q con il
          % valore del primo riflettore
          if i == 1
              Q=P;
          else
              Q=Q*P;
          end
    end
    R = A;
end