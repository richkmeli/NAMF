function [ x_new, num_iter ] = newton_mat_trid_sim( M , toll )
%RAD_MAT_TRID_SIM calcola l'autovalore massimo(destro) della matrice
%------------------INPUT
% M : matrice tridiagonale simmetrica
% toll : tolleranza
%-----------------OUTPUT
% x_new : autovalore massimo della matrice
% num_iter: numero iterazione
%-----------------------
    % troviamo il punto iniziale dove applicare newton con il metodo di
    % Gerschgorin. dal secondo teorema di Gerschgorin, in ogni cerchio si 
    % trova esattamente un autovalore.
    n= size(M,1);
    % R: vettore di lunghezza n, che contiene i raggi dei cerchi di Gerschgorin
    % R(i) = somma dei valori assoluti degli elementi non diagonali della riga
    % i-esima
    R = sum(abs(M-diag(diag(M))));
    % punto estremo del cerchio di Gerschgorin, punto del centro + punto del raggio
    x0 = max(diag(M)'+R);
    x_new = x0;
    % valore che varia per ogni iterazione
    x_old = x_new + toll +1;
    num_iter=0;
    while (abs(x_new - x_old) > toll)
        % VALUTAZIONE DEL POLINOMIO CARATTERISTICO
        % y_old : valore del polinomio carattatteristico in x_old, funzione 
        % dove valutiamo la x_old
        y_old = 0;
        D=zeros(n+1,1);
        % determinante della matrice di dimensione 0
        D(1)=1;
        % determinante della matrice di dimensione 1
        D(2)=M(1,1)-x_old;
        % determinante della matrice di dimensione da 2 a n
        for i = 2:n
            %D_n = det(J-lambda*I)
            D(i+1)=det(M(1:i,1:i)-x_new*eye(i));
        end
        y_old = D(n+1);

        % VALUTAZIONE DERIVATA DEL POLINOMIO CARATTERISTICO
        % Dy_old : valore della derivata prima del polinomio caratteristico.
        Dy_old = 0;
        D_deriv = zeros(n+1,1);
        % derivata determinante della matrice di dimensione 0
        D_deriv(1) = 0;
        % derivata determinante della matrice di dimensione 1
        D_deriv(2) = -1;
        % derivata determinante della matrice di dimensione da 2 a n
        for i = 2:n
           % D'_n=D_(n-1)(lambda) + (alfa-lambda)D'_(n-1)-(beta_n)^n*D'_(n-2)(lambda)
           % alpha(M(i,i)): elem diag principale, betaM(i-1,i): elem diag sup e ing
           D_deriv(i+1)=-D(i) + (M(i,i)-x_new)*D_deriv(i) - M(i,i-1)^2 *D_deriv(i-1);
        end
        Dy_old = D_deriv(n+1);
        % Calcolo del nuovo punto
        x_old=x_new;
        x_new=x_new - y_old/Dy_old;
        num_iter=num_iter+1;
    end
    
end

