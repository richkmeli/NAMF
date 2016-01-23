function [ L, U ] = Fatt_LU( A )
% Fattorizzazione LU
% ------------INPUT----------------
% A: matrice da fattorizzare
% ---------------------------------
    % dimensioni matrice n righe, m colonne
    [n,m] = size(A);
    % iterazioni sulle n righe
    L=zeros(n,m);
    U=A;
    % inizializziamo il pivot con valori di default, e quando U(i,i)=0 non
    % vengono scambiate le righe
    pivot=1:length(A);
    for i = 1:n
        % la condizione serve per evitare pivoting non necessari
        if U(i,i) ==0
            % PIVOTING
            % inizializziamo max con la prima riga i
            max = i;
            % impostiamo max uguale alla riga con il massimo valore , cioe
            % scegliamo la riga del pivot
            for j = i+1 : n
                if (abs(U(j,i)) > abs(U(max,i))) 
                    max = j;
                end
            end
            % scambio della riga i con la riga max(del pivot)
            L([i max],:) = L([max i],:);
            % memorizziamo il pivot per la permutazione
            pivot(i) = max;
            U([i max],:) = U([max i],:);
        end
        % SCALING
        for j=i+1 :n
            if U(i,i) ~=0
            % Costruice matrice valori coefficienti di gauss, -U(j,i) e tmp
            % per ottenere il segno corretto rispettivamente sulla matrice L e U
            tmp=-U(j,i)/U(i,i);
            % costruisce matrice riduzione di gauss
            U(j,:)=U(i,:)*tmp + U(j,:);
            L(j,i)=-tmp;
            end
        end
    end
    % Permutazione per ordinare la matrice L, cioè la matrice P è la
    % permutazione che applico su L
     P = diag(ones(size(U,1),1));
     for i = 1:n-1
         P([i pivot(i)],:) = P([pivot(i) i],:);
     end
    L = L + diag(ones(size(U,1),1));
     L = P\L;
end


