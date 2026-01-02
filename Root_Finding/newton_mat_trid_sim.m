function [x_new, num_iter] = newton_mat_trid_sim(M, toll)
% NEWTON METHOD FOR SYMMETRIC TRIDIAGONAL MATRIX EIGENVALUES
% Calculates the maximum (rightmost) eigenvalue of the matrix
% Input:
%        M     symmetric tridiagonal matrix
%        toll  tolerance
% Output:
%        x_new     maximum eigenvalue of the matrix
%        num_iter  number of iterations

    % Find the initial point where to apply Newton with Gerschgorin's method.
    % From Gerschgorin's second theorem, exactly one eigenvalue is found in each circle.
    n= size(M,1);
    % R: vector of length n, containing the radii of Gerschgorin circles
    % R(i) = sum of absolute values of non-diagonal elements of row i
    R = sum(abs(M-diag(diag(M))),2)';
    % extreme point of Gerschgorin circle, center point + radius point
    x0 = max(diag(M)'+R);
    x_new = x0;
    % value that varies for each iteration
    x_old = x_new + toll +1;
    num_iter=0;
    
    while (abs(x_new - x_old) > toll)
        % CHARACTERISTIC POLYNOMIAL EVALUATION
        % y_old: value of characteristic polynomial at x_old
        D=zeros(n+1,1);
        % determinant of matrix of dimension 0
        D(1)=1;
        % determinant of matrix of dimension 1
        D(2)=M(1,1)-x_old;
        % determinant of matrix of dimension from 2 to n
        for i = 2:n
            D(i+1)=(M(i,i)-x_old)*D(i) - M(i-1,i)*M(i,i-1)*D(i-1);
        end
        y_old = D(n+1);

        % CHARACTERISTIC POLYNOMIAL DERIVATIVE EVALUATION
        % Dy_old: value of first derivative of characteristic polynomial
        D_deriv = zeros(n+1,1);
        % derivative of determinant of matrix of dimension 0
        D_deriv(1) = 0;
        % derivative of determinant of matrix of dimension 1
        D_deriv(2) = -1;
        % derivative of determinant of matrix of dimension from 2 to n
        for i = 2:n
           % D'_n=D_(n-1)(lambda) + (alpha-lambda)D'_(n-1)-(beta_n)^2*D'_(n-2)(lambda)
           % alpha(M(i,i)): main diagonal element, beta M(i-1,i): super/sub diagonal element
           D_deriv(i+1)=-D(i) + (M(i,i)-x_old)*D_deriv(i) - M(i,i-1)^2 *D_deriv(i-1);
        end
        Dy_old = D_deriv(n+1);
        
        % Calculate new point
        x_old=x_new;
        if Dy_old ~= 0
            x_new=x_new - y_old/Dy_old;
        else
            break;
        end
        num_iter=num_iter+1;
    end
end