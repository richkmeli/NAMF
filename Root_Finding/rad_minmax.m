function [rad_min, rad_max] = rad_minmax(coeffp, toll, max_iter)
% MINIMUM AND MAXIMUM ROOTS
% Calculates the minimum and maximum real roots with tolerance 'toll' within
% a certain number of iterations of a polynomial whose coefficients are 'coeffp'
% using Newton's method.
% Input:
%        coeffp    polynomial coefficients
%        toll      tolerance for Newton's method
%        max_iter  maximum number of iterations
% Output:
%        rad_min   minimum real root of the polynomial
%        rad_max   maximum real root of the polynomial

    n=length(coeffp)-1;
    % HIRSCH METHOD
    % infinity norm, maximum of the sum of elements on rows
    % vector containing the sum of rows, first element treated
    % separately because there are no 1s on the first row of Frobenius matrix
    sum_rig=[abs((coeffp(n+1)/coeffp(1))) 1+abs(coeffp(n:-1:2)/coeffp(1))];
    max_cerc_rig = max(sum_rig);
    % 1-norm, maximum of the sum of elements on columns
    % vector containing the sum of columns
    sum_col=[ 1 sum(abs(coeffp(2:n+1)/coeffp(1)))];
    max_cerc_col = max(sum_col);
    V = min(max_cerc_rig, max_cerc_col);
    
    rad_min = rad_pol_newton(coeffp,-V, toll,max_iter);
    rad_max = rad_pol_newton(coeffp,V, toll,max_iter);
end