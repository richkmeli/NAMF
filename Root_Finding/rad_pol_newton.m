function [x_new] = rad_pol_newton(coeffp, x0, toll, max_iter)
% POLYNOMIAL ROOT FINDING WITH NEWTON'S METHOD
% Calculates the real root with a certain tolerance 'toll' within 'max_iter'
% iterations, starting from 'x0' with Newton's method
% Input:
%        coeffp    polynomial coefficients
%        x0        initial point
%        toll      tolerance, difference between root approximation between
%                  one iteration (x_old) and the next (x_new)
%        max_iter  maximum number of iterations
% Output:
%        x_new     real root of the polynomial

    x_old=x0;
    test=toll+1;
    num_iter = 0;
    while (test > toll && num_iter < max_iter)
        num_iter=num_iter+1;
        [y_x, dy_x] = Horner(coeffp,x_old);
        % check on first derivative, if it's a maximum or minimum point, stop
        if dy_x == 0
            test = 0;
        else
            x_new = x_old -(y_x/dy_x);
            test = abs(x_old - x_new);
        end
        x_old = x_new;
    end
end