function [max_iter] = max_iter_trapezi(der2_funz, a, b, toll)
% MAXIMUM ITERATIONS FOR TRAPEZOIDAL RULE
% Calculation of the maximum iterations necessary to approximate with a
% given tolerance 'toll' the integral of the function whose second derivative
% is 'der2_funz', using the trapezoidal rule.
% The admissible error must be bounded by the tolerance
% Input:
%        der2_funz  second derivative of the integrand function
%        a, b       integration limits
%        toll       required tolerance to approximate the integral
% Output:
%        max_iter   maximum iterations necessary

    % points at which to evaluate the second derivative, 101 for midpoint
    x= linspace(a,b,101);
    % step size, specific case of trapezoids n=1
    h=(b-a)/1;
    % from the numerical integration error formula we can derive:
    % upper bound of the second derivative value
    maggioraz = norm(der2_funz(x),inf);
    H = sqrt(12*toll/(maggioraz * h));
    max_iter = h / H;
end