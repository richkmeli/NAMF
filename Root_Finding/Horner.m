function [y, dy] = Horner(coeffp, x)
% HORNER'S ALGORITHM
% Algorithm that allows evaluating a polynomial x^n+a_1*x^(n-1)+...+a_n
% performing n additions and n multiplications. This allows evaluating a
% polynomial with fewer multiplications than traditional evaluation (n(n+1)/2).
% Ideal algorithm when searching for real roots with iterative methods.
% Input:
%        coeffp    polynomial coefficients
%        x         point at which to calculate polynomial and its derivative
% Output:
%        y         polynomial value at point x
%        dy        derivative value at point x

    % Calculate polynomial and derivative values at point x,
    % using the recursive definition
    n=length(coeffp)-1;
    y= coeffp(1);
    for i=2:n+1
        y= y*x + coeffp(i);
    end
    % Horner's method allows calculating the first and second derivatives
    % of p_n(x) with fewer operations. The first derivative will be of the
    % form x^(n-1)+b_1*x^(n-2)+...+b_(n-1)
    % number of output arguments required, to avoid derivative calculation
    % if not requested
    if nargout > 1
        % calculate polynomial derivative
        dy= coeffp(1)*n;
        for i= 2: n
            dy=dy*x+coeffp(i)*(n-i+1);
        end
    end
end