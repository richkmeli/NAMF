function Q=NewCot_aper_semp(fun,a,b,n)
% SIMPLE OPEN NEWTON-COTES (on one interval)
% Input:
%        fun   integrand function
%        a,b   integration limits
%        n+1   number of nodes
% Output:
%        Q     resulting quadrature

    h=(b-a)/(n+2);
    % interval on which to build the quadrature
    x=linspace(a+h,b-h,n+1)';

    switch n % degree
        case 0  % midpoint rule
            alpha=2;
        case 1
            alpha=[3/2 3/2];
    end
    Q=h*alpha*fun(x);
end