function Q=NewCot_chiu_semp(fun,a,b,n)
% SIMPLE CLOSED NEWTON-COTES (on one interval)
% Input:
%        fun   integrand function
%        a,b   integration limits
%        n+1   number of nodes
% Output:
%        Q     resulting quadrature

    h=(b-a)/n;
    % interval on which to build the quadrature
    x=linspace(a,b,n+1)';
    
    switch n % degree
        case 1  % trapezoidal rule
            alpha=[1/2 1/2];
        case 2  % Cavalieri-Simpson rule
            alpha=[1/3 4/3 1/3];
        case 3  % three-eighths rule
            alpha=[3/8 9/8 9/8 3/8];
        case 4
            alpha=[14/45 64/45 24/45 64/45 14/45];
        case 5
            alpha=[95/288 375/288 250/288 250/288 375/288 95/288];
    end
    Q=h*alpha*fun(x);
end