function Q=NewCot_aper_comp(fun,a,b,n,N)
% COMPOSITE OPEN NEWTON-COTES
% Input:
%        fun   integrand function
%        a,b   integration limits
%        n+1   number of nodes
%        N     number of subdivisions of [a,b]
% Output:
%        Q     resulting quadrature

    H=(b-a)/N;
    % interval on which to build the quadrature
    X=linspace(a,b,N+1)';
    Q=0;
    for i=1:N
        Q=Q + NewCot_aper_semp(fun,X(i),X(i+1),n);
    end
end