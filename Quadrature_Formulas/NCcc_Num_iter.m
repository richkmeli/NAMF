function [quadr, num_iter] = NCcc_Num_iter(funz,a,b,n,toll,max_iter)
% COMPOSITE CLOSED NEWTON-COTES WITH ITERATION COUNT
% Number of iterations to approximate the integral from a to b of function funz
% on n+1 nodes, to a tolerance of toll within a maximum of max_iter iterations
% Input:
%        funz      integrand function
%        a,b       integration limits
%        n         quadrature degree, n+1 nodes
%        toll      tolerance, difference between exact value and quadrature
%        max_iter  maximum number of iterations
% Output:
%        quadr     quadrature value
%        num_iter  number of iterations performed to reach tolerance < toll

val_esa=integral(funz,a,b);
quadr = NewCot_chiu_comp(funz,a,b,n,1);
% initial error
test=abs(val_esa-quadr);
num_iter=1;
while (test>toll && num_iter < max_iter)
    num_iter=num_iter+1;
    % calculate quadrature with (num_iter) subintervals
    quadr = NewCot_chiu_comp(funz,a,b,n,num_iter);
    test=abs(val_esa-quadr);
end