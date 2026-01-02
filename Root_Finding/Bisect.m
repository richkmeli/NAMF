function [x,fx,n]=Bisect(f,x1,x2,toll)
% Bisection method for finding roots of a function
% Input: 
%        f       string containing the function in x
%        x1,x2   left and right endpoints of the interval
%        toll    tolerance on the interval
% Output: 
%        n       number of iterations 
%        x       approximation of the zero
%        fx      function value at x

x=x1;
f1=eval(f);
if f1==0
    x=x1;
    fx=f1;
    n=0;
    return;
end
x=x2;
f2=eval(f);
if f2==0
    x=x2;
    fx=f2;
    n=0;
    return;
end
if sign(f1)*sign(f2)> 0
    error('** ERROR ** f(x1)*f(x2) > 0 - No sign change in interval');
end
% n = number of iterations necessary for tolerance toll
n=fix(log(abs(x2-x1)/toll)/log(2)+1); % +1 because fix truncates
for i=1:n
    x=x1+(x2-x1)/2; % more accurate than (x1+x2)/2
    fx=eval(f);
    if sign(f1)*sign(fx)>0
        x1=x;
        f1=fx;
    else
        x2=x;
    end
end
return