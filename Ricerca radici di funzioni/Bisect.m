function [x,fx,n]=Bisect(f,x1,x2,toll)
% Input: 
%        f    macro contenente la funzione in x
%        x1,x2   estremi sx e dx dell'intervallo
%        toll    tolleranza sull'intervallo
% Output: 
%        n   numero iterazioni 
%        x   approssimazione dello zero
%        fx  valore della funzione in x
x=x1;
f1=eval(f);
if f1==0;
    x=x1;
    fx=f1;
    n=0;
    return;
end
x=x2;
f2=eval(f);
if f2==0;
    x=x2;
    fx=f2;
    n=0;
    return;
end
if sign(f1)*sign(f2)> 0
    disp('** ERRDR ** f(x1)*f(x2) > 0 '),
    return,
end
% n= numero iter. neces. per prec. toll
n=fix(log(abs(x2-x1)/toll)/log(2)+1); % +1 perche' int tronca
for i=1:n
    x=x1+(x2-x1)/2; % piu' accurato di (x1+x2)/2
    fx=eval (f);
    if sign(f1)*sign(fx)>0;
        xn=x2;
        x1=x;
    else
        xn=x1;
        x2=x;
    end
end
return