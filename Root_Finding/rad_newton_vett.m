function [rad, num_iter] = rad_newton_vett(funz,der_funz,x0, toll, max_iter)
% NEWTON'S METHOD FOR ROOT FINDING
% Calculates the root and the number of iterations needed to obtain a root
% according to a given tolerance within a maximum number of iterations,
% starting from point x0
% Input:
%        funz      function for which to find the root
%        der_funz  derivative of the function
%        x0        initial point
%        toll      minimum acceptable variation between root approximation
%                  and that of the previous iteration
%        max_iter  maximum number of iterations
% Output:
%        rad       root approximation
%        num_iter  number of iterations in which the root is found

    num_iter=0;
    rad = zeros(max_iter+1,1);
    rad(1) = x0;
    test = toll+1;
    while (test > toll && num_iter < max_iter && der_funz(rad(num_iter+1))~=0)
        num_iter=num_iter+1;
        rad(num_iter+1)=rad(num_iter)-(funz(rad(num_iter))/der_funz(rad(num_iter)));
        % evaluate tolerance as the difference between current and previous root
        test = abs(rad(num_iter+1)-rad(num_iter));
    end
    rad = rad(num_iter+1);
end

