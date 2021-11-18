function [x, resnorm, residual] = levmarq(func, x0)
%LEVMARQ Summary of this function goes here
%   Detailed explanation goes here

x = x0;
r = 100;
delta_k = 1;
lambda = 0;
mu = 1/4;
eta = 3/4;

while norm(r) > 0.35
    [r, J] = func(x);
    lambda = 0;
    p = (J'*J+lambda*eye(length(x)))\(J'*r);
    disp(inv(J'*J+lambda*eye(length(x))))
    if norm(p) > delta_k
        lambdafun = @(l)(norm((J'*J+l*eye(length(x)))\(J'*r))-delta_k);
        lambda = fzero(lambdafun, 100);
        p = (J'*J+lambda*eye(length(x)))\J'*r;
    end
    disp(lambda)
        
    q = norm(f + J*p);
    [r_new, J_new] = func(x+p);
    rho = norm(r - r_new)/norm(-J*p);
    
    if rho <= mu
        delta_k = 1/2*delta_k;
    elseif rho < eta
        x = x + p;
    else
        delta_k = 2*delta_k;
        x = x + p;
    end
end

[r, J] = func(x);
resnorm = norm(r);
residual = r;
end