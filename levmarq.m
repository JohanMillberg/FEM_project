function [x, resnorm, residual] = levmarq(func, x0)
%LEVMARQ Summary of this function goes here
%   Detailed explanation goes here

x = x0; 
r = 100;
delta_k = 10;
lambda = 1/4;
mu = 1/4;
eta = 3/4;

for i = 1:5
    [f, r, J] = func(x);
    lambdafun = @(l)norm((J'*J+l*eye(length(x)))\J'*r)-delta_k;
    lambda = fzero(lambdafun, 1);
    p = (J'*J+lambda*eye(length(x)))\J'*r;
    
    q = norm(f + J*p);
    [f_new, r_new, J_new] = func(x+p);
    rho = norm(f - f_new)/norm(-J*p);
    
    if rho <= mu
        delta_k = 1/2*delta_k
    elseif rho < eta
        x = x + p
    else
        delta_k = 2*delta_k
        x = x + p
end

[r, J] = func(x);
resnorm = norm(r);
residual = r;
end