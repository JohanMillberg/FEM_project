function [x, resnorm, residual] = levmarq(func, x0)
%LEVMARQ Summary of this function goes here
%   Detailed explanation goes here

x = x0; 
r = 100;
delta_k = 0.01;
lambda = 1/4;
mu = 1/4;
eta = 3/4;

while norm(r) < 0.001
    [r, J] = func(x);
    p = inv(J'*J+lambda*eye(length(x)))*J'*r;
    x = x + p';
end

[r, J] = func(x);
resnorm = norm(r);
residual = r;
end