function [x, resnorm, residual] = levmarq(func, x0)
%LEVMARQ Summary of this function goes here
%   Detailed explanation goes here

x = x0; 
r = 100;
while norm(r) < 0.001
    [r, J] = func(x);
    mu = 0
    p = inv(J'*J+mu*eye(length(x)))*J'*r;
    x = x + p';
end

[r, J] = func(x);
resnorm = norm(r);
residual = r;
end