function [x, resnorm, residual] = levmarq(func, x0)
%LEVMARQ Summary of this function goes here
%   Detailed explanation goes here

x = x0;
for i = 1:10
    [r, J] = func(x);
    p = inv(J'*J)*J'*r;
    
    x = x + p';
end

[r, J] = func(x);
resnorm = norm(r);
residual = r;
end