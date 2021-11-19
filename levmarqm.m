function [x, resnorm, residual] = levmarqm(func, x0)
%LEVMARQ Summary of this function goes here
%   Detailed explanation goes here

x = x0;
r = 100;
lambda = 10;
nu = 1.5;

while norm(r) > 0.4
    [r, J] = func(x);
    p1 = (J'*J+lambda*eye(length(x)))\(J'*r);
    p2 = (J'*J+lambda/nu*eye(length(x)))\(J'*r);
    [r1, J1] = func(x+p1);
    [r2, J2] = func(x+p2);
    
    if norm(r1) < norm(r2)
        if norm(r1) < norm(r)
            x = x + p1;
        else
            lambda = lambda*nu;
        end
    else
        if norm(r2) < norm(r)
            x = x + p2;
            lambda = lambda/nu;
        else
            lambda = lambda*nu;
        end
    end
end

[r, J] = func(x);
resnorm = norm(r);
residual = r;
end