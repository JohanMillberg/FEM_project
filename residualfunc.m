function [fun,r,rgrad] = residualfunc(x)
%RESIDUALFUNC Summary of this function goes here
%   Detailed explanation goes here

f = @(x, t) x(1)*exp(x(2)*t);

t = [0.5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0]';
y = [6.8, 3.0, 1.5, 0.75, 0.48, 0.25, 0.2, 0.15]';
fun = f(x, t)
r = (y-f(x, t));
%rgrad = jacobian_fd(f, x, t, 0.025);
rgrad = [exp(x(2)*t), (t).*x(1).*exp(x(2)*t)];
end

