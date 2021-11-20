function [r,rgrad,f] = residualfunc(x)
%Function computing the residual 
%x: x vector

t = [0.5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0]';
y = [6.8, 3.0, 1.5, 0.75, 0.48, 0.25, 0.2, 0.15]';

f = @(x, t) x(1)*exp(x(2)*t);

%Optional residual gradient variable
%res_grad = exp(x(2)*t), (t).*x(1).*exp(x(2)*t);

r = (y-f(x, t));

if ~exist('grad', 'var')
    rgrad = jacobian_fd(f, x, t, 0.005);
else
    rgrad = [res_grad(x)];
end

