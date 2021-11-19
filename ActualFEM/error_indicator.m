function [sum_eta2, eta2] = error_indicator(x, func, alpha, xi, A)

M = mass_laplacian(x);

laplacian = -M\(A*xi);
f = func(x);
N = length(x) - 1;
eta2 = zeros(N,1);
for i = 1:N
    h = x(i+1) - x(i);
    eta2(i) = (h^2*trapezoidal([x(i) x(i+1)]', ...
        [(f(i) + alpha*laplacian(i)).^2 (f(i+1) + alpha*laplacian(i+1))^2]'));
end

sum_eta2 = sum(eta2);