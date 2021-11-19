function sum_eta2 = error_indicator(x, func, alpha, xi, A)

M = mass_laplacian(x);

laplacian = -M\(A*xi);

N = length(x) - 1;
eta2 = zeros(N,1);
for i = 1:N
    f = @(x)func(x) + alpha*(laplacian(i));
    h = x(i+1) - x(i);
    sub_interval = linspace(x(i),x(i+1),50);
    eta2(i) = (h^2*trapezoidal(sub_interval, f));
end

sum_eta2 = sum(eta2);