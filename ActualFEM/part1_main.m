a = -1;
b = 1;
N = 11;

alpha = 0.01;
rho = 10;
R = 0.5;
r = 0.3;
lambda = 0.9;

h = (b-a)/N;
x = (a:h:b)';

f = @(x) rho*(abs(R-abs(x)) <= r);

A = stiffness(x);
B = load_vector(x, rho, R, r, alpha);
xi = A\B;
[error, eta2] = error_indicator(x, f, alpha, xi, A);

TOL = 1e-3;

while  error > TOL && length(x) < 1e4
    for i = 1:length(eta2)
        if eta2(i) > lambda*max(eta2)
            x = [x; (x(i+1)+x(i))/2];
        end
    end
    x = sort(x);
    A = stiffness(x);
    B = load_vector(x, rho, R, r, alpha);
    xi = A\B;
    [error, eta2] = error_indicator(x, f, alpha, xi, A);
end
subplot(2,2,1)
plot(x,xi)
title('FEM-Solution')

subplot(2,2,2)
plot(x(2:end),[1./diff(x)])
title('Mesh distribution')

subplot(2,2,3)
M = mass_laplacian(x);
laplacian = -M\(A*xi);
plot(x, f(x)+alpha*laplacian)
title('Residuals')

subplot(2,2,4)
plot(x(1:end-1),eta2)
title('\eta(u_{h})')
