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

TOL = 1e-3;
error = 1;
while  error > TOL || N > 1e4
    A = stiffness(x, alpha);
    B = load_vector(x, rho, R, r);
    xi = A\B;
    [error, eta2] = error_indicator(x, f, alpha, xi, A);
    for i = 1:length(eta2)
        if eta2(i) > lambda*max(eta2)
            x = [x; (x(i+1)+x(i))/2];
        end
    end
    x = sort(x);
end

A = stiffness(x, alpha);
B = load_vector(x, rho, R, r);
xi = A\B;
plot(x,xi)