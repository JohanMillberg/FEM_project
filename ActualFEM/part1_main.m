a = -1;
b = 1;
N = 11;

alpha = 0.01;
rho = 10;
R = 0.5;
r = 0.3;

h = (b-a)/N;
x = (a:h:b)';

f = @(x) 2.*x.^0;

TOL = 1e-3;
error = 1;
while  error > TOL || N > 1e4
    A = stiffness(x, alpha);
    B = load_vector(x, rho, R, r);
    xi = A\B;
    error = error_indicator(x, f, alpha, xi, A);
end