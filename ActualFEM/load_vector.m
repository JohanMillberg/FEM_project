function B = load_vector(x, rho, R, r)

N = length(x) - 1;
B = zeros(N+1, 1);

for i = 1:N
    f = @(x) rho*(abs(R-abs(x)) <= r);
    h = x(i+1) - x(i);
    n = [i i+1];
    B(n) = B(n) + [f(x(i)); f(x(i+1))]*h/2;
end
% If nonhomogoenous Dirichlet B(end) = 7e6;