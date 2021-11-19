function B = load_vector(x, rho, R, r)

N = length(x) - 1;
B = zeros(N+1, 1);

for i = 1:N
    if abs(R-abs(x)) <= r
        f = @(x) rho*x^0;
    else
        f = @(x) 0*x^0;
    end
    h = x(i+1) - x(i);
    n = [i i+1];
    B(n) = B(n) + [f(x(i)); f(x(i+1))]*h/2;
end
% If nonhomogoenous Dirichlet B(end) = 7e6;