function B = load_vector(x, rho, R, r, alpha)

N = length(x) - 1;
B = zeros(N+1, 1);

f = @(x) rho*(abs(R-abs(x)) <= r);

for i = 1:N
    h = x(i+1) - x(i);
    n = [i i+1];
    B(n) = B(n) + [f(x(i)); f(x(i+1))]*h/2;
end
B = B.*1/alpha;