function M = mass_laplacian(x)

N = length(x) - 1;
M = zeros(N+1,N+1);

for i = 1:N
    h = x(i+1) - x(i);
    n = [i i+1];
    M(n,n) = h*[1/3 1/6; 1/6 1/3];
end