function result = trapezoidal(x, func)

ba = x(end)-x(1);
N = length(x) - 1;

trapz_vec = 2*ones(1,N+1);
trapz_vec(1,1) = 1;
trapz_vec(1,N+1) = 1;

result = ba/(2*N) * trapz_vec*func;
