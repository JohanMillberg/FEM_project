fun = @(b)residualfunc(b);
x0 = [1,1];

[x_sol, rn, r] = levmarq(fun, x0)