fun = @(b)residualfunc(b);
x0 = [10;1];

[x_sol, rn, r] = levmarqm(fun, x0)


t = [0.5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0]';
y = [6.8, 3.0, 1.5, 0.75, 0.48, 0.25, 0.2, 0.15]';
approx_y = fun(x_sol);

%legend(['Given points', 'Approximated solution\n'+string(x_sol(1))+'
plot(t,y,'x',t, approx_y)