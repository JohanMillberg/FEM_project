fun = @(b)residualfunc(b);
x0 = [1;1];

[x_sol, rn, r] = levmarqm(fun, x0)


t = [0.5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0]';
y = [6.8, 3.0, 1.5, 0.75, 0.48, 0.25, 0.2, 0.15]';
[a,b,approx_y] = fun(x_sol);

hold on
plot(t,y,'x',t, approx_y(x_sol, t))
xlabel('t');
ylabel('y');
legend(['Given points', 'Approximated solution '+string(round(x_sol(1),2))+'*exp('+string(round(x_sol(2),2))+'*t)']);
