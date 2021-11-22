% Script used to demonstrate the Levenberg-Marquardt algorithm

fun = @(b)residualfunc(b);
x0 = [1;1];

options = optimset();

[x_sol, rn, r, lambda, iter] = levmarqm(fun, x0);
[x_sol_mat, rn_mat] = lsqnonlin(fun, x0);

t = [0.5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0]';
y = [6.8, 3.0, 1.5, 0.75, 0.48, 0.25, 0.2, 0.15]';

[a,b,approx_y] = fun(x_sol);
[a,b,approx_y_mat] = fun(x_sol_mat);

disp(lambda(end))
figure()
hold on
plot(t,y,'x',t, approx_y(x_sol, t))
xlabel('t');
ylabel('y');
legend(['Data points', 'Approximated solution '+string(round(x_sol(1),2))+'*exp('+string(round(x_sol(2),2))+'*t)']);
hold off

figure()
disp(lambda(end))
hold on
plot(t, approx_y_mat(x_sol_mat, t), 'x-', t, approx_y(x_sol, t))
xlabel('t');
ylabel('y');
legend(['Lsqnonlin', 'Approximated solution '+string(round(x_sol(1),2))+'*exp('+string(round(x_sol(2),2))+'*t)']);
hold off

figure()
hold on
plot(iter, lambda)
xlabel('iteration');
ylabel('Damping factor \lambda')