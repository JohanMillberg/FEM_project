function J = jacobian_fd(f, x, t, h)
% Function used to numerically calculate the Jacobian using finite
% difference
%f: function head
%x: x vector, t: t vector, h: step length

m = length(t);
n = length(x);
J = zeros(m,n);

for i = 1:n
    for j = 1:m
        y = f(x, t(j));
        new_x = x;
        new_x(i) = new_x(i) + h;
        new_y = f(new_x, t(j));
        J(j,i) = new_y - y;
    end
end

J = J / h;
