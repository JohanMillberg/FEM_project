function J = jacobian_fd(f, x, h)

y = f(x);
m = length(y);
n = length(x);
J = zeros(m,n);

for i = 1:n
    new_x = x;
    new_x(i) = new_x(i) + h;
    new_y = f(new_x);
    J(:,i) = new_y - y;
end

J = J / h;
