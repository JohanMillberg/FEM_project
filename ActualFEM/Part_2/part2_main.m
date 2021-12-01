u_exact = @(x1, x2) 8.*pi^2.*sin(2.*pi.*x1).*sin(2*pi.*x2);

geometry = @circleg;
hmax = 1/5;
[p, e, t] = initmesh(geometry, 'hmax', hmax);
I = eye(length(p));
[A,b] = assemble(p,e,t);

x = p(1, e(1,:));
y = p(2, e(1,:));

A(e(1,:),:) = I(e(1,:),:);
b(e(1,:)) = u_exact(x,y);

xi = A\b;