geom = unitsquare();

[p,e,t] = initmesh(geom,'hmax',0.25);

[A,R,b,r] = assemble(p,e,t);
I = eye(length(p));
A(e(1,:),:) = I(e(1,:),:);

b(e(1,:)) = 0;

xi = (A+R)\(b+r);