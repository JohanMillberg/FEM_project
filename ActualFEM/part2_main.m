u_exact = @(x1, x2) sin(2.*pi.*x1).*sin(2*pi.*x2);

geometry = @circleg;

mesh_sizes = [1/2, 1/4, 1/8, 1/16, 1/32];
energy_norms = zeros(size(mesh_sizes));
convergence_rate = zeros(1,length(mesh_sizes)-1);
old_error = 0;

for k = 1:length(mesh_sizes)
    hmax = mesh_sizes(k);
    [p, e, t] = initmesh(geometry, 'hmax', hmax);
    I = eye(length(p));
    [A,b] = assemble(p,e,t);

    x = p(1, e(1,:));
    y = p(2, e(1,:));

    A(e(1,:),:) = I(e(1,:),:);
    b(e(1,:)) = u_exact(x,y);

    x = p(1,:);
    y = p(2,:);

    xi = A\b;

    exact_sol = u_exact(x,y);

    err = exact_sol' - xi;
    
    EnE = sqrt(err'*A*err);
        
    if k >= 2
        convergence_rate(k-1) = log(EnE/energy_norms(k-1)) / ...
            log(mesh_sizes(k)/mesh_sizes(k-1));
    end
    
    energy_norms(k) = EnE;

end
figure(1)
pdemesh(p,e,t,xi)
figure(2)
loglog(mesh_sizes(2:end),convergence_rate, ...
    mesh_sizes(2:end),energy_norms(2:end))
legend('\gamma','Energy norm')
