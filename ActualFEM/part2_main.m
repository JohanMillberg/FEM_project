u_exact = @(x1, x2) sin(2.*pi.*x1).*sin(2*pi.*x2);

geometry = @circleg;

mesh_sizes = [1/2, 1/4, 1/8, 1/16, 1/32];
energy_norms = zeros(size(mesh_sizes));
convergence_rate = zeros(1,length(mesh_sizes)-1);
old_error = 0;

% Calculates and plots the solution using the available mesh sizes
for k = 1:length(mesh_sizes)
    hmax = mesh_sizes(k);
    [p, e, t] = initmesh(geometry, 'hmax', hmax);
    I = eye(length(p));
    [A,b] = assemble1(p,e,t);

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
    figure(k)
    pdemesh(p,e,t,xi)

end

%Calculates h_max^gamma
hmax_gamma = zeros(1,length(mesh_sizes) - 1);
for j = 1:length(mesh_sizes)
    hmax_gamma(j) = mesh_sizes(j)^convergence_rate(end);
end

figure()
hold on
loglog(mesh_sizes,hmax_gamma, ...
    mesh_sizes,energy_norms)
legend('h_{max}^{\gamma}','Energy norm')
xlabel('h_{max}')
hold off


