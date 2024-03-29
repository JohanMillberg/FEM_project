%Sets all parameters
rho = 10;
alpha = 0.01;
R = 0.5;
r = 0.3;
T = 30;
time = 0:0.1:T;
mesh_sizes = [1/5, 1/20];

geometry = @circleg;

mass_loss = zeros(2,length(time));

%Calculates the solution for both mesh sizes
for k = 1:length(mesh_sizes)
    hmax = mesh_sizes(k);
    [p, e, t] = initmesh(geometry, 'hmax', hmax);
    I = eye(length(p));
    [A,b, M] = assemble2(p,e,t);
    A = alpha.*A;
    x = p(1, e(1,:));
    y = p(2, e(1,:));
    
    A(e(1,:),:) = I(e(1,:),:);
    b(e(1,:)) = 0;

    xi_0 = zeros(length(b), 1);
    
    for i = 1:length(p)
        x = p(1,i);
        y = p(2,i);
        if (x^2 + y^2)^(1/2) >= r && (x^2 + y^2)^(1/2) < R
            xi_0(i,1) = rho;
        end
    end
    %Uses Crank-Nicholson to make the time discretisation
    xi = zeros(length(b),length(time));
    xi(:,1) = xi_0;
    for i = 2:length(time)
        kn = time(i)- time(i-1);
        xi(:,i) = (1/kn.*M+(1/2).*A)\((1/kn.*M-(1/2) .* A)* ...
        xi(:,i-1)+b);
      
    end
    %Calculates the mass loss
    for i = 1:length(time)
    g = @(index) xi(index,1) - xi(index,i);
    mass_loss(k,i) = integration_2D(g,p,t);
    end
end

%Plots the mass loss

figure(1)
hold on
xlabel('t')
ylabel('Mass loss')
title('Mass loss over time')
plot(time,mass_loss(1,:),time,mass_loss(2,:))
legend('h_{max} = 1/5', 'h_{max} = 1/20')
hold off

