%Sets all parameters

rho = 10;
alpha = 0.01;
R = 0.5;
r = 0.3;
T = 30;
time = 0:0.1:T;
mesh_sizes = [1/5, 1/20];
beta = 1;
gamma = 0.2;

geometry = @circleg;

mass_loss = zeros(2,length(time));

%Calculates the solution for both mesh sizes
for k = 1:length(mesh_sizes)
    hmax = mesh_sizes(k);
    [p, e, t] = initmesh(geometry, 'hmax', hmax);
    I = eye(length(p));
    [A, M] = assemble_heat(p,e,t);
    A = alpha.*A;
    x = p(1, e(1,:));
    y = p(2, e(1,:));
    
    A(e(1,:),:) = I(e(1,:),:);

    xi_0 = zeros(length(A), 1);
    
    for i = 1:length(p)
        x = p(1,i);
        y = p(2,i);
        if (x^2 + y^2)^(1/2) >= r && (x^2 + y^2)^(1/2) < R
            xi_0(i,1) = rho;
        end
    end
    
    xi = zeros(length(A),length(time));
    xi(:,1) = xi_0;
    mass_loss(k,1) = 0;
    b_vec = zeros(length(A),length(time));
    
    %Uses Crank-Nicholson to make the time discretisation
    for i = 2:length(time)
        b = heat_source(beta,gamma,M,xi(:,i-1));
        b(e(1,:)) = 0;
        kn = time(i)- time(i-1);
        xi(:,i) = (1/kn.*M+(1/2).*A)\((1/kn.*M-(1/2) .* A)* ...
        xi(:,i-1)+b);
        b_vec(:,i) = b;
        
        %Calculates the mass loss
        g = @(index) xi(index,1) - xi(index,i);
        mass_loss(k,i) = integration_2D(g,p,t);
        
    end
   
    
end

figure(1)
hold on
xlabel('t')
ylabel('Mass loss')
title('Mass loss over time, \gamma = 0.2')
plot(time,mass_loss(1,:),time,mass_loss(2,:))
legend('h_{max} = 1/5', 'h_{max} = 1/20')
hold off

%for i = 1:length(time)
%    pdemesh(p,e,t,xi(:,i))
%    pause(0.001)
%end