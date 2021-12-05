function result = integration_2D(g,p,t)
result = 0;
for K = 1:size(t,2)
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    [area_K, ~, ~] = hat_gradients(x,y);

    sum_function = sum(g(nodes)); 
    result = result + area_K / 3 * sum_function;
end
