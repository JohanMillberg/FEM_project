function [A,b_vec] = assemble1(p,e,t)

N = size(p,2);
A = sparse(N,N);
b_vec = zeros(N,1);

for K = 1:size(t,2)
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    [area_K, b, c] = hat_gradients(x,y);
    xc = mean(x); yc = mean(y);
    AK = (b*b' + c*c') * area_K;
    bK = source_term1(xc,yc)*area_K / 3;
    A(nodes,nodes) = A(nodes, nodes) + AK;
    b_vec(nodes) = b_vec(nodes) + bK;
end
    
   