function [A, M] = assemble_heat(p,e,t)

N = size(p,2);
A = sparse(N,N);
M = sparse(N,N);

for K = 1:size(t,2)
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    [area_K, b, c] = hat_gradients(x,y);
    AK = (b*b' + c*c') * area_K;
    MK = [2, 1, 1; 1 2 1; 1 1 2;]/12 * area_K;
    A(nodes,nodes) = A(nodes, nodes) + AK;
    M(nodes, nodes) = M(nodes, nodes) + MK;
end
    