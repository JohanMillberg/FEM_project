function [A,b_vec] = assemble_2(p,e,t)

N = size(p,2);
A = sparse(N,N);
b_vec = zeros(N,1);

for K = 1:size(t,2)
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    [area_K, b, c] = hat_gradients(x,y);
    xc = mean(x); yc = mean(y);
    coeffs = cat(2, [1;1;1], x', y') \ eye(3);
    
    AK = zeros(3,3);
    for i=1:3
        for j=1:3
            AK(i,j) = coeffs(2,i)*coeffs(2,j) + coeffs(3,i)*coeffs(3,j); 
        end
    end
    AK = AK * area_K;
    bK = repmat(source_term(xc,yc) * area_K / 3, 3, 1);
    A(nodes,nodes) = A(nodes, nodes) + AK;
    b_vec(nodes) = b_vec(nodes) + bK;
end
    