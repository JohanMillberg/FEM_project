function b_vec = source_c(M,xi_u,xi_v)

b_vec = M*(xi_u.*(xi_v).^2);