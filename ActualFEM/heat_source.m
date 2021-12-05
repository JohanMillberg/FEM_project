function f = heat_source(beta,gamma,M,xi)

f = beta .* M * (xi-gamma.*xi.^2);