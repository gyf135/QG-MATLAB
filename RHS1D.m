function dahatdt = RHS1D(at_hat,kx)
    dahatdt = -kx.*kx.*at_hat;
end