function dahatdt = RHS2D(at_hat,kx,ky,I,J,dump, IJ2, tau_d, tau_f, nu, beta, psiR_hat)
    nx = length(kx);
    ny = length(ky);
    
    % Reshape to get the q1_hat and q2_hat matrices
    q1_hat_vec = at_hat(1:nx*ny);
    q2_hat_vec = at_hat(nx*ny+1:end);
    
    q1_hat = reshape(q1_hat_vec,nx,ny);
    q2_hat = reshape(q2_hat_vec,nx,ny);

    % Calculate psi from q
    psi1_hat = 0.5*(-(q1_hat+q2_hat)./(dump) - (q2_hat-q1_hat)./(-dump-2));
    psi2_hat = 0.5*(-(q1_hat+q2_hat)./(dump) + (q2_hat-q1_hat)./(-dump-2));
    
    % Calculate x derivative
    dpsi1dx_hat = 1i*I.*psi1_hat;
    dpsi2dx_hat = 1i*I.*psi2_hat;
    
    dq1dx_hat = 1i*I.*q1_hat;
    dq2dx_hat = 1i*I.*q2_hat;
    
    % Calculate y derivative
    dpsi1dy_hat = 1i*J.*psi1_hat;
    dpsi2dy_hat = 1i*J.*psi2_hat;
    
    dq1dy_hat = 1i*J.*q1_hat;
    dq2dy_hat = 1i*J.*q2_hat;
    
    % Non-linear terms
    NL1 = fft2((ifft2(dpsi1dx_hat)).*(ifft2(dq1dy_hat)) - (ifft2(dpsi1dy_hat)).*(ifft2(dq1dx_hat)));
    NL2 = fft2((ifft2(dpsi2dx_hat)).*(ifft2(dq2dy_hat)) - (ifft2(dpsi2dy_hat)).*(ifft2(dq2dx_hat)));
    
    % Time stepping
    RHS1 = -NL1 + (psi1_hat - psi2_hat - psiR_hat)/tau_d - nu*IJ2.*IJ2.*q1_hat;
    RHS2 = -NL2 - (psi1_hat - psi2_hat - psiR_hat)/tau_d + IJ2.*(psi2_hat)/tau_f - nu*IJ2.*IJ2.*q2_hat;

    RHS1_vec = reshape(RHS1, nx*ny, 1);
    RHS2_vec = reshape(RHS2, nx*ny, 1);
    
    dahatdt = [RHS1_vec;RHS2_vec];

end