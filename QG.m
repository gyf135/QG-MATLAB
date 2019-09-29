clear variables; close all; clc
% Physical params
tau_f = 15;
tau_d = 100;
nu    = 0.1;
beta  = 0.196;
sigma = 3.5;
tf    = 8000;

% Dimensions
Lx    = 46;
Ly    = 23;
M     = 40;
% dt    = 10/M;


% Wave numbers and viscosity
nx    = 50;
ny    = 50;

p     = 8.0;
% 
% nux   = 8.0*((Lx/nx/pi)^p)/dt/8.0;
% 
% nuy   = 8.0*((Ly/ny/pi)^p)/dt/8.0;

dx    = Lx/nx;
dy    = Ly/ny;
x     = linspace(0, Lx-dx, nx);
y     = linspace(0, Ly-dy, ny);

dt = 0.0002*dx^4/nu;

kx = 2*pi/Lx*[0:(nx/2) (-nx/2+1):-1];
ky = 2*pi/Ly*[0:(ny/2) (-ny/2+1):-1];
[J,I]  = meshgrid(ky,kx);

I2 = I.^2;
J2 = J.^2;

psiR = zeros(nx,ny);
for j = 1:length(y)
    psiR(:,j) = -sigma * tanh(y(j)/sigma);
end

psiR_hat = fft2(psiR);

psiInt1 = zeros(nx,ny);
psiInt2 = zeros(nx,ny);
for j = 1:length(y)
    for i = 1:length(x)
%         psiInt1(i,j) = 1e-2*sin(2*pi*x(i)/Lx)*sin(2*pi*y(j)/Ly);%*exp(-y(j)^2/(2*sigma^2));
%         psiInt2(i,j) = 1e-3*sin(2*pi*x(i)/Lx)*sin(2*pi*y(j)/Ly);%*exp(-y(j)^2/(2*sigma^2));
        psiInt1(i,j) = 1e-2*sin(2*pi*x(i)/Lx)*exp(-y(j)^2/(2*sigma^2));
        psiInt2(i,j) = 1e-3*sin(2*pi*x(i)/Lx)*exp(-y(j)^2/(2*sigma^2));

    end
end

psiInt1_hat = fft2(psiInt1);
psiInt2_hat = fft2(psiInt2);

q1_hat = -(I2+J2).*psiInt1_hat - (psiInt1_hat - psiInt2_hat);
q2_hat = -(I2+J2).*psiInt2_hat + (psiInt1_hat - psiInt2_hat);

q1 = real(ifft2(q1_hat));
q2 = real(ifft2(q2_hat));

for (j=1:ny)
    q1(:,j) = q1(:,j) + beta*y(j);
    q2(:,j) = q2(:,j) + beta*y(j);
end

q1_hat = fft2(q1);
q2_hat = fft2(q2);

plot(q1(length(y)/2,:),'k','LineWidth',2);
hold on


NSTEP = 10000;
    
IJ2 = I2+J2;
dump = IJ2;
dump(1,1) = 1;

% Stiff ODE solver
q1_hat_vec = reshape(q1_hat,nx*ny,1);
q2_hat_vec = reshape(q2_hat,nx*ny,1);

a_hat = [q1_hat_vec; q2_hat_vec];
t = linspace(0,10,11);
[t,at_hat] = ode15s(@(t,at_hat) RHS2D(at_hat,kx,ky,I,J,dump, IJ2, tau_d, tau_f, nu, beta, psiR_hat), t, a_hat);


for k = 1:size(at_hat,1)
    q1_hat_vec = at_hat(k, 1:nx*ny);
    q2_hat_vec = at_hat(k, nx*ny+1:end);
    
    q1 = real(ifft2(reshape(q1_hat_vec,nx,ny)));
    q2 = real(ifft2(reshape(q2_hat_vec,nx,ny)));
    
    plot(x,q1(length(y)/2,:));
    
    hold on;
end


%  Start time marching (Explicit method)
for n = 1:NSTEP
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
    
%     RHS1 =  -nu*(IJ2).*(IJ2).*q1_hat;
%     RHS2 =  -nu*(IJ2).*(IJ2).*q2_hat;
%         
    q1_hat = q1_hat + dt*RHS1;
    q2_hat = q2_hat + dt*RHS2;
    
    q1 = real(ifft2(q1_hat));
    q2 = real(ifft2(q2_hat));
%     figure
%     surf(q1);
    if (mod(n,NSTEP/10)==0)
        plot(q1(length(y)/2,:));
        hold on
    end
    
%     shading interp
end





% Test in 1D
% a = sin(2*pi/Lx*x);
% a_hat = fft(a);
% % for n = 1:5000
% % dadxx_hat = -kx.*a_hat;
% % % dadxx = real(ifft(dadxx_hat));
% % 
% % a_hat = a_hat-0.01*dt*dadxx_hat;
% % 
% % end
% 
% 
% t = linspace(0,1000,11);
% [t,at_hat] = ode15s(@(t,at_hat) -(kx.').*(kx.').*at_hat, t, a_hat);
% 
% 
% for k = 1:length(t)
%     an = real(ifft(at_hat(k,:)));
%     plot(x,an);
%     hold on
%     axis([0 Lx, -1 1]);
% end




