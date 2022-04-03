function Heat_1D_period
% set up
Lx = 1;
Nx = 64; dx = 2*Lx/Nx; x0 = -Lx:dx:Lx;

sigma = 1/4.1;

% stiffness matrix construction
e = ones(Nx,1)/dx^2;
StiffA = spdiags([e -2*e e],[-1:1],Nx,Nx);
StiffA(1,end) = 1/dx^2; StiffA(end,1) = 1/dx^2;
StiffA = sigma*StiffA;

% evolution
t = 0; t_Final = 0.1; dt = dx/10; time_counter = 1;
rho_ini = cos(pi*x0(1:end-1)/2); rho_ini = rho_ini(:); figure(1); plot(rho_ini);pause;
rho = zeros(Nx,floor(t_Final/dt)); rho(:,1) = rho_ini; rho_old = rho_ini;
Sigma_inv = eye(Nx) - dt*StiffA/3;
while t < t_Final
    rho_new = Sigma_inv\rho_old;
    figure(1); plot(x0(1:end-1),rho_new);pause;
    rho(:,time_counter+1) = rho_new; rho_old = rho_new;
    t = t+dt; time_counter = time_counter+1;
end

save(['1D/conv_delta/heat_hom_64.mat']);
end