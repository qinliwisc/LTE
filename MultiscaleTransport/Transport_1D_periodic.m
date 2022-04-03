function Transport_1D_periodic
global Flux_A Sigma_x_mat_s Sigma_mat_s Phi_mat_s epsilon dt Nx_H Np

% set up
Lx = 1;
Nx_H = 100; dx_H = 2*Lx/Nx_H; x0_H = -Lx:dx_H:Lx;
Nx_h = 20; dx_h = dx_H/Nx_h; x0_h = 0:dx_h:Lx;
Nv = 40; [v0,W0] = legendre_quad(Nv-1); Np = 20; [poly_l] = legendre_poly(v0,[0:1:Np-1]);

epsilon_index = 6; epsilon = 2^(-epsilon_index);
sigma = @(x)cos(20*pi*x)+4.1;

% stiffness matrix in x
Sigma_mat = zeros(Nx_H+1,Nx_H+1);
Phi_mat = zeros(Nx_H+1,Nx_H+1);
Sigma_x_mat = zeros(Nx_H+1,Nx_H+1);
for km = 1:Nx_H
    segment = x0_H(km) + [0:dx_h:dx_H]; segment = segment(:);
    [Sigma_mat_temp,Sigma_x_mat_temp,Phi_mat_temp,~] = BasisConstruct(sigma,segment);
    Sigma_mat([km,km+1],[km,km+1]) = Sigma_mat([km,km+1],[km,km+1]) + Sigma_mat_temp;
    Phi_mat([km,km+1],[km,km+1]) = Phi_mat([km,km+1],[km,km+1]) + Phi_mat_temp;
    Sigma_x_mat([km,km+1],[km,km+1]) = Sigma_x_mat([km,km+1],[km,km+1]) + Sigma_x_mat_temp;
end
Sigma_mat_s = Sigma_mat(1:end-1,1:end-1);
Sigma_mat_s(1,:) = Sigma_mat(1,1:end-1) + Sigma_mat(end,1:end-1);
Sigma_mat_s(:,1) = Sigma_mat(1:end-1,1) + Sigma_mat(1:end-1,end);
Sigma_mat_s(1,1) = Sigma_mat(1,1) + Sigma_mat(end,end);
Phi_mat_s = Phi_mat(1:end-1,1:end-1);
Phi_mat_s(1,:) = Phi_mat(1,1:end-1) + Phi_mat(end,1:end-1);
Phi_mat_s(:,1) = Phi_mat(1:end-1,1) + Phi_mat(1:end-1,end);
Phi_mat_s(1,1) = Phi_mat(1,1) + Phi_mat(end,end);
Sigma_x_mat_s = Sigma_x_mat(1:end-1,1:end-1);
Sigma_x_mat_s(1,:) = Sigma_x_mat(1,1:end-1) + Sigma_x_mat(end,1:end-1);
Sigma_x_mat_s(:,1) = Sigma_x_mat(1:end-1,1) + Sigma_x_mat(1:end-1,end);
Sigma_x_mat_s(1,1) = Sigma_x_mat(1,1) + Sigma_x_mat(end,end);


t = 0; tFinal = 0.1; dt = min(dx_H/5,epsilon/2);

% f_v = exp(-20*v0.^2); f_v = f_v(:); rho = ones(1,Nv)*diag(W0)*f_v; f_v = 2*f_v/rho;
% f_ini = f_v*cos(pi*x0_H(1:end-1)/2); proj_a = poly_l'*diag(W0)*f_ini;
% f_ini = f_v*cos(pi*x0_H(1:end-1)/2);  proj_a = poly_l'*diag(W0)*f_ini;
f_ini = ones(Nv,1)*cos(pi*x0_H(1:end-1)/2);  proj_a = poly_l'*diag(W0)*f_ini;

% flux in v
Flux_A = poly_l'*diag(v0.*W0)*poly_l; Flux_A(abs(Flux_A)<1e-14) = 0;
% evolution
time_counter = 1;

proj_a_total = zeros(Np,Nx_H,floor(tFinal/dt)); f_total = zeros(Nv,Nx_H,floor(tFinal/dt));
proj_a_total(:,:,time_counter) = proj_a; f_total(:,:,time_counter) = f_ini;

while t < tFinal
    t = t + dt; time_counter = time_counter + 1;

    data_pre = (Sigma_mat_s*proj_a')'; data_pre = data_pre(:);
    proj_a_new = gmres(@MAB_multi,data_pre,50,1e-8); proj_a_new = reshape(proj_a_new,Np,Nx_H);
    
    proj_a = proj_a_new;
    
    f = poly_l*proj_a; rho = proj_a(1,:)/sqrt(2);
    figure(4)
    plot(rho);
%     subplot(2,1,1); plot(rho);
%     subplot(2,1,2); mesh(f);
    pause;
    proj_a_total(:,:,time_counter) = proj_a;
    f_total(:,:,time_counter) = f;
end

save(['1D/20pi_transport_conv/ep_',num2str(epsilon_index),'_100_20.mat']);
end

function proj = MAB_multi(p)
global Flux_A Sigma_x_mat_s Sigma_mat_s Phi_mat_s epsilon dt Nx_H Np
p = reshape(p,Np,Nx_H);
term1 = (Sigma_mat_s*p')';
term2 = Flux_A*p; term2 = dt*(Sigma_x_mat_s*term2')'/epsilon;
term3 = zeros(size(p)); term3(2:end,:) = p(2:end,:); term3 = (Phi_mat_s*term3')'*dt/epsilon^2;
proj = term1 + term2 + term3;
proj = proj(:);
end

function [Sigma_mat,Sigma_x_mat,Phi_mat,Phi] = BasisConstruct(sigma,segment)
    Nx = length(segment); dx_h = segment(2)-segment(1);
    sigma_piece = sigma(segment); sigma_piece = sigma_piece(:);
    sigma_piece2 = sigma(segment(1:end-1)+dx_h/2); sigma_piece2 = sigma_piece2(:);
    sigg = [sigma_piece2(2:end),-sigma_piece2(1:end-1)-sigma_piece2(2:end),[0;sigma_piece2(2:end-1)]];
    StiffA = spdiags(sigg,[-1:1],Nx-2,Nx-2);
    Load = zeros(Nx-2,2); Load(1,1) = -sigma_piece2(1); Load(end,end) = -sigma_piece2(end);
    Phi = StiffA\Load; Phi = [1,0;Phi;0,1]; Phi_c = (Phi(1:end-1,:)+Phi(2:end,:))/2;

    Phi_end = Phi([1,end],:);
    Sigma_mat = Phi'*diag(sigma_piece)*Phi - Phi_end'*diag(sigma_piece(1,end))*Phi_end/2;
    Sigma_mat = dx_h*Sigma_mat;
    
    Phi_mat = Phi'*Phi - Phi_end'*Phi_end/2;
    Phi_mat = dx_h*Phi_mat;
    
    pxPhi = (Phi(2:end,:) - Phi(1:end-1,:))/dx_h;
    Sigma_x_mat = dx_h*pxPhi'*diag(sigma_piece2)*Phi_c;
% 
%     
%     pxPhi = (Phi(2:end,:) - Phi(1:end-1,:))/dx_h;
%     SigmaPhi = (sigma_piece*ones(1,2)).*Phi; pxSigmaPhi = (SigmaPhi(2:end,:) - SigmaPhi(1:end-1,:))/dx_h;
%     Sigma_x_mat = -dx_h*pxSigmaPhi'*diag(sigma_piece2)*pxPhi;

end