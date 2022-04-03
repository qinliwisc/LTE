function MsHeat_1D_period(delta)
% set up
Lx = 1;
Nx_H = 64; dx_H = 2*Lx/Nx_H; x0_H = -Lx:dx_H:Lx;
Nx_h = 40; dx_h = dx_H/Nx_h; x0_h = 0:dx_h:Lx;

sigma = @(x) 1./(cos(2*pi*delta*x)+4.1);
% sigma = @(x) 1./(cos(20*pi*x)+4.1);

% stiffness matrix construction
Sigma_mat = zeros(Nx_H+1,Nx_H+1);
Sigma_x_mat = zeros(Nx_H+1,Nx_H+1);
for km = 1:Nx_H
    segment = x0_H(km) + [0:dx_h:dx_H]; segment = segment(:);
    [Sigma_mat_temp,Sigma_x_mat_temp,~] = BasisConstruct(sigma,segment);
    Sigma_mat([km,km+1],[km,km+1]) = Sigma_mat([km,km+1],[km,km+1]) + Sigma_mat_temp;
    Sigma_x_mat([km,km+1],[km,km+1]) = Sigma_x_mat([km,km+1],[km,km+1]) + Sigma_x_mat_temp;
end
Sigma_mat_s = Sigma_mat(1:end-1,1:end-1);
Sigma_mat_s(1,:) = Sigma_mat(1,1:end-1)+Sigma_mat(end,1:end-1);
Sigma_mat_s(:,1) = Sigma_mat(1:end-1,1)+Sigma_mat(1:end-1,end);
Sigma_mat_s(1,1) = Sigma_mat(1,1) + Sigma_mat(end,end);
% Sigma_x_mat = -Sigma_x_mat^2;
Sigma_x_mat_s = Sigma_x_mat(1:end-1,1:end-1);
Sigma_x_mat_s(1,:) = Sigma_x_mat(1,1:end-1)+Sigma_x_mat(end,1:end-1);
Sigma_x_mat_s(:,1) = Sigma_x_mat(1:end-1,1)+Sigma_x_mat(1:end-1,end);
Sigma_x_mat_s(1,1) = Sigma_x_mat(1,1) + Sigma_x_mat(end,end);

% evolution
t = 0; t_Final = 0.1; dt = dx_H/5; time_counter = 1;
rho_ini = cos(pi*x0_H(1:end-1)/2); rho_ini = rho_ini(:); figure(2); plot(rho_ini);pause;
rho = zeros(Nx_H,floor(t_Final/dt)); rho(:,1) = rho_ini; rho_old = rho_ini;
Sigma_inv = Sigma_mat_s - dt*Sigma_x_mat_s/3;
while t < t_Final
    rho_new = Sigma_inv\(Sigma_mat_s*rho_old);
    figure(2); plot(x0_H(1:end-1),rho_new);pause;
    rho(:,time_counter+1) = rho_new; rho_old = rho_new;
    t = t+dt; time_counter = time_counter+1;
end

save(['1D/conv_delta/heat_delta_',num2str(delta),'_64_40.mat']);
end

function [Sigma_mat,Sigma_x_mat,Phi] = BasisConstruct(sigma,segment)
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
    
%     pxPhi = (Phi(2:end,:) - Phi(1:end-1,:))/dx_h;
%     SigmaPhi = (sigma_piece*ones(1,2)).*Phi; pxSigmaPhi = (SigmaPhi(2:end,:) - SigmaPhi(1:end-1,:))/dx_h;
%     Sigma_x_mat = -dx_h*pxSigmaPhi'*diag(sigma_piece2)*pxPhi;
    
    pxPhi = (Phi(2:end,:) - Phi(1:end-1,:))/dx_h;
    Sigma_x_mat = dx_h*pxPhi'*diag(sigma_piece2)*Phi_c;
%     
    pxPhi = (Phi(2:end,:) - Phi(1:end-1,:))/dx_h;
    SigmaPhi = (sigma_piece*ones(1,2)).*Phi; pxSigmaPhi = (SigmaPhi(2:end,:) - SigmaPhi(1:end-1,:))/dx_h;
    Sigma_x_mat = -dx_h*pxPhi'*diag(sigma_piece2)*pxSigmaPhi;

end