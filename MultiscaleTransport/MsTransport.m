function MsTransport
Lx = 1;
Nx_H = 200; dx_H = 2*Lx/Nx_H; x0_H = -Lx:dx_H:Lx;
Nx_h = 10; dx_h = dx_H/Nx_h; x0_h = 0:dx_h:Lx;

sigma = @(x)cos(50*pi*x)+1.1;

Phi = cell(Nx_H,1);
basis = cell(Nx_H+1,1);% for km = 1:Nx_H basis{km} = zeros(1,2);end
for km = 1:Nx_H
    segment = x0_H(km) + [0:dx_h:dx_H]; segment = segment(:);
    Phi_temp = BasisConstruct(sigma,segment);
    Phi{km} = Phi_temp;
    basis_temp = basis{km}; basis_segment = [Phi_temp(:,1),segment];
    basis_temp=[basis_temp;basis_segment];
    basis{km} = basis_temp;
    
    basis_temp = basis{km+1}; basis_segment = [Phi_temp(:,2),segment];
    basis_temp=[basis_temp;basis_segment];
    basis{km+1} = basis_temp;
end
% for km = 1:Nx_H+1
%     basis_temp = basis{km};
%     plot(basis_temp(:,2),basis_temp(:,1),'.-.');pause;
% end
Sigma_mat = zeros(Nx_H+1,Nx_H+1);
Sigma_x_mat = zeros(Nx_H+1,Nx_H+1);
for km = 1:Nx_H
    segment = x0_H(km) + [0:dx_h:dx_H];
    sigma_piece = sigma(segment); sigma_piece = sigma_piece(:);
    sigma_piece2 = sigma(segment(1:end-1)+dx_h/2); sigma_piece2 = sigma_piece2(:);

    Phi_temp = Phi{km}; Phi_temp_end = Phi_temp([1,end],:);
    Sigma_mat_temp = Phi_temp'*diag(sigma_piece)*Phi_temp - Phi_temp_end'*diag(sigma_piece(1,end))*Phi_temp_end/2;
    Sigma_mat([km,km+1],[km,km+1]) = Sigma_mat([km,km+1],[km,km+1]) + dx_h*Sigma_mat_temp;
    
    pxPhi = (Phi_temp(2:end,:) - Phi_temp(1:end-1,:))/dx_h;
    SigmaPhi = (sigma_piece*ones(1,2)).*Phi_temp;
    pxSigmaPhi = (SigmaPhi(2:end,:) - SigmaPhi(1:end-1,:))/dx_h;
    Sigma_x_mat_temp = pxSigmaPhi'*diag(sigma_piece2)*pxPhi;
    Sigma_x_mat([km,km+1],[km,km+1]) = Sigma_x_mat([km,km+1],[km,km+1]) - dx_h*Sigma_x_mat_temp/3;
end

t = 0; t_Final = 1; dt = dx_H/5; kt = 1;
rho_ini = cos(pi*x0_H(2:end-1)/2); rho_ini = rho_ini(:); plot(rho_ini);pause;
rho = zeros(Nx_H-1,t_Final/dt); rho(:,1) = rho_ini; rho_old = rho_ini;
Sigma_mat_s = Sigma_mat(2:end-1,2:end-1); Sigma_x_mat_s = Sigma_x_mat(2:end-1,2:end-1);
Sigma_inv = Sigma_mat_s - dt*Sigma_x_mat_s;
while t < t_Final
%     rho_new = Sigma_mat_s*rho_old + dt*Sigma_x_mat_s*rho_old/3;
%     rho_new = Sigma_mat_s\rho_new;
    rho_new = Sigma_inv\(Sigma_mat_s*rho_old);
%     plot(x0_H,[0;rho_new;0],'.-.',x0_H,sigma(x0_H));pause;
    rho(:,kt+1) = rho_new; rho_old = rho_new;
    t = t+dt; kt = kt+1;
end

save(['rho_',num2str(Nx_H),'_',num2str(Nx_h),'.mat']);
end

function Phi = BasisConstruct(sigma,segment)
Nx = length(segment); dx_h = segment(2)-segment(1);
sigma_value = sigma(segment(1:end-1)+dx_h); sigma_value = sigma_value(:);
sigg = [sigma_value(2:end),-sigma_value(1:end-1)-sigma_value(2:end),[0;sigma_value(2:end-1)]];
StiffA = spdiags(sigg,[-1:1],Nx-2,Nx-2);
Load = zeros(Nx-2,2); Load(1,1) = -sigma_value(1); Load(end,end) = -sigma_value(end);
Phi = StiffA\Load;
Phi = [1,0;Phi;0,1];
end