function MsHeat_2D_period
% set up
Lx = 1;
Nx_H = 100; dx_H = 2*Lx/Nx_H; x0_H = -Lx:dx_H:Lx;
[x0,y0] = meshgrid(x0_H,x0_H);
Nx_h = 20; dx_h = dx_H/Nx_h; x0_h = 0:dx_h:Lx;
[x0h,y0h] = meshgrid(x0_h,x0_h);

%sigma_fun = @(x,y)cos(10*pi*x).*cos(2*pi*y)+1;
sigma_fun = @(x,y)(2+1.8*sin(10*pi*x))./(2+1.8*cos(10*pi*y))+(2+sin(10*pi*y))./(2+1.8*sin(10*pi*x));

% stiffness matrix construction
Sigma_mat = zeros((Nx_H)^2,(Nx_H)^2);
Sigma_x_mat = zeros((Nx_H)^2,(Nx_H)^2);
for km = 1:Nx_H % for y
    for kn = 1:Nx_H % for x
        [indices] = index_periodic_map(km,kn,Nx_H);
        segment_x = x0_H(kn) + [0:dx_h:dx_H];
        segment_y = x0_H(km) + [0:dx_h:dx_H];
        [Sigma_mat_temp,Sigma_x_mat_temp,~] = BasisConstruct(sigma_fun,segment_x,segment_y);
        Sigma_mat(indices,indices) = Sigma_mat(indices,indices) + Sigma_mat_temp;
        Sigma_x_mat(indices,indices) = Sigma_x_mat(indices,indices) - Sigma_x_mat_temp/2;
    end
end

% evolution
t = 0; t_Final = 0.1; dt = dx_H/8; time_counter = 1;
rho_ini = cos(pi*x0(1:end-1,1:end-1)).*cos(pi*y0(1:end-1,1:end-1))+1.1;

% rho_ini = sqrt(2)*cos(pi*x0(1:end-1,1:end-1)/2).*cos(pi*y0(1:end-1,1:end-1)/2);

rho = cell(1,t_Final/dt); rho{1} = rho_ini; rho_old = rho_ini(:);
Sigma_inv = Sigma_mat - dt*Sigma_x_mat;
while t < t_Final
    rho_old = Sigma_mat*rho_old;
    rho_new = Sigma_inv\rho_old;
    figure(1); mesh(x0(1:end-1,1:end-1),y0(1:end-1,1:end-1),reshape(rho_new,Nx_H,Nx_H));title(num2str(time_counter));pause(0.01);
    rho{time_counter+1} = reshape(rho_new,Nx_H,Nx_H); rho_old = rho_new;
    t = t+dt; time_counter = time_counter+1
end

save(['2D/ex3_heat/rho_',num2str(Nx_H),'_',num2str(Nx_h),'.mat']);

end

function [Sigma_mat,Sigma_x_mat,Load] = BasisConstruct(sigma_fun,segment_x,segment_y)
    Nx = length(segment_x)-1; dx_h = segment_x(2)-segment_x(1);
    
    [xx,yy] = meshgrid(segment_x,segment_y);
    StiffA = zeros((Nx+1)^2,(Nx+1)^2); Load = zeros((Nx+1)^2,4);
    
    stiff_local_l = [1/2,0,0,-1/2;0,0,0,0;0,0,1/2,-1/2;-1/2,0,-1/2,1];
    stiff_local_r = [1/2,-1/2,0,0;-1/2,1,-1/2,0;0,-1/2,1/2,0;0,0,0,0];
    linear_increase = [0:Nx]/Nx; linear_decrease = -[0:Nx]/Nx+1;
    
    index_matrix = [1:(Nx+1)^2]; index_matrix = reshape(index_matrix,Nx+1,Nx+1);
    edges = [index_matrix(1,:)',index_matrix(:,end),index_matrix(end,:)',index_matrix(:,1)];
    remainder = index_matrix(2:end-1,2:end-1); remainder = remainder(:);
    edges_total = [edges(:,1);edges(2:end,2);edges(1:end-1,3);edges(2:end-1,4)];
    
    sigma_l = sigma_fun(xx(1:end-1,1:end-1) + dx_h/3,yy(1:end-1,1:end-1) + 2*dx_h/3);
    sigma_r = sigma_fun(xx(1:end-1,1:end-1) + 2*dx_h/3,yy(1:end-1,1:end-1) + dx_h/3);
    sigma = sigma_fun(xx,yy);
    sigma_c = sigma_fun(xx(1:end-1,1:end-1) + dx_h/2,yy(1:end-1,1:end-1) + dx_h/2);
    
    %constructing basis
    for km = 1:Nx
        for kn = 1:Nx
            indices = index_map(km,kn,Nx);
            StiffA(indices,indices) = StiffA(indices,indices) + sigma_l(km,kn)*stiff_local_l;
            StiffA(indices,indices) = StiffA(indices,indices) + sigma_r(km,kn)*stiff_local_r;
        end
    end
    Load(edges(:,1),1) = linear_decrease; Load(edges(:,2),1) = zeros(Nx+1,1);
    Load(edges(:,3),1) = zeros(Nx+1,1); Load(edges(:,4),1) = linear_decrease;
    
    Load(edges(:,1),2) = linear_increase; Load(edges(:,2),2) = linear_decrease;
    Load(edges(:,3),2) = zeros(Nx+1,1); Load(edges(:,4),2) = zeros(Nx+1,1);
    
    Load(edges(:,1),3) = zeros(Nx+1,1); Load(edges(:,2),3) = linear_increase;
    Load(edges(:,3),3) = linear_increase; Load(edges(:,4),3) = zeros(Nx+1,1);
    
    Load(edges(:,1),4) = zeros(Nx+1,1); Load(edges(:,2),4) = zeros(Nx+1,1);
    Load(edges(:,3),4) = linear_decrease; Load(edges(:,4),4) = linear_increase;
    
    Load_right = StiffA(remainder,edges_total)*Load(edges_total,:);
    Load(remainder,:) = - StiffA(remainder,remainder)\Load_right;

    % storing basis
    soln = cell(1,4); Sigma_mat = zeros(4,4); Sigma_x_mat = zeros(4,4);
    px_soln = cell(1,4); py_soln = cell(1,4);
    px_sigma_soln = cell(1,4); py_sigma_soln = cell(1,4);
    for kn = 1:4
        soln_temp = reshape(Load(:,kn),Nx+1,Nx+1);
%         mesh(xx,yy,soln_temp);pause;
        soln{kn} = soln_temp;
        px_soln{kn} = (soln_temp(2:end,2:end)+soln_temp(2:end,1:end-1)...
            - soln_temp(1:end-1,2:end) - soln_temp(1:end-1,1:end-1))/dx_h/2;
        py_soln{kn} = (soln_temp(1:end-1,2:end) + soln_temp(2:end,2:end)...
            - soln_temp(1:end-1,1:end-1) - soln_temp(2:end,1:end-1))/dx_h/2;
        
        px_sigma_soln{kn} = (sigma(2:end,2:end).*soln_temp(2:end,2:end)...
            + sigma(2:end,1:end-1).*soln_temp(2:end,1:end-1)...
            - sigma(1:end-1,2:end).*soln_temp(1:end-1,2:end)...
            - sigma(1:end-1,1:end-1).*soln_temp(1:end-1,1:end-1))/dx_h/2;
        py_sigma_soln{kn} = (sigma(1:end-1,2:end).*soln_temp(1:end-1,2:end)...
            + sigma(2:end,2:end).*soln_temp(2:end,2:end)...
            - sigma(1:end-1,1:end-1).*soln_temp(1:end-1,1:end-1)...
            - sigma(2:end,1:end-1).*soln_temp(2:end,1:end-1))/dx_h/2;
    end
    
    for km = 1:4
%         for kn = 1: km-1
%             Sigma_mat(km,kn) = Sigma_mat(kn,km);
%             Sigma_x_mat(km,kn) = Sigma_x_mat(kn,km);
%         end
        for kn = 1:4
            PhiSigmaPhi = soln{km}.*soln{kn}.*sigma;
            Sigma_mat(km,kn) = sum(sum(PhiSigmaPhi))*dx_h*dx_h;
            pSPhiSpPhi = px_sigma_soln{km}.*sigma_c.*px_soln{kn};
            pSPhiSpPhi = pSPhiSpPhi + py_sigma_soln{km}.*sigma_c.*py_soln{kn};
            Sigma_x_mat(km,kn) = sum(sum(pSPhiSpPhi))*dx_h*dx_h;
%             PhiSigmaPhi = soln{km}.*soln{kn};
%             Sigma_mat(km,kn) = sum(sum(PhiSigmaPhi))*dx_h*dx_h;
%             pSPhiSpPhi = px_soln{km}.*sigma_c.*px_soln{kn};
%             pSPhiSpPhi = pSPhiSpPhi + py_soln{km}.*sigma_c.*py_soln{kn};
%             Sigma_x_mat(km,kn) = sum(sum(pSPhiSpPhi))*dx_h*dx_h;
        end
    end
end

function indices = index_periodic_map(km,kn,Nx_H)    
    indices(1) = (km-1)*Nx_H + kn;
    indices(2) = (km-1)*Nx_H + kn+1;
    indices(3) = (km)*Nx_H + kn+1;
    indices(4) = (km)*Nx_H +kn;
    if kn == Nx_H
        indices(2) = (km-1)*Nx_H + 1;
        indices(3) = (km)*Nx_H + 1;
    end
    if km == Nx_H
        indices(3) = kn+1; if kn == Nx_H indices(3) = 1; end
        indices(4) = kn;
    end
end

function indices = index_map(km,kn,Nx_H)
    indices(1) = (km-1)*(Nx_H+1)+kn;
    indices(2) = (km-1)*(Nx_H+1)+kn+1;
    indices(3) = (km)*(Nx_H+1)+kn+1;
    indices(4) = (km)*(Nx_H+1)+kn;
end