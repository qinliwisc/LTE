function Transport_2D_periodic_symm
global Flux_A_x Flux_A_y Sigma_x_mat Sigma_y_mat Partial_x_mat Partial_y_mat 
global Sigma_mat Sigma_inv_mat Phi_mat epsilon dt Nx_H Np

% set up
Lx = 1;
Nx_H = 50; dx_H = 2*Lx/Nx_H; x0_H = -Lx:dx_H:Lx; [x0,y0] = meshgrid(x0_H,x0_H);
Nx_h = 20; dx_h = dx_H/Nx_h; x0_h = 0:dx_h:Lx; [x0h,y0h] = meshgrid(x0_h,x0_h);
Nv = 20; [v0,W0] = legendre_quad(Nv-1); Np = 4; [poly_l] = legendre_poly(v0,[0:1:Np-1]);
weight_adjust = ones(1,Nv)*diag(W0)*poly_l(:,1)/poly_l(1,1);
poly_l = sqrt(weight_adjust)*poly_l; W0 = W0/weight_adjust;


% sigma_fun = @(x,y)cos(10*pi*x).*cos(4*pi*y)+2;
sigma_fun = @(x,y)(2+1.8*sin(10*pi*x))./(2+1.8*cos(10*pi*y))+(2+sin(10*pi*y))./(2+1.8*sin(10*pi*x));
epsilon_index = 1.5; epsilon = 10^(-epsilon_index);

% flux in v

Flux_A_x = poly_l'*diag(cos(pi*v0).*W0)*poly_l; Flux_A_x(abs(Flux_A_x)<1e-14) = 0;
Flux_A_y = poly_l'*diag(sin(pi*v0).*W0)*poly_l; Flux_A_y(abs(Flux_A_y)<1e-14) = 0;

% stiffness matrix construction
Phi_mat = zeros((Nx_H)^2,(Nx_H)^2);
Sigma_mat = zeros((Nx_H)^2,(Nx_H)^2); Sigma_inv_mat = zeros((Nx_H)^2,(Nx_H)^2);
Sigma_x_mat = zeros((Nx_H)^2,(Nx_H)^2); Sigma_y_mat = zeros((Nx_H)^2,(Nx_H)^2);
Partial_x_mat = zeros((Nx_H)^2,(Nx_H)^2); Partial_y_mat = zeros((Nx_H)^2,(Nx_H)^2);
for km = 1:Nx_H % for y
    for kn = 1:Nx_H % for x
        [indices] = index_periodic_map(km,kn,Nx_H);
        segment_x = x0_H(kn) + [0:dx_h:dx_H];
        segment_y = x0_H(km) + [0:dx_h:dx_H];
        [Sigma_mat_temp,Sigma_inv_mat_temp,Sigma_x_mat_temp,Sigma_y_mat_temp,...
            Partial_x_mat_temp,Partial_y_mat_temp,Phi_mat_temp,~] = BasisConstruct(sigma_fun,segment_x,segment_y);
        Sigma_mat(indices,indices) = Sigma_mat(indices,indices) + Sigma_mat_temp;
        Sigma_inv_mat(indices,indices) = Sigma_inv_mat(indices,indices) + Sigma_inv_mat_temp;
        Phi_mat(indices,indices) = Phi_mat(indices,indices) + Phi_mat_temp;
        Sigma_x_mat(indices,indices) = Sigma_x_mat(indices,indices) + Sigma_x_mat_temp;
        Sigma_y_mat(indices,indices) = Sigma_y_mat(indices,indices) + Sigma_y_mat_temp;
        Partial_x_mat(indices,indices) = Partial_x_mat(indices,indices) + Partial_x_mat_temp;
        Partial_y_mat(indices,indices) = Partial_y_mat(indices,indices) + Partial_y_mat_temp;
    end
end

t = 0; tFinal = 0.1; dt = min(dx_H/2);
% rho_ini = sqrt(2)*cos(pi*x0(1:end-1,1:end-1)/2).*cos(pi*y0(1:end-1,1:end-1)/2);
rho_ini = cos(pi*x0(1:end-1,1:end-1)+dx_H/2).*cos(pi*y0(1:end-1,1:end-1)+dx_H/2)+1.1;
f_total = zeros(Nx_H,Nx_H,2*Nv);
proj_a_E = zeros(Nx_H,Nx_H,Np);proj_a_O = zeros(Nx_H,Nx_H,Np);
f_plane = rho_ini;
for km = 1:Nv
    f_total(:,:,km) = f_plane;
end
for km = 1:Nx_H
    for kn = 1:Nx_H
        f_temp = f_total(km,kn,1:Nv); f_temp = f_temp(:);
        proj_a_tempE = poly_l'*diag(W0)*f_temp;
        f_temp = f_total(km,kn,Nv+1:2*Nv); f_temp = f_temp(:);
        proj_a_tempO = poly_l'*diag(W0)*f_temp;
        for kp = 1:Np
            proj_a_E(km,kn,kp) = proj_a_tempE(kp);
            proj_a_O(km,kn,kp) = proj_a_tempO(kp);
        end
    end
end

% evolution
time_counter = 1;

while t < tFinal
    t = t + dt; time_counter = time_counter + 1;

    data_pre = zeros(Nx_H,Nx_H,2*Np);
    for km = 1:Np
        data_temp = proj_a_E(:,:,km); data_temp = data_temp(:);
        data_temp = Phi_mat*data_temp;
        data_pre(:,:,km) = reshape(data_temp,Nx_H,Nx_H);
        data_temp = proj_a_O(:,:,km); data_temp = data_temp(:);
        data_temp = Sigma_mat*data_temp;
        data_pre(:,:,km+Np) = reshape(data_temp,Nx_H,Nx_H);
    end
    data_pre = data_pre(:);
    proj_a_new = gmres(@MAB_multi,data_pre,50,1e-8);
    proj_a_new = reshape(proj_a_new,Nx_H,Nx_H,2*Np);
    proj_a_new_E = proj_a_new(:,:,1:Np);
    proj_a_new_O = proj_a_new(:,:,Np+1:2*Np);
    
    proj_a_E = proj_a_new_E; rho = proj_a_E(:,:,1);
    proj_a_O = proj_a_new_O; 

    figure(4)
    mesh(rho);
end

filename = ['2D/ex3_transport_symm/transport_',num2str(epsilon_index),'.mat'];
save(filename,'proj_a_E','proj_a_O','rho','tFinal','x0_H','x0_h','Np','Nv','sigma_fun','epsilon','dt','rho_ini');
end

function proj = MAB_multi(p)
%need sigma_inv_mat partial_x_mat partial_y_mat
global Flux_A_x Flux_A_y Sigma_x_mat Sigma_y_mat Sigma_mat Phi_mat epsilon dt Nx_H Np
global Partial_x_mat Partial_y_mat Sigma_inv_mat

p = reshape(p,Nx_H,Nx_H,2*Np);

term1E = zeros(Nx_H,Nx_H,Np); term1O = zeros(Nx_H,Nx_H,Np);
term2_1E = zeros(Nx_H,Nx_H,Np); term2_2E = zeros(Nx_H,Nx_H,Np);
term2_1O = zeros(Nx_H,Nx_H,Np); term2_2O = zeros(Nx_H,Nx_H,Np);
term3E = zeros(Nx_H,Nx_H,Np); term3O = zeros(Nx_H,Nx_H,Np);
for km = 1:Nx_H
    for kn = 1:Nx_H
        data_temp = p(km,kn,1:Np); data_temp = data_temp(:);
        term2_1E(km,kn,:) = Flux_A_x*data_temp;
        term2_2E(km,kn,:) = Flux_A_y*data_temp;
        term3E(km,kn,2:end) = data_temp(2:end);
        data_temp = p(km,kn,Np+1:2*Np); data_temp = data_temp(:);
        term2_1O(km,kn,:) = Flux_A_x*data_temp;
        term2_2O(km,kn,:) = Flux_A_y*data_temp;
        term3O(km,kn,:) = data_temp;
    end
end
for km = 1:Np
    data_temp = p(:,:,km); data_temp = data_temp(:);
    data_temp = Phi_mat*data_temp;
    term1E(:,:,km) = reshape(data_temp,Nx_H,Nx_H);
    data_temp = p(:,:,km+Np); data_temp = data_temp(:);
    data_temp = Sigma_mat*data_temp;
    term1O(:,:,km) = reshape(data_temp,Nx_H,Nx_H);
    
    data_temp = term2_1E(:,:,km); data_temp = data_temp(:);
    data_temp = Sigma_x_mat*data_temp;
    term2_1E(:,:,km) = dt*reshape(data_temp,Nx_H,Nx_H)/epsilon;
    data_temp = term2_2E(:,:,km); data_temp = data_temp(:);
    data_temp = Sigma_y_mat*data_temp;
    term2_2E(:,:,km) = dt*reshape(data_temp,Nx_H,Nx_H)/epsilon;
    
    data_temp = term2_1O(:,:,km); data_temp = data_temp(:);
    data_temp = Partial_x_mat*data_temp;
    term2_1O(:,:,km) = dt*reshape(data_temp,Nx_H,Nx_H)/epsilon;
    data_temp = term2_2O(:,:,km); data_temp = data_temp(:);
    data_temp = Partial_y_mat*data_temp;
    term2_2O(:,:,km) = dt*reshape(data_temp,Nx_H,Nx_H)/epsilon;
    
    data_temp = term3E(:,:,km); data_temp = data_temp(:);
    data_temp = Sigma_inv_mat*data_temp;
    term3E(:,:,km) = dt*reshape(data_temp,Nx_H,Nx_H)/epsilon^2;
    data_temp = term3O(:,:,km); data_temp = data_temp(:);
    data_temp = Phi_mat*data_temp;
    term3O(:,:,km) = dt*reshape(data_temp,Nx_H,Nx_H)/epsilon^2;
end

proj1 = term1E + term2_1O + term2_2O + term3E;
proj1 = proj1(:);
proj2 = term1O + term2_1E + term2_2E + term3O;
proj2 = proj2(:);

proj = [proj1;proj2];
end

function [Sigma_mat,Sigma_inv_mat,Sigma_x_mat,Sigma_y_mat,Partial_x_mat,Partial_y_mat,Phi_mat,Load] = BasisConstruct(sigma_fun,segment_x,segment_y)
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
    soln = cell(1,4);
    Sigma_inv_mat = zeros(4,4); Partial_x_mat = zeros(4,4); Partial_y_mat = zeros(4,4); 
    Sigma_mat = zeros(4,4); Sigma_x_mat = zeros(4,4); Sigma_y_mat = zeros(4,4); Phi_mat = zeros(4,4);
    px_soln = cell(1,4); py_soln = cell(1,4); soln_c = cell(1,4);
    for kn = 1:4
        soln_temp = reshape(Load(:,kn),Nx+1,Nx+1);
        soln{kn} = soln_temp;
        py_soln{kn} = (soln_temp(2:end,2:end)+soln_temp(2:end,1:end-1)...
            - soln_temp(1:end-1,2:end) - soln_temp(1:end-1,1:end-1))/dx_h/2;
        px_soln{kn} = (soln_temp(1:end-1,2:end) + soln_temp(2:end,2:end)...
            - soln_temp(1:end-1,1:end-1) - soln_temp(2:end,1:end-1))/dx_h/2;
        
        soln_temp = soln_temp(1:end-1,1:end-1) + soln_temp(1:end-1,2:end) + soln_temp(2:end,1:end-1) + soln_temp(2:end,2:end);
        soln_c{kn} = soln_temp/4;
    end
    
    for km = 1:4
        for kn = 1:4
            PhiSigmaPhi = soln_c{km}.*soln_c{kn}.*sigma_c;
            Sigma_mat(km,kn) = sum(sum(PhiSigmaPhi))*dx_h*dx_h;
            PhiSigmaPhi = soln_c{km}.*soln_c{kn}./sigma_c;
            Sigma_inv_mat(km,kn) = sum(sum(PhiSigmaPhi))*dx_h*dx_h;
            PhiPhi = soln_c{km}.*soln_c{kn};
            Phi_mat(km,kn) = sum(sum(PhiPhi))*dx_h*dx_h;
            
            pxPhiPhi = px_soln{km}.*sigma_c.*soln_c{kn};
            pyPhiPhi = py_soln{km}.*sigma_c.*soln_c{kn};
            Sigma_x_mat(kn,km) = sum(sum(pxPhiPhi))*dx_h*dx_h;
            Sigma_y_mat(kn,km) = sum(sum(pyPhiPhi))*dx_h*dx_h;
            pxPhiPhi = px_soln{km}.*soln_c{kn};
            pyPhiPhi = py_soln{km}.*soln_c{kn};
            Partial_x_mat(kn,km) = sum(sum(pxPhiPhi))*dx_h*dx_h;
            Partial_y_mat(kn,km) = sum(sum(pyPhiPhi))*dx_h*dx_h;
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