function CheckingConvergence(delta_v)

% set up
Lx = 1;
Nx_H = 1; dx_H = 2*Lx/Nx_H; x0_H = -Lx:dx_H:Lx; [x0,y0] = meshgrid(x0_H,x0_H);
Nx_h = 100; dx_h = dx_H/Nx_h; x0_h = 0:dx_h:Lx; [x0h,y0h] = meshgrid(x0_h,x0_h);

% sigma_fun = @(x,y)cos(10*pi*x).*cos(4*pi*y)+2;
% sigma_fun = @(x,y)(2+1.8*sin(10*pi*x))./(2+1.8*cos(10*pi*y))+(2+sin(10*pi*y))./(2+1.8*sin(10*pi*x));
% sigma_fun = @(x,y,delta)cos(delta*pi*x).*cos(delta*pi*y)+2;
sigma_fun = @(x,y,delta)1./(2+1.8*sin(delta*pi*y))./(2+1.8*sin(delta*pi*x));
sigma_fun_star = @(x,y) 1/2/sqrt(4-1.8^2)-x+x-y+y;

% stiffness matrix construction
Phi_mat = zeros((Nx_H)^2,(Nx_H)^2);
Sigma_mat = zeros((Nx_H)^2,(Nx_H)^2);
Sigma_star_mat = zeros((Nx_H)^2,(Nx_H)^2);
Sigma_x_mat = zeros((Nx_H)^2,(Nx_H)^2);
Sigma_y_mat = zeros((Nx_H)^2,(Nx_H)^2);
Xi_x_mat = zeros((Nx_H)^2,(Nx_H)^2);
Xi_y_mat = zeros((Nx_H)^2,(Nx_H)^2);
for km = 1:1 % for y
    for kn = 1:1 % for x
        [indices] = index_periodic_map(km,kn,Nx_H);
        segment_x = x0_H(kn) + [0:dx_h:dx_H];
        segment_y = x0_H(km) + [0:dx_h:dx_H];
        [Sigma_mat_temp,Sigma_star_mat_temp,Sigma_x_mat_temp,Sigma_y_mat_temp,Xi_x_mat_temp,Xi_y_mat_temp,Phi_mat_temp,~] = BasisConstruct(sigma_fun,sigma_fun_star,segment_x,segment_y,delta_v);
        Sigma_mat(indices,indices) = Sigma_mat(indices,indices) + Sigma_mat_temp;
        Sigma_star_mat(indices,indices) = Sigma_star_mat(indices,indices) + Sigma_star_mat_temp;
        Phi_mat(indices,indices) = Phi_mat(indices,indices) + Phi_mat_temp;
        Sigma_x_mat(indices,indices) = Sigma_x_mat(indices,indices) + Sigma_x_mat_temp;
        Sigma_y_mat(indices,indices) = Sigma_y_mat(indices,indices) + Sigma_y_mat_temp;
        Xi_x_mat(indices,indices) = Xi_x_mat(indices,indices) + Xi_x_mat_temp;
        Xi_y_mat(indices,indices) = Xi_y_mat(indices,indices) + Xi_y_mat_temp;
    end
end

save(['2D/delta_conv/known_A_star/',num2str(delta_v)],'Phi_mat_temp','Sigma_mat_temp','Sigma_star_mat_temp','Sigma_x_mat_temp','Sigma_y_mat_temp','Xi_x_mat_temp','Xi_y_mat_temp');
end

function [Sigma_mat,Sigma_star_mat,Sigma_x_mat,Sigma_y_mat,Xi_x_mat,Xi_y_mat,Phi_mat,Load] = BasisConstruct(sigma_fun,sigma_fun_star,segment_x,segment_y,delta)
    Nx = length(segment_x)-1; dx_h = segment_x(2)-segment_x(1);
    
    [xx,yy] = meshgrid(segment_x,segment_y);
    
    bilinear_compute = [1,-1,-1,1;1,1,-1,-1;1,1,1,1;1,-1,1,-1];
    bilinear_coeff = inv(bilinear_compute);
    soln_star = cell(1,4);
    for km = 1:4
        soln_temp = bilinear_coeff(1,km) + bilinear_coeff(2,km)*xx(1:end-1,1:end-1)...
            + bilinear_coeff(3,km)*yy(1:end-1,1:end-1) +...
            bilinear_coeff(4,km)*xx(1:end-1,1:end-1).*yy(1:end-1,1:end-1);
        soln_star{km} = soln_temp;
    end

    StiffA = zeros((Nx+1)^2,(Nx+1)^2); Load = zeros((Nx+1)^2,4);
    
    stiff_local_l = [1/2,0,0,-1/2;0,0,0,0;0,0,1/2,-1/2;-1/2,0,-1/2,1];
    stiff_local_r = [1/2,-1/2,0,0;-1/2,1,-1/2,0;0,-1/2,1/2,0;0,0,0,0];
    linear_increase = [0:Nx]/Nx; linear_decrease = -[0:Nx]/Nx+1;
    
    index_matrix = [1:(Nx+1)^2]; index_matrix = reshape(index_matrix,Nx+1,Nx+1);
    edges = [index_matrix(1,:)',index_matrix(:,end),index_matrix(end,:)',index_matrix(:,1)];
    remainder = index_matrix(2:end-1,2:end-1); remainder = remainder(:);
    edges_total = [edges(:,1);edges(2:end,2);edges(1:end-1,3);edges(2:end-1,4)];
    
    sigma_l = sigma_fun(xx(1:end-1,1:end-1) + dx_h/3,yy(1:end-1,1:end-1) + 2*dx_h/3,delta);
    sigma_r = sigma_fun(xx(1:end-1,1:end-1) + 2*dx_h/3,yy(1:end-1,1:end-1) + dx_h/3,delta);
    sigma_c = sigma_fun(xx(1:end-1,1:end-1) + dx_h/2,yy(1:end-1,1:end-1) + dx_h/2,delta);
    sigma_c_star = sigma_fun_star(xx(1:end-1,1:end-1) + dx_h/2,yy(1:end-1,1:end-1) + dx_h/2);
    
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
    soln = cell(1,4); Sigma_star_mat = zeros(4,4); 
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
            Phi_starSigmaPhi_star = soln_star{km}.*soln_star{kn}.*sigma_c_star;
            Sigma_star_mat(km,kn) = sum(sum(Phi_starSigmaPhi_star))*dx_h*dx_h;
            
            PhiSigmaPhi = soln_c{km}.*soln_c{kn}.*sigma_c;
            Sigma_mat(km,kn) = sum(sum(PhiSigmaPhi))*dx_h*dx_h;
            
            PhiPhi = soln_c{km}.*soln_c{kn};
            Phi_mat(km,kn) = sum(sum(PhiPhi))*dx_h*dx_h;
            
            pxPhiSigmaPhi = px_soln{km}.*sigma_c.*soln_c{kn};
            pyPhiSigmaPhi = py_soln{km}.*sigma_c.*soln_c{kn};
            Sigma_x_mat(kn,km) = sum(sum(pxPhiSigmaPhi))*dx_h*dx_h;
            Sigma_y_mat(kn,km) = sum(sum(pyPhiSigmaPhi))*dx_h*dx_h;
            
            pxPhiPhi = px_soln{km}.*soln_c{kn};
            pyPhiPhi = py_soln{km}.*soln_c{kn};
            Xi_x_mat(kn,km) = sum(sum(pxPhiPhi))*dx_h*dx_h;
            Xi_y_mat(kn,km) = sum(sum(pyPhiPhi))*dx_h*dx_h;
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