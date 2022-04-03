function CheckingConvergencePhi(delta)

if nargin == 0
    delta = 2;
end

% set up
Lx = 1;
Nx = 100; dx_h = 2*Lx/Nx; x0_h = -Lx:dx_h:Lx; [xx,yy] = meshgrid(x0_h,x0_h);

sigma_fun = @(x,y,delta)1./(2-1.8*sin(delta*pi*y))./(2-1.8*sin(delta*pi*x));
sigma_star = 1/2/sqrt(4-1.8^2);

% reference
bilinear_compute = [1,-1,-1,1;1,1,-1,-1;1,1,1,1;1,-1,1,-1];
bilinear_coeff = inv(bilinear_compute);
soln_star = cell(1,4); soln_star_c = cell(1,4);
py_soln_star = cell(1,4); px_soln_star = cell(1,4);
for km = 1:4
    soln_temp = bilinear_coeff(1,km) + bilinear_coeff(2,km)*xx(1:end,1:end)...
        + bilinear_coeff(3,km)*yy(1:end,1:end) ...
        + bilinear_coeff(4,km)*xx(1:end,1:end).*yy(1:end,1:end);
    soln_star{km} = soln_temp;
    
    soln_temp2 = soln_temp(1:end-1,1:end-1) + soln_temp(1:end-1,2:end) + soln_temp(2:end,1:end-1) + soln_temp(2:end,2:end);
    soln_star_c{km} = soln_temp2/4;
    
    py_soln_star{km} = (soln_temp(2:end,2:end)+soln_temp(2:end,1:end-1)...
        - soln_temp(1:end-1,2:end) - soln_temp(1:end-1,1:end-1))/dx_h/2;
    
    px_soln_star{km} = (soln_temp(1:end-1,2:end) + soln_temp(2:end,2:end)...
            - soln_temp(1:end-1,1:end-1) - soln_temp(2:end,1:end-1))/dx_h/2;
end

% stiffness matrix construction
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
soln = cell(1,4); soln_c = cell(1,4); py_soln = cell(1,4);
error_soln = cell(1,4);
for kn = 1:4
    soln_temp = reshape(Load(:,kn),Nx+1,Nx+1);
    soln{kn} = soln_temp;
    
    soln_temp2 = soln_temp(1:end-1,1:end-1) + soln_temp(1:end-1,2:end) + soln_temp(2:end,1:end-1) + soln_temp(2:end,2:end);
    soln_c{kn} = soln_temp2/4;
    
    py_soln{kn} = (soln_temp(2:end,2:end)+soln_temp(2:end,1:end-1)...
        - soln_temp(1:end-1,2:end) - soln_temp(1:end-1,1:end-1))/dx_h/2;

    px_soln{kn} = (soln_temp(1:end-1,2:end) + soln_temp(2:end,2:end)...
            - soln_temp(1:end-1,1:end-1) - soln_temp(2:end,1:end-1))/dx_h/2;
        
    error_soln{kn} = soln{kn} - soln_star{kn};
    error_norm(kn) = norm(reshape(error_soln{kn},(Nx+1)^2,1),2);
    
    error_py_soln{kn} = py_soln{kn} - py_soln_star{kn};
    error_px_soln{kn} = px_soln{kn} - px_soln_star{kn};
    error_py_norm(kn) = norm(reshape(error_py_soln{kn},Nx^2,1),2)^2 + ...
        norm(reshape(error_px_soln{kn},Nx^2,1),2)^2;
    error_py_norm(kn) = sqrt(error_py_norm(kn));
end

% reconstruct macro matrices
Sigma_star_mat = zeros(4,4); Sigma_mat = zeros(4,4);
for km = 1:4
    for kn = 1:4
        Phi_starSigmaPhi_star = soln_star_c{km}.*soln_star_c{kn}*sigma_star;
        Sigma_star_mat(km,kn) = sum(sum(Phi_starSigmaPhi_star))*dx_h*dx_h;

        PhiSigmaPhi = soln_c{km}.*soln_c{kn}.*sigma_c;
        Sigma_mat(km,kn) = sum(sum(PhiSigmaPhi))*dx_h*dx_h;
        
        
        pyPhiSigmaPhi_star = py_soln_star{km}.*soln_star_c{kn}*sigma_star;
        Sigma_y_star_mat(kn,km) = sum(sum(pyPhiSigmaPhi_star))*dx_h*dx_h;
        
        pyPhiSigmaPhi = py_soln{km}.*sigma_c.*soln_c{kn};
        Sigma_y_mat(kn,km) = sum(sum(pyPhiSigmaPhi))*dx_h*dx_h;
        
        
        pyPhiPhi_star = py_soln_star{km}.*soln_star_c{kn};
        Phi_y_star_mat(kn,km) = sum(sum(pyPhiPhi_star))*dx_h*dx_h;
        
        pyPhiPhi = py_soln{km}.*sigma_c;
        Phi_y_mat(kn,km) = sum(sum(pyPhiPhi))*dx_h*dx_h;
    end
end

save(['2D/delta_conv/delta_',num2str(delta),'.mat'],'error_py_soln','error_py_norm','Phi_y_star_mat','Phi_y_mat','Sigma_y_star_mat','Sigma_y_mat','Sigma_star_mat','Sigma_mat','soln','soln_star','error_soln','error_norm');

end

function indices = index_map(km,kn,Nx_H)
    indices(1) = (km-1)*(Nx_H+1)+kn;
    indices(2) = (km-1)*(Nx_H+1)+kn+1;
    indices(3) = (km)*(Nx_H+1)+kn+1;
    indices(4) = (km)*(Nx_H+1)+kn;
end