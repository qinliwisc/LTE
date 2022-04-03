
% This is a fully implicit scheme for even-odd formulation
% use periodic boundary condition

function implicit_asymmetric_vector
global Nx Nv epsilon v0 dt dx

Lx = 1;
Nx = 50; 
dx = 1/Nx;
x0 = dx/2:dx:Lx-dx/2;
Nv = 10;
dv = 2/Nv;
v0 = -1+dv/2:dv:1-dv/2;
[xx,vv] = meshgrid(x0,v0);
dt      = dx/3;
tFinal  = 0.5;

epsilon = 0.001;

e = ones(Nv,1);

%% initial 
%f_ini = (1e-9)^(1/4)*ones(1,Nx); f_ini((abs(x0-0.5)<0.2)) = 2;
f_ini = 1/2*(tanh(20*(x0-Lx/2+0.2))+tanh((Lx/2+0.2-x0)*20));
f_ini = ones(Nv,1)*f_ini;

f = f_ini;
%% evolution
t = 0;
while t < tFinal
    data_pre = f;
    data_pre = data_pre*epsilon^2/dt;
 
    ee_data = e*e'*data_pre/Nv;
    data_pre = dt*ee_data/epsilon^2 + (data_pre-ee_data)/(1+epsilon^2/dt);
    
	data_pre = data_pre(:);
    
    f = gmres(@MAB_multi,data_pre);
    f = reshape(f,Nv,Nx);
        
    figure(1)
    mesh(vv,xx,f); title('f');
    %pause;
    
    t = t + dt;
end

save('ep3.mat');
end

function Mb = MAB_multi(p)
global epsilon Nv Nx v0 dt dx
e = ones(Nv,1);
v = v0(:);

p = reshape(p,Nv,Nx);

Ap = ([p(:,2:end),p(:,1)] - [p(:,end),p(:,1:end-1)])/2/dx;
VV = v*ones(1,Nx);
Ap = epsilon*VV.*Ap;

eeAp = e*e'*Ap/Nv;
Mb = dt*eeAp/epsilon^2 + (Ap-eeAp)/(1+epsilon^2/dt);

Mb = Mb + p;

Mb = Mb(:);
end

