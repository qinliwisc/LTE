function Reference_compatible_initial(sigma_value,epsilon_value)
% %=============
% reference solution
% linear transport
% full domain, uncorrected boundary
% %=============


if nargin == 0
    clear;
    sigma_value = 400; epsilon_value = 1/400;
end
% discretization x
% dx = 0.0001;
dx = 0.001/floor(0.001/(epsilon_value/20)+1);
dx = min(dx,5e-4);
x = -1+dx/2:dx:1-dx/2; Nx = length(x);

% discretization v, polynomials
Np = 16; N_modes = 8;
[Vp, Wp] = half_legendre_quad(Np-1);
[Vp,ind] = sort(Vp,'descend'); Wp = Wp(ind)/2;
[Hp] = half_legendre_poly(Vp,[0:1:N_modes-1]);

Vm = -Vp; Wm = Wp;
[Vm,ind] = sort(Vm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

V_full = [Vp; Vm]; Weight = [Wp; Wm];

H_poly = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

Np = 2*Np;

Lambda = 1 + V_full*V_full'/2;


% evolution set up
sigma = sigma_value*(x<=0)+sigma_value*(x>0);
epsilon = epsilon_value;
% epsilon = 1e-2*ones(1,Nx);
sigma_all = ones(Np,1)*sigma;

Time = 0.003; dt = min(dx*epsilon/2/max(V_full)); t = 0;

% initial data
phi = zeros(Np,Nx);

% boundary data
fboundary_l.fun =  inline('v-v+1+t-t','v','t');
fboundary_r.fun =  inline('v-v+1+t-t','v','t');
boundary_l = feval(fboundary_l.fun,Vp,0);
boundary_r = feval(fboundary_r.fun,Vm,0);

phi(Np/2+1:Np,end) = boundary_r;
phi(1:Np/2,1) = boundary_l;

while t<Time
    %=== transport term ====
    phi_n = phi;
    for k = 1:Np/2 % v> 0
        U = phi(k,:); U = [U,U(end)];
        Flux = VanLeerFlux(U,dt,dx,V_full(k));
        phi_n(k,2:end) = phi(k,2:end) + (- Flux(2:end) + Flux(1:end-1))/epsilon;%./epsilon(2:end-1);
    end
    
    for k = 1:Np/2 % v < 0
        U = phi(k+Np/2,:); U = [U(1),U];
        Flux = VanLeerFlux(U,dt,dx,V_full(k+Np/2));
        phi_n(k+Np/2,1:end-1) = phi(k+Np/2,1:end-1) + (- Flux(2:end) + Flux(1:end-1))/epsilon;%./epsilon(2:end-1);
    end
    
    %=== collision term ====
    collision = Lambda*diag(Weight)*phi_n;
    collision = collision - phi_n;% collision = collision/2;
    collision = dt*collision.*sigma_all/epsilon;

    phi_n = phi_n + collision;
    phi = phi_n;

   
    phi(Np/2+1:Np,end) = boundary_r;
    phi(1:Np/2,1) = boundary_l;
    
    theta_all = Weight'*phi;
%     plot(x,theta_all,'.');pause(0.0005);

   t = t+dt
end
filename = ['data/heat_compatible_initial/transport_',num2str(sigma_value)];
save(filename);
return