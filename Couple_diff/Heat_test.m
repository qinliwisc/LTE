function Heat_test
%%=============
% heat equation, kinetic boundary
%%=============

% discretization x, odd number of grid points
% interface takes one grid
% dx = 0.01; x = -1:dx:1; Nx = length(x);
dx = 0.00005; x = -1+dx/2:dx:1-dx/2; Nx = length(x);
theta = x-x;

% discretization v, polynomials
Np = 16;
[Vp, Wp] = half_legendre_quad(Np-1);
[Vp,ind] = sort(Vp,'descend'); Wp = Wp(ind);% Wp = Wp(ind)/2;
[Hp] = half_legendre_poly(Vp,[0:1:Np-1]); Hp = Hp/2;

Vm = -Vp; Wm = Wp;
[Vm,ind] = sort(Vm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

V = [Vp; Vm]; W = [Wp; Wm];

H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

Np = 2*Np;

% time
Time = 0.002; dt = min(dx^2/2); t = 0;

% boundary data, kinetic
% fboundary.fun =  inline('v.^2');%-v+1');
% fboundary.fun =  inline('abs(v)');%-v+1');
fboundary.fun =  inline('100*abs(v)*t+1.5','v','t');
boundary_l = feval(fboundary.fun,Vp,0);
boundary_r = feval(fboundary.fun,Vm,0);
% boundary_r = boundary_r - boundary_r;

% === left boundary ========
f_outgoing = boundary_l;
[boundary_l_adjust,eta_l] = boundary_computing(V,W,H,f_outgoing);
% [boundary_l_adjust,eta_l] = boundary_computing(V,W,H,fboundary);
% === right boundary ==========
f_outgoing = boundary_r;
f_outgoing = flipud(f_outgoing);
[boundary_r_adjust,eta_r] = boundary_computing(V,W,H,f_outgoing);
boundary_r_adjust = flipud(boundary_r_adjust);

% === left boundary ========
f_outgoing = feval(fboundary.fun,Vp,Time);
[boundary_l_adjust,eta_l_f] = boundary_computing(V,W,H,f_outgoing);
% [boundary_l_adjust,eta_l] = boundary_computing(V,W,H,fboundary);
% === right boundary ==========
f_outgoing = feval(fboundary.fun,Vm,Time);
f_outgoing = flipud(f_outgoing);
[boundary_r_adjust,eta_r_f] = boundary_computing(V,W,H,f_outgoing);

d_eta_l = (eta_l_f - eta_l)*dt/Time;
d_eta_r = (eta_r_f - eta_r)*dt/Time;

% initial data, kinetic
% initial data
% phi_full = (V.^2)*ones(1,Nx);
% phi_full = (abs(V))*ones(1,Nx);
phi_full = ones(Np,1)*(sin(x*pi)+1.5);
% phi_full(:,(Nx-1)/2+1:Nx) = phi_full(:,(Nx-1)/2+1:Nx)+ones(Np,1)*x((Nx-1)/2+1:end);

% initial data, fluid
theta = W'*phi_full/2; theta = theta(:);
theta(1) = eta_l; theta(end) = eta_r;

% laplace operator
e = ones(Nx,1);
A = spdiags([e -2*e e],[-1:1],Nx-2,Nx-2);
A = [zeros(Nx-2,1),A,zeros(Nx-2,1)];
A(1,1) = 1;
A(end,end) = 1;

while t<Time
    
    theta(2:end-1) = theta(2:end-1) + dt/dx^2*A*theta/3;

    eta_l = eta_l + d_eta_l;
    eta_r = eta_r + d_eta_r;
	
    theta(1) =  eta_l;
    theta(end) = eta_r;

%     plot(x, theta,'.-.');pause(0.005);
    
    t = t+dt
end
save('data/heat_resolve_5n5_test');
return