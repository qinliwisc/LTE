function Heat_compatible_test
%%=============
% heat equation, kinetic boundary
%%=============

% discretization x, odd number of grid points
dx = 1e-3; x = -1+dx/2:dx:1-dx/2; Nx = length(x);
theta = x-x;

% discretization v, polynomials
Np = 16; N_modes = 8;
[Vp, Wp] = half_legendre_quad(Np-1);
[Vp,ind] = sort(Vp,'descend'); Wp = Wp(ind)/2;% Wp = Wp(ind)/2;
[Hp] = half_legendre_poly(Vp,[0:1:N_modes-1]);% Hp = Hp/2;

Vm = -Vp; Wm = Wp;
[Vm,ind] = sort(Vm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

V = [Vp; Vm]; W = [Wp; Wm];

H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

Np = 2*Np;

% time
Time = 0.001; dt = min(dx^2/2); t = 0;

% boundary data, kinetic
fboundary_l.fun =  inline('abs(v)+t-t','v','t');
fboundary_r.fun =  inline('abs(v)+t-t','v','t');
boundary_l = feval(fboundary_l.fun,Vp,0);
boundary_r = feval(fboundary_r.fun,Vm,0);

fboundary_fund.fun =  inline('abs(v)+t-t','v','t');
boundary_fund = feval(fboundary_fund.fun,Vp,0);
[boundary_fund_adjust,eta_fund] = boundary_computing(V,W,H,boundary_fund);


% === left boundary ========
f_outgoing = boundary_l;
[boundary_l_adjust,eta_l] = boundary_computing(V,W,H,f_outgoing);
% === right boundary ==========
f_outgoing = boundary_r;
f_outgoing = flipud(f_outgoing);
[boundary_r_adjust,eta_r] = boundary_computing(V,W,H,f_outgoing);

% initial data, kinetic
% initial data
phi_full = eta_fund*ones(Np,Nx);
% initial data, fluid
theta = W'*phi_full; theta = theta(:);
theta(1) = eta_l; theta(end) = eta_r;

% laplace operator
e = ones(Nx,1);
A = spdiags([e -2*e e],[-1:1],Nx-2,Nx-2);
A = [zeros(Nx-2,1),A,zeros(Nx-2,1)];
A(1,1) = 1;
A(end,end) = 1;

while t<Time
    
    theta(2:end-1) = theta(2:end-1) + dt/dx^2*2*A*theta/5;
    
    theta(1) =  eta_fund;
    theta(end) = eta_fund;

    plot(x, theta,'.-.');pause(0.005);
    
    t = t+dt
end
save('data/heat_compatible_test/heat_1n3');
return