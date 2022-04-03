function Heat_compatible_CN_to_play(cases,epsilon,fboundary_l,fboundary_r,phi_x,phi_v,ChangeInTime)
%%=============
% heat equation, kinetic boundary
%%=============

% discretization x, odd number of grid points
dx = 5e-4; x = -1+dx/2:dx:1-dx/2; Nx = length(x);
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
Time = 0.1; dt = min(dx/4); t = 0;

% boundary data, kinetic
boundary_l = feval(fboundary_l.fun,Vp,0);
boundary_r = feval(fboundary_r.fun,Vm,0);

% === extrapolation ========
v_func = inline('v');
f_outgoing = v_func(Vp);
[v_func_adjust,extrapolation] = boundary_computing(V,W,H,f_outgoing);


% === left boundary ========
f_outgoing = boundary_l;
[boundary_l_adjust,eta_l] = boundary_computing(V,W,H,f_outgoing);
% === right boundary ==========
f_outgoing = boundary_r;
f_outgoing = flipud(f_outgoing);
[boundary_r_adjust,eta_r] = boundary_computing(V,W,H,f_outgoing);
eta_l_o = eta_l; eta_r_o = eta_r;

% initial data, kinetic
% initial data
phi_full = zeros(Np,Nx);
phix = feval(phi_x.fun,x); phix = phix(:)';
phiv = feval(phi_v.fun,V,eta_l); phiv = phiv(:);
phi_full = phiv*phix;

% initial data, fluid
theta = W'*phi_full; theta = theta(:);
theta(1) = eta_l; theta(end) = eta_r;

% laplace operator
alpha = 2/5;
beta = -1.2*epsilon*extrapolation;
e = ones(Nx,1);
A = alpha*dt/2/dx^2*spdiags([e -2*e e],[-1:1],Nx,Nx);
A(1,1) = A(1,1) + alpha*dt/beta/dx;
A(1,2) = A(1,2) + alpha*dt/2/dx^2;
A(end,end) = A(1,1); A(end,end-1) = A(1,2);

t_rec = [0.001:0.001:0.4];
theta_rec = zeros(length(t_rec),length(x));
km = 1;

while t<Time
    t = t + dt
    
    if ChangeInTime == 1
        eta_l_o = eta_l;
        eta_r_o = eta_r;
        % boundary data, kinetic
        boundary_l = feval(fboundary_l.fun,Vp,t);
        boundary_r = feval(fboundary_r.fun,Vm,t);

        % === left boundary ========
        f_outgoing = boundary_l;
        [boundary_l_adjust,eta_l] = boundary_computing(V,W,H,f_outgoing);
        % === right boundary ==========
        f_outgoing = boundary_r;
        f_outgoing = flipud(f_outgoing);
        [boundary_r_adjust,eta_r] = boundary_computing(V,W,H,f_outgoing);
    end
    
    load = (eye(Nx)+A)*theta;
    load(1) = load(1) - alpha*dt/dx/beta*(eta_l_o+eta_l);
    load(end) = load(end) - alpha*dt/dx/beta*(eta_r_o+eta_r);

    theta = (eye(Nx)-A)\load;
    

    if (abs(t - t_rec(km))<dt)
        theta_rec(km,:) = theta;
        km = km+1;
    end
%     plot(x, theta,'.-.');pause;
    
end
save(['data/heat_compatible_to_play/case',num2str(cases),'/heat_CN_5n4_',num2str(1/epsilon)]);
return