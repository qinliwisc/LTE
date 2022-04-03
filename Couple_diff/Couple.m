function Couple(sigma)
%%=============
% coupled system
%%=============

if nargin == 0
    sigma = 100;
end

% discretization x, odd number of grid points
% interface takes one grid
dx = 0.01; x = -1+dx/2:dx:1-dx/2; Nx = length(x);

% discretization v, polynomials
Np = 16;
[Vp, Wp] = half_legendre_quad(Np-1);
[Vp,ind] = sort(Vp,'descend'); Wp = Wp(ind)/2;
[Hp] = half_legendre_poly(Vp,[0:1:Np-1]); Hp = Hp/2;

Vm = -Vp; Wm = Wp;
[Vm,ind] = sort(Vm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

V_full = [Vp; Vm]; Weight = [Wp; Wm];

H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

Np = 2*Np;

% evolution set up
epsilon = 1/sigma;

Lambda = 1 + V_full*V_full'/2;

Time = 0.001; dt = min(min(epsilon/2,dx*epsilon/2),dx^2/2); t = 0;

% boundary data, kinetic
f_boundary_l.fun =  inline('100*abs(v)*t+1','v','t');
f_boundary_r.fun =  inline('100*abs(v)*t+0.5','v','t');
    
% initial data, kinetic
[phi,theta] = Compatible_ini(dx,Np);
theta = theta*(0.25*cos(pi*x)+0.75);
theta = theta((x>0)); theta = theta(:);

boundary_m_old = phi(:,end);

% laplace operator
e = ones(Nx/2-1,1);
A = spdiags([e -2*e e],[-1:1],Nx/2-2,Nx/2-2);
A = [zeros(Nx/2-2,1),A,zeros(Nx/2-2,1)];
A(1,1) = 1;
A(end,end) = 1;
while t < Time
    %% left domain, kinetic
    % === transport term ====
    boundary_l_old = feval(f_boundary_l.fun,Vp,t);
    phi_n = phi;
    for k = 1:Np/2
        U = phi(k,:); U = [boundary_l_old(k),U,U(end)];
        Flux = VanLeerFlux(U,dt,dx,V_full(k));
        phi_n(k,:) = phi(k,:) + (- Flux(2:end) + Flux(1:end-1))/epsilon;
    end
    
    for k = 1:Np/2        
        U = phi(k+Np/2,:); U = [U(1),U,boundary_m_old(k)];
        Flux = VanLeerFlux(U,dt,dx,V_full(k+Np/2));
        phi_n(k+Np/2,:) = phi(k+Np/2,:) + (- Flux(2:end) + Flux(1:end-1))/epsilon;
    end

    phi = phi_n;
    
    %=== collision term ====
    collision = Lambda*diag(Weight)*phi_n;
    collision = collision - phi_n;    
    collision = dt*collision/epsilon;
    
    phi = phi + collision;

    %% right domain, diffusive
    f_outgoing = phi(1:Np/2,end);
    [boundary_m_new,eta_m_new] = boundary_computing(V_full,Weight,H,f_outgoing);
    boundary_r_new = feval(f_boundary_r.fun,Vm,t);
    f_outgoing = boundary_r_new;
    f_outgoing = flipud(f_outgoing);
    [boundary_r_adjust,eta_r_new] = boundary_computing(V_full,Weight,H,f_outgoing);

    theta(1) =  eta_m_new;    theta(end) = eta_r_new;
    theta(2:end-1) = theta(2:end-1) + dt/dx^2*2*A*theta/5;
    
    %% combine
    theta_all = Weight'*phi;
    theta_all = [theta_all,theta'];

    plot(x, theta_all);pause(0.005);
    
    %% update
    boundary_m_old = boundary_m_new;
    
    t = t+dt
end
filename = ['data/couple_',num2str(sigma)];
save(filename);
return