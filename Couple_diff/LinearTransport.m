function LinearTransport
%%=============
% linear transport equation, left domain
%%=============

coupling = 1;

% discretization x
dx = 0.01; x = -1:dx:0; Nx = 1/dx+1;

% discretization v, polynomials
Np = 16;
[Vp, Wp] = half_legendre_quad(Np-1);
[Vp,ind] = sort(Vp,'descend'); Wp = Wp(ind)/2;
[Hp] = half_legendre_poly(Vp,[0:1:Np-1]); Hp = Hp/2;

Vm = -Vp; Wm = Wp;
[Vm,ind] = sort(Vm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

V = [Vp; Vm]; W = [Wp; Wm];

H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

Np = 2*Np;

% evolution set up
epsilon = 1;

Time = 30; dt = min(epsilon/2,0.005); t = 0;

% initial data
phi = V*ones(1,Nx);

% boundary data
fboundary.fun =  inline('v');%-v+1');
boundary_l = feval(fboundary.fun,Vp);

while t<Time
    %=== boundary data =====
    f_outgoing = phi(1:Np/2,end);
    if (coupling == 1)
        boundary_r = boundary_computing(V,W,H,f_outgoing);
    else
        fboundary.fun =  inline('v-v');
        boundary_r = feval(fboundary.fun,Vm);
    end
    
    %=== transport term ====
    Flux_all=zeros(Np,Nx-1);
    phi_n = phi;
    for k = 1:Np/2
        U = phi(k,:);
        Flux = VanLeerFlux(U,dt,dx,V(k)); Flux_all(k,:) = Flux;
        phi_n(k,2:end-1) = phi(k,2:end-1) - Flux(2:end) + Flux(1:end-1);
        phi_n(k,1) = boundary_l(k);
%         phi_n(k,end) = phi_n(k,end) + Flux(end);
    end
    
    for k = 1:Np/2        
        U = phi(k+Np/2,:);
        Flux = VanLeerFlux(U,dt,dx,V(k+Np/2)); Flux_all(k+Np/2,:) = Flux;
        phi_n(k+Np/2,2:end-1) = phi(k+Np/2,2:end-1) - Flux(2:end) + Flux(1:end-1);
        phi_n(k+Np/2,end) = boundary_r(k);
%         phi_n(k+Np/2,1) = phi_n(k+Np/2,1) - Flux(1);
    end
    phi = phi_n;
    
    mesh(x(2:end-1),V,phi(:, 2:end-1));pause(0.01);

    %=== collision term ====
    collision = W'*phi;
    collision = ones(Np,1)*collision - phi;
    phi = phi + dt*collision/epsilon;
    
    t = t+dt;
end
return