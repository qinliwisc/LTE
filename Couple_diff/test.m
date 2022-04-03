% boundary test
clear;
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

boundary_l = ones(Np,1);

[boundary_l_adjust,eta_l] = boundary_computing(V,W,H,boundary_l(1:Np/2));