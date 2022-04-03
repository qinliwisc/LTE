function [F_distribution,theta] = Compatible_ini(dx,Np)

if nargin == 0
    clear;
    dx = 0.01; Np = 16;
end

x = [-1+dx/2:dx:-dx/2]; Nx = length(x);

Np = Np/2;
modes = 4;
[Vp, Wp] = half_legendre_quad(Np-1);
[Vp,ind] = sort(Vp,'descend'); Wp = Wp(ind)/2;
[Hp] = half_legendre_poly(Vp,[0:1:modes-1]); Hp = Hp/2;

Vm = -Vp; Wm = Wp;
[Vm,ind] = sort(Vm,'descend'); Wm = Wm(ind);
Hm = Hp; Hm = Hm([end:-1:1],:);

V_full = [Vp; Vm]; Weight = [Wp; Wm];

H_poly = [Hp(:,1:modes-1),Hp;Hm(:,1:modes-1),-Hm];

Np = 2*Np;

[V,D,Amat,Bmat] = GeneralizedEigen(Np/2,modes,0.01,V_full,Weight,H_poly);

D_inv = 1./D;
index = find((abs(D_inv)>1e-10)&(abs(D_inv)<1e10));
V_nonsingular = V(:,index); D_nonsingular = D(index);

f_boundary.fun = inline('v + 3');
f_boundary_value = feval(f_boundary.fun,V_full);
f_boundary_value(1:Np/2) = 3*ones(Np/2,1);
% f_boundary_value(1:Np) = 3*ones(Np,1);


coeff = H_poly' * diag(Weight) * f_boundary_value;
coeff(find(abs(coeff)<1e-14)) = 0;

x_diff = x(1:Nx) + dx/2;
% decay_speed = D;
decay_speed = 1./D_nonsingular;
% decay_speed(find((decay_speed)>0))=0;
decay = decay_speed*x_diff; decay = exp(-decay);

projection_boundary = V_nonsingular' * Bmat * coeff;
% projection_boundary = V' * Bmat * coeff;
projection_boundary(find(abs(projection_boundary)<1e-14)) = 0;
% projection_boundary(find(decay_speed==0)) = 0;
projection = repmat(projection_boundary,1,Nx);
projection = projection.*decay;

for k=1:Nx
%     coeff_x(:,k) = (V'*Amat) \ projection(:,k);
    coeff_x(:,k) = (V_nonsingular'*Bmat) \ projection(:,k);
end

F_distribution = H_poly*coeff_x;

% h = figure(2);
% set(gca,'fontsize',20);
% mesh(x(1:Nx),V_full,F_distribution);
% xlim([-1 0]);
% print(gcf,'-depsc2', 'data/compatible.eps');

F_distribution = ones(Np,Nx);

[f_incoming,theta] = ...
    boundary_computing(V_full,Weight,H_poly,F_distribution(1:Np/2,end));

return