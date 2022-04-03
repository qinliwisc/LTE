function [f_incoming,eta] = boundary_computing(X,W,H,f_outgoing)

if (nargin == 0)
    clear;
    N = 16;

    [Xp, Wp] = half_legendre_quad(N-1);
    [Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind); Wp = Wp/2;
    Hp = half_legendre_poly(Xp,[0:1:N-1]);% Hp = Hp/2;% Hp = Hp(ind,:)/2;

    Xm = -Xp; Wm = Wp;
    [Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
    Hm = Hp; Hm = Hm([end:-1:1],:);

    X = [Xp; Xm]; 
    W = [Wp; Wm]; 

    H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

    f_outgoing.fun = inline('v');
end

% recovry
chi_data = ones(length(X),1); N = length(X);

[fb_c] = BoundaryData_c_computing(X,W,H,0,f_outgoing); fb = H*fb_c;% plot(X,fb);pause;
[gb_c] = BoundaryData_c_computing(X,W,H,1,f_outgoing); gb = H*gb_c;% plot(X,gb);pause;

C = X'*diag(W)*gb;
D = X'*diag(W)*fb;
eta = D/C;

g_tilde = gb*eta;
phi = chi_data*eta;

f_incoming = fb - g_tilde + phi;
f_incoming = f_incoming(N/2+1:N);

return