function [V,D,A,B]=GeneralizedEigen(N,N_modes,dampcoef,X,W,H)

if nargin == 0
    clear;
    N = 64; dampcoef = 0.01;
    N_modes = 16;
    [Xp, Wp] = half_legendre_quad(N-1); Wp = Wp/2;
    [Xp,ind] = sort(Xp,'descend');
    [Hp] = half_legendre_poly(Xp,[0:1:N_modes-1]);

    Xm = -Xp; Wm = Wp;
    [Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
    [Hm] = half_legendre_poly(Xm, [0:1:N-1]);
    Hm = Hp; Hm = Hm(ind,:);

    X = [Xp; Xm]; W = [Wp; Wm]; 
    
    H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];
end

% Calculating A
[alpha,beta] = half_legendre_recurrence(N_modes-1);
Ap = diag(alpha) + diag(beta(1:end-1), 1) + diag(beta(1:end-1), -1);

A = [zeros(N_modes-1,N_modes-1),Ap(1:end-1,:)];
A = [A;Ap(:,1:end-1),zeros(N_modes,N_modes)]/2;

% Calculating B
Lambda = ones(length(X),length(X));
Lambda = Lambda + X*X'/2;% Lambda_norm = Lambda * W;
% Lambda_norm = Lambda_norm*ones(1,2*N);
% Lambda = Lambda./Lambda_norm;
% Lambda = Lambda/(W'*Lambda*W);

chi0X = ones(length(X),1);
matX0 = chi0X(:);

vu = X;

matvuX0 = (vu * ones(1, size(matX0, 2))) .* matX0;
matvu_Linv_vuX0 = (vu * ones(1, size(matX0, 2))) .* matvuX0;

LH_coef_vuX0 = matvuX0.'*diag(W)*H;
LH_coef_vu_Linv_vuX0 = matvu_Linv_vuX0.'*diag(W)*H;

LH = Lambda*diag(W)*H - H;
LH = LH - dampcoef * matvuX0 * LH_coef_vuX0;
LH = LH - dampcoef * matvu_Linv_vuX0 * LH_coef_vu_Linv_vuX0;

B = H.'*diag(W)*LH;
% B = (B + B.')/2;

[V,D]=eig(A',B');
D = diag(D);

% fprintf('Number of negative eigenvalue = %ld (should equal to %ld)\n', ...
%         sum(D < -eps), N-1);
% 
% % number of negative eigenvalues equals to number of even basis
% % functions
% 
% fprintf('Number of zero eigenvalue = %ld (should equal to 1)\n', ...
%         sum(abs(D) < 1e-13));

return