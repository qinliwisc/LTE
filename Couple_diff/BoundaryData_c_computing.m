function coeff = BoundaryData_c_computing(X,W,H,index,fb_ini_s)

if (nargin == 0)
    clear;
    N = 16;

    [Xp, Wp] = half_legendre_quad(N-1);
    [Xp,ind] = sort(Xp,'descend'); Wp = Wp(ind); Wp = Wp/2;
    Hp = half_legendre_poly(Xp,[0:1:N-1]);% Hp = Hp(ind,:)/2;

    Xm = -Xp; Wm = Wp;
    [Xm,ind] = sort(Xm,'descend'); Wm = Wm(ind);
    Hm = Hp; Hm = Hm([end:-1:1],:);

    X = [Xp; Xm]; 
    W = [Wp; Wm]; 

    H = [Hp(:,1:end-1),Hp;Hm(:,1:end-1),-Hm];

    fb_ini_s.fun = inline('v');
    index = 0;
end

dampcoef = 0.01;

%% initial data
modes = size(H,2); modes = (modes+1)/2;
N = length(X); N = N/2; Hp = H(1:N,modes:end);
Xp = X(1:N); Wp = W(1:N);

fb_ini(:,1) = fb_ini_s;%feval(fb_ini_s.fun,Xp);
fb_ini(:,2) = ones(N,1);

[V,D,Amat,Bmat] = GeneralizedEigen(N,modes,dampcoef,X,W,H);

if (sum(D>-1e-14) ~= modes)
    error('Wrong number of negative eigenvalues!');
end
Vcons = V(:,find(D>-1e-14));

% consistency check for the eignevectors
if (norm(Vcons.' * Amat - diag(D(D>-1e-14)) * Vcons.' * Bmat)>1e-5)
    error('Something wrong in the generalized eigenvalue');
end

A = zeros(2*modes-1, 2*modes-1);
b = zeros(2*modes-1, 1);

% The first N rows are constraints from the growing and zero modes
A(1:modes, :) = Vcons.'*Bmat'; 
b(1:modes, :) = 0;

% The rest rows are the Galerkin condition from the boundary values
A(modes+1:end, :) = Hp(:, 1:modes-1)' * diag(Wp.*Xp) * [Hp(:, 1:modes-1), Hp];
b(modes+1:end, :) = Hp(:, 1:modes-1)' * diag(Wp.*Xp) * fb_ini(:,index+1);
% A(N+1:end, :) = Hp(:, 1:N-1).' * diag(Wp.*Xp) * [Hp(:, 1:N-1), Hp];
% b(N+1:end, :) = Hp(1:N, 1:N-1).' * diag(Wp.*Xp) * fb_ini(:,index+1);

coeff = A \ b;

thresholding = zeros(2*modes-1,1);
thresholding(1:modes-1) = [1:1:modes-1];
thresholding = thresholding/modes;
thresholding = (1+cos(pi*thresholding))/2;

index = ones(length(coeff),1);
index = (abs(coeff)>2.5*1e-3);

coeff = coeff.*thresholding;%index;
return