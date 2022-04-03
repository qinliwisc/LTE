function [X,W] = legendre_quad(n) % the method presented in Golub's paper
if nargin==0
    n = 15;
end

[beta] = legendre_recurrence(n);

A = diag(beta(1:end-1), 1) + diag(beta(1:end-1), -1);

[V,X]=eig(A);
X = diag(X);
W = 2*vpa(V(1, :)'.^2);

W = double(W);

return