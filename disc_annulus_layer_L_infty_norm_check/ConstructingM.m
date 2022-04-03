clear;

Nr = 50; Nt = 12;

N_total = (Nr+1)*Nt;

Mb = zeros(N_total,N_total);

PP = eye(N_total);

for km = 1:N_total
    km
    b_temp = AB_multi_Guo(PP(:,km),Nr,Nt);
    Mb(:,km) = b_temp;
end

[V,D] = eig(Mb);

D = diag(D); D_real = real(D); D_real = sort(D_real);