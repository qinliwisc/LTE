function Mb = AB_multi(p,Nr,Na,Nt)

if nargin == 0
    kkk = load('ep0_d_2.mat');

    Nr = 20; Nt = 12; Na = 24;

    if nargin == 0
        p = kkk.f;
        p = p(:);
    end
end

epsilon = 0.01;

Lr = 1; dr = Lr/Nr; r0 = 0:dr:Lr; r0 = r0+dr;
dt = 2*pi/Nt; theta0 = -pi+dt/2:dt:pi-dt/2;
da = 2*pi/Na; alpha0 = da/2:da:2*pi-da/2;
[rr,aa,tt] = meshgrid(r0,alpha0,theta0);

e = ones(Nt,1);
v = theta0(:);

p = reshape(p,Na,Nr+1,Nt);
Ap = zeros(size(p));
Rp = zeros(size(p));
BCp = zeros(size(p));
boundary_D = zeros(size(p));
boundary_C = zeros(size(p));

for kt = 1:Nt
    p_temp = p(:,:,kt);
    Ap(:,:,kt) = ([p_temp(2:end,:);p_temp(1,:)] - [p_temp(end,:);p_temp(1:end-1,:)])/2/da;
    Ap(:,:,kt) = sin(v(kt))*Ap(:,:,kt);
    
    cosT = cos(v(kt));
    if cosT<0
        Rp_pre = (p_temp(:,2:end) - p_temp(:,1:end-1))/dr;
        Rp(:,1:end-1,kt) = cos(v(kt))*Rp_pre;
        boundary_D(:,end,kt) = p(:,end,kt);
    else
        Rp_pre = (p_temp(:,2:end) - p_temp(:,1:end-1))/dr;
        boundary_D(:,1,kt) = p(:,1,kt);
        Rp(:,2:end,kt) = cos(v(kt))*Rp_pre;
        if v(kt)<0
            boundary_C(1:Na,1,kt) = p(1:Na,1,kt) - p(end:-1:1,1,Nt/2+1-kt);
        else
            boundary_C(1:Na,1,kt) = p(1:Na,1,kt) - p(end:-1:1,1,3*Nt/2+1-kt);
        end
    end
end

for kr = 1:Nr
    for ka = 1:Na
        p_temp = p(ka,kr,:); p_temp = p_temp(:);
        BCp(ka,kr,:) = p_temp - e*e'*p_temp/Nt;
    end
end

Mb = Ap./rr + Rp + BCp/epsilon + boundary_D + boundary_C;

Mb = Mb(:);
% 
% for kr = 1:Nr+1
%     f_r0 = Mb(:,kr,:);
%     f_r0 = reshape(f_r0,Na,Nt);
%     mesh(theta0,alpha0,f_r0); title([num2str(r0(kr))]);
%     xlabel('\theta'); ylabel('\alpha');
%     pause;
% end
