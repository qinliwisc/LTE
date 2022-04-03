function Mb = AB_multi_Guo(p,Nr,Nt)

epsilon = 2^(-8);% Nr = 100; Nt = 60;

Lr = 8; dr = Lr/Nr; r0 = 0:dr:Lr;
dt = 2*pi/Nt; theta0 = -pi+dt/2:dt:pi-dt/2;

[rr,tt] = meshgrid(r0,theta0);

e = ones(Nt,1); v = theta0(:);

p = reshape(p,Nt,Nr+1);
Rp = zeros(size(p));
Tp = zeros(size(p));
BCp = zeros(size(p));
boundary_D = zeros(size(p));

for kt = 1:Nt
    p_temp = p(kt,:);
    
    Rp_pre = (p_temp(2:end) - p_temp(1:end-1))/dr;
    
    sinT = sin(v(kt));
    if sinT>0
        Rp(kt,2:end) = sin(v(kt))*Rp_pre;
        boundary_D(kt,1) = p(kt,1);
    else
        Rp(kt,1:end-1) = sin(v(kt))*Rp_pre;
        boundary_D(kt,end) = p(kt,end);
    end
end

cos_v = cos(v)*ones(1,Nr+1);
Tp(:,2:end-1) = ([p(2:end,2:end-1);p(1,2:end-1)] - [p(end,2:end-1);p(1:end-1,2:end-1)])/2/dt;
Tp(:,2:end-1) = cos_v(:,2:end-1).*Tp(:,2:end-1);
% Tp = ([p(2:end,:);p(1,:)] - [p(end,:);p(1:end-1,:)])/2/dt;
% Tp = cos_v.*Tp;

for kr = 2:Nr
    BCp(:,kr) = p(:,kr) - e*e'*p(:,kr)/Nt;
end

rrr = (epsilon)./(1-epsilon*rr);
% Mb = Rp + BCp + boundary_D;
Mb = Rp - rrr.*Tp + BCp + boundary_D;

Mb = Mb(:);

end
