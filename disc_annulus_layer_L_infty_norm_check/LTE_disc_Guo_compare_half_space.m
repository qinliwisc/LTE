function LTE_disc_Guo_compare_half_space(ep_value,N_value,origin_dis)
global Nr Nt epsilon dt dr theta0 rr

if nargin == 0
    ep_value = 15;
    N_value = 10;
    origin_dis = 0;
end

Lr = 15; Nr = Lr*N_value; dr = Lr/Nr; r0 = 0:dr:Lr;
Nt = 240; dt = 2/Nt; theta0 = -1+dt/2:dt:1-dt/2;
[rr,tt] = meshgrid(r0,theta0);

epsilon = 0;

%% initial 
f = ones(size(rr));
data_pre = zeros(size(f));
for kt = 1:Nt
    if theta0(kt)>=0
%         data_pre(kt,1) = theta0(kt);
        data_pre(kt,1) = 0;
    else
        data_pre(kt,end) = 1;
%         data_pre(kt,end) = 0.814;%0.5;%0.814;%0.7104;%0.7832;
    end
end
data_pre = data_pre(:);

f = gmres(@MAB_multi,data_pre,150,1e-10);

%MAB_multi * f = data_pre

f = reshape(f,Nt,Nr+1);

save(['Guo_ep_compare_half_space/ep',num2str(ep_value),'_01_R',num2str(Lr),'.mat']);

end


function Mb = MAB_multi(p)
global epsilon Nr Nt dr dt rr theta0
e = ones(Nt,1);
v = theta0(:);

p = reshape(p,Nt,Nr+1);
Tp = zeros(size(p));
Rp = zeros(size(p));
BCp = zeros(size(p));
boundary_D = zeros(size(p));

for kt = 1:Nt
    p_temp = p(kt,:);
    
    Rp_pre = (p_temp(2:end) - p_temp(1:end-1))/dr;

    %upwinding
    if v(kt)>=0
        Rp(kt,2:end) = v(kt)*Rp_pre;
        boundary_D(kt,1) = p(kt,1);
    else
        Rp(kt,1:end-1) = v(kt)*Rp_pre;
        boundary_D(kt,end) = p(kt,end);
    end
end


cos_v = cos(v)*ones(1,Nr+1);
Tp(:,2:end-1) = ([p(2:end,2:end-1);p(1,2:end-1)] - [p(end,2:end-1);p(1:end-1,2:end-1)])/2/dt;
Tp(:,2:end-1) = cos_v(:,2:end-1).*Tp(:,2:end-1);
% Tp = ([p(2:end,:);p(1,:)] - [p(end,:);p(1:end-1,:)])/2/dt;
% Tp = sin_v.*Tp;

BCp(:,2:Nr) = p(:,2:Nr) - e*e'*p(:,2:Nr)/Nt;

rrr = (epsilon)./(1-epsilon*rr);
Mb = Rp - rrr.*Tp + BCp + boundary_D;

Mb = Mb(:);

end

