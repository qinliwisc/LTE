function LTE_1D_oppo_ep(ep_value,N_value)
global Nr Nv epsilon dt dr v0 rr sigma_x

if nargin == 0
    ep_value = 4;
    N_value = 20;
end

Lr = 4; Nr = Lr*N_value; dr = Lr/Nr; r0 = 0:dr:Lr;
Nv = 120; dt = 2/Nv; v0 = -1+dt/2:dt:1-dt/2;
[rr,vv] = meshgrid(r0,v0);

epsilon = 2^(-ep_value);
% sigma = inline('1./(sin(2*pi*x)+1.5)+1');
sigma = inline('x-x+1');
sigma_x = sigma(r0);

%% initial 
f = ones(size(rr));
data_pre = zeros(size(f));
for kv = 1:Nv
    v = v0(kv);
    if v<=0
%         data_pre(kv,1) = 1;
        data_pre(kv,1) = 0;
%         data_pre(kt,1) = sinT;
%         data_pre(kt,1) = 0;
    else
%         data_pre(kv,end) = 0;
        data_pre(kv,end) = 1;
%         data_pre(kt,end) = sinT;
%         data_pre(kt,end) = 0.814;%0.5;%0.814;%0.7104;%0.7832;
    end
end
data_pre = data_pre(:);
% 
% Matrix = MatrixConstruct;
% f = Matrix\data_pre;
f = gmres(@MAB_multi,data_pre,150,1e-10);

f = reshape(f,Nv,Nr+1);

% mesh(f);
save(['test/test_g1_ep_',num2str(ep_value),'.mat']);


end


function Mb = MAB_multi(p)
global epsilon Nr Nv dr v0 sigma_x
e = ones(Nv,1);
v = v0(:);

p = reshape(p,Nv,Nr+1);
Rp = zeros(size(p));
BCp = zeros(size(p));
boundary_D = zeros(size(p));

for kv = 1:Nv
    p_temp = p(kv,:);
    
    Rp_pre = (p_temp(2:end) - p_temp(1:end-1))/dr;
    
    if v(kv)>=0
        Rp(kv,1:end-1) = -v(kv)*Rp_pre;
        boundary_D(kv,end) = p(kv,end);
    else
        Rp(kv,2:end) = -v(kv)*Rp_pre;
        boundary_D(kv,1) = p(kv,1);
    end
end

BCp(:,2:Nr) = p(:,2:Nr) - e*e'*p(:,2:Nr)/Nv;
BCp = (ones(Nv,1)*sigma_x).*BCp;

Mb = epsilon*Rp + BCp + boundary_D;

Mb = Mb(:);

end

