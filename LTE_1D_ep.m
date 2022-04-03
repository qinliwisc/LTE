function LTE_1D_ep(ep_value,N_value)
global Nr Nv epsilon dt dr v0 rr sigma_x

if nargin == 0
    ep_value = 4;
    N_value = 20;
end

Lr = 4; Nr = Lr*N_value; dr = Lr/Nr; r0 = 0:dr:Lr;
Nv = 120; dt = 2/Nv; v0 = -1+dt/2:dt:1-dt/2;
[rr,vv] = meshgrid(r0,v0);

epsilon = 2^(-ep_value);
sigma = inline('1./(sin(2*pi*x)+1.5)+1');
% sigma = inline('x-x+1');
sigma_x = sigma(r0);

%% initial 
f = ones(size(rr));
data_pre = zeros(size(f));
for kv = 1:Nv
    v = v0(kv);
    if v>=0
        data_pre(kv,1) = 1;
%         data_pre(kv,1) = 0;
%         data_pre(kt,1) = sinT;
%         data_pre(kt,1) = 0;
    else
        data_pre(kv,end) = 0;
%         data_pre(kv,end) = 1;
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
save(['test_sin/test_ep_',num2str(ep_value),'.mat']);

% handle_f = figure(1);
% set(gca,'fontsize',20);
% mesh(r0,theta0,f); title(['\epsilon = 2^{',num2str(-ep_value),'}'],'fontsize',20);
% xlabel('r','fontsize',20);ylabel('\theta','fontsize',20);zlabel('f','fontsize',20);
% print(gcf,'-depsc2',['Guo/ep',num2str(ep_value),'_071_R',num2str(Lr),'.eps']);
% close(handle_f);
% 
% average = mean(f(:,end)); dis = sum((f - average).^2)/Nt; dis = sqrt(dis);
% handle_f = figure(1);
% set(gca,'fontsize',20);
% semilogy(r0,dis,'.-.'); title(['\epsilon = 2^{',num2str(-ep_value),'}, |f - f_\infty|_2'],'fontsize',20);
% xlabel('r','fontsize',20);ylabel('error_2','fontsize',20);
% print(gcf,'-depsc2',['Guo/ep',num2str(ep_value),'_071_R',num2str(Lr),'_error.eps']);
% close(handle_f);
% 
% handle_f = figure(1);
% set(gca,'fontsize',20);
% plot(theta0,f(:,1),'.-.'); title(['\epsilon = 2^{',num2str(-ep_value),'}'],'fontsize',20);
% xlabel('r','fontsize',20);ylabel('f(r=0)','fontsize',20);
% print(gcf,'-depsc2',['Guo/ep',num2str(ep_value),'_071_R',num2str(Lr),'_origin.eps']);
% close(handle_f);
% 
% handle_f = figure(1);
% set(gca,'fontsize',20);
% plot(theta0,f(:,end),'.-.'); title(['\epsilon = 2^{',num2str(-ep_value),'}'],'fontsize',20);
% xlabel('r','fontsize',20);ylabel('f(r=0)','fontsize',20);
% print(gcf,'-depsc2',['Guo/ep',num2str(ep_value),'_071_R',num2str(Lr),'_infinity.eps']);
% close(handle_f);

end

function Mb = MatrixConstruct
global Nr Nt

N_total = (Nr+1)*Nt;

Mb = zeros(N_total,N_total);

PP = eye(N_total);

for km = 1:N_total
    km
    b_temp = AB_multi_Guo(PP(:,km),Nr,Nt);
    Mb(:,km) = b_temp;
end
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
        Rp(kv,2:end) = v(kv)*Rp_pre;
        boundary_D(kv,1) = p(kv,1);
    else
        Rp(kv,1:end-1) = v(kv)*Rp_pre;
        boundary_D(kv,end) = p(kv,end);
    end
end

BCp(:,2:Nr) = p(:,2:Nr) - e*e'*p(:,2:Nr)/Nv;
BCp = (ones(Nv,1)*sigma_x).*BCp;

Mb = epsilon*Rp + BCp + boundary_D;

Mb = Mb(:);

end

