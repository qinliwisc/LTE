function LTE_2D_theta_ep(ep_value,N_value)
global Nr Nt epsilon dt dr theta0 rr sigma_x

if nargin == 0
    ep_value = 1;
    N_value = 20;
end

Lr = 4; Nr = Lr*N_value; dr = Lr/Nr; r0 = 0:dr:Lr;
Nt = 120; dt = 2*pi/Nt; theta0 = -pi+dt/2:dt:pi-dt/2;
[rr,tt] = meshgrid(r0,theta0);

epsilon = 2^(-ep_value);
sigma = inline('sin(2*pi*x)+1');
sigma_x = sigma(r0);

%% initial 
f = ones(size(rr));
data_pre = zeros(size(f));
for kt = 1:Nt
    sinT = sin(theta0(kt));
    if sinT>=0
        data_pre(kt,1) = 1;
%         data_pre(kt,1) = sinT;
%         data_pre(kt,1) = 0;
    else
        data_pre(kt,end) = 0;
%         data_pre(kt,end) = sinT;
%         data_pre(kt,end) = 0.814;%0.5;%0.814;%0.7104;%0.7832;
    end
end
data_pre = data_pre(:);
% 
% Matrix = MatrixConstruct;
% f = Matrix\data_pre;
f = gmres(@MAB_multi,data_pre,150,1e-10);

f = reshape(f,Nt,Nr+1);


save(['test_ep_',num2str(ep_value),'.mat']);
% save(['Guo_ep/ep',num2str(ep_value),'_v0_R',num2str(Lr),'.mat']);

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
global epsilon Nr Nt dr dt rr theta0 sigma_x
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
    
    sinT = sin(v(kt));
    if sinT>=0
        Rp(kt,2:end) = sin(v(kt))*Rp_pre;
        boundary_D(kt,1) = p(kt,1);
    else
        Rp(kt,1:end-1) = sin(v(kt))*Rp_pre;
        boundary_D(kt,end) = p(kt,end);
    end
end

% cos_v = cos(v)*ones(1,Nr+1);
% Tp(:,2:end-1) = ([p(2:end,2:end-1);p(1,2:end-1)] - [p(end,2:end-1);p(1:end-1,2:end-1)])/2/dt;
% Tp(:,2:end-1) = cos_v(:,2:end-1).*Tp(:,2:end-1);
% Tp = ([p(2:end,:);p(1,:)] - [p(end,:);p(1:end-1,:)])/2/dt;
% Tp = sin_v.*Tp;

BCp(:,2:Nr) = p(:,2:Nr) - e*e'*p(:,2:Nr)/Nt;
BCp = (ones(Nt,1)*sigma_x).*BCp;

% rrr = (epsilon)./(1-epsilon*rr);
% Mb = epsilon*Rp - rrr.*Tp + BCp + boundary_D;
Mb = epsilon*Rp + BCp + boundary_D;

Mb = Mb(:);

end

