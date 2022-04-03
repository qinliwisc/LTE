function LTE_disc_dirichlet(ep_value,N_value,origin_dis)
global Nr Nt Na epsilon dt dr da theta0 rr

if nargin == 0
    ep_value = 10;
    N_value = 4;
    origin_dis = 0;
end

Lr = 1; Nr = 10*N_value; dr = 1/Nr; r0 = 0:dr:Lr; r0 = r0 + (origin_dis*10+1)*dr;
Nt = 24; dt = 2*pi/Nt; theta0 = -pi+dt/2:dt:pi-dt/2;
Na = 12*N_value; da = 2*pi/Na; alpha0 = da/2:da:2*pi-da/2;
[rr,aa,tt] = meshgrid(r0,alpha0,theta0);

epsilon = 2^(-ep_value);

%% initial 
f = ones(size(rr));
data_pre = zeros(size(f));
for kt = 1:Nt
    cosT = cos(theta0(kt));
    if cosT<0
        data_pre(:,end,kt) = -cosT*ones(Na,1);
    else
        data_pre(:,1,kt) = 0.7832*ones(Na,1);
%         data_pre(:,end,kt) = ones(Na,1);
    end
end
data_pre = data_pre(:);

f = gmres(@MAB_multi,data_pre,250,1e-10);

f = reshape(f,Na,Nr+1,Nt);

extrapo = f(:,end,:); extrapo = reshape(extrapo,Na,Nt); extrapo = sum(sum(extrapo))/Na/Nt;
save(['disc/ep',num2str(ep_value),'_',num2str(origin_dis),'.mat']);

handle_f = figure(1);
set(gca,'fontsize',20);
f_alpha = f(1,:,:); f_alpha = reshape(f_alpha,Nr+1,Nt);
mesh(f_alpha); title(['\epsilon = 1'],'fontsize',20);
print(gcf,'-depsc2',['disc/ep',num2str(ep_value),'_',num2str(origin_dis),'.eps']);
close(handle_f);

% for kr = 1:Nr+1
%     f_r0 = f(:,kr,:);
%     f_r0 = reshape(f_r0,Na,Nt);
%     mesh(theta0,alpha0,f_r0); title([num2str(r0(kr))]);
%     xlabel('\theta'); ylabel('\alpha');
%     pause;
% end
end

function Mb = MAB_multi(p)
global epsilon Nr Nt Na da dr rr theta0
e = ones(Nt,1);
v = theta0(:);

p = reshape(p,Na,Nr+1,Nt);
Ap = zeros(size(p));
Rp = zeros(size(p));
BCp = zeros(size(p));
boundary_D = zeros(size(p));

for kt = 1:Nt
    p_temp = p(:,:,kt);
    
    Rp_pre = (p_temp(:,2:end) - p_temp(:,1:end-1))/dr;
    
    cosT = cos(v(kt));
    if cosT<0
        Ap(:,1:Nr,kt) = ([p_temp(2:end,1:Nr);p_temp(1,1:Nr)] - [p_temp(end,1:Nr);p_temp(1:end-1,1:Nr)])/2/da;
        Ap(:,1:Nr,kt) = sin(v(kt))*Ap(:,1:Nr,kt);
        
        Rp(:,1:end-1,kt) = cos(v(kt))*Rp_pre;
        boundary_D(:,end,kt) = p(:,end,kt);
    else
        Ap(:,2:Nr+1,kt) = ([p_temp(2:end,2:Nr+1);p_temp(1,2:Nr+1)] - [p_temp(end,2:Nr+1);p_temp(1:end-1,2:Nr+1)])/2/da;
        Ap(:,2:Nr+1,kt) = sin(v(kt))*Ap(:,2:Nr+1,kt);
        
        Rp(:,2:end,kt) = cos(v(kt))*Rp_pre;
        boundary_D(:,1,kt) = p(:,1,kt);
    end
    
end

for kr = 1:Nr
    for ka = 1:Na
        p_temp = p(ka,kr,:); p_temp = p_temp(:);
        BCp(ka,kr,:) = p_temp - e*e'*p_temp/Nt;
    end
end

Mb = Ap + rr.*Rp + rr.*BCp/epsilon + boundary_D;

Mb = Mb(:);

end

