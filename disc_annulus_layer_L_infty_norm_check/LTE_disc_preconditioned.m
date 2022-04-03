function LTE_disc_preconditioned
global Nr Nt Na epsilon dt dr da theta0 rr r0

Lr = 1; Nr = 50; dr = 1/Nr; r0 = 0:dr:Lr;
Nt = 12; dt = 2*pi/Nt; theta0 = -pi+dt/2:dt:pi-dt/2;
Na = 60; da = 2*pi/Na; alpha0 = da/2:da:2*pi-da/2;
[rr,aa,tt] = meshgrid(r0,alpha0,theta0);

epsilon = 1;

%% initial 
f = ones(size(rr));
data_pre = zeros(size(f));
for kt = 1:Nt
    cosT = cos(theta0(kt));
    if cosT<0
        data_pre(:,end,kt) = -cosT*ones(Na,1);
    end
end
data_pre = data_pre(:);

f = gmres(@MAB_multi,data_pre,150,1e-10);

f = reshape(f,Na,Nr+1,Nt);

save('ep0.mat');

f_r0 = f(1,:,:);
f_r0 = reshape(f_r0,Nr+1,Nt);
handle_f = figure(1);
set(gca,'fontsize',20);
mesh(theta0,r0,f_r0); title(['f(\theta,r)'],'fontsize',20);
xlabel('\theta','fontsize',20); ylabel('\alpha','fontsize',20);
print(gcf,'-depsc2','f');
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
global epsilon Nr Nt Na da dr rr r0 theta0
e = ones(Nt,1);
v = theta0(:);

p = reshape(p,Na,Nr+1,Nt);
Ap = zeros(size(p));
Rp = zeros(size(p));
BCp = zeros(size(p));
Mb = zeros(size(p));
boundary_D = zeros(size(p));
boundary_C = zeros(size(p));

for kt = 1:Nt
    p_temp = p(:,:,kt);
    Ap(:,:,kt) = ([p_temp(2:end,:);p_temp(1,:)] - [p_temp(end,:);p_temp(1:end-1,:)])/2/da;
    Ap(:,:,kt) = sin(v(kt))*Ap(:,:,kt);
    
    cosT = cos(v(kt));
    if cosT<0
        Rp_pre = (p_temp(:,2:end) - p_temp(:,1:end-1))/dr;
        Rp(:,1:end-1,kt) = min(cos(v(kt)),0)*Rp_pre;
        boundary_D(:,end,kt) = p(:,end,kt);
    else
        Rp_pre = (p_temp(:,2:end) - p_temp(:,1:end-1))/dr;
        Rp(:,2:end,kt) = cos(v(kt))*Rp_pre;
        if v(kt)<0
            boundary_C(1:Na/2,1,kt) = p(1:Na/2,1,kt) - p(end:-1:Na/2+1,1,Nt/2+1-kt);
        else
            boundary_C(1:Na/2,1,kt) = p(1:Na/2,1,kt) - p(end:-1:Na/2+1,1,Nt+1-kt);
        end
    end
    
end

Mb = Ap + rr.*Rp;

for kr = 1:Nr
    for ka = 1:Na
        p_temp = p(ka,kr,:); p_temp = p_temp(:);
        BCp(ka,kr,:) = (p_temp - e*e'*p_temp/Nt)*r0(kr);
        
        Mb_local = Mb(ka,kr,:); Mb_local = Mb_local(:);
        Mb(ka,kr,:) = epsilon*Mb_local + e*e'*Mb_local*(1-epsilon)/Nt;
        Mb(ka,kr,:) = BCp(ka,kr,:) + Mb(ka,kr,:);
    end
end

Mb = Mb + boundary_D + boundary_C;

Mb = Mb(:);

end