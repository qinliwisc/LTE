function LTE_2D_ep(ep_value,Nr_value,Nv_value,freq)
global Nr Nt epsilon dt dr theta0 r0 sigma_x

if nargin == 0
    ep_value = 6;
    Nr_value = 40;
    Nv_value = 28;
    freq = 1;
end

Lr = 1; Nr = Lr*Nr_value; dr = Lr/Nr; r0 = 0:dr:Lr;
Nt = Nv_value; dt = 2*pi/Nt; theta0 = -pi+dt/2:dt:pi-dt/2;
[xx,yy,theta] = meshgrid(r0,r0,theta0);
[xx_s,yy_s] = meshgrid(r0,r0);

epsilon = 2^(-ep_value);
sigma = @(x,y,km)(1./(sin(km*pi*x)+sin(km*pi*y)+3.5)+1);
% sigma = inline('x-x+1');
sigma_x = sigma(xx_s,yy_s,freq)-sigma(xx_s,yy_s,freq)+1;

boundary_left = ones(Nr_value+1,1);
boundary_right = zeros(Nr_value+1,1);
boundary_bottom = zeros(Nr_value+1,1);
boundary_top = zeros(Nr_value+1,1);

%% initial 
f = ones(size(xx));
data_pre = zeros(size(f));
for kt = 1:Nt
    sinT = sin(theta0(kt));
    cosT = cos(theta0(kt));
    if cosT>=0
        data_pre(:,1,kt) = boundary_left;
    else
        data_pre(:,end,kt) = boundary_right;
    end
    if sinT>=0
        data_pre(1,:,kt) = boundary_bottom';
    else
        data_pre(end,:,kt) = boundary_top';
    end
end
data_pre = data_pre(:);

f = gmres(@MAB_multi,data_pre,150,1e-10);

f = reshape(f,Nr+1,Nr+1,Nt);

rho = f(:,:,1)/Nt;
for kt = 2:Nt
    rho = rho + f(:,:,kt)/Nt;
end
mesh(rho);

save(['test_2D_left/test_ep_',num2str(ep_value),'.mat']);


end


function Mb = MAB_multi(p)
global epsilon Nr Nt dr theta0 sigma_x
e = ones(Nt,1);

p = reshape(p,Nr+1,Nr+1,Nt);
Rp_x = zeros(size(p));
Rp_y = zeros(size(p));
BCp = zeros(size(p));
boundary_D = zeros(size(p));

for kt = 1:Nt
    p_temp = p(:,:,kt);
    
    Rp_pre_x = (p_temp(:,2:end) - p_temp(:,1:end-1))/dr;
    Rp_pre_y = (p_temp(2:end,:) - p_temp(1:end-1,:))/dr;
    
    sinT = sin(theta0(kt));
    cosT = cos(theta0(kt));
    if sinT<=0
        Rp_y(1:end-1,:,kt) = sinT*Rp_pre_y;
        boundary_D(end,:,kt) = p_temp(end,:);
    else
        Rp_y(2:end,:,kt) = sinT*Rp_pre_y;
        boundary_D(1,:,kt) = p_temp(1,:);
    end
    if cosT>=0
        Rp_x(:,2:end,kt) = cosT*Rp_pre_x;
        boundary_D(:,1,kt) = p_temp(:,1);
    else
        Rp_x(:,1:end-1,kt) = cosT*Rp_pre_x;
        boundary_D(:,end,kt) = p_temp(:,end);
    end
end

for kx = 2:Nr
    for ky = 2:Nr
        p_temp = p(kx,ky,:); p_temp = p_temp(:);
        p_temp = p_temp - e*e'*p_temp/Nt;
        BCp(kx,ky,:) = sigma_x(kx,ky)*p_temp;
    end
end

Mb = epsilon*Rp_x + epsilon* Rp_y + BCp + boundary_D;

Mb = Mb(:);

end

