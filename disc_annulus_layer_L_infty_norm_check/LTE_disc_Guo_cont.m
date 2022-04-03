function LTE_disc_Guo_cont
global Nr Nt epsilon dt dr theta0 rr

if nargin == 0
	epsilon = 0;%1/45;
    dr = 1/40;%epsilon;
    dt = 1/40;%epsilon;
end

Lr = 6; r0 = 0:dr:Lr; Nr = length(r0);
theta0 = -pi+dt/2:dt:pi-dt/2; Nt = length(theta0);
[rr,tt] = meshgrid(r0,theta0);

%% initial 
f = ones(size(rr));
data_pre1 = zeros(size(f));
data_pre2 = zeros(size(f));

for kt = 1:Nt
    sinT = sin(theta0(kt));
    if sinT>=0
%         data_pre1(kt,1) = sinT;
%         if (theta0(kt)<pi/2)
%             data_pre1(kt,1) = (theta0(kt)-pi/4)*4/pi;
%             if (theta0(kt)<pi/4)
%                 data_pre1(kt,1) = 0;
%             end
%         end
%         data_pre2(kt,1) = 0;
        data_pre1(kt,1) = 1 - abs(theta0(kt)-pi/2)*20/9/pi;%1-abs(theta0(kt)-pi/2)/pi*4;%sin(2*theta0(kt)-pi/2);%cos(2*theta0(kt)-pi/2)^2;
        if (theta0(kt)<pi/20)
            data_pre1(kt,1) = 0;
        end
        if (theta0(kt)>19*pi/20)
            data_pre1(kt,1) = 0;
        end
        data_pre2(kt,1) = 0;
    else
        data_pre1(kt,end) = 0;
        data_pre2(kt,end) = 1;
    end
end
data_pre1 = data_pre1(:);
data_pre2 = data_pre2(:);

f = gmres(@MAB_multi,data_pre1,250,1e-6);
g = gmres(@MAB_multi,data_pre2,250,1e-6);

f = reshape(f,Nt,Nr);
g = reshape(g,Nt,Nr);


alpha_limit = f(:,end)./(1-g(:,end)); alpha_limit = alpha_limit(floor(end/2):end);

alpha = mean(alpha_limit);

distance = norm(alpha_limit - alpha, inf);

if distance > 0.001
    fprintf('computing domain not large enough. \n');
end

soln = f + alpha*g;

sin_T = sin(theta0);
flux = sin_T*soln*dt;
flux2 = (sin_T.^2)*soln*dt;

save(['Guo_cont/ep_00_40_bump3.mat']);
% save(['Guo_cont/ep_00_40_dent.mat']);

end


function Mb = MAB_multi(p)
global epsilon Nr Nt dr dt rr theta0
e = ones(Nt,1);
v = theta0(:);

p = reshape(p,Nt,Nr);
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

cos_v = cos(v)*ones(1,Nr);
Tp(:,2:end-1) = ([p(2:end,2:end-1);p(1,2:end-1)] - [p(end,2:end-1);p(1:end-1,2:end-1)])/2/dt;
Tp(:,2:end-1) = cos_v(:,2:end-1).*Tp(:,2:end-1);
% Tp = ([p(2:end,:);p(1,:)] - [p(end,:);p(1:end-1,:)])/2/dt;
% Tp = sin_v.*Tp;

BCp(:,2:Nr) = p(:,2:Nr) - e*e'*p(:,2:Nr)/Nt;

rrr = (epsilon)./(1-epsilon*rr);
Mb = Rp - rrr.*Tp + BCp + boundary_D;

Mb = Mb(:);

end

