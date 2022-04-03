epsilon_num = [-4:-1:-7];
epsilon = 2.^epsilon_num;

for km_num = 1:length(epsilon_num)
    km = -epsilon_num(km_num);
    ep = load(['test_2D_left/test_ep_',num2str(km),'.mat']);
    ep_f = ep.f;
    
    ep = load(['test_2D_left/test_g_right10_ep_',num2str(km),'.mat']);
    g4 = ep.f;
    
    r = ep.r0; r = r(:); Nr = length(r);
    theta = ep.theta0; theta = theta(:); Nv = length(theta);
   
    Lf = zeros(size(ep_f));
    for kx = 1:Nr
        for ky = 1:Nr
            p_temp = ep_f(kx,ky,:);
            p_temp = p_temp(:);
            p_temp = ones(Nv,1)*sum(p_temp)/Nv - p_temp;
            Lf(kx,ky,:) = p_temp;
        end
    end
    CC(km_num) = abs(sum(sum(sum(g4.*Lf)))/Nr^2/Nv);
end

pp = polyfit(log10(epsilon),log10(CC),1);


figure(1)
plot(log10(epsilon),log10(CC));
hold on;
plot(log10(epsilon),pp(1)*(log10(epsilon))+pp(2),'.-.');
xlabel('log(\epsilon)'); ylabel('|G(x_0)|');
title(['Slope = ',num2str(pp(1))]);
hold off;


% for km = 1:1:8
%     ep = load(['test_sin/test_ep_',num2str(km-1),'.mat']);
%     ep_f = ep.f;
%     
%     ep = load(['test_sin/test_g0_ep_',num2str(km-1),'.mat']);
%     g0 = ep.f;
%     ep = load(['test_sin/test_g1_ep_',num2str(km-1),'.mat']);
%     g1 = ep.f;
%     
%     r = ep.r0; r = r(:);
%     v = ep.v0; v = v(:);
%    
%     f
%     Lf = ones(length(v),1)*sum(ep_f)/length(v)-ep_f;
%     aa(km) = sum(sum(g0.*Lf))/length(r)/length(v);
%     bb(km) = sum(sum(g1.*Lf))/length(r)/length(v);
% end
% 
% CC = min(abs(aa),abs(bb));
% % CC = sqrt((aa.^2+bb.^2)/2);
% epsilon = 2.^[0:-1:-7];
% pp = polyfit(log10(epsilon),log10(CC),1);
% 
% figure(1)
% plot(log10(epsilon),log10(CC));
% hold on;
% plot(log10(epsilon),pp(1)*(log10(epsilon))+pp(2),'.-.');
% hold off;
% 
% figure(1)
% for km = 1:1:4
%     kn = 2*km-2;
%     ep = load(['test_sin/test_ep_',num2str(kn),'.mat']);
%     ep_f = ep.f;
%     subplot(2,2,km); plot(ep.r0,sum(ep_f)/ep.Nv); xlabel('x');ylabel('\rho'); title(['\epsilon = 2^{-',num2str(kn),'}']);
% %     mesh(ep.r0,ep.v0,ep_f); xlabel('x'); ylabel('v'); title(['\epsilon = 2^(-kn)']);
%     hold on;
% end
% hold off;
% 
% 
% for km = 1:1:8
%     ep = load(['test/test_ep_',num2str(km-1),'.mat']);
%     ep_f = ep.f;
%     f_zero = 1-ep.rr/4;
%     error_f = ep_f-f_zero;
%     l2_error(km) = sqrt(sum(sum(error_f(:,10:71).^2))/120/81);
% end
% epsilon = 2.^[0:-1:-7];
% pp = polyfit(log10(epsilon),log10(l2_error),1);
% 
% figure(1)
% plot(log10(epsilon),log10(l2_error));
% hold on;
% plot(log10(epsilon),pp(1)*(log10(epsilon))+pp(2),'.-.');
% hold off;
% xlabel('\epsilon');ylabel('||f_0-\rho_0||_2');title(['Slope = ',num2str(pp(1))]);
% 
% for km = 1:1:9
%     ep = load(['ep_',num2str(km-1),'.mat']);
%     ep_f = load(['ep_',num2str(km-1),'_f.mat']);
%     g = ep.f;
%     lf = ones(ep_f.Nt,1)*sum(ep_f.f)/ep_f.Nt - ep_f.f;
%     varies(km,:) = sum(g.*lf)/ep_f.Nt;
%     condition(km) = max(abs(varies(km,:)))/min(abs(varies(km,:)));
% end
% 
% figure(1)
% loglog(2.^[0:1:8],condition,'.-.');
% 
% figure(2)
% subplot(2,2,1); plot(varies(1,:));title('\epsilon = 1');
% subplot(2,2,2); plot(varies(2,:));title('\epsilon = 2^{-2}');
% subplot(2,2,3); plot(varies(3,:));title('\epsilon = 2^{-4}');
% subplot(2,2,4); plot(varies(4,:));title('\epsilon = 2^{-6}');