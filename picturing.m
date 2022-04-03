bump_fun = load('Guo_cont/ep_00_40_bump3.mat');
dent_fun = load('Guo_cont/ep_00_40_dent.mat');

bump = bump_fun.soln(:,1);
dent = dent_fun.soln(:,1);
theta0 = bump_fun.theta0;

lambda = (1-mean(dent))/(1+mean(bump)-mean(dent));
soln = lambda * bump + (1-lambda) * dent;

handle_f = figure(1);
subplot(3,1,1);plot(bump_fun.theta0,bump,'.-');ylabel('f_1','fontsize',20);
subplot(3,1,2);plot(bump_fun.theta0,dent,'.-');ylabel('f_2','fontsize',20);
subplot(3,1,3);plot(bump_fun.theta0,soln,'.-');ylabel('F','fontsize',20);
xlabel('\theta','fontsize',20);
print(gcf,'-depsc2','Guo_cont/cont_dent_bump3');
close(handle_f);


handle_f = figure(1);
subplot(3,1,1);plot(bump_fun.theta0,bump,'.-');ylabel('f_1','fontsize',20);xlim([-1 1]);
subplot(3,1,2);plot(bump_fun.theta0,dent,'.-');ylabel('f_2','fontsize',20);xlim([-1 1]);
subplot(3,1,3);plot(bump_fun.theta0,soln,'.-');ylabel('F','fontsize',20);xlim([-1 1]);
xlabel('\theta','fontsize',20);
print(gcf,'-depsc2','Guo_cont/cont_dent_bump3_zoom');
close(handle_f);


handle_f = figure(1);
subplot(3,1,1);plot(bump_fun.theta0,(bump-mean(bump))./sin(theta0)','.-');ylabel('f_1','fontsize',20);
title('(f-<f>)/sin\theta','fontsize',20)
subplot(3,1,2);plot(bump_fun.theta0,(dent-mean(dent))./sin(theta0)','.-');ylabel('f_2','fontsize',20);
subplot(3,1,3);plot(bump_fun.theta0,(soln-mean(soln))./sin(theta0)','.-');ylabel('F','fontsize',20);
xlabel('\theta','fontsize',20);
print(gcf,'-depsc2','Guo_cont/cont_dent_bump3_lip');
close(handle_f);


% ep0 = load('const_incoming/ep0_cond.mat');
% ep1 = load('const_incoming/ep1_cond.mat');
% ep2 = load('const_incoming/ep2_cond.mat');
% 
% 
% f_r0 = ep0.f(1,:,:); theta0 = ep0.theta0; r0 = ep0.r0;
% Nr = ep0.Nr; Nt = ep0.Nt;
% f_r0 = reshape(f_r0,Nr+1,Nt);
% handle_f = figure(1);
% set(gca,'fontsize',20);
% mesh(theta0,r0,f_r0); title(['f(\theta,r)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('r','fontsize',20);
% print(gcf,'-depsc2','const_incoming/ep0_f');
% close(handle_f);
% handle_f = figure(3);
% set(gca,'fontsize',20);
% plot(theta0,f_r0(Nr+1,:)); title(['f(\theta,r=1)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('f','fontsize',20);
% print(gcf,'-depsc2','const_incoming/ep0_f_r1');
% close(handle_f);
% handle_f = figure(3);
% set(gca,'fontsize',20);
% plot(theta0,f_r0(Nr/2+1,:)); title(['f(\theta,r=0.5)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('f','fontsize',20);
% print(gcf,'-depsc2','const_incoming/ep0_f_r05');
% close(handle_f);
% 
% f_r0 = ep1.f(1,:,:);
% theta0 = ep1.theta0; r0 = ep1.r0;
% Nr = ep1.Nr; Nt = ep1.Nt;
% f_r0 = reshape(f_r0,Nr+1,Nt);
% handle_f = figure(1);
% set(gca,'fontsize',20);
% mesh(theta0,r0,f_r0); title(['f(\theta,r)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('r','fontsize',20);
% print(gcf,'-depsc2','const_incoming/ep1_f');
% close(handle_f);
% handle_f = figure(3);
% set(gca,'fontsize',20);
% plot(theta0,f_r0(Nr+1,:)); title(['f(\theta,r=1)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('f','fontsize',20);
% print(gcf,'-depsc2','const_incoming/ep1_f_r1');
% close(handle_f);
% handle_f = figure(3);
% set(gca,'fontsize',20);
% plot(theta0,f_r0(Nr/2+1,:)); title(['f(\theta,r=0.5)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('f','fontsize',20);
% print(gcf,'-depsc2','const_incoming/ep1_f_r05');
% close(handle_f);
% 
% 
% f_r0 = ep2.f(1,:,:); theta0 = ep2.theta0; r0 = ep2.r0;
% Nr = ep2.Nr; Nt = ep2.Nt;
% f_r0 = reshape(f_r0,Nr+1,Nt);
% handle_f = figure(1);
% set(gca,'fontsize',20);
% mesh(theta0,r0,f_r0); title(['f(\theta,r)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('r','fontsize',20);
% print(gcf,'-depsc2','const_incoming/ep2_f');
% close(handle_f);
% handle_f = figure(3);
% set(gca,'fontsize',20);
% plot(theta0,f_r0(Nr+1,:)); title(['f(\theta,r=1)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('f','fontsize',20);
% print(gcf,'-depsc2','const_incoming/ep2_f_r1');
% close(handle_f);
% handle_f = figure(3);
% set(gca,'fontsize',20);
% plot(theta0,f_r0(Nr/2+1,:)); title(['f(\theta,r=0.5)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('f','fontsize',20);
% print(gcf,'-depsc2','const_incoming/ep2_f_r05');
% close(handle_f);
% 
% 
% return
% 
% ep0 = load('v_incoming/ep0_cond_v.mat');
% ep1 = load('v_incoming/ep1_cond_v.mat');
% ep2 = load('v_incoming/ep2_cond_v.mat');
% 
% 
% f_r0 = ep0.f(1,:,:); theta0 = ep0.theta0; r0 = ep0.r0;
% Nr = ep0.Nr; Nt = ep0.Nt;
% f_r0 = reshape(f_r0,Nr+1,Nt);
% handle_f = figure(1);
% set(gca,'fontsize',20);
% mesh(theta0,r0,f_r0); title(['f(\theta,r)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('r','fontsize',20);
% print(gcf,'-depsc2','v_incoming/ep0_f');
% close(handle_f);
% handle_f = figure(3);
% set(gca,'fontsize',20);
% plot(theta0,f_r0(Nr+1,:)); title(['f(\theta,r=1)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('f','fontsize',20);
% print(gcf,'-depsc2','v_incoming/ep0_f_r1');
% close(handle_f);
% handle_f = figure(3);
% set(gca,'fontsize',20);
% plot(theta0,f_r0(Nr/2+1,:)); title(['f(\theta,r=0.5)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('f','fontsize',20);
% print(gcf,'-depsc2','v_incoming/ep0_f_r05');
% close(handle_f);
% 
% f_r0 = ep1.f(1,:,:);
% theta0 = ep1.theta0; r0 = ep1.r0;
% Nr = ep1.Nr; Nt = ep1.Nt;
% f_r0 = reshape(f_r0,Nr+1,Nt);
% handle_f = figure(1);
% set(gca,'fontsize',20);
% mesh(theta0,r0,f_r0); title(['f(\theta,r)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('r','fontsize',20);
% print(gcf,'-depsc2','v_incoming/ep1_f');
% close(handle_f);
% handle_f = figure(3);
% set(gca,'fontsize',20);
% plot(theta0,f_r0(Nr+1,:)); title(['f(\theta,r=1)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('f','fontsize',20);
% print(gcf,'-depsc2','v_incoming/ep1_f_r1');
% close(handle_f);
% handle_f = figure(3);
% set(gca,'fontsize',20);
% plot(theta0,f_r0(Nr/2+1,:)); title(['f(\theta,r=0.5)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('f','fontsize',20);
% print(gcf,'-depsc2','v_incoming/ep1_f_r05');
% close(handle_f);
% 
% 
% f_r0 = ep2.f(1,:,:); theta0 = ep2.theta0; r0 = ep2.r0;
% Nr = ep2.Nr; Nt = ep2.Nt;
% f_r0 = reshape(f_r0,Nr+1,Nt);
% handle_f = figure(1);
% set(gca,'fontsize',20);
% mesh(theta0,r0,f_r0); title(['f(\theta,r)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('r','fontsize',20);
% print(gcf,'-depsc2','v_incoming/ep2_f');
% close(handle_f);
% handle_f = figure(3);
% set(gca,'fontsize',20);
% plot(theta0,f_r0(Nr+1,:)); title(['f(\theta,r=1)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('f','fontsize',20);
% print(gcf,'-depsc2','v_incoming/ep2_f_r1');
% close(handle_f);
% handle_f = figure(3);
% set(gca,'fontsize',20);
% plot(theta0,f_r0(Nr/2+1,:)); title(['f(\theta,r=0.5)'],'fontsize',20);
% xlabel('\theta','fontsize',20); ylabel('f','fontsize',20);
% print(gcf,'-depsc2','v_incoming/ep2_f_r05');
% close(handle_f);