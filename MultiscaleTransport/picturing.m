
% foldername = '1D/20pi_heat/';
% heat200_10 = load([foldername,'heat_200_10.mat']);
% heat100_20 = load([foldername,'heat_100_20.mat']);
% heat50_40 = load([foldername,'heat_50_40.mat']);
% 
% figure_h = figure(2);
% set(gca,'fontsize',20);
% plot(heat200_10.x0_H(1:end-1),heat200_10.rho(:,end),...
%     heat100_20.x0_H(1:end-1),heat100_20.rho(:,end),'.-.',...
%     heat50_40.x0_H(1:end-1),heat50_40.rho(:,end),'o--');
% xlabel('x','fontsize',20);ylabel('\rho','fontsize',20);
% legend('N_x = 200','N_x = 100','N_x = 50');
% print(gcf,'-depsc',[foldername,'heat_resolution.eps']);
% close(figure_h);

% 
% foldername = '1D/20pi_transport/';
% trans200_0 = load([foldername,'ep_0_200_10.mat']);
% trans100_0 = load([foldername,'ep_0_100_20.mat']);
% trans200_1 = load([foldername,'ep_1_200_10.mat']);
% trans100_1 = load([foldername,'ep_1_100_20.mat']);
% trans200_2 = load([foldername,'ep_2_200_10.mat']);
% trans100_2 = load([foldername,'ep_2_100_20.mat']);
% heat = load(['1D/20pi_heat/heat_200_10.mat']);
% 
% figure_h = figure(1);
% set(gca,'fontsize',20);
% plot(trans200_0.x0_H(1:end-1),trans200_0.rho,trans100_0.x0_H(1:end-1),trans100_0.rho,'o--');
% xlabel('x','fontsize',20);ylabel('\rho','fontsize',20);
% legend('N_x = 200, \epsilon = 1','N_x = 100, \epsilon = 1');
% print(gcf,'-depsc',[foldername,'transport_resolution.eps']);
% close(figure_h);
% 
% figure_h = figure(2);
% set(gca,'fontsize',20);
% plot(trans100_0.x0_H(1:end-1),trans100_0.rho,'--',...
%     trans100_1.x0_H(1:end-1),trans100_1.rho,'*-.',...
%     trans100_2.x0_H(1:end-1),trans100_2.rho,'o--',...
%     heat.x0_H(1:end-1),heat.rho(:,end));
% xlabel('x','fontsize',20);ylabel('\rho','fontsize',20);
% legend('\epsilon = 1','\epsilon = 0.1','\epsilon = 0.01','heat');
% print(gcf,'-depsc',[foldername,'transport_epsilon.eps']);
% close(figure_h);


foldername = '2D/ex3_heat_symm/';
foldername2 = '2D/ex3_transport_symm/';
heat50 = load([foldername,'rho_50_40.mat']);
heat100 = load([foldername,'rho_100_20.mat']);
trans0 = load([foldername2,'transport_0.mat']);
trans1 = load([foldername2,'transport_1.mat']);
trans2 = load([foldername2,'transport_2.mat']);

rho_heat1 = heat50.rho{end};
rho_heat2 = heat100.rho{end};

figure_h = figure(2);
set(gca,'fontsize',20);
subplot(1,3,1);mesh(heat50.x0_H(1:end-1),heat50.x0_H(1:end-1),rho_heat1);
xlabel('x');ylabel('y');zlabel('\rho');zlim([0 2]);title('N_x = 50');
subplot(1,3,2);mesh(heat100.x0_H(1:end-1),heat100.x0_H(1:end-1),rho_heat2);
xlabel('x');ylabel('y');zlabel('\rho');zlim([0 2]);title('N_x = 100');
subplot(1,3,3);mesh(heat50.x0_H(1:end-1),heat50.x0_H(1:end-1),rho_heat1 - rho_heat2(1:2:end,1:2:end));
xlabel('x');ylabel('y'); zlim([-0.2 0.2]);title('error');
print(gcf,'-depsc',[foldername,'heat_resolution.eps']);
close(figure_h);

% figure_h = figure(1);
% plot(trans0.x0_H(1:end-1),trans0.rho(:,end/2),'--',...
%     trans1.x0_H(1:end-1),trans1.rho(:,end/2),'*-.',...
%     trans2.x0_H(1:end-1),trans2.rho(:,(end)/2),'o--',...
%     heat100.x0_H(1:end-1),rho_heat2(:,end/2),'.-.');
% xlabel('(x=0, y)','fontsize',20);ylabel('\rho','fontsize',20);
% legend('\epsilon = 1','\epsilon = 0.1','\epsilon = 0.01','heat');
% print(gcf,'-depsc',[foldername,'comparison_y.eps']);
% close(figure_h);
% 
% figure_h = figure(1);
% plot(trans0.x0_H(1:end-1),trans0.rho(end/2,:),'--',...
%     trans1.x0_H(1:end-1),trans1.rho(end/2,:),'*-.',...
%     trans2.x0_H(1:end-1),trans2.rho((end)/2,:),'o--',...
%     heat100.x0_H(1:end-1),rho_heat2((end)/2,:),'.-.');
% xlabel('(x, y=0)','fontsize',20);ylabel('\rho','fontsize',20);
% legend('\epsilon = 1','\epsilon = 0.1','\epsilon = 0.01','heat');
% print(gcf,'-depsc',[foldername,'comparison_x.eps']);
% close(figure_h);

heat = heat50;

figure_h = figure(1);
set(gca,'fontsize',20);
subplot(2,2,1);mesh(trans0.x0_H(1:end-1),trans0.x0_H(1:end-1),trans0.rho)
xlabel('x');ylabel('y');zlabel('\rho');zlim([0 2]);title('\epsilon = 1');
subplot(2,2,2);mesh(trans1.x0_H(1:end-1),trans1.x0_H(1:end-1),trans1.rho);
xlabel('x');ylabel('y');zlabel('\rho');zlim([0 2]);title('\epsilon = 0.1');
subplot(2,2,3);mesh(trans2.x0_H(1:end-1),trans2.x0_H(1:end-1),trans2.rho);
xlabel('x');ylabel('y');zlabel('\rho');zlim([0 2]);title('\epsilon = 0.01');
subplot(2,2,4);mesh(heat.x0_H(1:end-1),heat.x0_H(1:end-1),heat.rho{end});
xlabel('x');ylabel('y');zlabel('\rho');zlim([0 2]);title('heat');
print(gcf,'-depsc',[foldername2,'comparison.eps']);
close(figure_h);

rho_heat = heat.rho{end};
figure_h = figure(2);
set(gca,'fontsize',20);
subplot(1,3,1);mesh(trans2.x0_H(1:end-1),trans2.x0_H(1:end-1),trans2.rho)
xlabel('x');ylabel('y');zlabel('\rho');zlim([0 2]);title('\epsilon = 0.01');
subplot(1,3,2);mesh(heat.x0_H(1:end-1),heat.x0_H(1:end-1),heat.rho{end});
xlabel('x');ylabel('y');zlabel('\rho');zlim([0 2]);title('heat');
subplot(1,3,3);mesh(heat.x0_H(1:2:end-1),heat.x0_H(1:2:end-1),trans2.rho - rho_heat(1:2:end,1:2:end));
xlabel('x');ylabel('y'); zlim([-0.2 0.2]);title('error');
print(gcf,'-depsc',[foldername2,'error.eps']);
close(figure_h);

% 
% figure_h = figure(2);
% set(gca,'fontsize',20);
% mesh(heat.x0_H(1:end-1),heat.x0_H(1:end-1),heat.rho{end});
% xlabel('x','fontsize',20);ylabel('y','fontsize',20); zlim([0 2]);title('heat','fontsize',20);
% print(gcf,'-depsc',[foldername,'heat.eps']);
% close(figure_h);
% 
% figure_h = figure(2);
% set(gca,'fontsize',20);
% mesh(heat.x0_H,heat.x0_H,heat.sigma_fun(heat.y0,heat.x0));
% xlabel('x','fontsize',20);ylabel('y','fontsize',20); zlim([0 2]);title('\sigma','fontsize',20);
% print(gcf,'-depsc',[foldername,'sigma.eps']);
% close(figure_h);

