threshold_2 = 0.90;
threshold_inf = 0.90;

dx_stan = 5e-3; x_stan = -1+dx_stan/2:dx_stan:1-dx_stan/2;
Nx = length(x_stan);


filename = ['data/heat_compatible_to_play/case',num2str(cases),'/heat_implicit_dirichlet_1n3'];
heat = load(filename);
x = heat.x; theta_heat = heat.theta; theta_heat = theta_heat(:);
ratio_h = length(x)/Nx;
    
theta_heat = reshape(theta_heat,ratio_h,Nx);
theta_heat = sum(theta_heat,1)/ratio_h;
theta_heat = theta_heat(:);

theta_heat_t_2 = theta_heat((x_stan>-threshold_2)&(x_stan<threshold_2));
theta_heat_t_inf = theta_heat((x_stan>-threshold_inf)&(x_stan<threshold_inf));

% % theta_heat_t_2 = theta_heat((x>-threshold_2)&(x<threshold_2));
% theta_heat_t_2 = theta_heat((x>-threshold_2)&(x<threshold_2));
% theta_heat_t_inf = theta_heat((x>-threshold_inf)&(x<threshold_inf));

km = 0;
for k = x_grid(1:end)
    km = km+1;
    sigma_value = 2^k;
    filename = ['data/heat_compatible_to_play/case',num2str(cases),'/transport_',num2str(sigma_value)];
    transport(km) = load(filename);
    theta = transport(km).theta_all;
    x_transport = transport(km).x;
    ratio = length(x_transport)/Nx;
    
    theta = reshape(theta,ratio,Nx);
    theta = sum(theta,1)/ratio; theta = theta(:);
    
    theta_t_2 = theta((x_stan>-threshold_2)&(x_stan<threshold_2));
    theta_t_inf = theta((x_stan>-threshold_inf)&(x_stan<threshold_inf));
    
    error_2(km) = norm(theta - theta_heat,2)*sqrt(x(2) - x(1));
    error_t_inf(km) = norm(theta_t_inf - theta_heat_t_inf,inf);
    error_t_2(km) = norm(theta_t_2 - theta_heat_t_2,2)*sqrt(x(2) - x(1));
end


h = figure(2);
set(gca,'fontsize',20);
plot(transport(1).x,transport(1).theta_all,'o-.',...
    transport(2).x,transport(2).theta_all,'.-.',...
    transport(3).x,transport(3).theta_all,'*-.',...
    heat.x,heat.theta);
%     heat.x,heat.theta);
legend('\sigma = 64','\sigma = 128','\sigma = 256','heat','location','northwest');
xlabel('x');
ylabel('\Theta');
print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/pure_heat','_case',num2str(cases),'_im.eps']);
close(h);

h = figure(2);
set(gca,'fontsize',20);
plot(transport(1).x,transport(1).theta_all,'o-.',...
    transport(2).x,transport(2).theta_all,'.-.',...
   transport(3).x,transport(3).theta_all,'*-.',...
    heat.x,heat.theta);
%     heat.x,heat.theta);
legend('\sigma = 64','\sigma = 128','\sigma = 256','heat','location','northeast');
xlabel('x');
xlim([-1 -0.9]);
ylabel('\Theta');
print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/pure_heat_zoom_l','_case',num2str(cases),'_im.eps']);
close(h);

h = figure(2);
set(gca,'fontsize',20);
plot(transport(1).x,transport(1).theta_all,'o-.',...
    transport(2).x,transport(2).theta_all,'.-.',...
    transport(3).x,transport(3).theta_all,'*-.',...
    heat.x,heat.theta);
%     heat.x,heat.theta);
% legend('\sigma = 64','\sigma = 128','\sigma = 256','heat','location','southwest');
xlabel('x');
xlim([0.9 1]);
ylabel('\Theta');
print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/pure_heat_zoom_r','_case',num2str(cases),'_im.eps']);
close(h);

regression = polyfit(x_grid,log(error_2),[1]);
error_reg = polyval(regression,x_grid);
error_reg = exp(error_reg);
h = figure(2);
set(gca,'fontsize',20);
loglog(2.^x_grid,error_2,'o-.',2.^x_grid,error_reg);
legend('Error','linear regression');
% gtext(['slope = ',num2str(regression(1))],'fontsize',20);
xlabel('log(\sigma)'); ylabel('log(Error_2)');
xlim([2^(x_grid(1)-1),2^(x_grid(end)+1)]);
ylim([0.5*error_2(end),2*error_2(1)]);
title(['slope = ',num2str(regression(1))]);%,', error = ',num2str(error_2(1))]);
print(gcf,'-depsc2',['data/heat_compatible_to_play/case',num2str(cases),'/error_2','_case',num2str(cases),'_im.eps']);
close(h);


regression = polyfit(x_grid,log(error_t_2),[1]);
error_reg = polyval(regression,x_grid);
error_reg = exp(error_reg);
h = figure(2);
set(gca,'fontsize',20);
loglog(2.^x_grid,error_t_2,'o-.',2.^x_grid,error_reg);
legend('Error','linear regression');
% gtext(['slope = ',num2str(regression(1))],'fontsize',20);
xlabel('log(\sigma)'); ylabel('log(Error_{2,inner})');
xlim([2^(x_grid(1)-1),2^(x_grid(end)+1)]);
ylim([0.5*error_t_2(end),2*error_t_2(1)]);
title(['slope = ',num2str(regression(1))]);%,', error = ',num2str(error_t_2(1))]);
print(gcf,'-depsc2',['data/heat_compatible_to_play/case',num2str(cases),'/error_2_',num2str(100*threshold_2),'_case',num2str(cases),'_im.eps']);
close(h);


regression = polyfit(x_grid,log(error_t_inf),[1]);
error_reg = polyval(regression,x_grid);
error_reg = exp(error_reg);
h = figure(2);
set(gca,'fontsize',20);
loglog(2.^x_grid,error_t_inf,'o-.',2.^x_grid,error_reg);
legend('Error','linear regression');
% gtext(['slope = ',num2str(regression(1))],'fontsize',20);
xlabel('log(\sigma)'); ylabel('log(Error_{\inf})');
xlim([2^(x_grid(1)-1),2^(x_grid(end)+1)]);
ylim([0.5*error_t_inf(end),2*error_t_inf(1)]);
title(['slope = ',num2str(regression(1))]);%,', error = ',num2str(error_t_inf(1))]);
print(gcf,'-depsc2',['data/heat_compatible_to_play/case',num2str(cases),'/error_inf',num2str(100*threshold_inf),'_case',num2str(cases),'_im.eps']);
close(h);