threshold_2 = 0.90;

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

km = 0;
for k = x_grid(1:end)
    km = km+1;
    sigma_value = 2^k;
    filename = ['data/heat_compatible_to_play/case',num2str(cases),'/transport_',num2str(sigma_value)];
    transport(km) = load(filename);

    v_num = length(transport(km).V_full);
    Weight = transport(km).Weight;
    theta_transport = transport(km).theta_all;
    phi_transport = transport(km).phi;
    x_transport = transport(km).x;
    ratio = length(x_transport)/Nx;
    
    
    
    theta_transport = reshape(theta_transport',ratio,Nx);
    theta_transport = sum(theta_transport,1)/ratio;
    theta_transport = theta_transport(:);
    
    theta_t_2 = theta_transport((x_stan>-threshold_2)&(x_stan<threshold_2));
    
    error_2(km) = norm(theta_transport - theta_heat,2)*sqrt(x(2) - x(1));
    error_t_2(km) = norm(theta_t_2 - theta_heat_t_2,2)*sqrt(x(2) - x(1));
   
    phi_trans_reshape = zeros(v_num,Nx);
    for k_x = 1:Nx
        phi_temp = phi_transport(1:v_num,(k_x-1)*ratio+[1:ratio]);
        phi_temp = sum(phi_temp,2)/ratio;
        phi_trans_reshape(:,k_x) = phi_temp;
    end
    phi_heat = ones(v_num,1)*theta_heat';
    
    error_phi_v = Weight'*(phi_trans_reshape - phi_heat).^2;
    error_phi_v_threshold = error_phi_v((x_stan>-threshold_2)&(x_stan<threshold_2));
    error_phi(km) = sqrt(sum(error_phi_v)*dx_stan);
    error_t_phi(km) = sqrt(sum(error_phi_v_threshold)*dx_stan);
end

% 
% h = figure(2);
% set(gca,'fontsize',20);
% plot(transport(1).x,transport(1).theta_all,'o-.',...
%     transport(2).x,transport(2).theta_all,'.-.',...
%     transport(3).x,transport(3).theta_all,'*-.',...
%    transport(4).x,transport(4).theta_all,'--',...
%     heat.x,heat.theta);
% %     heat.x,heat.theta);
% legend('\epsilon = 1/32','\epsilon = 1/64','\epsilon = 1/128',...
%     '\epsilon = 1/256','heat','location','northwest');
% % legend('\sigma = 32','\sigma = 64','\sigma = 128','\sigma = 256','heat','location','northwest');
% xlabel('x');
% ylabel('\Theta');
% print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/pure_heat','_case',num2str(cases),'_im.eps']);
% close(h);
% 
% h = figure(2);
% set(gca,'fontsize',20);
% plot(transport(1).x,transport(1).theta_all,'o-.',...
%     transport(2).x,transport(2).theta_all,'.-.',...
%    transport(3).x,transport(3).theta_all,'*-.',...
%    transport(4).x,transport(4).theta_all,'--',...
%     heat.x,heat.theta);
% %     heat.x,heat.theta);
% legend('\epsilon = 1/32','\epsilon = 1/64','\epsilon = 1/128',...
%     '\epsilon = 1/256','heat','location','northeast');
% % legend('\sigma = 32','\sigma = 64','\sigma = 128','\sigma = 256','heat','location','northeast');
% xlabel('x');
% xlim([-1 -0.9]);
% ylabel('\Theta');
% print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/pure_heat_zoom_l','_case',num2str(cases),'_im.eps']);
% close(h);
% 
% h = figure(2);
% set(gca,'fontsize',20);
% plot(transport(1).x,transport(1).theta_all,'o-.',...
%     transport(2).x,transport(2).theta_all,'.-.',...
%     transport(3).x,transport(3).theta_all,'*-.',...
%    transport(4).x,transport(4).theta_all,'--',...
%     heat.x,heat.theta);
% %     heat.x,heat.theta);
% legend('\epsilon = 1/32','\epsilon = 1/64','\epsilon = 1/128',...
%     '\epsilon = 1/256','heat','location','northwest');
% xlabel('x');
% xlim([0.9 1]);
% ylabel('\Theta');
% print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/pure_heat_zoom_r','_case',num2str(cases),'_im.eps']);
% close(h);
% % 
% regression = polyfit(x_grid,log(error_2),[1]);
% regression_phi = polyfit(x_grid,log(error_phi),[1]);
% h = figure(2);
% set(gca,'fontsize',20);
% loglog(2.^x_grid,error_2,'o-.',2.^x_grid,error_phi,'*-.');
% legend('E_\theta','E_f');
% gtext(['slope = ',num2str(regression(1),'%6.3f')],'fontsize',20);
% gtext(['slope = ',num2str(regression_phi(1),'%6.3f')],'fontsize',20);
% xlabel('-log(\epsilon)');% ylabel('log(Error_2)');
% xlim([2^(x_grid(1)-1),2^(x_grid(end)+1)]);
% ylim([0.5*min(error_2(end),error_phi(end)),2*max(error_2(1),error_phi(1))]);
% % title(['slope = ',num2str(regression(1))]);%,', error = ',num2str(error_2(1))]);
% print(gcf,'-depsc2',['data/heat_compatible_to_play/case',num2str(cases),'/error','_case',num2str(cases),'.eps']);
% close(h);

regression = polyfit(x_grid,log(error_t_2),[1]);
regression_phi = polyfit(x_grid,log(error_t_phi),[1]);
h = figure(2);
set(gca,'fontsize',20);
loglog(2.^x_grid,error_t_2,'o-.',2.^x_grid,error_t_phi,'*-.');
legend('E_{\theta,inner}','E_{f,inner}');
gtext(['slope = ',num2str(regression(1),'%6.3f')],'fontsize',20);
gtext(['slope = ',num2str(regression_phi(1),'%6.3f')],'fontsize',20);
xlabel('-log(\epsilon)');% ylabel('log(Error_{2,inner})');
xlim([2^(x_grid(1)-1),2^(x_grid(end)+1)]);
ylim([0.5*min(error_t_2(end),error_t_phi(end)),2*max(error_t_2(1),error_t_phi(1))]);
% title(['slope = ',num2str(regression(1))]);%,', error = ',num2str(error_t_2(1))]);
print(gcf,'-depsc2',['data/heat_compatible_to_play/case',num2str(cases),'/error_',num2str(100*threshold_2),'_case',num2str(cases),'.eps']);
close(h);