% clear;

% fboundary_l.fun =  inline('abs(v)','v','t');
% fboundary_r.fun =  inline('abs(v)','v','t');
% 
% phi_x.fun = inline('1 + 0.25 * sin(pi*x)','x');
% phi_v.fun = inline('eta + v - v','v','eta');

x_grid = [5:8];

if(0)  % switch for running the simulation; change to 0 if only
       % replotting figures / change domains / etc.
%     Heat_compatible_to_play(cases,fboundary_l,fboundary_r,phi_x,phi_v,ChangeInTime);
    Heat_compatible_implicit_to_play_Dirichlet(cases,fboundary_l,fboundary_r,phi_x,phi_v,ChangeInTime);

    for k=5:7%x_grid(1:end)
        sigma_value = 2^k;
        epsilon_value = 2^(-k);
%         Heat_compatible_implicit_to_play(cases,epsilon_value,fboundary_l,fboundary_r,phi_x,phi_v,ChangeInTime);
%         Heat_compatible_CN_to_play(cases,epsilon_value,fboundary_l,fboundary_r,phi_x,phi_v,ChangeInTime);
%         cases
%         Reference_compatible_to_play...
%             (cases,sigma_value,epsilon_value,fboundary_l,fboundary_r,phi_x,phi_v,ChangeInTime);
%         cases
    end
end

% picturing_robin;
% picturing_time_dirichlet;
% picturing;
% picturing_compare;
% picturing_Robin_phi;
% picturing_Dirichlet_phi;
picturing_paper;

% threshold_2 = 1;
% threshold_inf = 0.80;
% 
% filename = ['data/heat_compatible_to_play/case',num2str(cases),'/heat_5n4'];
% heat = load(filename);
% x = heat.x; theta_heat = heat.theta; theta_heat = theta_heat(:);
% % 
% % x_adjust = reshape(x,2,length(x)/2); x = sum(x_adjust,1)/2;
% % theta_heat = reshape(theta_heat,2,length(theta_heat)/2);
% % theta_heat = sum(theta_heat,1)/2; theta_heat = theta_heat(:);
% theta_heat_t_2 = theta_heat((x>-threshold_2)&(x<threshold_2));
% theta_heat_t_2_2 = theta_heat((x>-threshold_inf)&(x<threshold_inf));
% theta_heat_t_inf = theta_heat((x>-threshold_inf)&(x<threshold_inf));
% 
% for k = x_grid(1:end)
%     km = k-5;
%     sigma_value = 2^k;
%     filename = ['data/heat_compatible_to_play/case',num2str(cases),'/transport_',num2str(sigma_value)];
%     transport(km) = load(filename);
%     theta = transport(km).theta_all;
%     x_transport = transport(km).x;
%     ratio = length(x_transport)/length(x);
%     
%     theta = reshape(theta,ratio,length(x));
%     theta = sum(theta,1)/ratio; theta = theta(:);
%     
%     theta_t_2 = theta((x>-threshold_2)&(x<threshold_2));
%     theta_t_inf = theta((x>-threshold_inf)&(x<threshold_inf));
%     
%     error_2(km) = norm(theta_t_2 - theta_heat_t_2,2)*sqrt(x(2) - x(1));
%     error_inf(km) = norm(theta_t_inf - theta_heat_t_inf,inf);
%     error_2_2(km) = norm(theta_t_inf - theta_heat_t_inf,2)*sqrt(x(2) - x(1));
% end
% 
% % 
% % h = figure(2);
% % set(gca,'fontsize',20);
% % plot(transport(1).x,transport(1).theta_all,'.-.',...
% %     transport(2).x,transport(2).theta_all,'.-.',...
% %     heat.x,heat.theta);
% % %     transport(3).x,transport(3).theta_all,'.-.',...
% % %     heat.x,heat.theta);
% % % legend('\sigma = 64','\sigma = 128','\sigma = 256','heat','location','southwest');
% % xlabel('x');
% % ylabel('\Theta');
% % print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/pure_heat.eps']);
% % close(h);
% % 
% % h = figure(2);
% % set(gca,'fontsize',20);
% % plot(transport(1).x,transport(1).theta_all,'.-.',...
% %     transport(2).x,transport(2).theta_all,'.-.',...
% %     heat.x,heat.theta);
% % %    transport(3).x,transport(3).theta_all,'.-.'
% % %     heat.x,heat.theta);
% % % legend('\sigma = 64','\sigma = 128','\sigma = 256','heat','location','southeast');
% % xlabel('x');
% % xlim([-1 -0.9]);
% % ylabel('\Theta');
% % print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/pure_heat_zoom_l.eps']);
% % close(h);
% % 
% % h = figure(2);
% % set(gca,'fontsize',20);
% % plot(transport(1).x,transport(1).theta_all,'.-.',...
% %     transport(2).x,transport(2).theta_all,'.-.',...
% %     heat.x,heat.theta);
% % %     transport(3).x,transport(3).theta_all,'.-.',...
% % %     heat.x,heat.theta);
% % % legend('\sigma = 64','\sigma = 128','\sigma = 256','heat','location','southwest');
% % xlabel('x');
% % xlim([0.9 1]);
% % ylabel('\Theta');
% % print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/pure_heat_zoom_r.eps']);
% % close(h);
% 
% regression = polyfit(x_grid,log(error_2),[1]);
% error_reg = polyval(regression,x_grid);
% error_reg = exp(error_reg);
% h = figure(2);
% set(gca,'fontsize',20);
% loglog(2.^x_grid,error_2,'.-.',2.^x_grid,error_reg);
% legend('Error','linear regression');
% % gtext(['slope = ',num2str(regression(1))],'fontsize',20);
% xlabel('log(\sigma)'); ylabel('log(Error_\inf)');
% title(['error v.s. sigma, slope = ',num2str(regression(1))]);
% print(gcf,'-depsc2',['data/heat_compatible_to_play/case',num2str(cases),'/error_2.eps']);
% close(h);
% 
% 
% regression = polyfit(x_grid,log(error_2_2),[1]);
% error_reg = polyval(regression,x_grid);
% error_reg = exp(error_reg);
% h = figure(2);
% set(gca,'fontsize',20);
% loglog(2.^x_grid,error_2_2,'.-.',2.^x_grid,error_reg);
% legend('Error','linear regression');
% % gtext(['slope = ',num2str(regression(1))],'fontsize',20);
% xlabel('log(\sigma)'); ylabel('log(Error_2_inner)');
% title(['error v.s. sigma, slope = ',num2str(regression(1))]);
% print(gcf,'-depsc2',['data/heat_compatible_to_play/case',num2str(cases),'/error_2_',num2str(100*threshold_inf),'.eps']);
% close(h);
% 
% 
% regression = polyfit(x_grid,log(error_inf),[1]);
% error_reg = polyval(regression,x_grid);
% error_reg = exp(error_reg);
% h = figure(2);
% set(gca,'fontsize',20);
% loglog(2.^x_grid,error_inf,'.-.',2.^x_grid,error_reg);
% legend('Error','linear regression');
% % gtext(['slope = ',num2str(regression(1))],'fontsize',20);
% xlabel('log(\sigma)'); ylabel('log(Error_\inf)');
% title(['error v.s. sigma, slope = ',num2str(regression(1))]);
% print(gcf,'-depsc2',['data/heat_compatible_to_play/case',num2str(cases),'/error_inf',num2str(100*threshold_inf),'.eps']);
% close(h);