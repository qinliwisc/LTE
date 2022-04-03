cases = 5;

x_grid = [5:7];

threshold_2 = 0.95;
threshold_inf = 0.95;

dx_stan = 1e-3; x_stan = -1+dx_stan/2:dx_stan:1-dx_stan/2;
Nx = length(x_stan);

km = 0;
for k = x_grid(1:end)
    km = km+1;
    sigma_value = 2^k;

    filename = ['data/heat_compatible_to_play/case',num2str(cases),'/heat_implicit_5n4_',num2str(sigma_value)];
    heat_c = load(filename);
    heat(km) = load(filename);
    x = heat_c.x; theta_heat = heat_c.theta; theta_heat = theta_heat(:);
    theta_rec_heat = heat(km).theta_rec;
    ratio_h = length(x)/Nx;
    
    filename = ['data/heat_compatible_to_play/case',num2str(cases),'/transport_',num2str(sigma_value)];
    transport(km) = load(filename);
    theta_refine = transport(km).theta_all;
    x_transport = transport(km).x;
    theta_rec_transport = transport(km).theta_rec;
    ratio = length(x_transport)/Nx;

    for time = 1:300
        theta_heat = theta_rec_heat(time,:);
        theta_heat = reshape(theta_heat,ratio_h,Nx);
        theta_heat = sum(theta_heat,1)/ratio_h;
        theta_heat = theta_heat(:);

        theta_heat_t_2 = theta_heat((x_stan>-threshold_2)&(x_stan<threshold_2));


        theta_refine = theta_rec_transport(time,:);
        theta_refine = reshape(theta_refine,ratio,Nx);
        theta_refine = sum(theta_refine,1)/ratio;
        theta_refine = theta_refine(:);

        theta_t_r_2 = theta_refine((x_stan>-threshold_2)&(x_stan<threshold_2));

%         plot(x_stan,theta_refine,'.-.',x_stan,theta_heat);pause(0.005);
        
        error_r_2(km,time) = norm(theta_t_r_2 - theta_heat_t_2,2)*sqrt(x_stan(2) - x_stan(1));
    end
    
    
end

plot([1:1:300],error_r_2(1,:),'.-.',[1:1:300],error_r_2(2,:),'o',[1:1:300],error_r_2(3,:));

% 
% 
% h = figure(2);
% set(gca,'fontsize',20);
% % plot(transport(1).x,transport(1).theta_all,'.-.',...
% %     transport(2).x,transport(2).theta_all,'.-.',...
% %     transport(3).x,transport(3).theta_all,'.-.',...
% %     heat(1).x,heat(1).theta,'.-.',...
% %     heat(2).x,heat(2).theta,'.-.',...
% %     heat(3).x,heat(3).theta,'.-.');
% plot(transport(1).x,transport(1).theta_all,'.-.',...
%     transport(2).x,transport(2).theta_all,'.-.',...
%     heat(1).x,heat(1).theta,'.-.',...
%     heat(2).x,heat(2).theta,'.-.');
% % legend('\sigma = 64, kinetic','\sigma = 128, kinetic','\sigma = 256, kinetic',...
% % '\sigma = 64, heat','\sigma = 128, heat','\sigma = 256, heat','location','southwest');
% % legend('\sigma = 32, kinetic','\sigma = 64, kinetic','\sigma = 128, kinetic',...
% % '\sigma = 32, heat','\sigma = 64, heat','\sigma = 128, heat','location','southwest');
% xlabel('x');
% ylabel('\Theta');
% print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/pure_heat.eps']);
% close(h);
% 
% h = figure(2);
% set(gca,'fontsize',20);
% plot(transport(1).x,transport(1).theta_all,'.-.',...
%     transport(2).x,transport(2).theta_all,'.-.',...
%     heat(1).x,heat(1).theta,'.-.',...
%     heat(2).x,heat(2).theta,'.-.');
% % plot(transport(1).x,transport(1).theta_all,'.-.',...
% %     transport(2).x,transport(2).theta_all,'.-.',...
% %     transport(3).x,transport(3).theta_all,'.-.',...
% %     heat(1).x,heat(1).theta,'.-.',...
% %     heat(2).x,heat(2).theta,'.-.',...
% %     heat(3).x,heat(3).theta,'.-.');
% % legend('\sigma = 64, kinetic','\sigma = 128, kinetic','\sigma = 256, kinetic',...
% % '\sigma = 64, heat','\sigma = 128, heat','\sigma = 256, heat','location','southwest');
% % legend('\sigma = 32, kinetic','\sigma = 64, kinetic','\sigma = 128, kinetic',...
% % '\sigma = 32, heat','\sigma = 64, heat','\sigma = 128, heat','location','northeast');
% xlabel('x');
% xlim([-1 -0.9]);
% ylabel('\Theta');
% print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/pure_heat_zoom_l.eps']);
% close(h);
% 
% h = figure(2);
% set(gca,'fontsize',20);
% plot(transport(1).x,transport(1).theta_all,'.-.',...
%     transport(2).x,transport(2).theta_all,'.-.',...
%     heat(1).x,heat(1).theta,'.-.',...
%     heat(2).x,heat(2).theta,'.-.');
% % plot(transport(1).x,transport(1).theta_all,'.-.',...
% %     transport(2).x,transport(2).theta_all,'.-.',...
% %     transport(3).x,transport(3).theta_all,'.-.',...
% %     heat(1).x,heat(1).theta,'.-.',...
% %     heat(2).x,heat(2).theta,'.-.',...
% %     heat(3).x,heat(3).theta,'.-.');
% % legend('\sigma = 64, kinetic','\sigma = 128, kinetic','\sigma = 256, kinetic',...
% % '\sigma = 64, heat','\sigma = 128, heat','\sigma = 256, heat','location','southwest');
% % legend('\sigma = 32, kinetic','\sigma = 64, kinetic','\sigma = 128, kinetic',...
% % '\sigma = 32, heat','\sigma = 64, heat','\sigma = 128, heat','location','southwest');
% xlabel('x');
% xlim([0.9 1]);
% ylabel('\Theta');
% print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/pure_heat_zoom_r.eps']);
% close(h);
% % 
% % regression = polyfit(x_grid,log(error_2),[1]);
% % error_reg = polyval(regression,x_grid);
% % error_reg = exp(error_reg);
% % h = figure(2);
% % set(gca,'fontsize',20);
% % loglog(2.^x_grid,error_2,'.-.',2.^x_grid,error_reg);
% % legend('Error','linear regression');
% % % gtext(['slope = ',num2str(regression(1))],'fontsize',20);
% % xlabel('log(\sigma)'); ylabel('log(Error_\inf)');
% % title(['error v.s. \sigma, slope = ',num2str(regression(1))]);
% % print(gcf,'-depsc2',['data/heat_compatible_to_play/case',num2str(cases),'/error_2_32.eps']);
% % close(h);
% % 
% % regression = polyfit(x_grid,log(error_inf),[1]);
% % error_reg = polyval(regression,x_grid);
% % error_reg = exp(error_reg);
% % h = figure(2);
% % set(gca,'fontsize',20);
% % loglog(2.^x_grid,error_inf,'.-.',2.^x_grid,error_reg);
% % legend('Error','linear regression');
% % % gtext(['slope = ',num2str(regression(1))],'fontsize',20);
% % xlabel('log(\sigma)'); ylabel('log(Error_\inf)');
% % title(['error v.s. \sigma, slope = ',num2str(regression(1))]);
% % print(gcf,'-depsc2',['data/heat_compatible_to_play/case',num2str(cases),'/error_inf_32_',num2str(100*threshold_inf),'.eps']);
% % close(h);
% 
% 
% regression = polyfit(x_grid,log(error_r_2),[1]);
% error_reg = polyval(regression,x_grid);
% error_reg = exp(error_reg);
% h = figure(2);
% set(gca,'fontsize',20);
% loglog(2.^x_grid,error_r_2,'.-.',2.^x_grid,error_reg);
% legend('Error','linear regression');
% % gtext(['slope = ',num2str(regression(1))],'fontsize',20);
% xlabel('log(\sigma)'); ylabel('log(Error_{inf})');
% title(['error v.s. \sigma, slope = ',num2str(regression(1))]);
% print(gcf,'-depsc2',['data/heat_compatible_to_play/case',num2str(cases),'/error_2_',num2str(100*threshold_2),'.eps']);
% close(h);
% 
% regression = polyfit(x_grid,log(error_r_inf),[1]);
% error_reg = polyval(regression,x_grid);
% error_reg = exp(error_reg);
% h = figure(2);
% set(gca,'fontsize',20);
% loglog(2.^x_grid,error_r_inf,'.-.',2.^x_grid,error_reg);
% legend('Error','linear regression');
% % gtext(['slope = ',num2str(regression(1))],'fontsize',20);
% xlabel('log(\sigma)'); ylabel('log(Error_{inf})');
% title(['error v.s. \sigma, slope = ',num2str(regression(1))]);
% print(gcf,'-depsc2',['data/heat_compatible_to_play/case',num2str(cases),'/error_inf_',num2str(100*threshold_inf),'.eps']);
% close(h);