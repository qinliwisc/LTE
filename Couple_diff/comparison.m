clear;

% 
% Heat;
% for k=6:9
%     sigma_value = 2^k;
%     epsilon_value = 2^(-k);
%     Reference(sigma_value,epsilon_value);
% end

threshold = 0.99;

filename = ['data/heat/heat_1n3'];
heat = load(filename);
x = heat.x; theta_heat = heat.theta; theta_heat = theta_heat(:);
% 
% x_adjust = reshape(x,2,length(x)/2); x = sum(x_adjust,1)/2;
% theta_heat = reshape(theta_heat,2,length(theta_heat)/2);
% theta_heat = sum(theta_heat,1)/2; theta_heat = theta_heat(:);
theta_heat_t = theta_heat((x>-threshold)&(x<threshold));

for k = 6:9
    km = k-5;
    sigma_value = 2^k;
    filename = ['data/heat/transport_',num2str(sigma_value)];
    transport(km) = load(filename);
    theta = transport(km).theta_all;
    x_transport = transport(km).x;
    ratio = length(x_transport)/length(x);
    
    theta = reshape(theta,ratio,length(x));
    theta = sum(theta,1)/ratio; theta = theta(:);
    
    theta_t = theta((x>-threshold)&(x<threshold));
    
%     error(km) = norm(theta_t - theta_heat_t,2)*sqrt(x(2) - x(1));
    error(km) = norm(theta_t - theta_heat_t,inf);
end


h = figure(2);
set(gca,'fontsize',20);
plot(transport(1).x,transport(1).theta_all,'.-.',...
    transport(2).x,transport(2).theta_all,'.-.',...
    transport(3).x,transport(3).theta_all,'.-.',...
    transport(4).x,transport(4).theta_all,'.-.',...
    heat.x,heat.theta);
legend('\sigma = 64','\sigma = 128','\sigma = 256','\sigma = 512','heat');
xlabel('x');
ylabel('\Theta');
print(gcf,'-depsc2', 'data/heat/pure_heat.eps');
close(h);

h = figure(2);
set(gca,'fontsize',20);
plot(transport(1).x,transport(1).theta_all,'.-.',...
    transport(2).x,transport(2).theta_all,'.-.',...
    transport(3).x,transport(3).theta_all,'.-.',...
    transport(4).x,transport(4).theta_all,'.-.',...
    heat.x,heat.theta);
legend('\sigma = 64','\sigma = 128','\sigma = 256','heat');
xlabel('x');
xlim([-1 -0.9]);
ylabel('\Theta');
print(gcf,'-depsc2', 'data/heat/pure_heat_zoom.eps');
close(h);

regression = polyfit([6:9],log(error),[1]);
error_reg = polyval(regression,[6:9]);
error_reg = exp(error_reg);
h = figure(2);
set(gca,'fontsize',20);
loglog(2.^[6:9],error,'.-.',2.^[6:9],error_reg);
legend('Error','linear regression');
gtext(['slope = ',num2str(regression(1))],'fontsize',20);
xlabel('log(\sigma)'); ylabel('log(Error_\inf)');
title('error v.s. sigma');
print(gcf,'-depsc2','data/heat/error.eps');
close(h);