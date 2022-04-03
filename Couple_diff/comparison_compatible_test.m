clear;


% Heat_compatible_test;
% for k=8:8
%     sigma_value = 2^k;
%     epsilon_value = 2^(-k);
%     Reference_compatible_test(sigma_value,epsilon_value);
% end

threshold = 0.98;

filename = ['data/heat_compatible_test/heat_1n3'];
heat = load(filename);
x = heat.x; theta_heat = heat.theta; theta_heat = theta_heat(:);
theta_heat_t = theta_heat((x>-threshold)&(x<threshold));


% filename = ['data/heat_compatible_test/heat_1n4'];
% heat = load(filename);
% x = heat.x; theta_heat = heat.theta; theta_heat = theta_heat(:);
% 
% x_adjust = reshape(x,10,length(x)/10); x = sum(x_adjust,1)/10;
% theta_heat = reshape(theta_heat,10,length(theta_heat)/10);
% theta_heat = sum(theta_heat,1)/10; theta_heat = theta_heat(:);
% theta_heat_t = theta_heat((x>-threshold)&(x<threshold));
% 
% 
% filename = ['data/heat_compatible_test/heat_5n4'];
% heat = load(filename);
% x = heat.x; theta_heat = heat.theta; theta_heat = theta_heat(:);
% 
% x_adjust = reshape(x,2,length(x)/2); x = sum(x_adjust,1)/2;
% theta_heat = reshape(theta_heat,2,length(theta_heat)/2);
% theta_heat = sum(theta_heat,1)/2; theta_heat = theta_heat(:);
% theta_heat_t = theta_heat((x>-threshold)&(x<threshold));

for k = 6:8
    km = k-5;
    sigma_value = 2^k;
    filename = ['data/heat_compatible_test/transport_',num2str(sigma_value)];
    transport(km) = load(filename);
    theta = transport(km).theta_all;
    x_transport = transport(km).x;
    ratio = length(x_transport)/length(x);
    
    theta = reshape(theta,ratio,length(x));
    theta = sum(theta,1)/ratio; theta = theta(:);
    
    theta_t = theta((x>-threshold)&(x<threshold));
    
    error(km) = norm(theta_t - theta_heat_t,2)*sqrt(x(2) - x(1));
%     error(km) = norm(theta_t - theta_heat_t,inf);
end


h = figure(2);
set(gca,'fontsize',20);
plot(transport(1).x,transport(1).theta_all,'.-.',...
    transport(2).x,transport(2).theta_all,'.-.',...
    transport(3).x,transport(3).theta_all,'.-.',...
    heat.x,heat.theta);
legend('\sigma = 64','\sigma = 128','\sigma = 256','heat');
xlabel('x');
ylabel('\Theta');
print(gcf,'-depsc2', 'data/heat_compatible_test/pure_heat.eps');
close(h);

h = figure(2);
set(gca,'fontsize',20);
plot(transport(1).x,transport(1).theta_all,'.-.',...
    transport(2).x,transport(2).theta_all,'.-.',...
    transport(3).x,transport(3).theta_all,'.-.',...
    heat.x,heat.theta);
legend('\sigma = 64','\sigma = 128','\sigma = 256','heat');
xlabel('x');
xlim([-1 -0.9]);
ylabel('\Theta');
print(gcf,'-depsc2', 'data/heat_compatible_test/pure_heat_zoom.eps');
close(h);

regression = polyfit([6:8],log(error),[1]);
error_reg = polyval(regression,[6:8]);
error_reg = exp(error_reg);
h = figure(2);
set(gca,'fontsize',20);
loglog(2.^[6:8],error,'.-.',2.^[6:8],error_reg);
legend('Error','linear regression');
gtext(['slope = ',num2str(regression(1))],'fontsize',20);
xlabel('log(\sigma)'); ylabel('log(Error_\inf)');
title('error v.s. sigma');
print(gcf,'-depsc2','data/heat_compatible_test/error.eps');
close(h);