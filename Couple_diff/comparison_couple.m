clear;


% for k=6:8
%     sigma_value = 2^k;
%     epsilon_value = 2^(-k);
%     Reference_couple(sigma_value,epsilon_value);
%     Couple(sigma_value);
% end

threshold = 0.5;

for k = 6:8
    sigma_value = 2^k;
    filename = ['data/transport_couple_',num2str(sigma_value)];
    transport(k-5) = load(filename);
    theta = transport(k-5).theta_all;
    filename = ['data/couple_',num2str(sigma_value)];
    couple(k-5) = load(filename);
    theta_couple = couple(k-5).theta_all;
    x = couple(k-5).x; x_transport = transport(k-5).x;
    ratio = length(x_transport)/length(x);
    
    theta = reshape(theta,ratio,length(x));
    theta = sum(theta,1)/20;
    
    theta_t = theta((x>-threshold)&(x<threshold));
    theta_couple_t = theta_couple((x>-threshold)&(x<threshold));
    error(k-5) = norm(theta_t - theta_couple_t,2)*sqrt(x(2) - x(1));
end

h = figure(2);
set(gca,'fontsize',20);
plot(transport(1).x,transport(1).theta_all,'.-.',couple(1).x,couple(1).theta_all,...
    transport(2).x,transport(2).theta_all,'.-.',couple(2).x,couple(2).theta_all,...
    transport(3).x,transport(3).theta_all,'.-.',couple(3).x,couple(3).theta_all);
legend('\sigma = 64','couple, 64','\sigma = 128','couple, 128',...
    '\sigma = 256','couple, 256');
xlabel('x');
ylabel('\Theta');
print(gcf,'-depsc2', 'data/couple_compatible.eps');
close(h);


h = figure(2);
set(gca,'fontsize',20);
plot(transport(1).x,transport(1).theta_all,'.-.',couple(1).x,couple(1).theta_all,...
    transport(2).x,transport(2).theta_all,'.-.',couple(2).x,couple(2).theta_all,...
    transport(3).x,transport(3).theta_all,'.-.',couple(3).x,couple(3).theta_all);
legend('\sigma = 64','couple, 64','\sigma = 128','couple, 128',...
    '\sigma = 256','couple, 256');
xlabel('x');
ylabel('\Theta');
xlim([-1 -0.9]);
print(gcf,'-depsc2', 'data/couple_compatible_zoom.eps');
close(h);

regression = polyfit([1:3],log(error),[1]);
error_reg = polyval(regression,[1:3]);
error_reg = exp(error_reg);
h = figure(2);
set(gca,'fontsize',20);
loglog(2.^[6:8],error,'.-.',2.^[6:8],error_reg);
legend('Error','linear regression');
gtext(['slope = ',num2str(regression(1))],'fontsize',20);
xlabel('log(\sigma)'); ylabel('log(Error_\inf)');
title('error v.s. sigma');
print(gcf,'-depsc2','data/couple_error.eps');
close(h);