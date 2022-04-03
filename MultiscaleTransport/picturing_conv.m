delta_v = [2:2:8]; delta_vr = [2:2:6];

% for km = 1:length(delta_v)
%     CheckingConvergencePhi(delta_v(km));
% end

for km = 1:4
    checkconv{km} = load(['2D/delta_conv/delta_',num2str(km*2),'.mat']);
    sigma(km) = checkconv{km}.Phi_y_mat(1,1);
    error(km) = checkconv{km}.error_py_norm(1);
%     error(km) = norm(reshape(checkconv{km}.error_py_soln{1},100^2,1),inf);
end
% error = zeros(1,3);
% for km = 1:3
%     error(km) = abs(sigma(km+1) - sigma(km))./abs(sigma(km));
% end
reg=polyfit(log(delta_v),log(error),1);
error_es = polyval(reg,log(delta_v)); error_es = exp(error_es);

figure(1)
plot(log(delta_v),log(error),'bo-',log(delta_v),log(error_es),'r--');
gtext(['slope = ',num2str(reg(1),'%1.4f')],'fontsize',20);
xlabel('log(\delta)','fontsize',20); ylabel('log(error)','fontsize',20);
legend('error: \Phi_y(1,1)','regression','fontsize',20);
print(['2D/delta_conv/phi_conv_H1.eps'],'-depsc');


delta = {};
sigma12 = zeros(1,8);
sigma12x = zeros(1,8);
sigma12y = zeros(1,8);
phi12 = zeros(1,8);

for km = 1:8
    delta{km} = load(['2D/delta_conv/',num2str(2*km),'.mat']);
end

for km = 1:8
    sigma12x(km) = delta{km}.Sigma_x_mat_temp(2,4);
    sigma12y(km) = delta{km}.Sigma_y_mat_temp(2,4);
    sigma12(km) = delta{km}.Sigma_mat_temp(1,2);
    phi12(km) = delta{km}.Phi_mat_temp(1,2);
end

delta_v = 1./(2*[1:8]);

phi_diff = (phi12(2:5) - phi12(1:4))./phi12(1:4);
sigma_diff = (sigma12(2:8) - sigma12(1:7))./sigma12(1:7);
sigmax_diff = (sigma12x(2:8) - sigma12x(1:7))./sigma12x(1:7);
sigmay_diff = (sigma12y(2:8) - sigma12y(1:7))./sigma12y(1:7);

reg=polyfit(log(delta_v(2:6)),log(-phi_diff(2:6)),1);
error_es = polyval(reg,log(delta_v(2:6))); error_es = exp(error_es);

figure(1)
plot(log(delta_v(2:6)),log(-phi_diff(2:6)),'bo-',log(delta_v(2:6)),log(error_es),'r--');
gtext(['slope = ',num2str(reg(1),'%1.4f')],'fontsize',20);
xlabel('log(\delta)','fontsize',20); ylabel('log(error)','fontsize',20);
print(['2D/delta_conv/phi12.eps'],'-depsc');

