% delta = [8,16,24,32,40,48,56,64];

delta = [8,24,40,56];

for km = 1:length(delta)
    trans{km} = load(['1D/conv_delta/trans_delta_',num2str(delta(km)),'_ep_10_64_40.mat']);
%     heat{km} = load(['1D/conv_delta/heat_delta_',num2str(delta(km)),'_64_40.mat']);
end

hom = load('1D/conv_delta/heat_hom_64.mat');


figure(2)
plot(trans{2}.x0_H(1:end-1),trans{2}.rho,'.-.',trans{3}.x0_H(1:end-1),trans{3}.rho,'o--',...
    trans{4}.x0_H(1:end-1),trans{4}.rho,'*--',hom.x0(1:end-1),hom.rho_new)
xlabel('x','fontsize',20); ylabel('\rho','fontsize',20);
legend('\delta = 1/24','\delta = 1/40','\delta = 1/56','homogenized');
print(['1D/conv_delta/compare.eps'],'-depsc');

for km = 1:length(delta)
    error(km) = norm(tran{km}.rho-hom.rho_new',2);
end

reg=polyfit(log(1./delta),log(error),1);
error_es = polyval(reg,log(1./delta)); error_es = exp(error_es);

figure(1)
plot(log(1./delta),log(error),'bo-',log(1./delta),log(error_es),'r--');
gtext(['slope = ',num2str(reg(1),'%1.4f')],'fontsize',20);
xlabel('log(\delta)','fontsize',20); ylabel('log(error)','fontsize',20);
legend('error','regression','fontsize',20);
print(['1D/conv_delta/rho_error.eps'],'-depsc');