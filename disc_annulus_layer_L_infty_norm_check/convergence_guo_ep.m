clear;
N_alpha = [20];
ep = [8:15];
% origin_dis = [1:5];


guo_ref_01 = load(['Guo_ep/ep00_01_R8.mat']);
guo_ref_v0 = load(['Guo_ep/ep00_v0_R8.mat']);
f_v0 = guo_ref_v0.f(:,end); f_01 = guo_ref_01.f(:,end);
alpha_limit = f_v0./(1-f_01); alpha_limit = alpha_limit(121:end);
alpha = mean(alpha_limit);
soln_ref = guo_ref_v0.f + alpha*guo_ref_01.f;

alpha = zeros(length(ep),1);

for km = 1:length(ep)
    guo_v0 = load(['Guo_ep/ep',num2str(ep(km)),'_v0_R',num2str(8),'.mat']);
    guo_01 = load(['Guo_ep/ep',num2str(ep(km)),'_01_R',num2str(8),'.mat']);
    f_v0 = guo_v0.f(:,end); f_01 = guo_01.f(:,end);
    alpha_limit = f_v0./(1-f_01);% plot(alpha_limit); pause;

    alpha_limit = alpha_limit(121:end);

    alpha(km) = mean(alpha_limit);
    soln = guo_v0.f + alpha(km)*guo_01.f;
    
    sin_T = sin(guo_v0.theta0);
    
    flux = sin_T*soln*guo_v0.dt;
    flux2 = (sin_T.^2)*soln*guo_v0.dt;
    
    handle_f = figure(1);
    subplot(2,1,1); plot(guo_v0.r0,flux); ylabel('flux','fontsize',20);
    subplot(2,1,2); plot(guo_v0.r0,flux2); ylabel('flux 2','fontsize',20);
    print(gcf,'-depsc2',['Guo_ep/flux_',num2str(ep(km)),'.eps']);
    close(handle_f);
        
    handle_f = figure(1);
    plot(guo_v0.theta0,soln(:,1) - soln_ref(:,1)); ylabel('outgoing','fontsize',20);
    print(gcf,'-depsc2',['Guo_ep/outgoing_',num2str(ep(km)),'.eps']);
    close(handle_f);
end

% handle_f = figure(1);
% set(gca,'fontsize',20);
% plot(2.^(-ep),diff_2);
% xlabel('\epsilon','fontsize',20);ylabel('2_error','fontsize',20);
% print(gcf,'-depsc2',['Guo_ep/diff_2.eps']);
% close(handle_f); 
% 
% 
% handle_f = figure(1);
% set(gca,'fontsize',20);
% plot(2.^(-ep),diff_inf);
% xlabel('\epsilon','fontsize',20);ylabel('inf_error','fontsize',20);
% print(gcf,'-depsc2',['Guo_ep/diff_inf.eps']);
% close(handle_f);