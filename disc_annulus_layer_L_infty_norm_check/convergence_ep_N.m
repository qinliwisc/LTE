clear;
N_alpha = [20];
ep = [12:15];
% origin_dis = [1:5];

% extrapo = zeros(length(ep),length(origin_dis));

for km = 1:length(ep)
    for kn = 1:length(N_alpha)
        LTE_disc_Guo(ep(km),N_alpha(kn),0);
%         extrapo(km,kn) = LTE_annulus(ep(km),N_alpha,origin_dis(kn));
    end
end

% for km = 1:length(ep)
%     for kn = 1:length(origin_dis)
%         annulus = load(['annulus_dis/ep',num2str(ep(km)),'_',num2str(origin_dis(kn)),'.mat']);
% 
% %         handle_f = figure(1);
% %         set(gca,'fontsize',20);
% %         f_alpha = annulus.f(1,:,:);
% %         f_alpha = reshape(f_alpha,N_alpha(km)*10+1,24);
% %         mesh(annulus.theta0,annulus.r0,f_alpha); title(['\epsilon = 2^{',num2str(-ep(kn)),'}, Nr = ',num2str(10*N_alpha(km))],'fontsize',20);
% %         xlabel('\theta','fontsize',20);ylabel('r','fontsize',20);zlabel('f','fontsize',20);
% %         print(gcf,'-depsc2',['annulus/ep',num2str(ep(kn)),'_',num2str(N_alpha(km)),'.eps']);
% %         close(handle_f);
%         
%         extrapo(km,kn) = annulus.extrapo;
%     end
% end
% 
% handle_f = figure(1); km = 1;
% for kn = 1:length(ep)
%     f_alpha = annulus(1,kn).f(1,:,:);
%     f_alpha = reshape(f_alpha,N_alpha(km)*10+1,24);
%     subplot(2,3,kn), mesh(annulus(1,km).theta0,annulus(1,km).r0,f_alpha);
%     title(['\epsilon = 2^{',num2str(-ep(kn)),'}'],'fontsize',20);
%     xlabel('\theta');ylabel('r');zlabel('f');
% end
% print(gcf,'-depsc2',['annulus/collection.eps']);
% close(handle_f);

% handle_f = figure(1);
% set(gca,'fontsize',20);
% mesh(N_alpha*10,2.^(-ep),extrapo); title(['extrapolation length'],'fontsize',20);
% xlabel('N_alpha');ylabel('N_epsilon');
% print(gcf,'-depsc2',['annulus/extrapolation.eps']);
% close(handle_f);
% % 
% 
% handle_f = figure(1);
% set(gca,'fontsize',20);
% plot(2.^(-ep),extrapo(:,end),'.-.'); title(['extrapolation length'],'fontsize',20);
% xlabel('\epsilon'); ylabel('extrapolation');
% print(gcf,'-depsc2',['annulus/extrapolation_150_plot.eps']);
% close(handle_f);

% handle_f = figure(1);
% set(gca,'fontsize',20);
% mesh((origin_dis*10+1)/40,ep,extrapo);
% title(['extrapolation length'],'fontsize',20);
% ylabel('\epsilon'); xlabel('inner radius');
% print(gcf,'-depsc2',['annulus_dis/extrapolation_40.eps']);
% close(handle_f);
% 
% mesh(extrapo);