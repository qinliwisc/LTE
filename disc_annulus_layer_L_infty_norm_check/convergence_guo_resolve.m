clear;
ep = [1:5];
epsilon = 20 + 5*ep;

flux = cell(5,1);
flux2 = cell(5,1);
ray = cell(5,1);
r0 = cell(5,1);


flux_ref = cell(5,1);
flux2_ref = cell(5,1);
ray_ref = cell(5,1);
r0_ref = cell(5,1);
for km = 1:5
    km
    guo = load(['Guo_resolve/ep_',num2str(epsilon(km)),'_1.mat']);
    ref = load(['Guo_resolve/ep_00_',num2str(epsilon(km)),'.mat']);
    soln_point(km) = guo.soln(1,2);
    for kn = 1:15%min(length(guo.soln(1,:)),length(guo.soln(:,1))/2)
        ray_temp(kn) = guo.soln(kn,kn);
        ray_temp_ref(kn) = ref.soln(kn,kn);
    end
    ray{km} = ray_temp;
    ray_ref{km} = ray_temp_ref;
    flux{km} = guo.flux;
    flux2{km} = guo.flux2;
    r0{km} = guo.r0;
    flux_ref{km} = ref.flux;
    flux2_ref{km} = ref.flux2;
    r0_ref{km} = ref.r0;
        
%     handle_f = figure(1);
%     mesh(guo.r0,guo.theta0,guo.soln);
%     print(gcf,'-depsc2',['Guo_resolve/soln_',num2str(epsilon(km)),'.eps']);
%     close(handle_f);
%     handle_f = figure(1);
%     plot(guo.theta0,guo.soln(:,1),ref.theta0,ref.soln(:,1));
%     print(gcf,'-depsc2',['Guo_resolve/soln_',num2str(epsilon(km)),'_boundary.eps']);
%     close(handle_f);

    diff = (guo.soln - ref.soln);
    diff_2(km) = sqrt(sum(sum(diff.^2))*ref.dr*ref.dt);
    diff_inf(km) = max(max(abs(diff)));
    diff_inf_end(km) = max((abs(diff(:,end))));
end

% handle_f = figure(1);
% mesh(ref.r0,ref.theta0,ref.soln);
% print(gcf,'-depsc2',['Guo_resolve/soln_ref.eps']);
% close(handle_f);
% handle_f = figure(1);
% plot(ref.theta0,ref.soln(:,end));
% print(gcf,'-depsc2',['Guo_resolve/soln_ref',num2str(epsilon(km)),'_end.eps']);
% close(handle_f);

handle_f = figure(1);
set(gca,'fontsize',20);
for km = 1:5
    plot(ray{km}-ray_ref{km},'.-.');
    hold on;
end
legend('epsilon = 1/25','epsilon = 1/30','epsilon = 1/35','epsilon = 1/40','epsilon = 1/45');
xlim([1 15]); ylabel('f_{\epsilon}-f','fontsize',20);
% print(gcf,'-depsc2',['Guo_resolve_revision/ray_diff.eps']);
% close(handle_f);


handle_f = figure(1);
set(gca,'fontsize',20);
plot(1./epsilon,diff_2,'o--');
xlabel('\epsilon','fontsize',20);title('||f-f_{\epsilon}||_2','fontsize',20,'Interpreter','tex');%ylabel('2 error','fontsize',20);
print(gcf,'-depsc2',['revision/diff_2.eps']);
%print(gcf,'-depsc2',['Guo_resolve/diff_2.eps']);
close(handle_f); 


handle_f = figure(1);
set(gca,'fontsize',20);
plot(1./epsilon,diff_inf_end,'o--');
xlabel('\epsilon','fontsize',20);title('|f_{\infty}-f_{\infty,\epsilon}|','fontsize',20,'Interpreter','tex');%ylabel('inf error','fontsize',20);
%print(gcf,'-depsc2',['Guo_resolve/diff_inf_end.eps']);
print(gcf,'-depsc2',['revision/diff_inf_end.eps']);
close(handle_f);