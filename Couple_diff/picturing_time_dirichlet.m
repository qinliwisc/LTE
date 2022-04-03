threshold_2 = 1;
threshold_inf = 1;

dx_stan = 5e-3; x_stan = -1+dx_stan/2:dx_stan:1-dx_stan/2;
Nx = length(x_stan);

filename = ['data/heat_compatible_to_play/case',num2str(cases),'/heat_5n3.mat'];
heat_c = load(filename);
x = heat_c.x; theta_heat = heat_c.theta; theta_heat = theta_heat(:);
theta_rec_heat = heat_c.theta_rec;
ratio_h = length(x)/Nx;


km = 0;
for k = x_grid(1:end)
    km = km+1;
    sigma_value = 2^k;

    filename = ['data/heat_compatible_to_play/case',num2str(cases),'/transport_',num2str(sigma_value)];
    transport(km) = load(filename);
    theta_refine = transport(km).theta_all;
    x_transport = transport(km).x;
    theta_rec_transport = transport(km).theta_rec;
    ratio = length(x_transport)/Nx;

    for time = 1:50
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

h = figure(2);
set(gca,'fontsize',20);
plot([1:1:50],error_r_2(1,:),'.-.',[1:1:50],error_r_2(2,:),'o',[1:1:50],error_r_2(3,:));
xlabel('x');
ylabel('\Theta');
print(gcf,'-depsc2', ['data/heat_compatible_to_play/case',num2str(cases),'/error_dirichlet_100.eps']);
close(h);