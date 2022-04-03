threshold_2 = 1;
threshold_inf = 0.95;
% 
% filename = ['data/heat_compatible_to_play/case',num2str(cases),'/heat_5n4'];
% heat_c = load(filename);
% % heat(km) = load(filename);
% x = heat_c.x; theta_heat = heat_c.theta; theta_heat = theta_heat(:);
% 
% theta_heat_t_2 = theta_heat((x>-threshold_2)&(x<threshold_2));
% theta_heat_t_inf = theta_heat((x>-threshold_inf)&(x<threshold_inf));

for k = x_grid(1:end)
    km = (k-4);
    sigma_value = 2^k;

    filename = ['data/heat_compatible_to_play/case',num2str(cases),'/heat_implicit_5n4_',num2str(sigma_value)];
    heat_c = load(filename);
    heat(km) = load(filename);
    x = heat_c.x; theta_heat = heat_c.theta; theta_heat = theta_heat(:);

    theta_heat_t_2 = theta_heat((x>-threshold_2)&(x<threshold_2));
    theta_heat_t_inf = theta_heat((x>-threshold_inf)&(x<threshold_inf));

    filename = ['data/heat_compatible_to_play/case',num2str(cases),'/transport_refine_75_',num2str(sigma_value)];
    transport_refine(km) = load(filename);
    theta_refine = transport_refine(km).theta_all;
    x_transport = transport_refine(km).x;
    ratio = length(x_transport)/length(x);
    
    theta_refine = reshape(theta_refine,ratio,length(x));
    theta_refine = sum(theta_refine,1)/ratio; theta_refine = theta_refine(:);
    
    theta_t_r_2 = theta_refine((x>-threshold_2)&(x<threshold_2));
    theta_t_r_inf = theta_refine((x>-threshold_inf)&(x<threshold_inf));

    filename = ['data/heat_compatible_to_play/case',num2str(cases),'/transport_refine_',num2str(sigma_value)];
    transport(km) = load(filename);
    theta = transport(km).theta_all;
    x_transport = transport(km).x;
    ratio = length(x_transport)/length(x);
    
    theta = reshape(theta,ratio,length(x));
    theta = sum(theta,1)/ratio; theta = theta(:);
    
    theta_t_2 = theta((x>-threshold_2)&(x<threshold_2));
    theta_t_inf = theta((x>-threshold_inf)&(x<threshold_inf));
    
    error_2(km) = norm(theta_t_2 - theta_heat_t_2,2)*sqrt(x(2) - x(1));
    error_inf(km) = norm(theta_t_inf - theta_heat_t_inf,inf);

    error_r_2(km) = norm(theta_t_r_2 - theta_heat_t_2,2)*sqrt(x(2) - x(1));
    error_r_inf(km) = norm(theta_t_r_inf - theta_heat_t_inf,inf);

    diff_2(km) = norm(theta_t_2 - theta_t_r_2,2)*sqrt(x(2) - x(1));
    diff_inf(km) = norm(theta_t_inf - theta_t_r_inf,inf);
end

diff_2
diff_inf