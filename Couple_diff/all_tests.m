clear;

cases = 1;

fboundary_l.fun =  inline('v-v+1+0*t','v','t');
fboundary_r.fun =  inline('v-v+1+0*t','v','t');

phi_x.fun = inline('x-x','x');
phi_v.fun = inline('v-v+eta-eta+1','v','eta');

ChangeInTime = 0;
comparison_compatible_to_play;


clear;
cases = 2;

fboundary_l.fun =  inline('1.5+100*t*abs(v)','v','t');
fboundary_r.fun =  inline('1.5+100*t*abs(v)','v','t');

phi_x.fun = inline('sin(pi*x) + 1.5','x');
phi_v.fun = inline('1 + v - v + eta - eta','v','eta');

ChangeInTime = 1;
comparison_compatible_to_play;
% 
% % 
% clear;
% cases = 3;
% 
% fboundary_l.fun =  inline('abs(v)+t-t','v','t');
% fboundary_r.fun =  inline('abs(v)+t-t','v','t');
% 
% phi_x.fun = inline('x-x+1','x');
% phi_v.fun = inline('eta+v-v','v','eta');
% 
% ChangeInTime = 0;
% comparison_compatible_to_play;
% 
% 
% clear;
% cases = 4;
% 
% fboundary_l.fun =  inline('abs(v)+t-t','v','t');
% fboundary_r.fun =  inline('abs(v)+t-t','v','t');
% 
% phi_x.fun = inline('1 + 0.25*sin(pi*x)','x');
% phi_v.fun = inline('eta + v - v','v','eta');
% 
% ChangeInTime = 0;
% comparison_compatible_to_play;
% 

clear;
cases = 5;

fboundary_l.fun =  inline('abs(v)*(1+100*t)','v','t');
fboundary_r.fun =  inline('abs(v)*(1+100*t)','v','t');

phi_x.fun = inline('x-x+1','x');
phi_v.fun = inline('abs(v)*eta + eta/2','v','eta');

ChangeInTime = 1;
comparison_compatible_to_play;


clear;
cases = 6;

fboundary_l.fun =  inline('abs(v)+t-t','v','t');
fboundary_r.fun =  inline('abs(v)+t-t','v','t');

phi_x.fun = inline('x-x+1','x');
phi_v.fun = inline('abs(v)','v','eta');

ChangeInTime = 0;
comparison_compatible_to_play;

clear;
cases = 7;

fboundary_l.fun =  inline('v - v + t-t','v','t');
fboundary_r.fun =  inline('v - v + t-t','v','t');

phi_x.fun = inline('sin(pi*x)','x');
phi_v.fun = inline('1+0.5*abs(v) + eta - eta','v','eta');

ChangeInTime = 0;
comparison_compatible_to_play;


clear;
cases = 8;

fboundary_l.fun =  inline('v - v + t-t','v','t');
fboundary_r.fun =  inline('v - v + t-t','v','t');

phi_x.fun = inline('sin(pi*x)','x');
phi_v.fun = inline('1+ v - v + eta - eta','v','eta');

ChangeInTime = 0;
comparison_compatible_to_play;

% 
% clear;
% cases = 9;
% 
% fboundary_l.fun =  inline('v - v + t-t','v','t');
% fboundary_r.fun =  inline('v - v + t-t','v','t');
% 
% phi_x.fun = inline('(x-1).^2.*(x+1).^2','x');
% phi_v.fun = inline('1+ v - v + eta - eta','v','eta');
% 
% ChangeInTime = 0;
% comparison_compatible_to_play;

% case 10-12 are with Robin boundary condition, just for fun
% clear;
% cases = 10;
% 
% fboundary_l.fun =  inline('v - v + t + 1','v','t');
% fboundary_r.fun =  inline('v - v + t + 1','v','t');
% 
% phi_x.fun = inline('x - x + 1','x');
% phi_v.fun = inline('v - v + eta','v','eta');
% 
% ChangeInTime = 1;
% comparison_compatible_to_play;
% % 
% clear;
% cases = 11;
% 
% fboundary_l.fun =  inline('v - v + t.^2 + 1','v','t');
% fboundary_r.fun =  inline('v - v + t.^2 + 1','v','t');
% 
% phi_x.fun = inline('x - x + 1','x');
% phi_v.fun = inline('v - v + eta','v','eta');
% 
% ChangeInTime = 1;
% comparison_compatible_to_play;
% 
% 
% clear;
% cases = 12;
% 
% fboundary_l.fun =  inline('v - v + t.^2 + 1','v','t');
% fboundary_r.fun =  inline('v - v + t.^2 + 1','v','t');
% 
% phi_x.fun = inline('x - x + 1','x');
% phi_v.fun = inline('v - v + 1','v','eta');
% 
% ChangeInTime = 1;
% comparison_compatible_to_play;

% 
% 
% clear;
% cases = 13;
% 
% fboundary_l.fun =  inline('v - v + t.^4 + 1','v','t');
% fboundary_r.fun =  inline('v - v + t.^4 + 1','v','t');
% 
% phi_x.fun = inline('x - x + 1','x');
% phi_v.fun = inline('v - v + 1','v','eta');
% 
% ChangeInTime = 1;
% comparison_compatible_to_play;

% 
% 
% clear;
% cases = 14;
% 
% fboundary_l.fun =  inline('v - v + t - t','v','t');
% fboundary_r.fun =  inline('v - v + t - t','v','t');
% 
% phi_x.fun = inline('sin(pi*x)','x');
% phi_v.fun = inline('abs(v)','v','eta');
% 
% ChangeInTime = 0;
% comparison_compatible_to_play;