clear all;

syms t;
a        = 1;
Ts       = 0.065;
fx(t)    = (a*sqrt(2)*cos(t)) / (sin(t)^2+1);
fy(t)    = (a*sqrt(2)*cos(t)*sin(t)) / (sin(t)^2+1);
dxdt(t)  = diff(fx,t);
dydt(t)  = diff(fy,t);

t_eval   = 0:0.01:2*pi;
dxdt_num = double(dxdt(t_eval))./Ts;
dydt_num = double(dydt(t_eval))./Ts;
vel_t    = sqrt(dxdt_num.^2 + dydt_num.^2);

figure;
plot(t_eval, vel_t)

