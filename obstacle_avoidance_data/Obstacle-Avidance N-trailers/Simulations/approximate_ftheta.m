clear all; clc;

t       = casadi.MX.sym('s');
p       = 5;
q       = 5;
Ts      = 0.65;
ds      = Ts;
%
traj = 'infty';

if strcmp(traj,'cuadrado')
    fx      = p.*(sqrt(cos(t)^2).*cos(t) + sqrt(sin(t)^2).*sin(t));
    fy      = q.*(sqrt(cos(t)^2).*cos(t) - sqrt(sin(t)^2).*sin(t));
elseif strcmp(traj,'circulo')
    fx      = p.*cos(t);
    fy      = q.*sin(t);
elseif strcmp(traj,'infty')
    a   = 4;
    c   = 10;
    b   = 0.5;
    fx  = (a*sqrt(2).*cos(t))./(sin(t).^2+1);
    fy  = (c*sqrt(2).*cos(t).*sin(t))./(sin(t).^2 + b);
end
%
fx_fun  = casadi.Function('fx',{t},{fx});
fy_fun  = casadi.Function('fy',{t},{fy});
%
dfxdt   = fx_fun.jacobian;
dfydt   = fy_fun.jacobian;
%
di      = 0.01;
I       = -pi:di:3*pi;
alpha0  = atan((fy_fun(di)-fy_fun(0))/(fx_fun(di)-fx_fun(0)));
theta   = alpha0;
integrando = [];
args    = [];

for i=I
    arg         = (dfydt(i,[])*dfxdt(i+di,[]) - dfxdt(i,[])*dfydt(i+di,[])) / (dfxdt(i,[])*dfxdt(i+di,[]) + dfydt(i,[])*dfydt(i+di,[]));
    args        = [args, arg];
%     integrando  = [integrando, atan(arg)];
    integrando  = [integrando, atan2((dfydt(i,[])*dfxdt(i+di,[]) - dfxdt(i,[])*dfydt(i+di,[])), (dfxdt(i,[])*dfxdt(i+di,[]) + dfydt(i,[])*dfydt(i+di,[])))];
    theta       = [theta, theta(end)+integrando(end)];
end


% Polyfit
p = polyfit(I,full(theta(1:end-1)),20);
theta_est = polyval(p,0:di:2.5*pi);

% fit
% fo = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,0],...
%                'Upper',[Inf,max(I)],...
%                'StartPoint',[1 1]);
% ft = fittype('a*(x-b)^n','problem','n','options',fo);
% % [theta_est,gof2] = fit(I',full(theta(1:end-1))',ft,'problem',2);
% [theta_est,gof2] = fit(I',full(theta(1:end-1))','poly9');
% theta_est = theta_est_prms.

figure; hold on; grid on;
plot(I, full(-theta(1:end-1)),'b')
% plot(I, theta_est,'r')
plot(I, -theta_est,'r')
% plot(theta_est,'r')

