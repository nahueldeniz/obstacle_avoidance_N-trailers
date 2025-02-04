clear all;

dt = 0.01;
t=0:dt:2*pi;
% x = 2*cos(t);
% y = 3*sin(t);
d = 0.5;
x = (d*sqrt(2).*cos(t))./(sin(t).^2+1);
y = (d*sqrt(2).*cos(t).*sin(t))./(sin(t).^2+1);
% dx = -2*sin(t);
% dy = 3*cos(t);
dx = (2^(1/2)*d*sin(t).*(sin(t).^2 - 3))./(sin(t).^2 + 1).^2;   %-(sqrt(2)*d*sin(t))./(sin(t).^2+1) - (2*sqrt(2)*d.*cos(t).^2.*sin(t))/(sin(t).^2+1).^2;
dy = -(2^(1/2)*d*(3*sin(t).^2 - 1))./(sin(t).^2 + 1).^2;        %(d*sqrt(2).*(cos(t).^2-sin(t).^2))./(sin(t).^2+1) - (2*sqrt(2)*d.*cos(t).^2.*sin(t).^2)/(sin(t).^2+1).^2;
c = sin(t);

alpha_aux = [];
alpha = 0;%pi/2;
beta = pi/2;
tic;
for i=1:length(t)-1
%     alpha_aux = [alpha_aux, atan2( (dy(i)*dx(i+1) - dx(i)*dy(i+1)), (dx(i)*dx(i+1) + dy(i)*dy(i+1)) )];
    alpha_aux = [alpha_aux, c(i)*dt];
    %     alpha_aux = [alpha_aux, full(ref_theta(t(i)))];
%     alpha_aux = [alpha_aux, dt*(x(i)*x(i+1)+y(i)*y(i+1)-1)*pi/2];    
%     alpha = [alpha, alpha(1) - sum(alpha_aux(1:end))];    
    alpha = [alpha, sum(alpha_aux(1:end))];    
    %
%     Dx = x(i+1)-x(i);
%     Dy = y(i+1)-y(i);
%     beta = [beta, atan2(Dy,Dx)];
%     if (beta(end)-beta(end-1)>pi)
%         beta(end) = beta(end)-2*pi;
%     elseif (beta(end)-beta(end-1)< -pi)
%         beta(end) = beta(end)+2*pi;
%     end
end
t1 = toc

figure; grid on; hold on;
plot(t,alpha,'b');
% plot(t(1:end),beta,'r-.');
% plot3(x,y,alpha,'b');plot3(x,y,beta,'r')
%
s      = casadi.MX.sym('s');
% ds     = 0.01;
% %
% % fx     = (d*sqrt(2)*cos(s))/(sin(s)^2+1);
% % fy     = (d*sqrt(2)*cos(s)*sin(s))/(sin(s)^2 + 1);
% % ftheta = atan2(  ((d*sqrt(2).*cos(s+ds).*sin(s+ds))./(sin(s+ds).^2+0.4) - (a*sqrt(2).*cos(s).*sin(s))./(sin(s).^2+0.4)),...
% %                         ((a*sqrt(2).*cos(s+ds))./(sin(s+ds).^2+1) - (a*sqrt(2).*cos(s))./(sin(s).^2+1)) );
% 
% dfx_s    = (d*sqrt(2)*sin(s)*(sin(s)^2-3))/(sin(s)^2+1)^2;
% dfy_s    = -(d*sqrt(2)*(3*sin(s)^2-1))/(sin(s)^2+1)^2;
% dfx_sp1  = (d*sqrt(2)*sin(s+ds)*(sin(s+ds)^2-3))/(sin(s+ds)^2+1)^2;
% dfy_sp1  = -(d*sqrt(2)*(3*sin(s+ds)^2-1))/(sin(s+ds)^2+1)^2;

ft  = sin(s);%atan((dfy_s*dfx_sp1 - dfx_s*dfy_sp1) / (dfx_s*dfx_sp1+dfy_s*dfy_sp1));
ode     = struct('x',s,'ode',ft);
INT     = casadi.integrator('T','idas',ode,struct('t0',0,'tf',dt));
%
Ref_theta = casadi.Function('th', {s}, {ft});
k1       = Ref_theta(s);
k2       = Ref_theta(s + dt / 2 * k1);
k3       = Ref_theta(s + dt / 2 * k2);
k4       = Ref_theta(s + dt * k3);
x_rk4    = s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
FTHETA   = casadi.Function('FTHETA', {s}, {x_rk4});
FTHETA2   = casadi.Function('FTHETA2', {s}, {s + dt*Ref_theta(s)});

alpha2 = [];
alpha0 = pi/2;
alpha3 = [];
alpha4 = [];
for i=0:dt:2*pi
    r       = INT('x0',alpha0);
    alpha2  = [alpha2, full(r.xf)];
    alpha0  = r.xf;
    %
    alpha3 = [alpha3, full(FTHETA(i))];
    %
    alpha4 = [alpha4, full(FTHETA2(i))];
end

plot(t,alpha2,'g')
plot(t,alpha3,'r')
plot(t,alpha4,'c-.')

