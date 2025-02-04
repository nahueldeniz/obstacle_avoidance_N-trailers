clear all;

dt      = 0.01;
x       = casadi.MX.sym('x');

xdot    = -sin(x);
f1      = casadi.Function('Cos_x',{x},{xdot});

k1      = f1(x);
k2      = f1(x + dt / 2 * k1);
k3      = f1(x + dt / 2 * k2);
k4      = f1(x + dt * k3);
rk4     = x + (dt/6 * (k1 + 2*k2 + 2*k3 + k4));
F       = casadi.Function('INT', {x}, {rk4});

K1      = casadi.Function('K1',{x},{k1});
K2      = casadi.Function('K2',{x},{k2});
K3      = casadi.Function('K3',{x},{k3});
K4      = casadi.Function('K4',{x},{k4});
%
f       = [];
F1      = [];
F2      = 0;
x0      = 0.5;
X       = x0;

for i=0:dt:pi
    X  = [X, full(F(X(end)))];
end


figure; hold on;
plot(X,'b')
% plot(full(F1),'g-.')
% plot(full(F2),'r-.')