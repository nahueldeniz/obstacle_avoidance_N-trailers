syms a b c t real;
dt = 0.3;

x           = a*cos(t);
y           = b*sin(t);
dx = -a*sin(t);
dy = b*cos(t);

% circumnference
% num = -a^2 * cos(t) * (sin(t)*cos(dt) + cos(t)*sin(dt)) + a^2 * sin(t) * (cos(t)*cos(dt) - sin(t)*sin(dt));
% den =  a^2 * sin(t) * (sin(t)*cos(dt) + cos(t)*sin(dt)) + a^2 * cos(t) * (cos(t)*cos(dt) - sin(t)*sin(dt));

% ellipse
num = -a*b * cos(t) * (sin(t)*cos(dt) + cos(t)*sin(dt)) + a*b * sin(t) * (cos(t)*cos(dt) - sin(t)*sin(dt));
den = a^2  * sin(t) * (sin(t)*cos(dt) + cos(t)*sin(dt)) + b^2 * cos(t) * (cos(t)*cos(dt) - sin(t)*sin(dt));

alpha = atan(num/den);
alpha0 = pi/2;

ALPHA = alpha0 - int(alpha,t,'IgnoreAnalyticConstraints',true,'IgnoreSpecialCases',true)/dt
% ALPHA = alpha0 + 100*t*atan(dt);

%% Now, verify this numerically
clear all;
dt = 0.3;
t = 0:dt:2*pi;
a = 4;
x = a*cos(t);
y = a*sin(t);
alpha0 = atan2( (a*sin(t(1)+dt) - a*sin(t(1))), (a*cos(t(1)+dt) - a*cos(t(1))) );
theta = alpha0 + (1/dt)*t.*atan(dt); 
 
figure; hold on;
plot3(x,y,theta)