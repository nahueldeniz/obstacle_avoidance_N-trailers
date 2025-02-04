% parametric Cartesian equations for free, based on the parametric 
% equations of the Lamé curve:

t = 0:0.01:2*pi;
p = 5
q = 5;

x = p.*(sqrt(cos(t).*cos(t)).*cos(t) + sqrt(sin(t).*sin(t)).*sin(t));
y = q.*(sqrt(cos(t).*cos(t)).*cos(t) - sqrt(sin(t).*sin(t)).*sin(t));

figure; hold on;
for i=1:length(t)
    plot(x(i), y(i),'ro')
    pause(10/length(t))
end
%%
syms x y t

x = p.*(sqrt(cos(t).*cos(t)).*cos(t) + sqrt(sin(t).*sin(t)).*sin(t));
y = q.*(sqrt(cos(t).*cos(t)).*cos(t) - sqrt(sin(t).*sin(t)).*sin(t));

