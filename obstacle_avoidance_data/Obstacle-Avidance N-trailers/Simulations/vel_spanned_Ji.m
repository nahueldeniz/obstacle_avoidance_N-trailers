clear all;

Lhi         = 0.35;
Li          = 0.5;
v_min       = -1;
v_max       = 1;
w_min       = -1;
w_max       = 1;
Ts          = 0.3;
Nt          = 6;

range_beta  = -pi/2:1*pi/180:pi/2;
P           = {};
clrs        = ['b','r'];

P0          = polytope([w_min v_max; w_max v_max; w_max v_min; w_min v_min]);
volPoly     = zeros(1,Nt+1);
volPoly(1)  = volume(P0);
I           = P0;
Pi          = {P0};
wv          = extreme(I);
w0          = wv(:,1);
v0          = wv(:,2);

figure; hold on;
plot(P0,'b')

for nt=1:Nt    
    for i=1:length(range_beta)          
        ui       = zeros(2,length(w0)*length(v0));
        uiminus1 = zeros(2,length(w0)*length(v0));
        for k=1:length(v0)
            for l=1:length(w0)
                uiminus1(:,(k-1)*length(v0)+l) = [w0(l); v0(k)];
                Ji          = [-Lhi*cos(range_beta(i))/Li, sin(range_beta(i))/Li; Lhi*sin(range_beta(i)), cos(range_beta(i)) ];                
                ui(:,(k-1)*length(v0)+l) = Ji * uiminus1(:,(k-1)*length(v0)+l);
            end
        end        
        ui      = unique(ui','rows');
        Paux    = reduce(polytope(ui));
        P       = [P(:)' {Paux}];
    end    
    I = P{1};
    for i=2:length(P)
        tic
        I = intersect(I,P{i});
        I = reduce(I);
        [length(P), i, toc]
    end
    wv      = extreme(I);
    wv      = wv(1:3:end,:);
    w0      = wv(:,1);
    v0      = wv(:,2);
    volPoly(nt+1) = volume(I);
    Pi      = [Pi(:)' {I}];
    %
    P       = {};    
    plot(I,clrs(mod(nt,2)+1))
end
xlabel('$\omega_i\,(rad/s)$','interpreter','latex','fontsize',15);
ylabel('$v_i\,(m/s)$','interpreter','latex','fontsize',15);
daspect([1 1 1]);

%%
Ts = 0.3;
aw = -1:0.1:1;
av = -1:0.1:1;
Ns = 1:20;

[AW, NS] = meshgrid(aw, Ns);
[AV, NS] = meshgrid(av, Ns);

alpha = AW.*(Ts^2).*NS;
dista = AV.*(Ts^2).*NS;


figure;
subplot(1,2,1); surf(AW,NS,alpha); 
xlabel('$a_{\omega}$','Interpreter','latex');
ylabel('$N_s$','Interpreter','latex');
zlabel('$\alpha$','Interpreter','latex');

subplot(1,2,2); surf(AV,NS,dista); 
xlabel('$a_{v}$','Interpreter','latex');
ylabel('$N_s$','Interpreter','latex');
zlabel('$d$','Interpreter','latex');


%%
qx      = 1;
qy      = 1;
qtheta  = 1;
x       = -10:0.5:10;
y       = -10:0.5:10;
theta   = -pi:5*pi/180:pi;

[X,Y,THETA] = ndgrid(x,y,theta);

J = (X.^2)*qx + (Y.^2)*qy + (THETA.^2)*qtheta;

figure; hold on;

for i=1:size(X,3)
    surf(X(:,:,i),Y(:,:,i),THETA(:,:,i),J(:,:,i))
end
% alpha(0.5)
view(45,45)



