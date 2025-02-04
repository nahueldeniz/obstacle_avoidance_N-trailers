clear all;
%
traj_field_exp = 'flinfit'; % 'infty', 'rect', 'flinfit'
obsMethod = 'gauss'; % 'barrier', 'gauss', 'hard', ''
%
dirBase = '/media/nahuel/DATA/LIDAR/';
%
if strcmp(traj_field_exp, 'infty')
    obs             = 0;
    if obs ==0
        num_exps    = [1:3]; % no obs
    elseif obs==1
        num_exps    = [2,4];   % obs    
    end    
elseif strcmp(traj_field_exp, 'flinfit')
    obs             = 2;
    if obs==1
        num_exps    = [5:7];   % obs    
    elseif obs==2
        num_exps    = 8:11;   % obs 
    else
        num_exps    = 1:4;
    end
elseif strcmp(traj_field_exp, 'rect')
    obs             = 2;
    if obs==1
        num_exps    = [1:3,5];   % obs    
    elseif obs==2
        num_exps    = 1:5;   % obs 
    end
end
%
pltNomPath          = true;
pltEst              = true; % plt estimated position or measured
%
LAT0                = -33.034115;
LON0                = -71.592205;
%
devRefGlobal_0      = [];
devRefGlobal_N      = [];
devPathGlobal_0     = [];
devPathGlobal_N     = [];
devRefPathGlobal_0  = [];
devRefPathGlobal_N  = [];
%
figure(1); grid on; hold on;
%
% Find osbatcles for plotting purposes
ii          = num_exps(end);
fileName    = strcat([dirBase,traj_field_exp,'-', num2str(ii),'-obs-',num2str(obs),'-',obsMethod,'.mat']);
load(fileName);
globalObs   = [];
for i=1:length(S.sensors.vlp16)
    S.obs.ptCloud                       = S.obs.ptCloud(end);
    xy0                                 = S.data.xest(2*S.config.N+2:2*S.config.N+3,i);
    theta0                              = S.data.xest(S.config.N+1,i);
    % Remove ground
    pc                                  = S.sensors.vlp16{i}.Location;
    pc                                  = pointCloud(pc);
    [~,~,outlierIndicesW,~]             = pcfitplane(pc,S.obs.filterGround.disDev,[0 0 1],S.obs.filterGround.angleDev);
    S.obs.ptCloudAux                    = select(pc,outlierIndicesW);
    % Select ROI with LiDAR at centre
    S.obs.xyz                           = S.obs.ptCloudAux.Location(abs(S.obs.ptCloudAux.Location(:,1))<=S.config.lidarLimits.X & S.obs.ptCloudAux.Location(:,2)<=S.config.lidarLimits.Y & S.obs.ptCloudAux.Location(:,3)<=S.config.lidarLimits.Z,:);
    S.obs.ptCloudAux                    = pointCloud(S.obs.xyz);
    [labels,numClusters]                = pcsegdist(S.obs.ptCloudAux,S.obs.pcSegdist.minDistance,'NumClusterPoints',[S.obs.pcSegdist.minNroPoints,S.obs.pcSegdist.maxNroPoints]);
    S.obs.localObstacles                = findObstacles(S.obs.ptCloudAux, xy0, labels, numClusters, 0); % last argument passed as a constant
    S.obs.globalObstacles               = localToGlobalObs(S,S.obs.localObstacles,theta0,xy0);
    %
    globalObs = [globalObs; S.obs.globalObstacles];
end

[idx,obstacles] = kmeans(globalObs,2);
% Plot the nominal path
clear alpha;
plot(S.path.coordinates(1,:),S.path.coordinates(2,:),'k-.','LineWidth',2);
% Plot the obstacles
if~isempty(obstacles)
    nroObs = size(obstacles,1);
    for i=1:nroObs
        r   = sqrt(obstacles(i,3)^2+obstacles(i,4)^2);
        if r<1
            rs  = r + S.system.b/2 + 0.3;
            circles(obstacles(i,1),obstacles(i,2),rs,'color','blue','edgecolor','none'); alpha(0.1);
        end
    end
    for i=1:nroObs
        r   = sqrt(obstacles(i,3)^2+obstacles(i,4)^2);
        if r<1
            circles(obstacles(i,1),obstacles(i,2),r,'color','red','edgecolor','red'); %alpha(1);
        end
    end
end
% Plot the trajectories
for ii=num_exps
    fileName    = strcat([traj_field_exp,'-', num2str(ii),'-obs-',num2str(obs),'-',obsMethod,'.mat']);
    load(fileName);
    % #####################################################################
    for j=1:S.config.N+1
        figure(1)
        if j==1
            patchline(S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j)*2,:),'linestyle','-','edgecolor','y','linewidth',1.5,'edgealpha',0.3);
        elseif j~=S.config.N+1
            patchline(S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j)*2,:),'linestyle','-','edgecolor','b','linewidth',1.5,'edgealpha',0.3);
        else
            patchline(S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j)*2,:),'linestyle','-','edgecolor','r','linewidth',1.5,'edgealpha',0.45);        
        end        
    end  
    if ii==4
        indx = [20, ceil(length(S.data.mhempc.performance.xest{1})*0.3), ceil(length(S.data.mhempc.performance.xest{1})*0.75), ceil(length(S.data.mhempc.performance.xest{1})*0.95)];
        for i=1:length(indx)
            j=1;
            patchline(S.data.mhempc.performance.qref{indx(i)}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.qref{indx(i)}(2*S.config.N+1+(j)*2,:),'linestyle','-','edgecolor','y','linewidth',5,'edgealpha',1);
            j=2;
            patchline(S.data.mhempc.performance.qref{indx(i)}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.qref{indx(i)}(2*S.config.N+1+(j)*2,:),'linestyle','-','edgecolor','b','linewidth',5,'edgealpha',1);
            j=3;
            patchline(S.data.mhempc.performance.qref{indx(i)}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.qref{indx(i)}(2*S.config.N+1+(j)*2,:),'linestyle','-','edgecolor','r','linewidth',5,'edgealpha',1);
        end
        for i=1:length(indx)
            plot_mono2(S, S.data.mhempc.performance.xest{1}(:,indx(i))); hold on;
        end
    end

    drawnow;
end
%
for ii=4
    fileName    = strcat([traj_field_exp,'-', num2str(ii),'-obs-',num2str(obs),'-',obsMethod,'.mat']);
    load(fileName);
    % #####################################################################
    if ii==4
        indx = [20, ceil(length(S.data.mhempc.performance.xest{1})*0.3), ceil(length(S.data.mhempc.performance.xest{1})*0.75), ceil(length(S.data.mhempc.performance.xest{1})*0.95)];
        for i=1:length(indx)
            j=1;
            patchline(S.data.mhempc.performance.qref{indx(i)}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.qref{indx(i)}(2*S.config.N+1+(j)*2,:),'linestyle','-','edgecolor','y','linewidth',5,'edgealpha',1);
            j=2;
            patchline(S.data.mhempc.performance.qref{indx(i)}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.qref{indx(i)}(2*S.config.N+1+(j)*2,:),'linestyle','-','edgecolor','b','linewidth',5,'edgealpha',1);
            j=3;
            patchline(S.data.mhempc.performance.qref{indx(i)}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.qref{indx(i)}(2*S.config.N+1+(j)*2,:),'linestyle','-','edgecolor','r','linewidth',5,'edgealpha',1);
        end
        for i=1:length(indx)
            plot_mono2(S, S.data.mhempc.performance.xest{1}(:,indx(i))); hold on;
        end
    end

    drawnow;
end
%
grid on;
ax              = gca;
ax.GridAlpha    = 0.4;
box on
ax.FontSize     = 15;
ax.GridAlpha    = 0.3;
daspect([1 1 1])
xlim([0 25])
ylim([4 8])
xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',19);
ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',19);
%
exportgraphics(figure(1),'rect-obs2-gauss.eps','contenttype','vector')


%% Visualise pointCloud Data
% 
fileName = 'rect-2-obs-1-gauss.mat';
load(fileName);
%
nroPCs = size(S.sensors.vlp16,2);
% ax                      = figure(1);
% ax.Color                = 'w';
figure; grid on;

pc                      = S.sensors.vlp16{1}.Location;
pc                      = pointCloud(pc);
[~,~,outlierIndicesW,~] = pcfitplane(pc,S.obs.filterGround.disDev,[0 0 1],S.obs.filterGround.angleDev);
pc                      = select(pc,outlierIndicesW);
theta                   = S.data.xest(S.config.N+1,2);
rot = [cos(theta) sin(theta) 0; ...
      -sin(theta) cos(theta) 0; ...
               0          0  1];
trans                   = [0 0 0];
tform                   = rigid3d(rot, trans);
pcMerged                = pctransform(pc, tform);
XX = [];
YY = [];
ZZ = [];
for i=2:length(S.sensors.vlp16)
    pc                      = S.sensors.vlp16{i}.Location;
    pc                      = pointCloud(pc);

    [~,~,outlierIndicesW,~] = pcfitplane(pc,S.obs.filterGround.disDev,[0 0 1],S.obs.filterGround.angleDev);
    pc                      = select(pc,outlierIndicesW);
    %
    theta                   = -S.data.xest(S.config.N+1,i+1);
    xy_h                    = S.data.xest(2*S.config.N+2:2*S.config.N+3,i+1);
    %
    X = [];
    Y = [];
    Z = []; 
    for j=1:length(pc.Location)
        if ~isnan(pc.Location(j,1))
            R   = [cos(-pi/2-theta) -sin(-pi/2-theta); sin(-pi/2-theta) cos(-pi/2-theta)];
            xy  = R * [pc.Location(j,1);pc.Location(j,2)];
            X   = [X; xy(1)+xy_h(1)];
            Y   = [Y; xy(2)+xy_h(2)];
            Z   = [Z; pc.Location(j,3)];
            if (Y(end)>4 && Y(end)<6 && X(end)>13.5 && X(end)<25) || (Y(end)>5 && Y(end)<8.5 && X(end)>20 && X(end)<25)
                X(end) = [];
                Y(end) = [];
                Z(end) = [];
            end
        end
    end   
    XX = [XX; X];
    YY = [YY; Y];
    ZZ = [ZZ; Z];
    view([-111.475,39])
    xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',15);
    ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',15);
    zlabel({'$z\,(m)$'},'interpreter','latex','fontsize',15);
%     pcMerged    = pcmerge(pcMerged, ptCloudOut, 0.1);
    %
%     view(player,pcMerged);
%     figure(1)
%     pcshow(ptCloudOut);
    xlim([0 25]);
    ylim([-2 10]);
    zlim([-0.5 2]);
    %
    drawnow
    i
end
% obs1 = 
XXYY      = [XX,YY];
indx_moving = find(XXYY(:,1)>5 & XXYY(:,1)<14 & XXYY(:,2)>3 & XXYY(:,2)<8.5);
indx_moving = [indx_moving; [find(XXYY(:,1)>15 & XXYY(:,1)<17 & XXYY(:,2)>5.5 & XXYY(:,2)<8.5)] ];
indx_static = find(XXYY(:,1)>17.25 & XXYY(:,1)<20 & XXYY(:,2)>5 & XXYY(:,2)<8.75);

XX_noObs = XX;
YY_noObs = YY;
ZZ_noObs = ZZ;

XX_noObs([indx_moving;indx_static],:) = [];
YY_noObs([indx_moving;indx_static],:) = [];
ZZ_noObs([indx_moving;indx_static],:) = [];
%%
figure
hold on;
% scatter3(XX_noObs,YY_noObs,ZZ_noObs,'b');
% scatter3(XX(indx_static,1),YY(indx_static,1),ZZ(indx_static,1),'r');
% scatter3(XX(indx_moving(1:200),1),YY(indx_moving(1:200),1),ZZ(indx_moving(1:200),1),'r');
% scatter3(smooth(XX(indx_moving,1),200),smooth(YY(indx_moving,1),200),smooth(ZZ(indx_moving,1),200),'r');

% hold on;
% plot(S.path.coordinates(1,:),S.path.coordinates(2,:),'k-.','LineWidth',3)
% plot(S.data.mhempc.performance.xest{1}(6,:),S.data.mhempc.performance.xest{1}(7,:),'y','LineWidth',3)
% plot(S.data.mhempc.performance.xest{1}(8,:),S.data.mhempc.performance.xest{1}(9,:),'b','LineWidth',3)
% plot(S.data.mhempc.performance.xest{1}(10,:),S.data.mhempc.performance.xest{1}(11,:),'r','LineWidth',3)
% figure;
% pcshow(ptCloudOut);
grid on;
ax              = gca;
ax.Color        = 'w';
ax.GridAlpha    = 0.4;
box on
ax.FontSize     = 18;
ax.GridAlpha    = 0.3;
view([-111.475,39])
daspect([1 1 1])
xlim([0 25]);
ylim([-2 10]);
zlim([-0.5 2]);
xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',22);
ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',22);
zlabel({'$z\,(m)$'},'interpreter','latex','fontsize',22);

% exportgraphics(ax,'pc-obs1-2.eps','contenttype','vector')

pCloud1 = pointCloud([XX_noObs,YY_noObs,ZZ_noObs]);
pcshow(pCloud); xlim([0 25]); ylim([-2 10]); zlim([-0.5 2],'color','b'); view([-1.752324218750000e+02,84.9921875])


%%
% Plot times
clear ALL;

% axF1        = figure; hold on;

meanHard    = [];
meanBarr    = [];
meanGauss   = [];
meanMheHard = [];
meanMheBarr = [];
meanMheGauss = [];
maxDeMax    = 1;
maxDeMaxMhe = 1;
%
fileName    = 'infty2-obs-3-hard-10trials.mat';
load(fileName);
for i=1:S.config.NUM_SIMS
%     subplot(2,3,1); hold on;
%     plot(S.data.mhempc.performance.exec_times{i}.t_tot,'c')    
%     plot(S.data.mhempc.performance.exec_times{i}.t_mpc,'m')
%     subplot(2,3,4); hold on;
%     plot(S.data.mhempc.performance.exec_times{i}.t_mhe,'y')
    meanHard    = [meanHard, mean(S.data.mhempc.performance.exec_times{i}.t_tot)];
    meanMheHard = [meanMheHard, mean(S.data.mhempc.performance.exec_times{i}.t_mhe)];
end
%
% clear ALL;
fileName    = 'infty2-obs-3-barrier-10trials.mat';
load(fileName);
for i=1:S.config.NUM_SIMS
%     subplot(2,3,2); hold on;
%     plot(S.data.mhempc.performance.exec_times{i}.t_tot,'c')    
%     plot(S.data.mhempc.performance.exec_times{i}.t_mpc,'m')
%     subplot(2,3,5); hold on;
%     plot(S.data.mhempc.performance.exec_times{i}.t_mhe,'y')
    meanBarr    = [meanBarr, mean(S.data.mhempc.performance.exec_times{i}.t_tot)];
    meanMheBarr = [meanMheBarr, mean(S.data.mhempc.performance.exec_times{i}.t_mhe)];
end
%
% clear ALL;
fileName    = 'infty2-obs-3-soft-10trials.mat';
load(fileName);
for i=1:S.config.NUM_SIMS
%     subplot(2,3,3); hold on;
%     plot(S.data.mhempc.performance.exec_times{i}.t_tot,'c')    
%     plot(S.data.mhempc.performance.exec_times{i}.t_mpc,'m')
%     subplot(2,3,6); hold on;
%     plot(S.data.mhempc.performance.exec_times{i}.t_mhe,'y')
    meanGauss   = [meanGauss, mean(S.data.mhempc.performance.exec_times{i}.t_tot)];
    meanMheGauss = [meanMheGauss, mean(S.data.mhempc.performance.exec_times{i}.t_mhe)];
end
%
% maxDeMax = max(meanHard);
% if max(meanBarr)> maxDeMax
%     maxDeMax = max(meanBarr);
% end
% if max(meanGauss)> maxDeMax
%     maxDeMax = max(meanGauss);
% end
% %
% maxDeMaxMhe = max(meanMheHard);
% if max(meanMheBarr)> maxDeMaxMhe
%     maxDeMaxMhe = max(meanMheBarr);
% end
% if max(meanMheGauss)> maxDeMaxMhe
%     maxDeMaxMhe = max(meanMheGauss);
% end
%
meanHard        = meanHard./maxDeMax;
meanBarr        = meanBarr./maxDeMax;
meanGauss       = meanGauss./maxDeMax;
%
meanMheHard     = meanMheHard./maxDeMaxMhe;
meanMheBarr     = meanMheBarr./maxDeMaxMhe;
meanMheGauss    = meanMheGauss./maxDeMaxMhe;
%
%
%
% Config plots
fontSize        = 19;
fontSizeLatex   = 25;
%
%
%
bp1             = figure; hold on;
boxplot([meanHard; meanBarr; meanGauss]','Labels',{'Hard','Barrier','Gauss'});
ylim([0 30])
grid on;
ax              = gca;
ax.GridAlpha    = 0.4;
box on
ax.FontSize     = fontSize;
ax.GridAlpha    = 0.3;
ylabel({'$t\,(s)$'},'interpreter','latex','fontsize',fontSizeLatex);
exportgraphics(bp1,'bxp1.eps','contenttype','vector')
% close all;
% 
%
%
bp2             = figure; hold on;
boxplot([meanMheHard; meanMheBarr; meanMheGauss]','Labels',{'Hard','Barrier','Gauss'});
ylim([0.1 0.2])
grid on;
ax              = gca;
ax.GridAlpha    = 0.4;
box on
ax.FontSize     = fontSize;
ax.GridAlpha    = 0.3;
ylabel({'$t\,(s)$'},'interpreter','latex','fontsize',fontSizeLatex);
exportgraphics(bp2,'bxp2.eps','contenttype','vector')
% close all;
%
%
% Static obstacle
statobs         = figure; grid on; hold on;
circles(0,0,10,'points',4,'rotation',45,'facecolor','y'); 
circles(0,0,0.5,'facecolor','r','edgecolor','r'); 
grid on;
ax              = gca;
ax.GridAlpha    = 0.4;
box on
ax.FontSize     = fontSize-1;
ax.GridAlpha    = 0.3;
daspect([1 1 1])
xlim([-8 8])
ylim([-8 8])
xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',fontSizeLatex-1);
ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',fontSizeLatex-1);
exportgraphics(statobs,'staticobs.eps','contenttype','vector')
% close all;
%
%
% Dynamic obstacle
dynamicobs      = figure; grid on; hold on; 
circles(0,0,10,'points',4,'rotation',45,'facecolor','y'); 
circles(0:1:4,((0:1:4).^2)/10,0.5,'facecolor','r','edgecolor','r'); 
daspect([1 1 1])
grid on;
ax              = gca;
ax.GridAlpha    = 0.4;
box on
ax.FontSize     = fontSize-1;
ax.GridAlpha    = 0.3;
daspect([1 1 1])
xlim([-8 8])
ylim([-8 8])
xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',fontSizeLatex-1);
ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',fontSizeLatex-1);
exportgraphics(dynamicobs,'dynamicobs.eps','contenttype','vector')
% close all;
%
%
% Obstacle modelling
p1 = [1 1.375];
p2 = [1.5 1.75 ];
p3 = [1.6 1.62];
p4 = [1.87 1.62];
p5 = [2 1.375];
p6 = [1.75 1.15];
p7 = [1.6 1.375];
p8 = [1.5 1.15];
p9 = [1.6 1.08];
p10 = [1.25 1];

Pxy = [p1' p2' p3' p4' p5' p6' p7' p8' p9' p10' p1'];
Px = Pxy(1,:);
Py = Pxy(2,:);

obstacle = figure; hold on; grid on;
circles(1.5, 1.375, 0.9,'color','g','edgecolor','none')
alpha(0.15)
circles(1.5, 1.375, 0.8,'color','b','edgecolor','none')
alpha(0.15)
circles(1.5, 1.375, 0.5,'color','r','edgecolor','none')
alpha(0.15)
line(Px,Py,'linewidth',2,'color','k')

grid on;
ax              = gca;
ax.GridAlpha    = 0.4;
box on
ax.FontSize     = fontSize;
ax.GridAlpha    = 0.3;
daspect([1 1 1])
xlim([0.5 2.5])
ylim([0.4 2.4])
xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',fontSizeLatex);
ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',fontSizeLatex);
exportgraphics(obstacle,'obs_modelling.eps','contenttype','vector')
% close all;




%%
% plot gaussian functions for the paper
x       = -1:0.05:1;
y       = -1:0.05:1;
x0      = 0;
y0      = 0;
ro      = 0.1;
rv      = 0.2;
rs      = 0.3;
[X,Y]   = meshgrid(x,y);
%
go1     = 2.*exp((-1/(2*ro^2)) .* ((X-x0).^2 + (Y-y0).^2) );
[r,c]   = size(go1);
COgo1   = zeros(r,c,3);
COgo1(:,:,1) = 1;
%
gv1     = 2.*exp((-1/(2*rv^2)) .* ((X-x0).^2 + (Y-y0).^2) );
COgv1   = zeros(r,c,3);
COgv1(:,:,2) = 1;
%
gs1     = 2.*exp((-1/(2*rs^2)) .* ((X-x0).^2 + (Y-y0).^2) );
COgs1   = zeros(r,c,3);
COgs1(:,:,3) = 1;
%
qd1     = 0.75.*(X+0.25).^2 + 0.75.*(Y+0.25).^2;
COqd1   = zeros(r,c,3);
COqd1(:,:,1) = 0.8;
COqd1(:,:,2) = 0.8;
COqd1(:,:,3) = 0.8;
% 0.80,0.80,0.80
%
zdiff   = gs1 - qd1;
C       = contours(X, Y, zdiff, [0 0]);
xL      = C(1, 2:end);
yL      = C(2, 2:end);
zL      = interp2(X, Y, gs1, xL, yL);
%
% ax1 = figure; hold on; grid on;
% clear alpha;
% surf(X,Y,go1,COgo1,'EdgeAlpha',1); alpha(0.7);
% fontsize(ax1,33,"pixels")
% legend({'$\rho_i$'},'interpreter','latex','FontSize',30,'Location','east')
% xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',39)
% ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',39)
% view([-49.35,19.90])
% exportgraphics(ax1,'gauss1.eps','contenttype','vector')
% 
% ax2 = figure; hold on; grid on;
% surf(X,Y,gv1,COgv1,'EdgeAlpha',0.4); alpha(0.7);
% fontsize(ax2,33,"pixels")
% legend({'$\rho_v$'},'interpreter','latex','FontSize',30,'Location','east')
% xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',39)
% ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',39)
% view([-49.35,19.90])
% exportgraphics(ax2,'gauss2.eps','contenttype','vector')
% 
% ax3 = figure; hold on; grid on;
% surf(X,Y,gs1,COgs1,'EdgeAlpha',0.2); alpha(0.7);
% fontsize(ax3,33,"pixels")
% legend({'$\rho_s$'},'interpreter','latex','FontSize',30,'Location','east')
% xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',39)
% ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',39)
% view([-49.35,19.90])
% exportgraphics(ax3,'gauss3.eps','contenttype','vector')

surf(X,Y,qd1,COqd1,'EdgeAlpha',0.1); alpha(0.06);
line(xL, yL, zL, 'color', 'k', 'linewidth', 3);
gsqd1   = gs1+qd1;
clrSur  = ones([size(gsqd1),3]);
clrSur(:,:,1) = 0.3.*clrSur(:,:,1);
clrSur(:,:,2) = 0.75.*clrSur(:,:,2);
clrSur(:,:,3) = 0.93.*clrSur(:,:,3);

surf(X,Y,gsqd1,COqd1,'EdgeAlpha',0.35); alpha(0.5);
indx0   = find(X==0);
indx    = X(indx0);
indy    = Y(indx0);
zL      = gsqd1(indx0);
line(indx, indy, zL, 'color', [1.00,0.07,0.65], 'linewidth', 10);
%
indx04   = find(X>=-0.452 & X<=-0.449);
indx_aux = X(indx04);
indx2    = [];
indy_aux = Y(indx04);
zL_aux   = gsqd1(indx04);
zL2      = zL;
for ii=1:length(zL_aux)
    l       = (cos((0:length(zL_aux)-1).*2*pi./length(zL_aux))+1)./2;
    zL2(ii)  = l(ii)*zL(ii) + (1.05-l(ii))*zL_aux(ii);
    indx2   = [indx2; round(l(ii)*indx0(ii) + (1-l(ii))*indx04(ii))];
end

line(smooth(X(indx2),5), smooth(indy,1), smooth(zL2,10), 'color', [0.00,0.70,1.00], 'linewidth', 9);

ax = gca;

% exportgraphics(ax,'gaussian_functions.eps','contenttype','vector')
exportgraphics(ax,'traj-costs.eps','contenttype','vector')


%%
function plot_mono3(S, qk, clr, N, transps)
    if ~isempty(clr)
        clrTractor          = clr;
        clrWheel            = clr;
        clrAxe              = clr;
        clrLongAxe          = clr;
        clrTrailerLoad      = clr;
        clrLastTrailerLoad  = clr;
    else
        clrTractor          = 'y';
        clrWheel            = 'k';
        clrAxe              = 'k';
        clrLongAxe          = 'k';
        clrTrailerLoad      = 'b';
        clrLastTrailerLoad  = 'r';
    end
    % Rotate the Husky and wheels according to their attitude
    thetas  = qk(N+1:2*N+1);
    betas   = qk(1:N);    
    xy0     = qk(2*N+2:2*N+3);
    xyi     = xy0;
    for i=1:N
%         xyi(:,i+1) = xyi(:,i) - [S.system.Lhi(i)*cos(thetas(i)) + S.system.Li(i)*cos(thetas(i+1)); S.system.Lhi(i)*sin(thetas(i)) + S.system.Li(i)*sin(thetas(i+1))]; 
        xyi = [xyi, qk(2*N+3+(i-1)*2+1:2*N+3+(i)*2)];
    end
    %
    R                       = [cos(thetas(1)-pi/2), -sin(thetas(1)-pi/2); sin(thetas(1)-pi/2), cos(thetas(1)-pi/2)];
    tractorBodyPlt          = R * S.system.XYtracBody + xyi(:,1);
    tractorWheelLeftPlt     = R * S.system.XYtracWheelLeft + xyi(:,1);
    tractorWheelRightPlt    = R * S.system.XYtracWheelRight + xyi(:,1);
    tractorAxePlt           = R * S.system.XYtracAxe + xyi(:,1);
    %
    patchline(tractorBodyPlt(1,:),tractorBodyPlt(2,:),'edgecolor',clrTractor,'facecolor',clrTractor,'linewidth',3,'edgealpha',transps(1));
    patchline(tractorWheelLeftPlt(1,:),tractorWheelLeftPlt(2,:),'edgecolor',clrWheel,'facecolor',clrWheel,'linewidth',4,'edgealpha',transps(1));
    patchline(tractorWheelRightPlt(1,:),tractorWheelRightPlt(2,:),'edgecolor',clrWheel,'facecolor',clrWheel,'linewidth',4,'edgealpha',transps(1));
    patchline(tractorAxePlt(1,:),tractorAxePlt(2,:),'edgecolor',clrAxe,'facecolor',clrAxe,'linewidth',2,'edgealpha',transps(1));
    %
    for i=1:N
        R                       = [cos(thetas(i+1)-pi/2), -sin(thetas(i+1)-pi/2); sin(thetas(i+1)-pi/2), cos(thetas(i+1)-pi/2)];
        trailerAxePlt           = R * S.system.XYtrailerAxe + xyi(:,i+1);
        trailerWheelLeftPlt     = R * S.system.XYtrailerWheelLeft + xyi(:,i+1);
        trailerWheelRightPlt    = R * S.system.XYtrailerWheelRight + xyi(:,i+1);
        trailerLongAxePlt       = R * S.system.XYtrailerLongAxe((i-1)*2+1:i*2,:) + xyi(:,i+1);
        trailerLoadPlt          = R * S.system.XYtrailerLoad((i-1)*2+1:i*2,:) + xyi(:,i+1);        
        if i~= S.config.N
            trailerClr = clrTrailerLoad;
        else
            trailerClr = clrLastTrailerLoad;
        end
        circles(trailerLongAxePlt(1,2),trailerLongAxePlt(2,2),S.system.r/2,'edgecolor',clrLongAxe,'facecolor','none')
        patchline([trailerLongAxePlt(1,2),xyi(1,i)],[trailerLongAxePlt(2,2),xyi(2,i)],'edgecolor',clrLongAxe,'facecolor',clrLongAxe,'linewidth',1,'edgealpha',transps(i+1));
        patchline(trailerAxePlt(1,:),trailerAxePlt(2,:),'edgecolor',clrAxe,'facecolor',clrAxe,'linewidth',2,'edgealpha',transps(i+1));
        patchline(trailerWheelLeftPlt(1,:),trailerWheelLeftPlt(2,:),'edgecolor',clrWheel,'facecolor',clrWheel,'linewidth',4,'edgealpha',transps(i+1));
        patchline(trailerWheelRightPlt(1,:),trailerWheelRightPlt(2,:),'edgecolor',clrWheel,'edgecolor',clrWheel,'linewidth',4,'edgealpha',transps(i+1));
        patchline(trailerLongAxePlt(1,:),trailerLongAxePlt(2,:),'edgecolor',clrLongAxe,'facecolor',clrLongAxe,'linewidth',2,'edgealpha',transps(i+1));
        patchline(trailerLoadPlt(1,:),trailerLoadPlt(2,:),'edgecolor',trailerClr,'facecolor',trailerClr,'linewidth',2,'edgealpha',transps(i+1));
    end   
    %
    drawnow limitrate
end

function plot_mono2(S, qk, clr)
    if nargin == 3
        clrTractor          = clr;
        clrWheel            = clr;
        clrAxe              = clr;
        clrLongAxe          = clr;
        clrTrailerLoad      = clr;
        clrLastTrailerLoad  = clr;
        noRef = true;
        transp = 0.5;
    else
        clrTractor          = 'y';
        clrWheel            = 'k';
        clrAxe              = 'k';
        clrLongAxe          = 'k';
        clrTrailerLoad      = 'b';
        clrLastTrailerLoad  = 'r';
        noRef = false;
        transp = 1;
    end
    %
    if ~noRef
%         plot(S.path.coordinates(1,:),S.path.coordinates(2,:),'color',[0.6 0.6 0.6],'LineWidth',8); grid on; daspect([1 1 1]); hold on;
%         xlim([min(S.path.coordinates(1,:))-1 max(S.path.coordinates(1,:))+1]); ylim([min(S.path.coordinates(2,:))-1 max(S.path.coordinates(2,:))+2]);
        %
%         for i=1:S.config.N+1
%             if i==1
%                 plot(S.algorithms.mpcCasadi.Qtraj(2*S.config.N+1+(i-1)*2+1,:),S.algorithms.mpcCasadi.Qtraj(2*S.config.N+1+i*2,:),'y','LineWidth',4);
%             elseif i~=S.config.N+1
%                 plot(S.algorithms.mpcCasadi.Qtraj(2*S.config.N+1+(i-1)*2+1,:),S.algorithms.mpcCasadi.Qtraj(2*S.config.N+1+i*2,:),'b','LineWidth',4);
%             else
%                 plot(S.algorithms.mpcCasadi.Qtraj(2*S.config.N+1+(i-1)*2+1,:),S.algorithms.mpcCasadi.Qtraj(2*S.config.N+1+i*2,:),'r','LineWidth',4);
%             end
%         end        
%         patch(S.config.initUncertaintySet.x,S.config.initUncertaintySet.y,'r','edgecolor','none','facealpha',0.2)
    end    
%     plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+2,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+3,end),'m+','linewidth',1.5,'markersize',20)
%     plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+2,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+3,end),'mo','linewidth',1.5,'markersize',20)
%     plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+4,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+5,end),'mx','linewidth',1.5,'markersize',20)
%     plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+4,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+5,end),'mo','linewidth',1.5,'markersize',20)
    % Rotate the Husky and wheels according to their attitude
    thetas  = qk(S.config.N+1:2*S.config.N+1);
    betas   = qk(1:S.config.N);    
    xy0     = qk(2*S.config.N+2:2*S.config.N+3);
    xyi     = xy0;
    for i=1:S.config.N
%         xyi(:,i+1) = xyi(:,i) - [S.system.Lhi(i)*cos(thetas(i)) + S.system.Li(i)*cos(thetas(i+1)); S.system.Lhi(i)*sin(thetas(i)) + S.system.Li(i)*sin(thetas(i+1))]; 
        xyi = [xyi, qk(2*S.config.N+3+(i-1)*2+1:2*S.config.N+3+(i)*2)];
    end

%     xyi(:,end) = qk(2*S.config.N+4:2*S.config.N+5);

    R                       = [cos(thetas(1)-pi/2), -sin(thetas(1)-pi/2); sin(thetas(1)-pi/2), cos(thetas(1)-pi/2)];
    tractorBodyPlt          = R * S.system.XYtracBody + xyi(:,1);
    tractorWheelLeftPlt     = R * S.system.XYtracWheelLeft + xyi(:,1);
    tractorWheelRightPlt    = R * S.system.XYtracWheelRight + xyi(:,1);
    tractorAxePlt           = R * S.system.XYtracAxe + xyi(:,1);
    %
    patchline(tractorBodyPlt(1,:),tractorBodyPlt(2,:),'edgecolor',clrTractor,'linewidth',3,'edgealpha',transp);
    line(tractorWheelLeftPlt(1,:),tractorWheelLeftPlt(2,:),'color',clrWheel,'linewidth',4);
    line(tractorWheelRightPlt(1,:),tractorWheelRightPlt(2,:),'color',clrWheel,'linewidth',4);
    line(tractorAxePlt(1,:),tractorAxePlt(2,:),'color',clrAxe,'linewidth',2);
    %
    for i=1:S.config.N
        R                       = [cos(thetas(i+1)-pi/2), -sin(thetas(i+1)-pi/2); sin(thetas(i+1)-pi/2), cos(thetas(i+1)-pi/2)];
        trailerAxePlt           = R * S.system.XYtrailerAxe + xyi(:,i+1);
        trailerWheelLeftPlt     = R * S.system.XYtrailerWheelLeft + xyi(:,i+1);
        trailerWheelRightPlt    = R * S.system.XYtrailerWheelRight + xyi(:,i+1);
        trailerLongAxePlt       = R * S.system.XYtrailerLongAxe((i-1)*2+1:i*2,:) + xyi(:,i+1);
        trailerLoadPlt          = R * S.system.XYtrailerLoad((i-1)*2+1:i*2,:) + xyi(:,i+1);        
        if i~= S.config.N
            trailerClr = clrTrailerLoad;
        else
            trailerClr = clrLastTrailerLoad;
        end
        circles(trailerLongAxePlt(1,2),trailerLongAxePlt(2,2),S.system.r/2,'edgecolor',clrLongAxe,'facecolor','none')
        line([trailerLongAxePlt(1,2),xyi(1,i)],[trailerLongAxePlt(2,2),xyi(2,i)],'color',clrLongAxe,'linewidth',1);
        line(trailerAxePlt(1,:),trailerAxePlt(2,:),'color',clrAxe,'linewidth',2);
        line(trailerWheelLeftPlt(1,:),trailerWheelLeftPlt(2,:),'color',clrWheel,'linewidth',4);
        line(trailerWheelRightPlt(1,:),trailerWheelRightPlt(2,:),'color',clrWheel,'linewidth',4);
        line(trailerLongAxePlt(1,:),trailerLongAxePlt(2,:),'color',clrLongAxe,'linewidth',2);
        patchline(trailerLoadPlt(1,:),trailerLoadPlt(2,:),'edgecolor',trailerClr,'linewidth',2,'edgealpha',transp);
    end
   
    %
%     plot(S.path.references_mhempc(2*S.config.N+2,end),S.path.references_mhempc(2*S.config.N+3,end),'y+','linewidth',3,'markersize',20)
%     plot(S.path.references_mhempc(2*S.config.N+2,end),S.path.references_mhempc(2*S.config.N+3,end),'yo','linewidth',3,'markersize',20)
%     plot(S.path.references_mhempc(2*S.config.N+4,end),S.path.references_mhempc(2*S.config.N+5,end),'rx','linewidth',3,'markersize',20)
%     plot(S.path.references_mhempc(2*S.config.N+4,end),S.path.references_mhempc(2*S.config.N+5,end),'ro','linewidth',3,'markersize',20)
    %
%     plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+2,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+3,end),'m+','linewidth',1.5,'markersize',20)
%     plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+2,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+3,end),'mo','linewidth',1.5,'markersize',20)
%     plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+4,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+5,end),'mx','linewidth',1.5,'markersize',20)
%     plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+4,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.N+5,end),'mo','linewidth',1.5,'markersize',20)
    %
%     if ~noRef
%         hold off;
% 
%         if~isempty(S.path.obstacles)
%             nroObs = size(S.path.obstacles,1);
%             for i=1:nroObs
%                 rs = S.path.obstacles(i,3) + S.system.b/2 + 0.3;
%                 circles(S.path.obstacles(i,1),S.path.obstacles(i,2),rs,'color','blue','edgecolor','none'); alpha(0.1);
%             end
%             for i=1:nroObs
%                 circles(S.path.obstacles(i,1),S.path.obstacles(i,2),S.path.obstacles(i,3),'color','red','edgecolor','red'); %alpha(1);
%             end
%         end
%     end
    %
    drawnow limitrate
end

function h = stacked_bar3(array)


    if any(array(:) < 0)
        error('Only positive values supported')
    end
    
    dims = size(array);
    if any(dims==0)
        error('Empty dimensions are not supported')
    end    

    switch length(dims)
        case 2
            ns = 1;
        case 3
            ns = dims(3);
        otherwise
            error('Must be a 3D array')
    end
    nr = dims(1);
    nc = dims(2);
    
    ax = newplot;
%     co = ax.ColorOrder;    
    co = [1 1 0;0 0 1;1 0 0];
    h = gobjects(1,ns);
    view(ax,3)
    xlim(ax,[.5 nc+.5])
    ylim(ax,[.5 nr+.5])
    
    bw = .4;
    offmat = [-bw, +bw, 0; ...
              -bw, -bw, 0; ...
              +bw, -bw, 0; ...
              +bw, +bw, 0];
    sidemat = [1, 2, 2, 1; ...
               2, 3, 3, 2; ...
               3, 4, 4, 3; ...
               4, 1, 1, 4] ...
            + repmat([0, 0, 4*nr*nc, 4*nr*nc],[4, 1]);
    topmat = (1:4) + 4*nr*nc;

    top = zeros(dims(1:2));
    for s = 1:ns
        bottom = top;
        top = bottom + array(:,:,s);

        verts = zeros(4*nr*nc*2, 3);
        faces = ones(5*nr*nc, 4);
        for r = 1:nr
            for c = 1:nc
                vindex = 4*(r-1 + nr*(c-1));
                lindex = 5*(r-1 + nr*(c-1));
                rindex = 4*(r-1 + nr*(c-1));
                verts(vindex +           (1:4)', :) = repmat([c,r,bottom(r,c)],[4,1]) + offmat;
                verts(vindex + 4*nr*nc + (1:4)', :) = repmat([c,r,   top(r,c)],[4,1]) + offmat;
                faces(lindex + (1:5)',:) = rindex + [sidemat; topmat];
            end
        end
        
        cix = 1+mod(s-1, size(co,1));
        h(s) = patch('Vertices', verts, ...
                     'Faces', faces, ...
                     'FaceColor', co(cix,:), ...
                     'Parent', ax);
                 
        bottom = top;
    end
end

function obs = findObstacles(ptCloudAux, posH, labels, numClusters, indx)
obs         = [];
obsToDel    = [];
    for i=1:numClusters
        obstacle_i      = find(labels == i);
        pt_obstacle     = select(ptCloudAux,obstacle_i);
        x_mean          = mean(pt_obstacle.Location(:,1));
        y_mean          = mean(pt_obstacle.Location(:,2));
        deltaX          = (pt_obstacle.XLimits(2)-pt_obstacle.XLimits(1));
        deltaY          = (pt_obstacle.YLimits(2)-pt_obstacle.YLimits(1));
        deltaZ          = (pt_obstacle.ZLimits(2)-pt_obstacle.ZLimits(1));        
        newObs          = [x_mean, y_mean, deltaX, deltaY];
        obs             = [obs; newObs];
    end
    obsToDel        = unique(obsToDel);
    obs(obsToDel,:) = [];
end

function Obstacles = localToGlobalObs(S,local,theta0,xy_husky)
Obstacles = local;
    if isempty(local)
        return;
    else 
    obsOutOfLims = [];
        for i=1:size(local,1)
            d = sqrt(local(i,1)^2+local(i,2)^2); 
            alpha = atan(local(i,1)/local(i,2));    % x and y components are swapped since the lidar read in this way.
            xol = d * cos(theta0-alpha);
            yol = d * sin(theta0-alpha);
            xog = xy_husky(1) + xol;
            yog = xy_husky(2) + yol;
            if (xog >= S.config.x_min) && (xog <= S.config.x_max) && (yog >= S.config.y_min) && (yog <= S.config.y_max)
                Obstacles(i,1) = xog;
                Obstacles(i,2) = yog;
            else
                obsOutOfLims = [obsOutOfLims; i];
            end
        end   
        Obstacles(obsOutOfLims,:) = [];
    end
end