% Script for plotting results from simulations
disObsPerSegment    = NaN(S.config.Nt+1,length(S.data.mhempc.performance.obsPos{1}));
disPathPerSegment   = NaN(S.config.Nt+1,length(S.data.mhempc.performance.obsPos{1}));

for i=1:length(S.data.mhempc.performance.obsPos{1})
    for j=1:S.config.Nt+1
        posIthSegment = S.data.mhempc.performance.xest{1}(2*S.config.Nt+1+(j-1)*2+1:2*S.config.Nt+1+j*2,i);
        disObsAux = zeros(size(S.data.mhempc.performance.obsPos{1}{1},1),1);
        for k=1:size(S.data.mhempc.performance.obsPos{1}{1},1)
            disObsAux(k) = norm(S.data.mhempc.performance.obsPos{1}{1}(k,1:2)'-posIthSegment)-S.data.mhempc.performance.obsPos{1}{1}(k,3);
        end
        disObsPerSegment(j,i) = min(disObsAux);        
        distances = sqrt((S.path.coordinates(1,:)-posIthSegment(1)).^2+(S.path.coordinates(2,:)-posIthSegment(2)).^2);
%         distances = norm(posIthSegment-S.data.mhempc.performance.mpcRefs{1}(1:2,i));
        disPathPerSegment(j,i) = min(distances);
    end
end

% Config. Figures:
labelFontSizeBP = 9;
labelFontSize   = 17;
axFontSize      = 15;
axFontSizeBP    = 16;
gridalp         = 0.4;
linW            = 1.5;
t               = linspace(0,length(S.data.mhempc.performance.obsPos{1})*S.config.Ts,length(S.data.mhempc.performance.obsPos{1}));

% Fig. 1: distance to path
figure; hold on; grid on;
for i=1:S.config.Nt+1
    if i==1
        clr = 'y';
    elseif i==S.config.Nt+1
        clr = 'r';
    else
        clr = 'b';
    end
    plot(t,disPathPerSegment(i,:),clr,'LineStyle','-')
end

ax = gca; ax.GridAlpha = gridalp; ax.Box='on'; ax.FontSize = axFontSize;
xlabel({'$t\,(s)$'},'interpreter','latex','fontsize',labelFontSize);
ylabel({'$\Pi$'},'interpreter','latex','fontsize',labelFontSize);
% xticks(0.1:0.1:1)
% xlim([0.1 1])
% ylim([0 0.2])
% exportgraphics(figure(1),strcat(['ce-',num2str(expers),'-sims.eps']),'contenttype','vector')


% Fig. 2: min. distance to obstacles
figure; hold on; grid on;
for i=1:S.config.Nt+1
    if i==1
        clr = 'y';
    elseif i==S.config.Nt+1
        clr = 'r';
    else
        clr = 'b';
    end
    plot(t,disObsPerSegment(i,:),clr,'LineStyle','-')
end
