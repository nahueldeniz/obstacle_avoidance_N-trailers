clear all;
% Config what and how to plot
flgPltMono  = true;
indxPlt     = NaN(12,2);
%
dataDirBase1    = '/media/nahuel/DATA/obstacle_avoidance_data/Obstacle-Avidance N-trailers/Simulations/workspaces_original_data/';
dataDirBase     = dataDirBase1;
fileNameBase    = 'obsAvoidanceGNT-Nt5-Nc25-Np65-numMovObs2-numStatObs2-mthd-proposed';
FILENAMES       = {};
NtVal           = 1:3;
stObs           = 6;
dynObs          = 6;
methods         = {'proposed'};
NcVal           = 35;
NpVal           = 65;
Dev             = NaN(numel(NtVal)*numel(stObs),2);
minDis          = NaN(numel(NtVal)*numel(stObs),2);
ctrlEff         = NaN(numel(NtVal)*numel(stObs),2);
execTime        = NaN(numel(NtVal)*numel(stObs),2);

row     = 1;
row2    = 1;
for j=stObs
    for i=NtVal
        for k=1:length(methods)
            if i==1
                NcVal = 25;
            else
                NcVal = 35;
            end
            fileName = strcat([dataDirBase,fileNameBase,'.mat']);
            FILENAMES{(row2-1)*3+i,k} = fileName;
            try
                clear S;
                load(fileName);
                vehicleDims = [S.system.long;S.system.Lh(:,2)];
                tractorDev  = zeros(length(S.data.mhempc.performance.xest{1})-1,1);
                disObs      = zeros(length(S.data.mhempc.performance.xest{1})-1,1);
                for l=1:length(S.data.mhempc.performance.xest{1})-1                    
                    posTractor      = S.data.mhempc.performance.xest{1}(2*S.config.Nt+2:2*S.config.Nt+3,l);
                    distances       = sqrt((S.path.coordinates(1,:)-posTractor(1)).^2+(S.path.coordinates(2,:)-posTractor(2)).^2);
                    tractorDev(l)   = min(distances);
                    %
                    disObsPerSegment    = NaN(S.config.Nt+1,length(S.data.mhempc.performance.obsPos{1}));
                    for m=1:S.config.Nt+1                        
                        posIthSegment   = S.data.mhempc.performance.xest{1}(2*S.config.Nt+1+(m-1)*2+1:2*S.config.Nt+1+m*2,l);
                        attIthSegment   = S.data.mhempc.performance.xest{1}(S.config.Nt+m,l);
                        for n=1:size(S.data.mhempc.performance.obsPos{1}{1},1)
                            v                       = S.data.mhempc.performance.obsPos{1}{l}(n,1:2)-posIthSegment';
                            u                       = [cos(attIthSegment) sin(attIthSegment)];
                            projUoverV              = (dot(u,v)/norm(v)^2).*v;
                            dim                     = (S.system.width/2) + norm(projUoverV)*max((vehicleDims(m)-S.system.width)/2,0);
                            disObsPerSegment(m,n)   = norm(S.data.mhempc.performance.obsPos{1}{l}(n,1:2)'-posIthSegment)-(S.path.radii(n)+dim);%S.data.mhempc.performance.obsPos{1}{1}(n,3);  
                        end                        
                    end
                    disObs(l) = min(min(disObsPerSegment));
                end
                %
                Dev((row-1)*1+1:row*1,k)        = mean(tractorDev);
                %
                [minDis((row-1)*1+1:row*1,k),indxPlt((row2-1)*3+i,k)]     = min(disObs);
                %
                ctrlEff((row-1)*1+1:row*1,k)    = sqrt(sum(S.data.mhempc.performance.ctrl{1}(1,:).^2+S.data.mhempc.performance.ctrl{1}(2,:).^2))/length(S.data.mhempc.performance.ctrl{1});
                %
                execTime((row-1)*1+1:row*1,k)   = mean(S.exec_time.t_mpc);
                %
            catch
                fileName
                a=1;
            end
        end
        row = row+1;
    end
    row2 = row2+1;
end
% Normalise exwecution-times
execTime = execTime./max(max(execTime));
%

% *************************************************************************
% Config fontsize for plots
% *************************************************************************
axFontSize      = 20;
fontSizeLatex   = 25;
gridalp         = 0.4;
% *************************************************************************
% plot the collisions
% *************************************************************************
if flgPltMono
    for k=1%:2
        for i=1:size(minDis,1)
            if minDis(i,k)<0
                load(FILENAMES{i,k});
                %
                qk = S.data.mhempc.performance.xest{1}(:,indxPlt(i,k));
                figure; hold on; grid on; daspect([1 1 1]);
                plot_mono2(S, qk, [], S.path.listOfObsStr{indxPlt(i,k)});
                % save file
                set(gca,'FontSize',axFontSize);
                ax = gca; ax.GridAlpha = gridalp; ax.Box='on'; ax.FontSize = axFontSize;
                figName = strcat(['collision_',methods{k},'_Nt',num2str(S.config.Nt),'_numObs',num2str(S.config.maxStaticObs),'.eps']);
                xlabel({'$x\textnormal{(m)}$'},'interpreter','latex','fontsize',fontSizeLatex);
                ylabel({'$y\textnormal{(m)}$'},'interpreter','latex','fontsize',fontSizeLatex);
                exportgraphics(ax,figName,'contenttype','vector')
            end
        end    
    end
end
% *************************************************************************
% plot min deviation
% *************************************************************************
figure; hold on; grid on;
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      0.55, 0.60, 0.79;...
      0.90, 0.70, 0.30]; 
group_names = {'NMPPFC'};
condition_names = {'10','11','12'};
group_inx = 1:2;
h = dabarplot(Dev','errorbars','SD','groups',group_inx,...
    'xtlabels', condition_names,'legend',group_names,'color',c,...
    'errorhats',0);
ylabel({'$\Delta\,\textnormal{(m)}$'},'interpreter','latex','fontsize',fontSizeLatex);
xl = xlim; xlim([xl(1), xl(2)]);%+1]);  % make more space for the legend
set(gca,'FontSize',axFontSize);
h.lg.Location='northwest';
ax = gca; ax.GridAlpha = gridalp; ax.Box='on'; ax.FontSize = axFontSize;
exportgraphics(ax,strcat('meanDev_12obs.eps'),'contenttype','vector')
% *************************************************************************
% plot min distance to obstacles
% *************************************************************************
figure; hold on; grid on;
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      0.55, 0.60, 0.79;...
      0.90, 0.70, 0.30]; 
group_names = {'NMPPFC'};
condition_names = {'10','11','12'};
group_inx = 1:2;
h = dabarplot(minDis','errorbars','SD','groups',group_inx,...
    'xtlabels', condition_names,'legend',group_names,'color',c,...
    'errorhats',0);
ylabel({'$\nabla\,\textnormal{(m)}$'},'interpreter','latex','fontsize',fontSizeLatex);
xl = xlim; xlim([xl(1), xl(2)]);%+1]);  % make more space for the legend
set(gca,'FontSize',axFontSize);
h.lg.Location='northwest';
ax = gca; ax.GridAlpha = gridalp; ax.Box='on'; ax.FontSize = axFontSize;
exportgraphics(ax,strcat('mindDis_12obs.eps'),'contenttype','vector')
% *************************************************************************
% plot control effort
% *************************************************************************
figure; hold on; grid on;
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      0.55, 0.60, 0.79;...
      0.90, 0.70, 0.30]; 
group_names = {'NMPPFC'};
condition_names = {'10','11','12'};
group_inx = 1:2;
h = dabarplot(ctrlEff','errorbars','SD','groups',group_inx,...
    'xtlabels', condition_names,'legend',group_names,'color',c,...
    'errorhats',0);
ylabel({'$\Psi\,\sqrt{(m^2+rad^2)/s^2}$'},'interpreter','latex','fontsize',fontSizeLatex);
xl = xlim; xlim([xl(1), xl(2)]);%+1]);  % make more space for the legend
set(gca,'FontSize',axFontSize);
h.lg.Location='northwest';
ax = gca; ax.GridAlpha = gridalp; ax.Box='on'; ax.FontSize = axFontSize;
exportgraphics(ax,strcat('ctrlEffort_12obs.eps'),'contenttype','vector')
% *************************************************************************
% plot mean execution time
% *************************************************************************
figure; hold on; grid on;
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      0.55, 0.60, 0.79;...
      0.90, 0.70, 0.30]; 
group_names = {'NMPPFC'};
condition_names = {'10','11','12'};
group_inx = 1:2;
h = dabarplot(execTime','errorbars','SD','groups',group_inx,...
    'xtlabels', condition_names,'legend',group_names,'color',c,...
    'errorhats',0);
ylabel({'$\eta\,\textnormal{(s)}$'},'interpreter','latex','fontsize',fontSizeLatex);
xl = xlim; xlim([xl(1), xl(2)]);%+1]);  % make more space for the legend
set(gca,'FontSize',axFontSize);
h.lg.Location='northwest';
ax = gca; ax.GridAlpha = gridalp; ax.Box='on'; ax.FontSize = axFontSize;
exportgraphics(ax,strcat('execTimev_12obs.eps'),'contenttype','vector')
% *************************************************************************








%
%
function plot_mono2(S, qk, clr,listOfObsAtK)
    if nargin == 3
        clrTractor          = clr;
        clrWheel            = clr;
        clrAxe              = clr;
        clrLongAxe          = clr;
        clrTrailerLoad      = clr;
        clrLastTrailerLoad  = clr;
        noRef               = true;
    else
        clrTractor          = 'y';
        clrWheel            = 'k';
        clrAxe              = 'k';
        clrLongAxe          = 'k';
        clrTrailerLoad      = 'b';
        clrLastTrailerLoad  = 'r';
        noRef               = false;
    end
    %
    coordinates         = S.path.coordinates;
    Nt                  = S.config.Nt;
    listOfObs           = listOfObsAtK;
    safeMargin          = S.path.safeMargin;
    vehicleDims         = S.path.vehicleDims;
    segmentTosteer      = S.config.segmentTosteer;
    ref                 = S.controller.ref;
    xReachable          = ref.xReachable;
    yReachable          = ref.yReachable;
    nearObs             = S.path.nearObs;
    flg                 = S.path.occluded.flg;
    p1                  = S.path.occluded.p1;
    p2                  = S.path.occluded.p2;
    p3                  = S.path.occluded.p3;
    XYtrailerAxe        = S.system.XYtrailerAxe;
    XYtrailerWheelLeft  = S.system.XYtrailerWheelLeft;
    XYtrailerWheelRight = S.system.XYtrailerWheelRight;
    XYtrailerLongAxe    = S.system.XYtrailerLongAxe;
    

    if ~noRef        
        % Plot the path
        plot(coordinates(1,:),coordinates(2,:),'color',[0.7 0.7 0.7],'LineWidth',8); grid on; daspect([1 1 1]); hold on;
        %
        if~isempty(listOfObs)
            nroObs = size(listOfObs,1);
            for i=1:nroObs/2
                outter1 = circles(listOfObs(i,1),listOfObs(i,2),listOfObs(i,3),'color',[6 18 131]./255);
                alpha(outter1,0.2)
                outter2 = circles(listOfObs(i,1),listOfObs(i,2),listOfObs(i,3)-safeMargin,'color',[19 141 144]./255);
                alpha(outter2,0.4)
                circles(listOfObs(i,1),listOfObs(i,2),listOfObs(i,3)-safeMargin-vehicleDims,'color',[253 60 60]./255);
            end
            for i=nroObs/2:nroObs
                outter1 = circles(listOfObs(i,1),listOfObs(i,2),listOfObs(i,3),'color',[0 207 250]./255);
                alpha(outter1,0.2)
                outter2 = circles(listOfObs(i,1),listOfObs(i,2),listOfObs(i,3)-safeMargin,'color',[76 63 84]./255);
                alpha(outter2,0.4)
                circles(listOfObs(i,1),listOfObs(i,2),listOfObs(i,3)-safeMargin-vehicleDims,'color',[255 0 56]./255);
            end
        end
        %
    end    
    % Plot the current reference
    if segmentTosteer == 0
        clr         = 'y';
        clrReach    = [0.95 0.95 0.1];
    elseif segmentTosteer == Nt
        clr         = 'r';
        clrReach    = [1 0.1 0.1];
    else
        clr = 'b';
        clrReach    = [0.1 0.1 1];
    end
    plot(ref.x,ref.y,clr,'marker','+','linewidth',2,'markersize',15)
    plot(ref.x,ref.y,clr,'marker','o','linewidth',2,'markersize',15)
    % plot the reachable reference
    plot(xReachable,yReachable,'color',clrReach,'marker','+','linewidth',2,'markersize',15)
    plot(xReachable,yReachable,'color',clrReach,'marker','o','linewidth',2,'markersize',15)
    %
    betas   = qk(1:Nt);    
    thetas  = qk(Nt+1:2*Nt+1); 
    xy0     = qk(2*Nt+2:2*Nt+3);
    xyi     = zeros(2,Nt+1);
    xyi(:,1) = xy0;
    for i=1:Nt
        xy_aux      = qk(2*Nt+3+(i-1)*2+1:2*Nt+3+(i)*2);
        xyi(:,i+1)  = xy_aux;
    end
    % plot occlusion's triangle
    if flg
        T = delaunayTriangulation([p1;p2;p3]);
        if ~isempty(T.ConnectivityList)
            hold on;
            triplot(T,'color',[72 151 216]./255);
            hold off;
        end
    end 
    %
    R                       = [cos(thetas(1)-pi/2), -sin(thetas(1)-pi/2); sin(thetas(1)-pi/2), cos(thetas(1)-pi/2)];
    Rh                      = [cos(-thetas(1)-pi/2), -sin(-thetas(1)-pi/2); sin(-thetas(1)-pi/2), cos(-thetas(1)-pi/2)];

    tractorBodyPlt          = S.system.XYtracBody*Rh + xyi(:,1);
    hold on;
    tractorBodyPlt.plot('color','y');
    hold off;

    tractorWheelLeftPlt     = R * S.system.XYtracWheelLeft + xyi(:,1);
    tractorWheelRightPlt    = R * S.system.XYtracWheelRight + xyi(:,1);
    tractorAxePlt           = R * S.system.XYtracAxe + xyi(:,1);
    %
    line(tractorWheelLeftPlt(1,:),tractorWheelLeftPlt(2,:),'color',clrWheel,'linewidth',4);
    line(tractorWheelRightPlt(1,:),tractorWheelRightPlt(2,:),'color',clrWheel,'linewidth',4);
    line(tractorAxePlt(1,:),tractorAxePlt(2,:),'color',clrAxe,'linewidth',2);
    %
    for i=1:Nt
        R                       = [cos(thetas(i+1)-pi/2), -sin(thetas(i+1)-pi/2); sin(thetas(i+1)-pi/2), cos(thetas(i+1)-pi/2)];
        Rtr                     = [cos(-thetas(i+1)+pi/2), -sin(-thetas(i+1)+pi/2); sin(-thetas(i+1)+pi/2), cos(-thetas(i+1)+pi/2)];
        trailerAxePlt           = R * XYtrailerAxe + xyi(:,i+1);
        trailerWheelLeftPlt     = R * XYtrailerWheelLeft + xyi(:,i+1);
        trailerWheelRightPlt    = R * XYtrailerWheelRight + xyi(:,i+1);
        trailerLongAxePlt       = R * XYtrailerLongAxe((i-1)*2+1:i*2,:) + xyi(:,i+1);
        trailerLoadPlt          = S.system.XYtrailerLoad{i}*Rtr + xyi(:,i+1);
        hold on;
        if i==Nt
            trailerLoadPlt.plot('color','r');
        else
            trailerLoadPlt.plot('color','b');
        end

        if i~= Nt
            trailerClr = clrTrailerLoad;
        else
            trailerClr = clrLastTrailerLoad;
        end
        circles(trailerLongAxePlt(1,2),trailerLongAxePlt(2,2),S.system.r/2,'edgecolor',clrLongAxe,'facecolor','none')
        line([trailerLongAxePlt(1,2),xyi(1,i)],[trailerLongAxePlt(2,2),xyi(2,i)],'color',clrLongAxe,'linewidth',1)    
        line(trailerAxePlt(1,:),trailerAxePlt(2,:),'color',clrAxe,'linewidth',2)
        line(trailerWheelLeftPlt(1,:),trailerWheelLeftPlt(2,:),'color',clrWheel,'linewidth',4)
        line(trailerWheelRightPlt(1,:),trailerWheelRightPlt(2,:),'color',clrWheel,'linewidth',4)
        line(trailerLongAxePlt(1,:),trailerLongAxePlt(2,:),'color',clrLongAxe,'linewidth',2)
    end
    %
    if ~noRef
        hold off;
    end
    xlim([min(coordinates(1,:))-2 max(coordinates(1,:))+2]); ylim([min(coordinates(2,:))-4 max(coordinates(2,:))+4]);
    drawnow limitrate    
end