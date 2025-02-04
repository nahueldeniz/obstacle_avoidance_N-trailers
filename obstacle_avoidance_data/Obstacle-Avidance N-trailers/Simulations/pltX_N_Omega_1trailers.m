plt3d       = true;
pltMono     = false;
%
coordX0     = 2*S.config.Nt+2;
coordY0     = 2*S.config.Nt+3;
coordTheta0 = S.config.Nt+1;
%
XXred       = Polyhedron('lb',[S.config.stab.XxMin;S.config.stab.XyMin;S.config.stab.XthetaMin],'ub',[S.config.stab.XxMax;S.config.stab.XyMax;S.config.stab.XthetaMax]);
xx          = Polyhedron('lb',[-0.01;-0.01;-3*pi/180],'ub',[0.01;0.01;3*pi/180]);
XX          = XXred+xx;
setTgt      = Polyhedron('lb',[S.config.stab.xMin;S.config.stab.yMin;S.config.stab.thetaMin],'ub',[S.config.stab.xMax;S.config.stab.yMax;S.config.stab.thetaMax]);

XXplt       = Polyhedron('lb',[S.config.stab.XxMin;S.config.stab.XyMin;S.config.stab.XthetaMin],'ub',[S.config.stab.XxMax;S.config.stab.XyMax;S.config.stab.XthetaMax]);
setTgtPlt   = Polyhedron('lb',[S.config.stab.xMin;S.config.stab.yMin;S.config.stab.thetaMin],'ub',[S.config.stab.xMax;S.config.stab.yMax;S.config.stab.thetaMax]);

labelFontSize = 40;
axFontSize    = 20;
gridalp       = 0.4;

for j=3%1:length(S.stability.theta1Range)
    iInindx     = [];
    iOutindx    = [];
    trajIn      = [];
    trajOut     = [];
    for i=j:length(S.stability.theta1Range):S.config.stab.nx0Stab    
        cond1 = sum(setTgt.contains([S.trajStab{i}([coordX0,coordY0,coordTheta0],:)]));        % Is the prediction Nc sampled ahed inside the Inner (target) set?
        cond2 = S.config.Nc+1;%sum(XX.contains([S.trajStab{(j-1)*S.config.stab.nx0Stab+i}([coordX0,coordY0,coordTheta0],1:end)]));     % Is every predicted point inside the Oouter set?
        cond3 = sum(sum(S.trajStab{i}(:,2:end)==0));                                        % If every predicted point is ZERO, the solver failed to solve the optimisation problem...
        % Constraints verification
        cond4 = sum(S.UStab{i}(1,:)>S.mpc.mpcCasadi.u_ub(1)*1.01);
        cond5 = sum(S.UStab{i}(1,:)<S.mpc.mpcCasadi.u_lb(1)*1.01);
        cond6 = sum(S.UStab{i}(2,:)>S.mpc.mpcCasadi.u_ub(2)*1.01);
        cond7 = sum(S.UStab{i}(2,:)<S.mpc.mpcCasadi.u_lb(2)*1.01);
        dUw   = diff([0 S.UStab{i}(1,:)]);
        cond8 = sum(dUw < S.mpc.mpcCasadi.du_lb(1)*1.01);
        cond9 = sum(dUw > S.mpc.mpcCasadi.du_ub(1)*1.01);
        dUv   = diff([0 S.UStab{i}(2,:)]);
        cond10 = sum(dUv < S.mpc.mpcCasadi.du_lb(2)*1.01);
        cond11 = sum(dUv > S.mpc.mpcCasadi.du_ub(2)*1.01);

        dUvx0 = diff(S.trajStab{i}(2*S.config.Nt+2,:));
        cond12 = sum(dUvx0 < S.mpc.mpcCasadi.du_lb(2)*1.01);
        cond13 = sum(dUvx0 > S.mpc.mpcCasadi.du_ub(2)*1.01);

        dUvy0 = diff(S.trajStab{i}(2*S.config.Nt+3,:));
        cond14 = sum(dUvy0 < S.mpc.mpcCasadi.du_lb(2)*1.01);
        cond15 = sum(dUvy0 > S.mpc.mpcCasadi.du_ub(2)*1.01);

        dUw0 = diff(S.trajStab{i}(S.config.Nt+1,:));
        cond16 = sum(dUw0 < S.mpc.mpcCasadi.du_lb(1)*1.01);
        cond17 = sum(dUw0 > S.mpc.mpcCasadi.du_ub(1)*1.01);

        condU = cond4||cond5||cond6||cond7||cond8||cond9||cond10||cond11||cond12||cond13||cond14||cond15||cond16||cond17;
%         condU = cond12||cond13||cond14||cond15||cond16||cond17;

        if (cond1>=1) && (cond2==(S.config.Nc+1)) && ~(cond3==S.system.nq*S.config.Nc) && ~condU
            iInindx = [iInindx, i];
            trajIn  = [trajIn; [S.trajStab{i}(:,1);1]' ];
        else
            iOutindx = [iOutindx, i];
            trajOut  = [trajOut; [S.trajStab{i}(:,1);0]' ];
        end

        i
    end
    %
    if ~isempty(trajIn)
        trajInOut               = [trajIn;trajOut];
        X_N_Omega_Iprojection   = Polyhedron(trajIn(:,[coordX0,coordY0,coordTheta0]));
        X_N_Omega_Iprojection   = X_N_Omega_Iprojection.minVRep;
        %
        figure; hold on; grid on; xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',labelFontSize);ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',labelFontSize);
        zlabel({'$\theta_{0,0}\,(rad)$'},'interpreter','latex','fontsize',labelFontSize); %title(['Nc= ',num2str(S.config.Nc),', Vol= ',num2str(X_N_Omega_Iprojection.volume),',','$\beta_1=$',num2str(S.stability.theta1Range(j))])                
        ax = gca; ax.GridAlpha = gridalp; ax.Box='on'; ax.FontSize = axFontSize;
        if plt3d
            xnHandle    = X_N_Omega_Iprojection.plot('color','b');
            alpha(xnHandle,0.15);
%             xxhandle    = XXplt.plot('color','r');
%             alpha(xxhandle, 0.05);            
%             for i=iInindx(1:100:end)
%                 sphHandle1 = plotSpheres(S.trajStab{i}(coordX0,1),S.trajStab{i}(coordY0,1),S.trajStab{i}(coordTheta0,1),3.0,'b');
%                 alpha(sphHandle1, 0.17);
% %                 sphHandle2 = plot3(S.trajStab{i}(coordX0,1),S.trajStab{i}(coordY0,1),S.trajStab{i}(coordTheta0,1),'bo','markersize',3);
%                 
%                 sphHandle3 = plot3(S.trajStab{i}(coordX0,:),S.trajStab{i}(coordY0,:),S.trajStab{i}(coordTheta0,:),'c');
% %                 sphInHandle = plotSpheres(S.trajStab{i}(coordX0,end),S.trajStab{i}(coordY0,end),S.trajStab{i}(coordTheta0,end),0.65,'r');
%             end
            tgthandle   = setTgtPlt.plot('color','r');
            alpha(tgthandle, 1);
        else
            for i=iInindx
                if ~pltMono
                    plot(S.trajStab{i}(coordX0,1),S.trajStab{i}(coordY0,1),'bo');
                    plot(1,S.trajStab{i}(coordTheta0,1),'ro');
                    plot(-1,S.trajStab{i}(1,1),'mo');
                    plot(-0.5,S.trajStab{i}(S.config.Nt+2,1),'co');
                    plot(S.trajStab{i}(2*S.config.Nt+4,1),S.trajStab{i}(2*S.config.Nt+5,1),'go');
                    %
                    xlim([S.config.stab.XxMin S.config.stab.XxMax]);ylim([S.config.stab.XyMin S.config.stab.XyMax]);
                else
                    plot_mono2(S, S.trajStab{i}(:,1));
                end
            end        
        end
        xlim([S.config.stab.XxMin-0.5 S.config.stab.XxMax+0.5]);ylim([S.config.stab.XyMin-0.5 S.config.stab.XyMax+0.5]); daspect([1 1 1])
        zlim([-pi pi]);zticks([-pi 0 pi])
        zticklabels([{'-\pi','0','\pi'},'interpreter','latex'])
    end
end



% figure; hold on; grid on; xlabel('x');ylabel('y');zlabel('theta');
% xlim([S.config.stab.XxMin S.config.stab.XxMax]);ylim([S.config.stab.XyMin S.config.stab.XyMax]); daspect([1 1 1])
% if plt3d    
%     XXplt.plot; alpha(0.1);
%     setTgtPlt.plot; alpha(0.1);
% end
% for i=iOutindx
%     if plt3d
%         plotSpheres(S.trajStab{i}(coordX0,1),S.trajStab{i}(coordY0,1),S.trajStab{i}(1,1),2,'b');alpha(0.1);
%         plot3(S.trajStab{i}(coordX0,:),S.trajStab{i}(coordY0,:),S.trajStab{i}(1,:),'c')
%         plotSpheres(S.trajStab{i}(coordX0,end),S.trajStab{i}(coordY0,end),S.trajStab{i}(1,end),2,'r');
%     else
%         if ~pltMono
%             plot(S.trajStab{i}(coordX0,1),S.trajStab{i}(coordY0,1),'bo');
%             plot(1,S.trajStab{i}(coordTheta0,1),'ro');
%             plot(-1,S.trajStab{i}(1,1),'mo');
%             plot(-0.5,S.trajStab{i}(S.config.Nt+2,1),'co');
%             plot(S.trajStab{i}(2*S.config.Nt+4,1),S.trajStab{i}(2*S.config.Nt+5,1),'go');
%             %
%             xlim([S.config.stab.XxMin S.config.stab.XxMax]);ylim([S.config.stab.XyMin S.config.stab.XyMax]);
%         else
%             plot_mono2(S, S.trajStab{i}(:,1));
%         end
%     end
% end




























function plot_mono2(S, qk, clr)
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
    if isempty(S.mpc.mpcCasadi.Qtraj)
        Qtraj = repmat(S.mpc.mpcCasadi.q0bar, 1, S.config.Nc);
    else
        Qtraj = S.mpc.mpcCasadi.Qtraj;
    end 
    thetas  = qk(S.config.Nt+1:2*S.config.Nt+1); 
    xy0     = qk(2*S.config.Nt+2:2*S.config.Nt+3);
    xyi     = xy0;
    for i=1:S.config.Nt
        xy_aux  = qk(2*S.config.Nt+3+(i-1)*2+1:2*S.config.Nt+3+(i)*2);
        xyi     = [xyi, xy_aux];
        %
        pltLine = [];
        if ~isempty(S.path.nearObs)%{i+1}
            for j=1:size(S.path.nearObs{i+1},1)
                xy      = S.path.nearObs{i+1}(j,1:2)';
                pltLine = [pltLine, [xy_aux, xy]];
            end
            if ~isempty(pltLine)
                if i==S.config.Nt
                    line(pltLine(1,:),pltLine(2,:),'color',[1 0 0]);
                else
                    line(pltLine(1,:),pltLine(2,:),'color',[0 0 1]);
                end
            end
        end        
    end
    %
    R                       = [cos(thetas(1)-pi/2), -sin(thetas(1)-pi/2); sin(thetas(1)-pi/2), cos(thetas(1)-pi/2)];
    tractorBodyPlt          = R * S.system.XYtracBody + xyi(:,1);
    tractorWheelLeftPlt     = R * S.system.XYtracWheelLeft + xyi(:,1);
    tractorWheelRightPlt    = R * S.system.XYtracWheelRight + xyi(:,1);
    tractorAxePlt           = R * S.system.XYtracAxe + xyi(:,1);
    %
    line(tractorBodyPlt(1,:),tractorBodyPlt(2,:),'color',clrTractor,'linewidth',3);
    line(tractorWheelLeftPlt(1,:),tractorWheelLeftPlt(2,:),'color',clrWheel,'linewidth',4);
    line(tractorWheelRightPlt(1,:),tractorWheelRightPlt(2,:),'color',clrWheel,'linewidth',4);
    line(tractorAxePlt(1,:),tractorAxePlt(2,:),'color',clrAxe,'linewidth',2);
    %
    for i=1:S.config.Nt
        R                       = [cos(thetas(i+1)-pi/2), -sin(thetas(i+1)-pi/2); sin(thetas(i+1)-pi/2), cos(thetas(i+1)-pi/2)];
        trailerAxePlt           = R * S.system.XYtrailerAxe + xyi(:,i+1);
        trailerWheelLeftPlt     = R * S.system.XYtrailerWheelLeft + xyi(:,i+1);
        trailerWheelRightPlt    = R * S.system.XYtrailerWheelRight + xyi(:,i+1);
        trailerLongAxePlt       = R * S.system.XYtrailerLongAxe((i-1)*2+1:i*2,:) + xyi(:,i+1);
        trailerLoadPlt          = R * S.system.XYtrailerLoad((i-1)*2+1:i*2,:) + xyi(:,i+1);        
        if i~= S.config.Nt
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
        line(trailerLoadPlt(1,:),trailerLoadPlt(2,:),'color',trailerClr,'linewidth',2)
    end
    %
    if ~noRef
        hold off;
    end    
    %
    drawnow limitrate    
end
