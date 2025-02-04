% #########################################################################
% 
% Estimation and Control for Generalised N-Trailer Systems
% Author: Nestor. N. Deniz - 2022
%
% #########################################################################notification icon


function S = husky_obsAvoidanceNew()
    clear all; clc;
    import casadi.*
    % ---------------------------------------------------------------------
    S = init();
    % Init ROS's things ---------------------------------------------------
    if S.config.dataOnline == true; figure('units','normalized','outerposition',[0 0 1 1]); hold on; 
%         dim = [.2 .5 .3 .3];
%         ann = annotation('textbox',dim,'String','','verticalalignment','top','horizontalalignment','left','FitBoxToText','on','Interpreter','Latex','FontSize',30);
    end
    if ~S.config.SIM == true; S = ROS(S); end
    for num_sims = 1:S.config.NUM_SIMS
        S = call_init_functions(S);
        %
        S.config.iters = 1;
        while (S.config.time(end) < S.config.tf)% && ~(S.path.reach_end_mhempc || ~S.config.mpc)
            if check_end_condition(S); break; else
                S = CALL_SIMULATOR(S);
            end
            % #############################################################        
            S.config.time   = [S.config.time, S.config.time(end)+S.config.Ts];
            S.config.iters  = S.config.iters+1;
            a               = [num_sims, S.config.iters, S.exec_time.t_tot(end)/S.config.Ts]
            % -------------------------------------------------------------
            if S.config.dataOnline == true
                plot_mono2(S, S.data.xest(:,end));
                % plot_mono2(S, S.algorithms.mpcCasadi.qRefOpt1(:,1),'m');
                % plot_mono2(S, S.algorithms.mpcCasadi.qRefOpt1(:,end),'c');
            end
            S.exec_time.t_mis = [S.exec_time.t_mis, toc];
            % -------------------------------------------------------------
            if S.config.mpc == true
                S.data.mhempc.performance.xfut{num_sims, S.config.iters-1}  = S.algorithms.mpcCasadi.Qtraj;
                S.data.mhempc.performance.J{num_sims, S.config.iters-1}     = S.algorithms.mpcCasadi.Jnum;
                S.data.mhempc.performance.qref{num_sims, S.config.iters-1}  = S.algorithms.mpcCasadi.qRefOpt1;
                S.data.mhempc.performance.XYR{num_sims, S.config.iters-1}   = S.algorithms.mpcCasadi.XYRobst;
            end
        end
        if ~S.config.SIM
            S = write_ctrl(S, 0, 0); % Stop the vehicle
        end
        % -----------------------------------------------------------------
        S.data.mhempc.performance.xsim{num_sims} = S.data.xsim;            
        S.data.mhempc.performance.ysim{num_sims} = S.data.ysim;
        S.data.mhempc.performance.xest{num_sims} = S.data.xest;
        S.data.mhempc.performance.ctrl{num_sims} = S.algorithms.Controls;
        S.data.mhempc.performance.umeas{num_sims} = S.sensors.velocities;
        S.data.mhempc.performance.mpcRefs{num_sims} = S.data.references;
        S.data.mhempc.performance.exec_times{num_sims} = S.exec_time;
        if S.config.save_workspace == true
            if S.config.lidarOn == true
                stop(S.vlp16);
                S = rmfield(S,'vlp16');
                % S = rmfield(S,'obs');
            end
            % clk = clock; save([date,'-',num2str(clk(1:5))]);
        end
%         if ~S.config.reference; ref = S.data.xsim; save(S.path.name,'ref'); end
    end
end

function S = build_setup(S)
    % SIMULATION PARAMETERS ===============================================
    % S.ROS.IMU_ZEROING_VEC1  = -2.68;
    % S.ROS.IMU_ZEROING_VEC2  = -0.46;
    % S.ROS.IMU_ZEROING_MIC1  = 2.11;
    S.ROS.ENCODERS_PULSE_TO_RAD = 0.6 * pi / 180;
    S.config.ZEROING_VEC1   = true;
    S.config.ZEROING_VEC2   = false;
    S.config.ZEROING_MIC1   = false;
    S.config.num_meas_zeroing = 100;
    % Mobile Sensors ______________________________________________________
    S.mobile.useMobileSensors = true;
    S.mobile.nroMobileSensors = 1;
    % Solver ______________________________________________________________
    S.config.solver         = 'casadi'; % options: 'casadi', 'acado'
    % Simulation or field experiment (Husky)_______________________________
    S.config.SIM            = true;
    S.config.N              = 2;
    S.config.verbose        = false;
    S.config.calcMtxConvCoord = true; if ~S.config.calcMtxConvCoord; warning('WARNING!!! The matrix for correcting x-y coordinates is not being computed...'); end;
    S.config.t0             = 0;
    S.config.tf             = 3500;
    S.config.Ts             = 0.65;
    S.config.same_seed      = false;
    S.config.EXPORT         = false; if ~S.config.EXPORT; warning('WARNING!!! MHE and MPC were not compiled...'); end
    S.config.iters          = 0;
    S.config.time           = 0;
    S.config.NUM_SIMS       = 1;
    S.config.updtAC         = false;
    S.config.outputs        = [1:S.config.N+1,2*S.config.N+2:2*S.config.N+3];   % 2*S.config.N+6:2*S.config.N+7];
    S.config.Nc             = 10;%6;                                               % Lenght of the control horizon
    S.config.Ne             = 6;                                                % Lenght of the estimation window
    % Estimator and Control algorithm to execute ________________________
    S.config.obs_strategy   = 'gauss';
    S.config.mpc            = true;
    S.config.mhe            = true;
    % Disturbances ______________________________________________________
    S.config.noise          = 'gaussian'; % 'gaussian' or 'uniform'
    S.config.noise_lvl      = 0.*[(0.2*pi/180).*ones(S.config.N, 1); 0.2*pi/180; [0.05; 0.05]];%; [0.05; 0.05]]; % measurement noise amplitude: [betas'; theta_0; x_0; y_0; w0; v0]
    S.config.initUncertainty = 10*pi/180;
    S.config.slip           = [1; 1];
    S.config.procDist.type  = 'normal';
    S.config.procDist.amp   = 0.0.*[(0.5*pi/180).*ones(S.config.N,1); zeros(S.config.N+1,1); zeros(2,1); zeros(2,1)];
    S.config.gpsOutlier     = false;
    S.config.gpsOutlierAmp  = 2;
    S.config.gpsOutlierProbab = 10;
    S.config.IMUbias        = false;
    S.config.IMUbiasDev     = 5/180;
    S.config.model.uncty    = false;
    S.config.model.dev      = 55;                                       % porcentual value of uncertainty
    S.config.iterCorrecOn   = 200;
    S.config.CtrlNoise_lvl  = [0; 0];
    % Reference velocities ________________________________________________
    S.config.vNr            = 0;
    S.config.vNTarget       = 0.3;
    S.config.timeToReachVNr = 5; % seconds
    S.config.deltavN        = S.config.vNTarget/(S.config.timeToReachVNr/S.config.Ts + 2*(S.config.Ne+1));
    S.config.wNr            = 0;
    S.config.uref           = [S.config.wNr; S.config.vNr];    
    % GPS's relative position to centre of vehicle ________________________
    S.config.gps_x          = -0.3;
    S.config.gps_y          = 0;
    S.config.gps_d          = sqrt(S.config.gps_x^2 + S.config.gps_y^2);
    S.config.gps_fcx        = 0;
    S.config.gps_fcy        = 0;
    S.config.SIMGPS         = false; % use to simulate date from GPS when it does not have signal. Just for indoor testing purposes.
    % Obstacles ___________________________________________________________
    S.config.lidarLimits.X  = 5;%1.0;       % distance to points from the lidar
    S.config.lidarLimits.Y  = 10;%2.3;   
    S.config.lidarLimits.Z  = 1;
    % boundaries for ignoring obstacles
    S.config.x_max          = 13;
    S.config.x_min          = 2;
    S.config.y_max          = 8;
    S.config.y_min          = 1.5;
    S.config.lidarOn        = false;
    S.config.obsDetection   = false;
    % Plot data online during experiments. USeful for debbuging purposes __
    S.config.dataOnline     = true;
    S.config.save_workspace = false;
end

function S = CALL_SIMULATOR(S)
    % Solve estimation problem ________________________________
    S = call_ESTIMATOR(S);
    % Obstacle detection ______________________________________
    S = call_OBSDETECTOR(S);
    % Path-tracking problem ___________________________________
    S = call_PFA(S);
    % Solve control problem ___________________________________
    S = call_CONTROLLER(S);                
    % Apply controls to the Husky _____________________________
    S = call_SYSTEM(S);                
    % Update measurements _____________________________________
    S = update_measurements(S);
    % Perform real-time iteration _____________________________
    S = call_RTF(S);
end

function S = call_RTF(S)
    tic;
    S.exec_time.t_tot   = [ S.exec_time.t_tot, S.exec_time.t_mhe(end)+S.exec_time.t_pfa(end)+S.exec_time.t_mpc(end)+...
                            S.exec_time.t_ctrl(end)+S.exec_time.t_mis(end)+S.exec_time.t_obsdetector(end)+S.exec_time.t_sensors(end)];
    S.exec_time.t_acum  = S.exec_time.t_tot(end);
    while ((S.exec_time.t_acum + toc) < S.config.Ts) && ~S.config.SIM; end
end

function plot_frames(S, qk, qref, clr)
    if nargin == 4
        clrTractor          = clr;
        clrWheel            = clr;
        clrAxe              = clr;
        clrLongAxe          = clr;
        clrTrailerLoad      = clr;
        clrLastTrailerLoad  = clr;
        noRef = true;
    else
        clrTractor          = 'y';
        clrWheel            = 'k';
        clrAxe              = 'k';
        clrLongAxe          = 'k';
        clrTrailerLoad      = 'b';
        clrLastTrailerLoad  = 'r';
        noRef = false;
    end
    %
    if ~noRef
        plot(S.path.coordinates(1,:),S.path.coordinates(2,:),'color',[0.6 0.6 0.6],'LineWidth',8); grid on; daspect([1 1 1]); hold on;
        xlim([min(S.path.coordinates(1,:))-1 max(S.path.coordinates(1,:))+1]); ylim([min(S.path.coordinates(2,:))-1 max(S.path.coordinates(2,:))+1]);
    end    
    % Rotate the Husky and wheels according to their attitude
    thetas  = qk(S.config.N+1:2*S.config.N+1);
    betas   = qk(1:S.config.N);    
    xy0     = qk(2*S.config.N+2:2*S.config.N+3);
    xy0ref  = qref(2*S.config.N+2:2*S.config.N+3);
    xyi     = xy0;
    xyiref  = xy0ref;
    for i=1:S.config.N
%         xyi(:,i+1) = xyi(:,i) - [S.system.Lhi(i)*cos(thetas(i)) + S.system.Li(i)*cos(thetas(i+1)); S.system.Lhi(i)*sin(thetas(i)) + S.system.Li(i)*sin(thetas(i+1))]; 
        xyi     = [xyi, qk(2*S.config.N+3+(i-1)*2+1:2*S.config.N+3+(i)*2)];
        xyiref  = [xyiref, qref(2*S.config.N+3+(i-1)*2+1:2*S.config.N+3+(i)*2)];
    end
%     xyi(:,end) = qk(2*S.config.N+4:2*S.config.N+5);
    % Plot the frames -----------------------------------------------------
    val = [];
    pos = [];
    R = [cos(-pi/2), -sin(-pi/2); sin(-pi/2) cos(-pi/2)];
    for i=1:S.config.N+1
        dis                 = sqrt( (S.path.coordinates(1,:)-xyi(1,i)).^2 + (S.path.coordinates(2,:)-xyi(2,i)).^2 );
        [~, pos_aux]        = min(dis);
        %
        xy_pathF            = S.path.coordinates(:,pos_aux);
        alphaF              = atan2( S.path.coordinates(2,pos_aux+1)-S.path.coordinates(2,pos_aux), S.path.coordinates(1,pos_aux+1)-S.path.coordinates(1,pos_aux) );
        uF1                 = cos(alphaF);
        vF1                 = sin(alphaF);
        uv2                 = R*[uF1;vF1];
        uF2                 = uv2(1);
        vF2                 = uv2(2);
        quiver(xy_pathF(1), xy_pathF(2), uF1, vF1, 'm', 'LineWidth',2);
        quiver(xy_pathF(1), xy_pathF(2), uF2, vF2, 'm', 'LineWidth',2);
        %
        xy_pathM            = qref(2*S.config.N+2+(i-1)*2:2*S.config.N+1+i*2);
        alphaM              = atan2( -xy_pathM(2)+xyi(2,i), -xy_pathM(1)+xyi(1,i) );
        uM1                 = cos(qref(S.config.N+i));
        vM1                 = sin(qref(S.config.N+i));
        uv2                 = R*[uM1;vM1];
        uM2                 = uv2(1);
        vM2                 = uv2(2);
%         quiver(xy_pathM(1), xy_pathM(2), uM1, vM1, 'c', 'LineWidth',2);
%         quiver(xy_pathM(1), xy_pathM(2), uM2, vM2, 'c', 'LineWidth',2);
    end
    % plot the vehicle's segments -----------------------------------------
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
        line([trailerLongAxePlt(1,2),xyi(1,i)],[trailerLongAxePlt(2,2),xyi(2,i)],'color',clrLongAxe,'linewidth',1)    
        line(trailerAxePlt(1,:),trailerAxePlt(2,:),'color',clrAxe,'linewidth',2)
        line(trailerWheelLeftPlt(1,:),trailerWheelLeftPlt(2,:),'color',clrWheel,'linewidth',4)
        line(trailerWheelRightPlt(1,:),trailerWheelRightPlt(2,:),'color',clrWheel,'linewidth',4)
        line(trailerLongAxePlt(1,:),trailerLongAxePlt(2,:),'color',clrLongAxe,'linewidth',2)
        line(trailerLoadPlt(1,:),trailerLoadPlt(2,:),'color',trailerClr,'linewidth',2)
    end
    % ---------------------------------------------------------------------    
    %
    if ~noRef
        hold off;
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
    else
        clrTractor          = 'y';
        clrWheel            = 'k';
        clrAxe              = 'k';
        clrLongAxe          = 'k';
        clrTrailerLoad      = 'b';
        clrLastTrailerLoad  = 'r';
        noRef = false;
    end
    %
    if isempty(S.algorithms.mpcCasadi.Qtraj)
        Qtraj = repmat(S.algorithms.mpcCasadi.q0bar, 1, S.config.Nc);
    else
        Qtraj = S.algorithms.mpcCasadi.Qtraj;
    end

    if ~noRef
        plot(S.path.coordinates(1,:),S.path.coordinates(2,:),'color',[0.6 0.6 0.6],'LineWidth',8); grid on; daspect([1 1 1]); hold on;
        xlim([min(S.path.coordinates(1,:))-2 max(S.path.coordinates(1,:))+2]); ylim([min(S.path.coordinates(2,:))-4 max(S.path.coordinates(2,:))+4]);
        %
        for i=1:S.config.N+1
            if i==1
                plot(Qtraj(2*S.config.N+1+(i-1)*2+1,:),Qtraj(2*S.config.N+1+i*2,:),'y','LineWidth',4);
            elseif i~=S.config.N+1
                plot(Qtraj(2*S.config.N+1+(i-1)*2+1,:),Qtraj(2*S.config.N+1+i*2,:),'b','LineWidth',4);
            else
                plot(Qtraj(2*S.config.N+1+(i-1)*2+1,:),Qtraj(2*S.config.N+1+i*2,:),'r','LineWidth',4);
            end
        end
        %
        if~isempty(S.path.obstacles)
            nroObs = size(S.path.obstacles,1);
            for i=1:nroObs
                circles(S.path.obstacles(i,1),S.path.obstacles(i,2),S.path.obstacles(i,3),'color','red')
            end
        end
%         patch(S.config.initUncertaintySet.x,S.config.initUncertaintySet.y,'r','edgecolor','none','facealpha',0.2)
    end    
    plot(S.algorithms.mpcCasadi.qRefOpt1(2*S.config.N+2,end),S.algorithms.mpcCasadi.qRefOpt1(2*S.config.N+3,end),'m+','linewidth',1.5,'markersize',20)
    plot(S.algorithms.mpcCasadi.qRefOpt1(2*S.config.N+2,end),S.algorithms.mpcCasadi.qRefOpt1(2*S.config.N+3,end),'mo','linewidth',1.5,'markersize',20)
    plot(S.algorithms.mpcCasadi.qRefOpt1(2*S.config.N+4,end),S.algorithms.mpcCasadi.qRefOpt1(2*S.config.N+5,end),'mx','linewidth',1.5,'markersize',20)
    plot(S.algorithms.mpcCasadi.qRefOpt1(2*S.config.N+4,end),S.algorithms.mpcCasadi.qRefOpt1(2*S.config.N+5,end),'mo','linewidth',1.5,'markersize',20)
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
    line(tractorBodyPlt(1,:),tractorBodyPlt(2,:),'color',clrTractor,'linewidth',3);
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
        line([trailerLongAxePlt(1,2),xyi(1,i)],[trailerLongAxePlt(2,2),xyi(2,i)],'color',clrLongAxe,'linewidth',1)    
        line(trailerAxePlt(1,:),trailerAxePlt(2,:),'color',clrAxe,'linewidth',2)
        line(trailerWheelLeftPlt(1,:),trailerWheelLeftPlt(2,:),'color',clrWheel,'linewidth',4)
        line(trailerWheelRightPlt(1,:),trailerWheelRightPlt(2,:),'color',clrWheel,'linewidth',4)
        line(trailerLongAxePlt(1,:),trailerLongAxePlt(2,:),'color',clrLongAxe,'linewidth',2)
        line(trailerLoadPlt(1,:),trailerLoadPlt(2,:),'color',trailerClr,'linewidth',2)
    end
    
%     C = [S.algorithms.mpcCasadi.OFCosts.Jobs,...
%         S.algorithms.mpcCasadi.OFCosts.Jxyref,...
%         S.algorithms.mpcCasadi.OFCosts.Jws0,...
%         S.algorithms.mpcCasadi.OFCosts.Jgeom,...
%         S.algorithms.mpcCasadi.OFCosts.Jq,...
%         S.algorithms.mpcCasadi.OFCosts.JqN,...
%         S.algorithms.mpcCasadi.OFCosts.Ju,...
%         S.algorithms.mpcCasadi.OFCosts.Jsigma,...
%         S.algorithms.mpcCasadi.OFCosts.Jwqref1,...
%         S.algorithms.mpcCasadi.OFCosts.Jwq0]%,...
%         sum([S.algorithms.mpcCasadi.OFCosts.Jobs+S.algorithms.mpcCasadi.OFCosts.Jxyref+S.algorithms.mpcCasadi.OFCosts.Jws0+...
%         S.algorithms.mpcCasadi.OFCosts.Jgeom, S.algorithms.mpcCasadi.OFCosts.Jq, S.algorithms.mpcCasadi.OFCosts.JqN,...
%         S.algorithms.mpcCasadi.OFCosts.Ju+S.algorithms.mpcCasadi.OFCosts.Jsigma+S.algorithms.mpcCasadi.OFCosts.Jwqref1+...
%         S.algorithms.mpcCasadi.OFCosts.Jwqref2+S.algorithms.mpcCasadi.OFCosts.Jwq0])-S.algorithms.mpcCasadi.Jnum]

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
    if ~noRef
        hold off;
    end
    %
    drawnow limitrate
end

function S = call_init_functions(S)
    S = reserve_temp_memory(S);
    S = gen_init_conditions(S);
    S = init_flags_and_counters(S);    

    % temporarily here:
    if S.config.lidarOn
        S = init_obsDetection(S);
        start(S.vlp16);
    end
end

function S = compute_tras_rot_mtx(S)
    sdpvar a11 a12 a21 a22;
    sdpvar x1bar y1bar x2bar y2bar;
    % Besides the reference point, two more are needed to obtaint the local
    % reference frame
    % AZOTEA AC3E *********************************************************
    % Coordinates of point (x, 0)
    S.ROS.local_coord_1.lat     = -33.034213;       % hand coded value, measurement from rtk
    S.ROS.local_coord_1.long    = -71.592168;   % hand coded value, measurement from rtk
    % CANCHA DE FUTBOL ****************************************************
%     S.ROS.local_coord_1.lat     = -33.03517;       % hand coded value, measurement from rtk
%     S.ROS.local_coord_1.long    = -71.594195;   % hand coded value, measurement from rtk
    % *********************************************************************
    [x1, y1]                    = latlon2xy(S.ROS.local_coord_1.lat, S.ROS.local_coord_1.long, S.ROS.LAT0, S.ROS.LON0);
    x1                          = x1*1000;
    y1                          = y1*1000;
    x                           = norm([x1 y1]);
    % AZOTEA AC3E *********************************************************
    % Coordinates of point (0, y)
    S.ROS.local_coord_2.lat     = -33.034088;        % hand coded value, measurement from rtk
    S.ROS.local_coord_2.long    = -71.59211333;         % hand coded value, measurement from rtk
    % CANCHA DE FUTBOL ****************************************************
%     S.ROS.local_coord_2.lat     = -33.034403333;        % hand coded value, measurement from rtk
%     S.ROS.local_coord_2.long    = -71.594075;         % hand coded value, measurement from rtk
    % *********************************************************************
    [x2, y2]                    = latlon2xy(S.ROS.local_coord_2.lat, S.ROS.local_coord_2.long, S.ROS.LAT0, S.ROS.LON0);
    x2                          = x2*1000;
    y2                          = y2*1000;
    y                           = norm([x2 y2]);
    % With the "origin" and a point alongside each axe, compute the mtx
    A                           = [a11 a12; a21 a22];
    v                           = [x1; y1; x2; y2];
    b                           = [x1bar; y1bar; x2bar; y2bar];
    Constraints                 = [[A zeros(2); zeros(2) A]*v - b == zeros(4,1); x1bar*x2bar + y1bar*y2bar == 0];
    %
    obj                         = (x1bar - x)^2 + (y2bar - y)^2;
    % 
    optimize(Constraints, obj);
    %
    S.ROS.Mtx                   = value(A);
end

function S = update_measurements(S)
    tic;
    if S.config.SIM == true
        if strcmp(S.config.noise,'gaussian')
            noise = randn(S.system.ny,1).*S.config.noise_lvl;
        else
            noise = (2*rand(S.system.ny,1)-1).*S.config.noise_lvl;
        end
        if S.config.gpsOutlier == true
            if randi(S.config.gpsOutlierProbab) == 1
                noise(S.config.N+2:S.config.N+3) = S.config.gpsOutlierAmp*randn(2,1);
            end
        end
        if S.config.IMUbias == true
            noise(S.config.N+1) = noise(S.config.N+1) + S.config.IMUbiasDev*S.config.iters*S.config.Ts/3600; % deviation  per hour     %*S.data.xsim(S.config.N+1,end);
        end

        S.data.measNoise = [S.data.measNoise, noise];
        S.data.ysim = [S.data.ysim, S.data.xsim(S.config.outputs,end) + noise];
        %
        unoise = randn(S.system.nu,1).*S.config.CtrlNoise_lvl;
        S.data.UmeasNoise = [S.data.UmeasNoise, unoise];
        S.sensors.velocities = S.algorithms.Controls(:,end) + unoise;
    else
        S           = read_sensors(S);
        S           = measurements_vector(S);
        S.data.ysim = [S.data.ysim, S.sensors.measurement]; % <-- MEASUREMENTS are stored in this vector in this order: [beta_1, beta_2, theta_0, x_0, y_0]
    end
    S = update_EstimatorMeasurement(S);
    S.exec_time.t_sensors = [S.exec_time.t_sensors, toc];
end

function flag = check_end_condition(S)
    flag        = false;
    condition   = S.algorithms.mpcCasadi.s_pos(1);
    % condition = norm(S.data.xest(2*S.config.N+2:2*S.config.N+3,end)-S.path.coordinates(:,end));
    if condition >= 2*pi
    % if condition <= 1
        flag = true;
    end
end

function S = ROS(S)
    S = init_ROS(S);
    S = create_obj_sens(S);
    S = create_obj_vel(S);    
    S = zeroing_imu(S);
end

function S = MOBILE(S)
    S = init_mobile(S);
end

function S = init_mobile(S)
    conn_devices = size(mobiledevlist,1);
    while conn_devices < S.mobile.nroMobileSensors
        conn_devices = size(mobiledevlist,1);
        fprintf('Waiting for mobiles to connect...\n')
    end
    S.mobile.devs    = {};
    S.mobile.devs{1} = 'POCO M4 Pro 5G';
    S.mobile.devs{2} = 'POCO M4 Pro 5G';
    S.mobile.devs{3} = 'POCO M4 Pro 5G';
    S.mobile.devs{4} = 'POCO M4 Pro 5G';
    %
    S.mobile.sensors = {};
    for i=1:S.config.nroMobileSensors
        S.mobile.sensors{i} = mobiledev(S.mobile.devs{i});
    end
      
end

function S = readMobile(S)
    for i=1:S.config.nroMobileSensors
        data = [S.mobile.sensors{i}.Acceleration'; S.mobile.sensors{i}.AngularVelocity';...
            S.mobile.sensors{i}.MagneticField'; S.mobile.sensors{i}.Orientation';...
            S.mobile.sensors{i}.Latitude; S.mobile.sensors{i}.Longitude; S.mobile.sensors{i}.Speed;...
            S.mobile.sensors{i}.Course; S.mobile.sensors{i}.Altitude; S.mobile.sensors{i}.HorizontalAccuracy ];
        S.mobile.data{i} = [S.mobile.data{i}, data];
    end
    
end

function S = init_ROS(S)
    % Init coordinates of my local 2D plane
    S = get_reference(S);
    
    if S.config.calcMtxConvCoord == true
        S = compute_tras_rot_mtx(S);
    else
        S.ROS.Mtx = eye(2);
    end
    
    S.sensors.measurements = [];
    % Create and init Husky's velocities
    S.algorithms.mpc.Husky_v0 = 0;
    S.algorithms.mpc.Husky_w0 = 0;
    % Here some instructions need to be followed from matlab's terminal
    rosshutdown;
    input('use el comando roscore desde terminal...\n')
    %
    rosshutdown;
    rosinit;
    %
    fprintf('\n---------------------------------------------------------\n');
    fprintf('RECORDAR RESETEAR LOS ENCODERS ANTES DEL EXPERIMENTOS...\n');
    fprintf('---------------------------------------------------------\n');
    %
    fprintf(['\n\nroslaunch husky_base base.launch...\n' ...
        'roslaunch parches rtk.launch...\n' ...
        'roslaunch vectornav vectornav.launch...\n' ...
        'roslaunch parches parche_imu_ins.launch...\n' ...
        'roslaunch vectornav vectornav2.launch...\n' ...
        'roslaunch imu_3dm_gx4 imu.launch...\n' ...
        'roslaunch parches encoders.launch...\n' ...
        'roslaunch parches parche_speed_holder.launch...\n'])
    aux1 = input('Presione enter...\n');
end

function S = create_obj_sens(S)
    sub_rtk         = rossubscriber('/fix');
    sub_vec1        = rossubscriber('/vectornav/IMU');
    sub_vec1_vel    = rossubscriber('/imu_INS');
    % sub_vec2        = rossubscriber('/vectornav2/IMU');
    sub_encoders    = rossubscriber('/enc');
    % sub_micro1      = rossubscriber('/imu/pose');
%     sub_lidar       = rossubscriber('/angles');
    %
    S.ROS.rtk            = sub_rtk;
    S.ROS.vectornav1     = sub_vec1;
    % S.ROS.vectornav2     = sub_vec2;
    S.ROS.vectornav1_vel = sub_vec1_vel;
    % S.ROS.microstrain    = sub_micro1;
    S.ROS.IncEncoders    = sub_encoders;
%     S.ROS.betasFromLidar = sub_lidar;

    % COrection factor for unwraping phase on real-time
    S.ROS.pose.correction_vec1              = 0;      % This value should be substracted to every pose measurement
    S.ROS.sensors.vectornav_euler_vec1      = [];
    S.ROS.sensors.vectornav_euler_vec1_Old  = 0;     % Value for computinng the difference. If it is bigger in absolute value than pi, a correction is needed
    %
    % S.ROS.pose.correction_vec2              = 0;
    % S.ROS.sensors.vectornav_euler_vec2      = [];
    % S.ROS.sensors.vectornav_euler_vec2_Old  = 0;
    %
    % S.ROS.pose.correction_micro1            = 0;      % This value should be substracted to every pose measurement
    % S.ROS.sensors.microstrain_euler_micro1    = [];
    % S.ROS.sensors.microstrain_euler_micro1_Old  = 0;
end

function S = create_obj_vel(S)
    [pub,msg]       = rospublisher('/cmd_vel_aux','geometry_msgs/Twist');
    %
    S.ROS.CTRL.pub       = pub;
    S.ROS.CTRL.msg       = msg;
end

function S = get_reference(S)
    % Coordiantes of the origin of our local reference 2D-axes x-y. Values
    % charged by hand before carry out the experiments.
    % AZOTEA AC3E *********************************************************
    S.ROS.LAT0   = -33.034115;
    S.ROS.LON0  = -71.592205;
    % CANCHA DE FUTBOL ****************************************************
%     S.ROS.LAT0   = -33.03453;
%     S.ROS.LON0  = -71.594498;
    %
end

function S = read_rtk(S)
    if ~S.config.SIMGPS
        S.ROS.sensors.gpsrtk = receive(S.ROS.rtk);
        
        if isnan(S.ROS.sensors.gpsrtk.Latitude) || isnan(S.ROS.sensors.gpsrtk.Longitude)
            S = write_ctrl(S, 0, 0);
            while isnan(S.ROS.sensors.gpsrtk.Latitude) || isnan(S.ROS.sensors.gpsrtk.Longitude)
                fprintf('NO GPS SIGNAL...\n')
                pause(0.2);
            end
        end
        % Convert lat &v lon to x-y coordinated in a tangential plane ro
        % earth's surface (tangent to my reference point)
        [x,y]   = latlon2xy(S.ROS.sensors.gpsrtk.Latitude,S.ROS.sensors.gpsrtk.Longitude,S.ROS.LAT0,S.ROS.LON0);
        % Convert Km to m
        x_raw   = x*1000;
        y_raw   = y*1000;
        % Adjust to my local 2D plane through the rotation matrix
        xy_cor  = S.ROS.Mtx * [x_raw; y_raw];
        % The, store values to be used later
        S.ROS.sensors.rtk.x0 = xy_cor(1) - S.config.gps_fcx;
        S.ROS.sensors.rtk.y0 = xy_cor(2) - S.config.gps_fcy;
% pos=[S.ROS.sensors.rtk.x0, S.ROS.sensors.rtk.y0]
    else
        S.ROS.sensors.rtk.x0 = S.data.xsim(2*S.config.N+2,end);
        S.ROS.sensors.rtk.y0 = S.data.xsim(2*S.config.N+3,end);
        %
        fprintf('Simulated data from GPS...\n')
    end
end

function S = write_ctrl(S, w0, v0)
    tic;
    S.ROS.CTRL.msg.Linear.X    = v0;
    S.ROS.CTRL.msg.Angular.Z   = w0;
    send(S.ROS.CTRL.pub,S.ROS.CTRL.msg);
    S.exec_time.t_ctrl = [S.exec_time.t_ctrl, toc];
end

function S = zeroing_imu(S)
%
    fprintf('\nZeroing imus...\n');
%
    if S.config.ZEROING_VEC1
        zero_vec1  = 0;
        for i=1:S.config.num_meas_zeroing
            vec1_orientation    = receive(S.ROS.vectornav1);
            quat_vec1           = vec1_orientation.Orientation;
            raw_vec1            = quat2eul([quat_vec1.X,quat_vec1.Y,quat_vec1.Z,quat_vec1.W]);
            zero_vec1           = zero_vec1 + raw_vec1(3);
        end
        S.ROS.IMU_ZEROING_VEC1 = zero_vec1 / S.config.num_meas_zeroing;
    else
        S.ROS.IMU_ZEROING_VEC1 = 2.48;
    end
    %
    if S.config.ZEROING_VEC2
        zero_vec2     = 0;
        for i=1:S.config.num_meas_zeroing
            vec2_orientation    = receive(S.ROS.vectornav2);
            quat_vec2           = vec2_orientation.Orientation;
            raw_vec2            = quat2eul([quat_vec2.X,quat_vec2.Y,quat_vec2.Z,quat_vec2.W]);
            zero_vec2           = zero_vec2 + raw_vec2(3);
        end
        S.ROS.IMU_ZEROING_VEC2 = zero_vec2 / S.config.num_meas_zeroing;
    else
        S.ROS.IMU_ZEROING_VEC2 = -2.93;
    end
    %
    if S.config.ZEROING_MIC1
        zero_mic1     = 0;
        for i=1:S.config.num_meas_zeroing
            micro1_orientation  = receive(S.ROS.microstrain);
            quat_mic1           = micro1_orientation.Pose.Orientation;
            raw_mic1            = quat2eul([quat_mic1.X,quat_mic1.Y,quat_mic1.Z,quat_mic1.W]);
            zero_mic1           = zero_mic1 + raw_mic1(3);
        end
        S.ROS.IMU_ZEROING_MIC1 = zero_mic1 / S.config.num_meas_zeroing;
    else
        S.ROS.IMU_ZEROING_MIC1 = -2.96;
    end
    %
    fprintf('Zeroing imus ready!\n\n');
end

function S = read_vectornav(S) % measure theta0 and the speeds
    % pose
    S.ROS.sensors.vectornav             = receive(S.ROS.vectornav1);
    quat                                = S.ROS.sensors.vectornav.Orientation;
    S.ROS.sensors.vectornav_euler_vec1  = quat2eul([quat.X,quat.Y,quat.Z,quat.W]);    
    % Unwrap phase from online data _______________________________________
    if (S.ROS.sensors.vectornav_euler_vec1(3)-S.ROS.sensors.vectornav_euler_vec1_Old) >= pi
        S.ROS.pose.correction_vec1 = S.ROS.pose.correction_vec1 + 2*pi;
    elseif (S.ROS.sensors.vectornav_euler_vec1(3)-S.ROS.sensors.vectornav_euler_vec1_Old) <= -pi
        S.ROS.pose.correction_vec1 = S.ROS.pose.correction_vec1 - 2*pi;         
    end         
    %
    S.ROS.sensors.vectornav_euler_vec1_Old = S.ROS.sensors.vectornav_euler_vec1(3);    
    % Do Not compute the attitude angle in my reference frame -------------
    S.ROS.sensors.vectornav_theta0 = -S.ROS.sensors.vectornav_euler_vec1(3) + S.ROS.IMU_ZEROING_VEC1 + S.ROS.pose.correction_vec1;
    % Measure the speed ---------------------------------------------------
    S.ROS.sensors.vectornav_vel     = receive(S.ROS.vectornav1_vel);
    S.ROS.sensors.vectornav_NedVelX = S.ROS.sensors.vectornav_vel.Data(1);
    S.ROS.sensors.vectornav_NedVelY = S.ROS.sensors.vectornav_vel.Data(2);
    S.ROS.sensors.vectornav_v0      = sqrt(S.ROS.sensors.vectornav_NedVelX^2 + S.ROS.sensors.vectornav_NedVelY^2);
    S.ROS.sensors.vectornav_w0      = -1 * S.ROS.sensors.vectornav.AngularVelocity.Z;
    % Compute correction factor ___________________________________________
    S.config.gps_fcx = S.config.gps_d * sin(S.ROS.sensors.vectornav_theta0);
    S.config.gps_fcy = S.config.gps_d * cos(S.ROS.sensors.vectornav_theta0);
end

function S = read_vectornav2(S) % measure theta2
    % pose
    S.ROS.sensors.vectornav2            = receive(S.ROS.vectornav2);
    quat                                = S.ROS.sensors.vectornav2.Orientation;
    S.ROS.sensors.vectornav_euler_vec2  = quat2eul([quat.X,quat.Y,quat.Z,quat.W]);    
    % Unwrap pahse from online data _______________________________________
    if (S.ROS.sensors.vectornav_euler_vec2(3)-S.ROS.sensors.vectornav_euler_vec2_Old) >= pi
        S.ROS.pose.correction_vec2 = S.ROS.pose.correction_vec2 + 2*pi;
    elseif (S.ROS.sensors.vectornav_euler_vec2(3)-S.ROS.sensors.vectornav_euler_vec2_Old) <= -pi
        S.ROS.pose.correction_vec2 = S.ROS.pose.correction_vec2 - 2*pi;         
    end         
    %
    S.ROS.sensors.vectornav_euler_vec2_Old = S.ROS.sensors.vectornav_euler_vec2(3);
    % No compute the attitude angle in my reference frame -----------------
    S.ROS.sensors.vectornav_theta2 = -S.ROS.sensors.vectornav_euler_vec2(3) + S.ROS.IMU_ZEROING_VEC2 + S.ROS.pose.correction_vec2;
end

function S = read_microstrain(S)
    S.ROS.sensors.microstrain   = receive(S.ROS.microstrain);
    quat                        = S.ROS.sensors.microstrain.Pose.Orientation;
    S.ROS.sensors.microstrain_euler_micro1 = quat2eul([quat.X,quat.Y,quat.Z,quat.W]);
    % Unwrap pahse from online data _______________________________________
    if (S.ROS.sensors.microstrain_euler_micro1(3)-S.ROS.sensors.microstrain_euler_micro1_Old) >= pi
        S.ROS.pose.correction_micro1 = S.ROS.pose.correction_micro1 + 2*pi;
    elseif (S.ROS.sensors.microstrain_euler_micro1(3)-S.ROS.sensors.microstrain_euler_micro1_Old) <= -pi
        S.ROS.pose.correction_micro1 = S.ROS.pose.correction_micro1 - 2*pi;         
    end         
    %
    S.ROS.sensors.microstrain_euler_micro1_Old = S.ROS.sensors.microstrain_euler_micro1(3);
    % No compute the attitude angle in my reference frame -----------------
    S.ROS.sensors.microstrain_theta1 = -S.ROS.sensors.microstrain_euler_micro1(3) + S.ROS.IMU_ZEROING_MIC1 + S.ROS.pose.correction_micro1;
end

function S = read_encoders(S)
    S.ROS.sensors.encoder.data      = receive(S.ROS.IncEncoders);    
    S.ROS.sensors.encoders.beta1    = double(S.ROS.sensors.encoder.data.Data(1)) * S.ROS.ENCODERS_PULSE_TO_RAD;
    S.ROS.sensors.encoders.beta2    = double(S.ROS.sensors.encoder.data.Data(2)) * S.ROS.ENCODERS_PULSE_TO_RAD;
end

function S = read_sensors(S)
    S = read_vectornav(S);
    % S = read_vectornav2(S);
    % S = read_microstrain(S);
    S = read_encoders(S);
    S = read_rtk(S);
end

function S = measurements_vector(S)
%     S.sensors.measurement   = [S.ROS.sensors.encoders.beta1; S.ROS.sensors.encoders.beta2; S.ROS.sensors.vectornav_theta0; S.ROS.sensors.rtk.x0; S.ROS.sensors.rtk.y0; S.ROS.sensors.vectornav_w0; S.ROS.sensors.vectornav_v0];
    S.sensors.measurement   = [S.ROS.sensors.encoders.beta1; S.ROS.sensors.encoders.beta2; S.ROS.sensors.vectornav_theta0; S.ROS.sensors.rtk.x0; S.ROS.sensors.rtk.y0];
    S.sensors.velocities    = [S.sensors.velocities, [S.ROS.sensors.vectornav_w0; S.ROS.sensors.vectornav_v0]];
    %
    S.sensors.measurements  = [S.sensors.measurements, S.sensors.measurement];   
    % Store the attitude of each trailer
    S.sensors.theta0        = [S.sensors.theta0, S.ROS.sensors.vectornav_theta0];
    % S.sensors.theta1        = [S.sensors.theta1, S.ROS.sensors.microstrain_theta1];
    % S.sensors.theta2        = [S.sensors.theta2, S.ROS.sensors.vectornav_theta2];
end

function S = update_EstimatorMeasurement(S)
    if S.config.mhe
        updateMeasurement(S.algorithms.mheCasadi, double(S.data.ysim(:,end)));    % in S.data.ysim(:,end-2) I save the velocities too
        updateInput(S.algorithms.mheCasadi, double(S.sensors.velocities(:,end)));
    else
        updateMeasurement(S.algorithms.ekfCasadi, double(S.data.ysim(1:end-2,end)));
        updateInput(S.algorithms.ekfCasadi, double(S.sensors.velocities(:,end)));
    end
end

function S = call_ESTIMATOR(S)
    tic;
%     if S.config.mhe
%         solve(S.algorithms.mheCasadi);
%         q_k = S.algorithms.mheCasadi.q_k;
%     else
%         solve(S.algorithms.ekfCasadi);
%         q_k = S.algorithms.ekfCasadi.x_k;
%     end
q_k = S.data.xsim(:,end);
    S.data.xest = [S.data.xest, q_k];
    S.exec_time.t_mhe  = [S.exec_time.t_mhe, toc];
end

function S = call_PFA(S)
    tic;
    q_k             = S.data.xest(:,end);
    pos0_k          = q_k(2*S.config.N+2:2*S.config.N+3);
    % =====================================================================
    % Search for obstacles in the neighbourhood and pass coordinates in
    % case of finding any
    % NEED TO BE REIMPLEMENTED FOR ACADO SOLVERS!!!! **********************
    if isempty(S.path.obstacles)
        setXYRobstacle(S.algorithms.mpcCasadi, [1e6, 1e6, 1]); % no obstacle found
    else
        dis             = pdist2(S.path.obstacles(:,1:2),pos0_k');
        [~, indx]       = sort(dis);
%         nroObs          = numel(indx);

        S.path.obstacles = S.path.obstacles(indx,:);
        S.path.obstacles = S.path.obstacles(1:min(3,size(S.path.obstacles,1)),:);

        nroObs          = size(S.path.obstacles,1);

        if nroObs<3
            setXYRobstacle(S.algorithms.mpcCasadi, [S.path.obstacles(indx,:);repmat([1e6, 1e6, 1],3-nroObs,1)])
        else
            setXYRobstacle(S.algorithms.mpcCasadi, S.path.obstacles);
        end
    end
    %
    S.controller.ref = S.algorithms.mpcCasadi.qRefOpt1;
    %
% if indx >= length(S.path.ref)-1
%     S.path.last_tgt             = true;
%     S.algorithms.mpc.last_tgt   = 1;
% end
    %
    S.exec_time.t_pfa = [S.exec_time.t_pfa, toc];
end

function S = call_SYSTEM(S)
    if S.config.SIM == true
        if S.config.mpc
            simulation_input.x = S.data.xsim(:,end);
            simulation_input.u = S.algorithms.Controls(:,end);%S.algorithms.mpcCasadi.u_k;            
        else
%             simulation_input.x = [S.data.xsim(1:end-2,end); S.Michalek2017.u0];
%             simulation_input.u = [0;0];
            simulation_input.x = S.data.xsim(:,end);
            simulation_input.u = S.Michalek2017.u0;
        end
        if S.config.slip == 1
            simulation_input.slip = zeros(2*S.config.N+1+2*(S.config.N+1),1); %+7
        else
            q1                      = S.dynamic.FNt(simulation_input.x,diag(S.config.slip-[1;1])*simulation_input.u);
            q2                      = S.dynamic.FNt(simulation_input.x,[0; 0]);
            simulation_input.slip   = q1-q2;
        end
        if sum(S.config.procDist.amp) > 0
            if strcmp(S.config.procDist.type,'gaussian') || strcmp(S.config.procDist.type,'normal')
                simulation_input.dist = S.config.procDist.amp.*randn(2*S.config.N+1+2*(S.config.N+1),1); % +7
            elseif strcmp(S.config.procDist.type,'uniform')
                simulation_input.dist = S.config.procDist.amp.*(2.*randn(2*S.config.N+1+2*(S.config.N+1),1)-1); % +7
            end        
        else
            simulation_input.dist = zeros(2*S.config.N+1+2*(S.config.N+1),1); %+7
        end
        S.data.procDist = [S.data.procDist, simulation_input.dist];
        S.data.slip     = [S.data.slip, full(simulation_input.slip)];
        %
%         states      = S.dynamic.FNt(simulation_input.x, simulation_input.u) + simulation_input.slip + simulation_input.dist;
        states      = FNt(full(simulation_input.x), full(simulation_input.u)) + simulation_input.slip + simulation_input.dist;   % generated mex file
        S.data.xsim = [S.data.xsim, full(states)];
    else
        S = write_ctrl(S, S.algorithms.mpc.Husky_w0, S.algorithms.mpc.Husky_v0);
    end
end

function S = call_CONTROLLER(S)
    if S.config.vNr < S.config.vNTarget
        S.config.vNr = S.config.vNr + S.config.deltavN;
        setVref(S.algorithms.mpcCasadi,S.config.vNr);
    end
    if S.config.mpc        
        S = call_MPC(S);
    else
        S = call_RSNTMichalek(S);
    end
end

function S = call_MPC(S)
    tic;
    % Solve control problem _______________________________________________
    setq0(S.algorithms.mpcCasadi,S.data.xest(:,end));
    solve(S.algorithms.mpcCasadi);
    if S.config.iters > S.config.Ne+1        
        S.algorithms.Controls  = [S.algorithms.Controls, S.algorithms.mpcCasadi.u_k];
    else
        S.algorithms.Controls  = [S.algorithms.Controls, [0;0]];
    end
    % To the last estimated velocties, add the last acceleration computed _
%     S.algorithms.Controls_w0 = [S.algorithms.Controls_w0, S.data.xest(2*S.config.N+6,end)+S.algorithms.Controls(1,end)];
%     S.algorithms.Controls_v0 = [S.algorithms.Controls_v0, S.data.xest(2*S.config.N+7,end)+S.algorithms.Controls(2,end)];
    S.algorithms.Controls_w0 = [S.algorithms.Controls_w0, S.algorithms.Controls(1,end)];
    S.algorithms.Controls_v0 = [S.algorithms.Controls_v0, S.algorithms.Controls(2,end)];
    % Velocities to apply to the Husky ____________________________________
    S.algorithms.mpc.Husky_w0       = S.algorithms.Controls_w0(end);
    S.algorithms.mpc.Husky_v0       = S.algorithms.Controls_v0(end);    
    %
    S.exec_time.t_mpc              = [S.exec_time.t_mpc, toc];
    %
    S.data.references = [S.data.references, S.algorithms.mpcCasadi.qRefOpt1];
end

function S = init_obsDetection(S)    
    %
    S.obs.pcMerge.gridSize        = 0.1;      % size of the voxels for averaging points
    S.obs.pcMerge.nroClouds       = 3;
    S.obs.pcSegdist.minNroPoints  = 15;       % minimum number of points of a cluster
    S.obs.pcSegdist.maxNroPoints  = inf;      % maximum number of points of a cluster
    S.obs.pcSegdist.minDistance   = 0.8;        % min disance between point of different clusters
    S.obs.filterWalls.area        = 2;        % area= deltaX+h. Clusters with area greater than this parameter are discarded
    S.obs.filterWalls.height      = 2;        % cluster with a height greater than this parameter are discarded
    S.obs.filterGround.disDev     = 0.2;      % deviation of points to the plane that is considered as ground
    S.obs.filterGround.angleDev   = 5;        % angular deviation of the normal vector to the plance considered as ground
    %
    S.obs.newObstacle.minDis      = 1;        % when an obstacle is found, if it is near to some other, they are merged.
    %
    S.obs.remMovingObs.searchRad  = 5;        % radius to find old obstacles
    S.obs.remMovingObs.disBetObs  = 1;        % when a new obstacle is found, and it is near to some obstacle from the list, this one is replaced by the new, assuming that it is a moving obstacle
    S.obs.remMovingObs.lossVotes  = 2;        % when osbtacles from the list are not found in the search radius, they loss votes
    %
    S.obs.validateObs.votes       = 10;       % minimum number of votes to consider a cluster as an obstacle, every time the same cluster is found, it gains a vote
    S.obs.validateObs.timeIndx    = 30;       % an obstacle with few votes and that it is not detected since a "long" time ago, is dicarded
    % variables to store pointclouds and data
    S.obs.ptCloud                 = [];
    S.obs.ptCloudAvg              = [];
    S.obs.ptCloudMap              = [];
    S.obs.Obstacles               = [];
    S.obs.ObstaclesAux            = [];
    %
    if ~exist('S.vlp16','var')
        S.vlp16 = velodynelidar('VLP16');              
    end
%     for i=1:S.obs.pcMerge.nroClouds
%         [S.obs.ptCloud, S.obs.timeStamps] = read(S.vlp16, 'latest');
%         S.obs.ptCloudAvg = [S.obs.ptCloud, S.obs.ptCloudAvg];
%     end
end

function S = call_OBSDETECTOR(S)
    tic;        
    if S.config.lidarOn
        if S.vlp16.NumPointCloudsAvailable > 0
            [S.obs.ptCloud, S.obs.timeStamps] = read(S.vlp16,'all');
            S.sensors.vlp16 = [S.sensors.vlp16, {S.obs.ptCloud}];
            if S.config.obsDetection
                S.obs.ptCloud = S.obs.ptCloud(end);
                if S.config.SIM
                    xy0     = S.data.xest(2*S.config.N+2:2*S.config.N+3,end);
                    theta0  = S.data.xest(S.config.N+1,end);
                else
                    xy0     = [S.ROS.sensors.rtk.x0; S.ROS.sensors.rtk.y0];
                    theta0  = S.ROS.sensors.vectornav_theta0;
                end
                % Remove ground
                [~,~,outlierIndicesW,~]             = pcfitplane(S.obs.ptCloud,S.obs.filterGround.disDev,[0 0 1],S.obs.filterGround.angleDev);
                S.obs.ptCloudAux                    = select(S.obs.ptCloud,outlierIndicesW);
                % Select ROI with LiDAR at centre
                S.obs.xyz                           = S.obs.ptCloudAux.Location(abs(S.obs.ptCloudAux.Location(:,1))<=S.config.lidarLimits.X & S.obs.ptCloudAux.Location(:,2)<=S.config.lidarLimits.Y & S.obs.ptCloudAux.Location(:,3)<=S.config.lidarLimits.Z,:);
                S.obs.ptCloudAux                    = pointCloud(S.obs.xyz);
                [labels,numClusters]                = pcsegdist(S.obs.ptCloudAux,S.obs.pcSegdist.minDistance,'NumClusterPoints',[S.obs.pcSegdist.minNroPoints,S.obs.pcSegdist.maxNroPoints]);
                S.obs.localObstacles                = findObstacles(S.obs.ptCloudAux, xy0, labels, numClusters, 0); % last argument passed as a constant
                S.obs.globalObstacles               = localToGlobalObs(S,S.obs.localObstacles,theta0,xy0);
                %
% S.path.obstacles = [];    
                for i=1:size(S.obs.globalObstacles,1)    
                    S = gen_obstacle(S,S.obs.globalObstacles(i,1),S.obs.globalObstacles(i,2),max(S.obs.globalObstacles(i,3:4)));
                end
            end
        end       
    end
    S.exec_time.t_obsdetector = [S.exec_time.t_obsdetector, toc];    
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

function obs = findObstacles(ptCloudAux, posH, labels, numClusters, indx)
obs         = [];
obsToDel    = [];
% ptCloudAvg = [];
    for i=1:numClusters
        obstacle_i      = find(labels == i);
        pt_obstacle     = select(ptCloudAux,obstacle_i);
% if i==1
%     ptCloudAvg = pt_obstacle;
% else
%     ptCloudAvg = pcmerge(ptCloudAvg, pt_obstacle);
% end

        x_mean          = mean(pt_obstacle.Location(:,1));
        y_mean          = mean(pt_obstacle.Location(:,2));
%         disGlo          = norm([x_mean, y_mean]);
        deltaX          = (pt_obstacle.XLimits(2)-pt_obstacle.XLimits(1));
        deltaY          = (pt_obstacle.YLimits(2)-pt_obstacle.YLimits(1));
        deltaZ          = (pt_obstacle.ZLimits(2)-pt_obstacle.ZLimits(1));        
%         h               = deltaZ;
%         area            = deltaX*h;
%         rho             = (pt_obstacle.Location(:,1)-mean(pt_obstacle.Location(:,1)))'*(pt_obstacle.Location(:,2)-mean(pt_obstacle.Location(:,2)))/length(pt_obstacle.Location)^2;% assumming here equiprobability of each point...
%         disH            = norm([x_mean-posH(1), y_mean-posH(2)]);
%         newObs          = [x_mean, y_mean, disGlo, deltaX, deltaY, rho, area, h, disH, indx, 0];
        newObs          = [x_mean, y_mean, deltaX, deltaY];
        obs             = [obs; newObs];
    end
    obsToDel        = unique(obsToDel);
    obs(obsToDel,:) = [];
% if ~isempty(ptCloudAvg)    
%     pcshow(ptCloudAvg)   
% end

end

function Obstacles = remSporaidcObs(globalObstacles, Obstacles, posH, theta, R, distanceBetObs, lossVotes)
    if isempty(Obstacles)
        return;
    end
    if ~isempty(globalObstacles)
        nroGloObs = size(globalObstacles,1);
    else
        nroGloObs = 1;
    end
    
    nroObs  = size(Obstacles,1);    
        
    for i=1:nroGloObs
        for j=1:nroObs
            posObs      = Obstacles(j,1:2)';
            if ~isempty(globalObstacles)
                posNewObs = globalObstacles(i,1:2)';
                disBetObs = sqrt((posNewObs(1)-posObs(1))^2 + (posNewObs(2)-posObs(2))^2);
            end                        
            vecObs      = posObs-posH;
            vecObs      = vecObs ./ norm(vecObs);
            dotProd     = cos(theta)*vecObs(1)+sin(theta)*vecObs(2);
            r           = sqrt((posObs(1)-posH(1))^2 + (posObs(2)-posH(2))^2);
            if r<=R && dotProd>0
                if isempty(globalObstacles)
                    Obstacles(j,11) = Obstacles(j,11) - lossVotes;
                    Obstacles(j,11) = min(Obstacles(j,11), 0);
                elseif disBetObs <= distanceBetObs
                    x_mean  = globalObstacles(i,1);
                    y_mean  = globalObstacles(i,2);
                    disGlo  = sqrt(x_mean^2+y_mean^2);
                    deltaX  = globalObstacles(i,4);
                    deltaY  = globalObstacles(i,5);
                    rho     = globalObstacles(i,6);
                    h       = globalObstacles(i,8);
                    area    = deltaX*h;
                    indx    = globalObstacles(i,10);
                    disH    = min([globalObstacles(i,9), Obstacles(j,9)]);
                    votes   = Obstacles(j,11);
                    newObs  = [x_mean, y_mean, disGlo, deltaX, deltaY, rho, area, h, disH, indx, votes];
                    Obstacles(j,:) = newObs;
                end
            end
        end
    end
end

function listObsAux = checkIfNewObs(listObsAux, recentFoundObs, minDis)
    if isempty(listObsAux)
        listObsAux = recentFoundObs;
        return;
    end
    if isempty(recentFoundObs)
        return;
    end    

    nroObsAux   = size(listObsAux,1);
    nroFoundAux = size(recentFoundObs,1);
    
    indxNewObs  = [];
    for i=1:nroFoundAux % recently found  
        countMatchs = 0;
        indxMatches = [];
        for j=1:nroObsAux % elements in the auxiliar list
            dij = sqrt((listObsAux(j,1)-recentFoundObs(i,1))^2 + (listObsAux(j,2)-recentFoundObs(i,2))^2);
            if dij <= minDis
                listObsAux(j,11) = listObsAux(j,11)+1;
                countMatchs      = countMatchs+1;
                indxMatches      = [indxMatches, [i;j]];
            else
                indxNewObs       = [indxNewObs; i];
            end
        end        
        if ~isempty(indxMatches)
            x_mean  = mean(listObsAux(indxMatches(2,:)',1));
            y_mean  = mean(listObsAux(indxMatches(2,:)',2));
            disGlo  = sqrt(x_mean^2+y_mean^2);
            deltaX  = max(listObsAux(indxMatches(2,:)',4));
            deltaY  = max(listObsAux(indxMatches(2,:)',5));
            rho     = mean(listObsAux(indxMatches(2,:)',6));
            h       = max(listObsAux(indxMatches(2,:)',8));
            area    = deltaX*h;
            indx    = max(listObsAux(indxMatches(2,:)',10));
            disH    = min(listObsAux(indxMatches(2,:)',9));
            votes   = max(listObsAux(indxMatches(2,:)',11));
            newObs  = [x_mean, y_mean, disGlo, deltaX, deltaY, rho, area, h, disH, indx, votes];
            listObsAux(indxMatches(2,:),:) = [];
            listObsAux = [listObsAux; newObs];
            nroObsAux  = size(listObsAux,1);
        end
    end
    if ~isempty(indxNewObs)
        indxNewObs = unique(indxNewObs);
        listObsAux = [listObsAux; recentFoundObs(indxNewObs,:)];  
    end
    % Sort list
    [~,indx] = sort(listObsAux(:,11));
    indx = flipud(indx);
    listObsAux = listObsAux(indx,:);
end

function StatObs = remWalls(Obstacles, areaThr, heightThr)
    StatObs  = Obstacles;
    if isempty(Obstacles)
        return;
    end    

    area     = StatObs(:,7);
    h        = StatObs(:,8);
    I        = length(area);

    obsToDel = [];
    for i=1:I
        if area(i) > areaThr || h(i) >= heightThr % then remove the farest obstacle
            obsToDel = [obsToDel; i];
        end
    end
% a= obsToDel    
    StatObs(obsToDel,:) = [];
end

function S = call_RSNTMichalek(S)
    tic;
    if S.config.iters >= S.config.iterCorrecOn
        setCorrection(S.Michalek2017, true);
    end
    % Solve control problem _______________________________________________
    S.path.references_mhempc    = [S.path.references_mhempc, S.controller.ref];
    update_q(S.Michalek2017, S.data.xest(:,end));
    set_thetaNr(S.Michalek2017, S.controller.ref(2*S.config.N+1));
    set_xNr(S.Michalek2017, S.controller.ref(2*S.config.N+4));
    set_yNr(S.Michalek2017, S.controller.ref(2*S.config.N+5));  

    solve(S.Michalek2017);

    S.algorithms.Controls  = [S.algorithms.Controls, full(S.Michalek2017.u0)];
%     % To the last estimated velocties, add the last acceleration computed _
%     S.algorithms.Controls_w0 = [S.algorithms.Controls_w0, S.data.xest(2*S.config.N+6,end)+S.algorithms.Controls(1,end)];
%     S.algorithms.Controls_v0 = [S.algorithms.Controls_v0, S.data.xest(2*S.config.N+7,end)+S.algorithms.Controls(2,end)];
%     % Velocities to apply to the Husky ____________________________________
%     S.algorithms.mpc.Husky_w0       = S.algorithms.Controls_w0(end);
%     S.algorithms.mpc.Husky_v0       = S.algorithms.Controls_v0(end);
    S.exec_time.t_mpc              = [S.exec_time.t_mpc, toc];
end

function S = compute_curvature(S)
    S.path.curvature = [];
    for i=3:length(S.path.coordinates)
        x1 = S.path.coordinates(1,i);
        x2 = S.path.coordinates(1,i-1);
        x3 = S.path.coordinates(1,i-2);
        y1 = S.path.coordinates(2,i);
        y2 = S.path.coordinates(2,i-1);
        y3 = S.path.coordinates(2,i-2);
        %
        a = sqrt((x1-x2)^2 + (y1-y2)^2);
        b = sqrt((x2-x3)^2 + (y2-y3)^2);
        c = sqrt((x3-x1)^2 + (y3-y1)^2);
        s = (a+b+c)/2;
        A = sqrt(s*(s-a)*(s-b)*(s-c));
        %
        S.path.curvature = [S.path.curvature, 4*A/(a*b*c)];
    end
%     figure;plot3(S.path.coordinates(1,3:end),S.path.coordinates(2,3:end),1./S.path.curvature,'r'); grid on; zlim([0 10]);
end

function S = gen_path(S,type)
    if strcmp(type,'infinity')
        a = 4;
        c = 8;
        b = 1;
        t = 0:0.0005:2*pi;
        x = (a*sqrt(2).*cos(t))./(sin(t).^2+1);
        y = (c*sqrt(2).*cos(t).*sin(t))./(sin(t).^2 + b);
        %
        S.path.coorection_x = 5.5;  % correction for the field experiemnts in order to fit the path in my locala reference frame
        S.path.coorection_y = 4.5;
        S.path.coordinates  = [x+S.path.coorection_x;y+S.path.coorection_y];
        %
        S.path.s            = casadi.MX.sym('s');
        S.path.ds           = S.config.Ts;
        %
        S.path.fx           = (a*sqrt(2).*cos(S.path.s))./(sin(S.path.s).^2+1) + S.path.coorection_x;
        S.path.fy           = (c*sqrt(2).*cos(S.path.s).*sin(S.path.s))./(sin(S.path.s).^2 + b) + S.path.coorection_y;
        %
        fx_fun              = casadi.Function('fx_fun',{S.path.s},{S.path.fx});
        fy_fun              = casadi.Function('fy_fun',{S.path.s},{S.path.fy});
        %
        dfxdt               = fx_fun.jacobian;
        dfydt               = fy_fun.jacobian;
        % Numerical approximation of the function that determines the att.
        % val.
        di                  = 0.01;
        I                   = -pi:di:3*pi;
        alpha0              = atan((fy_fun(di)-fy_fun(0))/(fx_fun(di)-fx_fun(0)));
        theta               = alpha0;
        integrando          = [];
        args                = [];
        
        for i=I
            arg         = (dfydt(i,[])*dfxdt(i+di,[]) - dfxdt(i,[])*dfydt(i+di,[])) / (dfxdt(i,[])*dfxdt(i+di,[]) + dfydt(i,[])*dfydt(i+di,[]));
            args        = [args, arg];
        %     integrando  = [integrando, atan(arg)];
            integrando  = [integrando, atan2((dfydt(i,[])*dfxdt(i+di,[]) - dfxdt(i,[])*dfydt(i+di,[])), (dfxdt(i,[])*dfxdt(i+di,[]) + dfydt(i,[])*dfydt(i+di,[])))];
            theta       = [theta, theta(end)+integrando(end)];
        end               
        % Polyfit
        p = polyfit(I,full(theta(1:end-1)),20);
        S.path.st            = casadi.MX.sym('st');
        S.path.ftheta        = [];%sum(S.path.s.*p);
        for i=1:length(p)
            if i==1
                S.path.ftheta = -p(i)*S.path.st^(length(p)-i);
            else
                S.path.ftheta = S.path.ftheta - p(i)*S.path.st^(length(p)-i);
            end            
        end
        %
        S.path.r      = [ S.path.fx;  
                          S.path.fy;
                          S.path.ftheta];
        %
        S = compute_curvature(S);
    elseif strcmp(type,'flat_infinity')
        a = 7;
        b = 7.5;
        t = 0:0.0005:2*pi;
        x = a*sin(t);
        y = b*sin(t).^2.*cos(t);
        %
        S.path.coorection_x = 6;  % correction for the field experiemnts in order to fit the path in my locala reference frame
        S.path.coorection_y = 4.5;
        S.path.coordinates  = [x+S.path.coorection_x;y+S.path.coorection_y];
        %
        S.path.s            = casadi.MX.sym('s');
        S.path.ds           = S.config.Ts;
        %
        S.path.fx           = a*sin(S.path.s) + S.path.coorection_x;
        S.path.fy           = b*sin(S.path.s).^2*cos(S.path.s) + S.path.coorection_y;
        %
        % fx_fun              = casadi.Function('fx_fun',{S.path.s},{S.path.fx});
        % fy_fun              = casadi.Function('fy_fun',{S.path.s},{S.path.fy});
        %
        % dfxdt               = fx_fun.jacobian;
        % dfydt               = fy_fun.jacobian;
% 
        % % Numerical approximation of the function that determines the att.
        % % val.
        % di                  = 0.01;
        % I                   = -pi:di:3*pi;
        % alpha0              = atan((fy_fun(di)-fy_fun(0))/(fx_fun(di)-fx_fun(0)));
        % theta               = alpha0;
        % integrando          = [];
        % args                = [];
        % 
        % for i=I
        %     arg         = (dfydt(i,[])*dfxdt(i+di,[]) - dfxdt(i,[])*dfydt(i+di,[])) / (dfxdt(i,[])*dfxdt(i+di,[]) + dfydt(i,[])*dfydt(i+di,[]));
        %     args        = [args, arg];
        % %     integrando  = [integrando, atan(arg)];
        %     integrando  = [integrando, atan2((dfydt(i,[])*dfxdt(i+di,[]) - dfxdt(i,[])*dfydt(i+di,[])), (dfxdt(i,[])*dfxdt(i+di,[]) + dfydt(i,[])*dfydt(i+di,[])))];
        %     theta       = [theta, theta(end)+integrando(end)];
        % end               
        % % Polyfit
        % p = polyfit(I,full(theta(1:end-1)),20);
        S.path.st            = casadi.MX.sym('st');
        % S.path.ftheta        = [];%sum(S.path.s.*p);
        % for i=1:length(p)
        %     if i==1
        %         S.path.ftheta = -p(i)*S.path.st^(length(p)-i);
        %     else
        %         S.path.ftheta = S.path.ftheta - p(i)*S.path.st^(length(p)-i);
        %     end            
        % end
        %
S.path.ftheta        = 0;%sum(S.path.s.*p);        
        S.path.r      = [ S.path.fx;  
                          S.path.fy;
                          S.path.ftheta];
        %
        S = compute_curvature(S);
    elseif strcmp(type ,'segment')
        x0 = 1;
        xf = 23;
        dt = 0.01;
        t = 0:dt:2*pi;
        x = x0 + (xf-x0) * t / (2*pi);
        y0 = 6;
        y = repmat(y0,1,length(x));
        %
        S.path.coordinates  = [x ; y];
        %
        S.path.s            = casadi.MX.sym('s');
        S.path.st           = casadi.MX.sym('st');
        S.path.ds           = S.config.Ts;
        %
        S.path.fx           = x0 + (xf-x0) * S.path.s / (2*pi);
        S.path.fy           = y0;
        %
        fx_fun              = casadi.Function('fx_fun',{S.path.s},{S.path.fx});
        fy_fun              = casadi.Function('fy_fun',{S.path.s},{S.path.fy});
        % Numerical approximation of the function that determines the att.
        % val.
        di                  = 0.01;
        alpha0              = atan((fy_fun(di)-fy_fun(0))/(fx_fun(di)-fx_fun(0)));
        theta               = alpha0;
        S.path.ftheta       = theta;
        %
        S.path.r            = [ S.path.fx;  
                                S.path.fy;
                                S.path.ftheta];
        %
        S = compute_curvature(S);
    elseif strcmp(type ,'invert_segment')
        x0 = 23;
        xf = 1;
        dt = 0.01;
        t = 0:dt:2*pi;
        x = x0 + (xf-x0) * t / (2*pi);
        y0 = 6;
        y = repmat(y0,1,length(x));
        %
        S.path.coordinates  = [x ; y];
        %
        S.path.s            = casadi.MX.sym('s');
        S.path.st           = casadi.MX.sym('st');
        S.path.ds           = S.config.Ts;
        %
        S.path.fx           = x0 + (xf-x0) * S.path.s / (2*pi);
        S.path.fy           = y0;
        %
        fx_fun              = casadi.Function('fx_fun',{S.path.s},{S.path.fx});
        fy_fun              = casadi.Function('fy_fun',{S.path.s},{S.path.fy});
        % Numerical approximation of the function that determines the att.
        % val.
        di                  = 0.01;
        alpha0              = atan((fy_fun(di)-fy_fun(0))/(fx_fun(di)-fx_fun(0)));
        theta               = alpha0;
        S.path.ftheta       = theta;
        %
        S.path.r            = [ S.path.fx;  
                                S.path.fy;
                                S.path.ftheta];
        %
        S = compute_curvature(S);
    elseif strcmp(type,'rectangular')
        t = 0:0.01:2*pi;
        p = 3;
        q = 2.5;

        x = p.*(sqrt(cos(t).*cos(t)).*cos(t) + sqrt(sin(t).*sin(t)).*sin(t));
        y = q.*(sqrt(cos(t).*cos(t)).*cos(t) - sqrt(sin(t).*sin(t)).*sin(t));

        S.path.coorection_x = p+4;  % correction for the field experiemnts in order to fit the path in my locala reference frame
        S.path.coorection_y = q+1.5;
        S.path.coordinates  = [x+S.path.coorection_x;y+S.path.coorection_y];
        %
        S.path.s            = casadi.MX.sym('s');
        S.path.ds           = S.config.Ts;
        %
        S.path.fx           = p.*(sqrt(cos(S.path.s).*cos(S.path.s)).*cos(S.path.s) + sqrt(sin(S.path.s).*sin(S.path.s)).*sin(S.path.s)) + S.path.coorection_x;
        S.path.fy           = q.*(sqrt(cos(S.path.s).*cos(S.path.s)).*cos(S.path.s) - sqrt(sin(S.path.s).*sin(S.path.s)).*sin(S.path.s)) + S.path.coorection_y;
        %
        fx_fun              = casadi.Function('fx_fun',{S.path.s},{S.path.fx});
        fy_fun              = casadi.Function('fy_fun',{S.path.s},{S.path.fy});
        %
        dfxdt               = fx_fun.jacobian;
        dfydt               = fy_fun.jacobian;

        % Numerical approximation of the function that determines the att.
        % val.
        di                  = 0.01;
        I                   = -pi:di:3*pi;
        alpha0              = atan((fy_fun(di)-fy_fun(0))/(fx_fun(di)-fx_fun(0)));
        theta               = alpha0;
        integrando          = [];
        args                = [];
        
        for i=I
            arg         = (dfydt(i,[])*dfxdt(i+di,[]) - dfxdt(i,[])*dfydt(i+di,[])) / (dfxdt(i,[])*dfxdt(i+di,[]) + dfydt(i,[])*dfydt(i+di,[]));
            args        = [args, arg];
        %     integrando  = [integrando, atan(arg)];
            integrando  = [integrando, atan2((dfydt(i,[])*dfxdt(i+di,[]) - dfxdt(i,[])*dfydt(i+di,[])), (dfxdt(i,[])*dfxdt(i+di,[]) + dfydt(i,[])*dfydt(i+di,[])))];
            theta       = [theta, theta(end)+integrando(end)];
        end               
        % Polyfit
        p = polyfit(I,full(theta(1:end-1)),20);
        S.path.st            = casadi.MX.sym('st');
        S.path.ftheta        = [];%sum(S.path.s.*p);
        for i=1:length(p)
            if i==1
                S.path.ftheta = -p(i)*S.path.s^(length(p)-i);
            else
                S.path.ftheta = S.path.ftheta - p(i)*S.path.s^(length(p)-i);
            end            
        end
        %
        S.path.r      = [ S.path.fx;  
                          S.path.fy;
                          S.path.ftheta];
        %
        S = compute_curvature(S);
    elseif strcmp(type,'test')
        % TEST TRAJECTORY #################################################
%         t = 0:3500;
%         x = 3*sin(2*pi.*t/t(end)) + 10;
%         y = 6*cos(2*pi.*t/t(end)) + 10;
%         S.path.coordinates = [x;y];
%         %
        N  = (S.config.tf-S.config.t0)/S.config.Ts;
        t           = linspace(0,2*pi,N);
        x           = 10.*cos(t)+5.*cos(5.*t) - 10;
        y           = 10.*sin(t)+5.*sin(5.*t) + 5;
        S.path.coordinates = [x;y];
    elseif strcmp(type,'monaco')
        load('monaco_xy');
        y = y.*0.7;
        S.path.coordinates  = [x(1:545); y(1:545)];
        
    elseif strcmp(type,'test_betas')
        x = 2:0.1:7.5;
        y = repmat(0.5,1,length(x));

        t = -pi/2:0.1:pi/2;
        x1 = 1.5.*cos(t) + 7.5;
        y1 = 1.5.*sin(t) + 2;

        x = [x, x1];
        y = [y, y1];

        x1 = x(end):-0.1:5;
        y1 = repmat(y(end),1,length(x1));

        x = [x, x1];
        y = [y, y1];

        t = -pi/2:-0.1:-3*pi/2;
        x1 = 1.5.*cos(t) + 4.5;
        y1 = 1.5.*sin(t) + 5;

        x = [x, x1];
        y = [y, y1];

        x1 = x(end):0.1:7.5;
        y1 = repmat(y(end),1,length(x1));

        x = [x, x1];
        y = [y, y1];

        t = -pi/2:0.1:pi/2;
        x1 = 1.5.*cos(t) + 7.5;
        y1 = 1.5.*sin(t) + 8;

        x = [x, x1];
        y = [y, y1];

        x1 = x(end):-0.1:5;
        y1 = repmat(y(end),1,length(x1));

        x = [x, x1];
        y = [y, y1];

        t = pi/2:0.1:3*pi/2;
        x1 = 4.5.*cos(t) + 5;
        y1 = 4.5.*sin(t) + 5;

        x = [x, x1];
        y = [y, y1];

        d = -max(S.system.Lh(:,2))*S.config.N*log(tan(0.1/2));
        y1 = y(1)-d:0.1:y(1);
        x1 = x(1).*ones(size(y1));

        y_aux = [y1,y];
        S.path.coordinates = [x1,x;y_aux];  
        
    elseif strcmp(type,'rectangular_terrace') || strcmp(type,'complex1_terrace')
        if strcmp(type,'rectangular_terrace')
%            bag = rosbag('/home/kernighan/2022-07-19-15-32-59.bag');             
            x = repmat(8.2,1,101);
            y = 20:0.1:30;
            x = [x, 8.2:-0.1:1.0];
            y = [y, repmat(30, 1, ceil((8.2-1.0)/0.1)+1)];
            x = [x, repmat(0.8,1,101)];
            y = [y, 30:-0.1:20];
            x = [x, 0.8:0.1:8.2];
            y = [y, repmat(20, 1, ceil((8.2-0.8)/0.1)+1)];
           
            d = -max(S.system.Lh(:,2))*S.config.N*log(tan(0.1/2));
            y1 = y(1)-d:0.1:y(1);
            x1 = x(1).*ones(size(y1));

            y_aux = [y1,y];
            y_aux = y_aux+1.25;
            S.path.coordinates = [x1,x;y_aux];                      
           
        elseif strcmp(type,'complex1_terrace')
           bag = rosbag('/home/nahuel/Dropbox/PosDoc/AC3E/NMHE for incremental encoders/2022-07-19-15-36-09.bag'); 
%            bag = rosbag('/home/kernighan/Documents/mhe-mpc-for-N-trailers/Simulaciones/ACADO/2022-08-10-17-00-40.bag');
%            bag = rosbag('/home/nahuel/Dropbox/PosDoc/AC3E/NMHE-NMPC-for-N-trailer/mhe-mpc-for-N-trailers/Simulaciones/ACADO/2022-07-19-15-36-09.bag');            
%             bag         = rosbag('/home/nahuel/Dropbox/PosDoc/AC3E/NMHE for incremental encoders//test_betas_3.bag');            
            bSel        = select(bag, 'Topic', 'fix');
            msgStructs  = readMessages(bSel);

            lat         = NaN(1,bSel.NumMessages);
            lon         = NaN(1,bSel.NumMessages);
            xraw        = NaN(1,bSel.NumMessages);
            yraw        = NaN(1,bSel.NumMessages);
            xcor        = NaN(1,bSel.NumMessages);
            ycor        = NaN(1,bSel.NumMessages);
            % ESQUINA DE LA AZOTEA TOMDA COMO ORIGEN DE COORDENADAS
            lat0        = -33.03422;%S.ROS.LAT0; %   = -33.03422;
            lon0        = -71.591885;%S.ROS.LON0; %  = -71.591885;
           
            %
            sdpvar a11 a12 a21 a22;
            sdpvar x1bar y1bar x2bar y2bar;
            % Besides reference point, two more are needed to obtaint the local
            % reference frame
            % Coordinates of point (x, 0)
            lat_local_coord_1 = -33.03421;%S.ROS.local_coord_1.lat;
            long_local_coord_1 = -71.591825;%S.ROS.local_coord_1.long;

            [x1, y1]    = latlon2xy(lat_local_coord_1, long_local_coord_1, lat0, lon0);
            x1          = x1*1000;
            y1          = y1*1000;
            x           = norm([x1 y1]);
            % Coordinates of point (0, y)
            lat_local_coord_2 = -33.034158333;%S.ROS.local_coord_2.lat;
            long_local_coord_2 = -71.591905;%S.ROS.local_coord_2.long;

            [x2, y2]    = latlon2xy(lat_local_coord_2, long_local_coord_2, lat0, lon0);
            x2          = x2*1000;
            y2          = y2*1000;
            y           = norm([x2 y2]);
            %
            A           = [a11 a12; a21 a22];
            v           = [x1; y1; x2; y2];
            b           = [x1bar; y1bar; x2bar; y2bar];
            Constraints = [[A zeros(2); zeros(2) A]*v - b == zeros(4,1); x1bar*x2bar + y1bar*y2bar == 0];
            %
            obj         = (x1bar - x)^2 + (y2bar - y)^2;
            % 
            optimize(Constraints, obj);
            %
            Mtx   = value(A);

            %

            for i=1:bSel.NumMessages
                lat(i)  = msgStructs{i}.Latitude;
                lon(i)  = msgStructs{i}.Longitude;
                [xm,ym] = latlon2xy(lat(i),lon(i),lat0,lon0);
                xraw(i) = xm*1000;
                yraw(i) = ym*1000;
                xy      = Mtx*[xraw(i);yraw(i)];
                xcor(i) = xy(1);
                ycor(i) = xy(2);
            end

            x = smooth(xcor,35)';
            y = smooth(ycor,35)';

            S.path.coordinates = [x; y];            
        end
%         for i=1:15
%             xcor = smooth(xcor,'lowess')';
%             ycor = smooth(ycor,'lowess')';
%         end
        
    elseif strcmp(type,'circular_terrace')           
        a = 4;
        b = 4;
%         b = 0.4;
        t = 0:0.0005:2*pi;        
        x = a*cos(t);
        y = b*sin(t);
        %
        S.path.coorection_x = 5.5;
        S.path.coorection_y = 4.5;
        S.path.coordinates  = [x+S.path.coorection_x;y+S.path.coorection_y];
        %
        S.path.s            = casadi.MX.sym('s');
        S.path.ds           = S.config.Ts;
        %
        S.path.fx           = a*cos(S.path.s)+S.path.coorection_x;
        S.path.fy           = b*sin(S.path.s)+S.path.coorection_y;
        %
        S.path.st           = casadi.MX.sym('st');
        S.path.ftheta       = pi/2 + (1/S.config.Ts)*S.path.st * atan(S.path.ds);
        S.path.r            = [ S.path.fx;  
                                S.path.fy;
                                S.path.ftheta];
        %
        S = compute_curvature(S);
    elseif strcmp(type,'ellipsoidal1_terrace')           
           y = 18:0.1:23.5;
           x = repmat(8,1,length(y));
           
           t = 0:0.01:0.75;
           x = [x, 4.5 + 3.5.*cos(2*pi*t)];
           y = [y, 23.5 + 2.5.*sin(2*pi*t)];
           
           x1 = x(end):0.05:5.5;
           y1 = repmat(y(end),1,length(x1));
           
           x = [x,x1];
           y = [y,y1];
           
           y1 = y(end):-0.05:y(1)-2;
           x1 = repmat(x(end),1,length(y1));
           
           x = [x,x1];
           y = [y,y1];
           
           t = pi:0.001:2*pi;
           
           x1 = 1.25*cos(t) + 6.75;
           y1 = 1.25*sin(t) + 16;
           
           x = [x,x1];
           y = [y,y1];
           
           y1 = y(end):0.05:y(1);
           x1 = repmat(x(end),1,length(y1));
           
           x = [x, x1];
           y = [y, y1];
           
           x = x-1;
           y = y+4;

           S.path.coordinates = [x; y];
    elseif strcmp(type,'ellipsoidal2_terrace')           
           y = 29:-0.1:23.5;
           x = repmat(8,1,length(y));
           
           t = 0:0.01:1;
           x = [x, 4.5 + 3.5.*cos(-2*pi*t)];
           y = [y, 23.5 + 2.5.*sin(-2*pi*t)];
           
           x = x;
           y = y+2;

           S.path.coordinates = [x; y];
    elseif strcmp(type,'agricultural')
            x0 = 5;
            y0 = 1;
            
            factor = 3;

            yinit = 0;
            
            d1 = 1;
            d2 = 0.5;
            d3 = 1;
            
            a1 = 2.5/factor;
            b1 = 2.5/factor;
            a2 = 3.5/factor;
            b2 = 2.5/factor;
            a3 = 2.5/factor;
            b4 = 2.5/factor;
            
            h  = 10/factor;
            
            delta = 0.00001;
            alpha1 = 25*pi/180;
            alpha2 = 35*pi/180;
            alpha3 = 25*pi/180;
            
            xc1 = x0+d1/2;
            yc1 = y0+h+b1;
            
            xc2 = x0+d2/2;
            yc2 = y0+b2;
            
            xc3 = x0+d1+d2+d3/2;
            yc3 = y0+h+b1;
            
            t1  = 3*pi/2-alpha1:-delta:-pi/2+alpha1;
            xe1 = xc1+a1*cos(t1);
            ye1 = yc1+b1*sin(t1);
            
            t2  = pi-(pi/2-alpha2):delta:2*pi+(pi/2-alpha2);
            xe2 = xc2+a2*cos(t2);
            xe2 = xe2 + xe1(end) - xe2(1);
            ye2 = yc2+b2*sin(t2);
            
            t3  = 3*pi/2-alpha3:-delta:-pi/2+alpha3;
            xe3 = xc3+a3*cos(t3);
            xe3 = xe3+xe2(end)-xe3(1); 
            ye3 = yc3+b4*sin(t3);
            
            y   = 2+h/2:delta:ye1(1);
            x   = repmat(xe1(1),1,length(y));
            
            path = [x, xe1;y, ye1];
            
            y = ye1(end):-delta:ye2(1);
            x = repmat(xe1(end),1,length(y));
            
            path = [path, [x,xe2;y, ye2]];
            
            y = ye2(end):delta:ye3(1);
            x = repmat(xe2(end),1,length(y));
            
            path = [path, [x,xe3;y, ye3]];
            
            y = ye3(end):-delta:yinit-delta;
            x = repmat(xe3(end),1,length(y));
            
            path = [path,[x;y]];
            
            x = path(1,end):-delta:path(1,1);
            y = repmat(path(2,end),1,length(x));
            
            path = [path,[x;y]];
            
            y = path(2,end):delta:path(2,1);
            x = repmat(path(1,end),1,length(y));
            
            path = [path,[x;y]];
            S.path.coordinates = path;
    else    
        % REAL TRAJECTORY #####################################################
        % First segment
        % Path ________________________________________________________________
        S.path.x0               = 5;
        S.path.y0               = 5;    
        S.path.R1               = 8;
        S.path.R2               = 6;
        S.path.R3               = 4;
        S.path.R4               = 3;
        S.path.R5               = 1;
        S.path.L1               = 20;
        S.path.L2               = 10;
        S.path.L3               = 2*(S.path.R1+S.path.R2+S.path.R3+S.path.R4+S.path.R5);
        S.path.deltaL           = 0.1;
        %
        S.path.coordinates = [repmat(S.path.x0, 1, S.path.L1/S.path.deltaL); S.path.y0:S.path.deltaL:S.path.L1+S.path.y0-S.path.deltaL];
        % First circunference
        N   = pi*S.path.R1/S.path.deltaL - S.path.deltaL;
        t   = pi:-1/N:0;
        xc1 = S.path.coordinates(1,end);
        yc1 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [xc1 + S.path.R1 + S.path.R1.*cos(t); yc1 + S.path.R1.*sin(t)]];
        % Second segment
        x2 = S.path.coordinates(1,end);
        y2 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [repmat(x2, 1, S.path.L2/S.path.deltaL); y2:-S.path.deltaL:y2-S.path.L2+S.path.deltaL] ];
        % Second circunference
        N   = pi*S.path.R2/S.path.deltaL - S.path.deltaL;
        t   = pi:1/N:2*pi;
        xc2 = S.path.coordinates(1,end);
        yc2 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [xc2 + S.path.R2 + S.path.R2.*cos(t); yc2 + S.path.R2.*sin(t)]];
        % Third segment
        x3 = S.path.coordinates(1,end);
        y3 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [repmat(x3, 1, S.path.L2/S.path.deltaL); y3:S.path.deltaL:y3+S.path.L2-S.path.deltaL] ];
        % Third circunference
        N   = pi*S.path.R3/S.path.deltaL - S.path.deltaL;
        t   = pi:-1/N:0;
        xc3 = S.path.coordinates(1,end);
        yc3 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [xc3 + S.path.R3 + S.path.R3.*cos(t); yc3 + S.path.R3.*sin(t)]];
        % Fourth segment
        x4 = S.path.coordinates(1,end);
        y4 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [repmat(x4, 1, S.path.L2/S.path.deltaL); y4:-S.path.deltaL:y4-S.path.L2+S.path.deltaL] ];
        % Fourth circunference
        N   = pi*S.path.R4/S.path.deltaL - S.path.deltaL;
        t   = pi:1/N:2*pi;
        xc4 = S.path.coordinates(1,end);
        yc4 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [xc4 + S.path.R4 + S.path.R4.*cos(t); yc4 + S.path.R4.*sin(t)]];
        % Fifth segment
        x5 = S.path.coordinates(1,end);
        y5 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [repmat(x5, 1, S.path.L2/S.path.deltaL); y5:S.path.deltaL:y5+S.path.L2-S.path.deltaL] ];
        % Fifth circunference
        N   = pi*S.path.R5/S.path.deltaL - S.path.deltaL;
        t   = pi:-1/N:0;
        xc5 = S.path.coordinates(1,end);
        yc5 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [xc5 + S.path.R5 + S.path.R5.*cos(t); yc5 + S.path.R5.*sin(t)]];
        % Sexth segment
        x6 = S.path.coordinates(1,end);
        y6 = S.path.coordinates(2,end);

        y6_aux = y6:-S.path.deltaL:S.path.y0-0.5;%+S.path.deltaL;
        %
        S.path.coordinates = [S.path.coordinates, [repmat(x6, 1, length(y6_aux)); y6_aux] ];
        % Merge with first point
        x8 = S.path.coordinates(1,end);
        y8 = S.path.coordinates(2,end);
        S.path.coordinates = [S.path.coordinates, [x8:-S.path.deltaL:S.path.x0; [repmat(y8, 1, (S.path.L3/S.path.deltaL)/2), repmat(S.path.y0, 1, (S.path.L3/S.path.deltaL)/2)]] ];
        %
        yval = S.path.coordinates(2,end);
        for i=length(S.path.coordinates):-1:2
            if S.path.coordinates(2,i) ~= yval
                S.path.coordinates(1,i+1) = S.path.coordinates(1,i);
                break;
            end
        end
        S.path.coordinates = [S.path.coordinates, S.path.coordinates(:,1)];
        %
        S.path.num_points_ttaj = length(S.path.coordinates);
%         % Gen path to intial point x0
%         x9 = S.path.coordinates(1,end);
%         y9 = S.path.coordinates(2,end);
%         trajX = x9-S.path.deltaL:-S.path.deltaL:S.init_condition.x0(2*S.config.N+2)-2*(1+S.config.N)*(S.system.Lh1+S.system.L1);
%         trajY = y9.*ones(size(trajX));
%         S.path.coordinates = [S.path.coordinates, [trajX; trajY]];
%         %
%         x10 = S.path.coordinates(1,end);
%         y10 = S.path.coordinates(2,end);    
%         trajY = y10-S.path.deltaL:-S.path.deltaL:S.init_condition.x0(2*S.config.N+3);
%         trajX = x10.*ones(size(trajY));
%         S.path.coordinates = [S.path.coordinates, [trajX; trajY]];
%         %
%         x11 = S.path.coordinates(1,end);
%         y12 = S.path.coordinates(2,end);
%         trajX = x11+S.path.deltaL:S.path.deltaL:S.init_condition.x0(2*S.config.N+2);
%         trajY = y12.*ones(size(trajX));
%         S.path.coordinates = [S.path.coordinates, [trajX; trajY]];
        % Remove inconsistencies
        pos_to_del = [];
        for i=2:length(S.path.coordinates)
            dx = S.path.coordinates(1,i)-S.path.coordinates(1,i-1);
            dy = S.path.coordinates(2,i)-S.path.coordinates(2,i-1);
            if dx==0 && dy ==0
                pos_to_del = [pos_to_del,i];                
            end
        end
        S.path.coordinates(:,pos_to_del) = [];
    end
    % Compute length of the path ------------------------------------------
    len = 0;
    for i=2:length(S.path.coordinates)
        len = len + norm(S.path.coordinates(:,i)-S.path.coordinates(:,i-1));
    end
    S.path.length = len;
    % Gen Slippage vector for each coordinate of the path -----------------
%     S.path.slippage = ones(S.system.np,length(S.path.coordinates));
end

function S = gen_obstacle(S,x,y,r)
    S.path.obstacles = [S.path.obstacles; [x, y, r]];
end

function S = gen_dynamic_and_solvers(S)
    % System's dimension __________________________________________________        
    S.system.nu      = 2;                                                              % input dimension: [w0,v0]
    S.system.ny      = S.config.N + 3;% + S.system.nu;      % output dimension: [beta_i,theta_0,xy_0]
    S.system.nv      = S.system.ny;  
    S.system.nq      = S.config.N + S.config.N+1 + 2*(S.config.N+1);
    % System's parameters -------------------------------------------------
    S.system.Lh1     = 0.342;
    S.system.L1      = 1.08;
    S.system.Lh2     = 0;
    S.system.L2      = 0.78;
    % ---------------------------------------------------------------------
    S.system.Lh3     = 0.15;%-0.342;
    S.system.L3      = 0.7;
    S.system.Lh4     = -0.15;%0.342;
    S.system.L4      = 0.7;
    S.system.Lh5     = 0;
    S.system.L5      = 0.7;
    S.system.Lh6     = -0.342;
    S.system.L6      = 0.7;
    S.system.Lh7     = 0.342;
    S.system.L7      = 0.7;
    S.system.Lh8     = 0;
    S.system.L8      = 0.7;
    S.system.Lh9     = 0;%0.3;%0.342;%0.048;
    S.system.L9      = 0.25;%0.78;%0.229;
    S.system.Lh10    = 0;%0.3;%0.048;
    S.system.L10     = 0.25;%0.78;%0.229;
    %
    S.system.Lhi     = [S.system.Lh1;S.system.Lh2;S.system.Lh3;S.system.Lh4;S.system.Lh5;S.system.Lh6;S.system.Lh7;S.system.Lh8;S.system.Lh9;S.system.Lh10]; 
    S.system.Li      = [S.system.L1;S.system.L2;S.system.L3;S.system.L4;S.system.L5;S.system.L6;S.system.L7;S.system.L8;S.system.L9;S.system.L10]; 
    %
    S.system.r       = 0.15;           % wheel radius
    S.system.b       = 0.55;%0.229;     % vehicles's width    
    %
    S.system.XYtracAxe           = [-S.system.b/2 S.system.b/2; 0 0];
    S.system.XYtracWheelLeft     = [-S.system.b/2 -S.system.b/2; S.system.r/2 -S.system.r/2];
    S.system.XYtracWheelRight    = [S.system.b/2 S.system.b/2; -S.system.r/2 S.system.r/2];
    S.system.XYtracBody          = [-S.system.b/2*0.75 0 S.system.b/2 -S.system.b/2*0.75; -S.system.b/2*0.75 S.system.b/2*2.75 -S.system.b/2*0.75 -S.system.b/2*0.75];
    %
    S.system.XYtrailerAxe        = [-S.system.b/2 S.system.b/2; 0 0];
    S.system.XYtrailerWheelLeft  = [-S.system.b/2 -S.system.b/2; S.system.r/2 -S.system.r/2];
    S.system.XYtrailerWheelRight = [S.system.b/2 S.system.b/2; -S.system.r/2 S.system.r/2];
    S.system.XYtrailerLoad       = [];
    S.system.XYtrailerLongAxe    = [];
    for i=1:S.config.N
        S.system.XYtrailerLongAxe = [S.system.XYtrailerLongAxe; [0 0; 0 S.system.Li(i)]];
        S.system.XYtrailerLoad = [S.system.XYtrailerLoad; [-S.system.b/2*0.9 -S.system.b/2*0.9 S.system.b/2*0.9 S.system.b/2*0.9 -S.system.b/2*0.9; -S.system.r/2*1.2 S.system.Li(i)*0.65 S.system.Li(i)*0.65 -S.system.r/2*1.2 -S.system.r/2*1.2]];
    end
    % Casadi variabes ----------------------------------------------------
    S.dynamic.q      = casadi.MX.sym('q',S.system.nq);
    S.dynamic.u      = casadi.MX.sym('u',S.system.nu);
%     S.dynamic.p      = casadi.MX.sym('beta_0',S.config.N);   % initial joint-angle values
    
    % Uncertain system's parameters ---------------------------------------
    if S.config.model.uncty
        S.system.Lh1_u   = S.system.Lh1 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.L1_u    = S.system.L1 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.Lh2_u   = S.system.Lh2 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.L2_u    = S.system.L2 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.Lh3_u   = S.system.Lh3 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.L3_u    = S.system.L3 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.Lh4_u   = S.system.Lh4 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.L4_u    = S.system.L4 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.Lh5_u   = S.system.Lh5 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.L5_u    = S.system.L5 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.Lh6_u   = S.system.Lh6 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.L6_u    = S.system.L6 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.Lh7_u   = S.system.Lh7 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.L7_u    = S.system.L7 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.Lh8_u   = S.system.Lh8 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.L8_u    = S.system.L8 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.Lh9_u   = S.system.Lh9 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.L9_u    = S.system.L9 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.Lh10_u  = S.system.Lh10 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
        S.system.L10_u   = S.system.L10 * (1 + (S.config.model.dev/100) * ((rand*2)-1));
    else
        S.system.Lh1_u   = S.system.Lh1;
        S.system.L1_u    = S.system.L1;
        S.system.Lh2_u   = S.system.Lh2;
        S.system.L2_u    = S.system.L2;
        S.system.Lh3_u   = S.system.Lh3;
        S.system.L3_u    = S.system.L3;
        S.system.Lh4_u   = S.system.Lh4;
        S.system.L4_u    = S.system.L4;
        S.system.Lh5_u   = S.system.Lh5;
        S.system.L5_u    = S.system.L5;
        S.system.Lh6_u   = S.system.Lh6;
        S.system.L6_u    = S.system.L6;
        S.system.Lh7_u   = S.system.Lh7;
        S.system.L7_u    = S.system.L7;
        S.system.Lh8_u   = S.system.Lh8;
        S.system.L8_u    = S.system.L8;
        S.system.Lh9_u   = S.system.Lh9;
        S.system.L9_u    = S.system.L9;
        S.system.Lh10_u  = S.system.Lh10;
        S.system.L10_u   = S.system.L10;
    end   
    % System's variables  -------------------------------------------------
    % Since ACADO does not allows to incorporate directly rate of change
    % constraints, I've add the controls w0 and v0 as additional states and
    % define as controls the variaton of w0 and v0 and add constraints on
    % these new inputs. Proper differential equations were added too.
    
    Control delta_w0 delta_v0;
    %
    switch(S.config.N)
        case 1            
            % CASADI FORMULATION ******************************************
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1, sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)), cos(S.dynamic.q(1))];
            S.system.J       = [S.system.J1];
            %
            S.dynamic.f_rhs  = [ [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                 %
                                 [1 0] * S.dynamic.u;...
                                 [1 0] * S.system.J1 * S.dynamic.u;...
                                 %
                                 [0 cos(S.dynamic.q(2))] * S.dynamic.u;...
                                 [0 sin(S.dynamic.q(2))] * S.dynamic.u;...
                                 [0 cos(S.dynamic.q(3))] * S.system.J1 * S.dynamic.u;...
                                 [0 sin(S.dynamic.q(3))] * S.system.J1 * S.dynamic.u  ];
            %
            S.system.J1_u       = [-S.system.Lh1_u * cos(S.dynamic.q(1)) / S.system.L1_u, sin(S.dynamic.q(1))/S.system.L1_u; S.system.Lh1_u*sin(S.dynamic.q(1)), cos(S.dynamic.q(1))];
            S.system.J_u        = [S.system.J1_u];
            S.dynamic.f_rhs_u   = [ [1 0]*(eye(2)-S.system.J1_u) * S.dynamic.u;...
                                 %
                                 [1 0] * S.dynamic.u;...
                                 [1 0] * S.system.J1_u * S.dynamic.u;...
                                 %
                                 [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                 [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                 [0 cos(S.dynamic.q(2*S.config.N+1))] * S.system.J1_u * S.dynamic.u;...
                                 [0 sin(S.dynamic.q(2*S.config.N+1))] * S.system.J1_u * S.dynamic.u];
                                 %
%                                  S.dynamic.u(1);...
%                                  S.dynamic.u(2) ];
            S.system.Lh      = [S.system.Lh1_u,S.system.L1_u];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 2
            % CASADI FORMULATION ******************************************
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2, sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)), cos(S.dynamic.q(2))];            
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1, sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)), cos(S.dynamic.q(1))];            
            S.system.J       = [S.system.J2, S.system.J1];

            S.dynamic.f_rhs  = [    [1 0] * (eye(2)-S.system.J1) * S.dynamic.u;
                                    [1 0] * (eye(2)-S.system.J2) * S.system.J1 * S.dynamic.u;
                                    [1 0] * S.dynamic.u;
                                    [1 0] * S.system.J1 * S.dynamic.u;
                                    [1 0] * S.system.J2 * S.system.J1 * S.dynamic.u;
                                    [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;
                                    [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;
                                    [0 cos(S.dynamic.q(S.config.N+2))] * S.system.J1 * S.dynamic.u;
                                    [0 sin(S.dynamic.q(S.config.N+2))] * S.system.J1 * S.dynamic.u;
                                    [0 cos(S.dynamic.q(S.config.N+3))] * S.system.J2 * S.system.J1 * S.dynamic.u;
                                    [0 sin(S.dynamic.q(S.config.N+3))] * S.system.J2 * S.system.J1 * S.dynamic.u];
            %
            S.system.J2_u    = [-S.system.Lh2_u * cos(S.dynamic.q(2)) / S.system.L2_u, sin(S.dynamic.q(2))/S.system.L2_u; S.system.Lh2_u*sin(S.dynamic.q(2)), cos(S.dynamic.q(2))];            
            S.system.J1_u    = [-S.system.Lh1_u * cos(S.dynamic.q(1)) / S.system.L1_u, sin(S.dynamic.q(1))/S.system.L1_u; S.system.Lh1_u*sin(S.dynamic.q(1)), cos(S.dynamic.q(1))];            
            S.system.J_u     = [S.system.J2, S.system.J1];

            S.dynamic.f_rhs_u  = [  [1 0] * (eye(2)-S.system.J1_u) * S.dynamic.u;
                                    [1 0] * (eye(2)-S.system.J2_u) * S.system.J1_u * S.dynamic.u;
                                    [1 0] * S.dynamic.u;
                                    [1 0] * S.system.J1_u * S.dynamic.u;
                                    [1 0] * S.system.J2_u * S.system.J1_u * S.dynamic.u;
                                    [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;
                                    [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;
                                    [0 cos(S.dynamic.q(S.config.N+2))] * S.system.J1_u * S.dynamic.u;
                                    [0 sin(S.dynamic.q(S.config.N+2))] * S.system.J1_u * S.dynamic.u;
                                    [0 cos(S.dynamic.q(S.config.N+3))] * S.system.J2_u * S.system.J1_u * S.dynamic.u;
                                    [0 sin(S.dynamic.q(S.config.N+3))] * S.system.J2_u * S.system.J1_u * S.dynamic.u;];
            %
            S.system.Lh      = [    S.system.Lh2_u,S.system.L2_u;...
                                    S.system.Lh1_u,S.system.L1_u];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 3
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.N+2))]*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+2))]*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+3))]*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+3))]*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(2*S.config.N+1))]*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.N+1))]*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u];
            %
            S.system.J3_u      = [-S.system.Lh3_u * cos(S.dynamic.q(3)) / S.system.L3_u sin(S.dynamic.q(3))/S.system.L3_u; S.system.Lh3_u*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2_u      = [-S.system.Lh2_u * cos(S.dynamic.q(2)) / S.system.L2_u sin(S.dynamic.q(2))/S.system.L2_u; S.system.Lh2_u*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1_u      = [-S.system.Lh1_u * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1_u*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J_u       = [S.system.J3_u, S.system.J2_u, S.system.J1_u];
            %
            S.dynamic.f_rhs_u  = [   [1 0]*(eye(2)-S.system.J1_u) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2_u)*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3_u)*S.system.J2_u * S.system.J1_u * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.N+2))]*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+2))]*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+3))]*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+3))]*S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(2*S.config.N+1))]*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.N+1))]*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u];
            %
            S.system.Lh      = [    S.system.Lh3_u,S.system.L3_u;...
                                    S.system.Lh2_u,S.system.L2_u;...
                                    S.system.Lh1_u,S.system.L1_u];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 4
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J4, S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.N+2))] * S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+2))] * S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+3))] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+3))] * S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+4))] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+4))] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(2*S.config.N+1))] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.N+1))] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u];
            S.system.J4_u      = [-S.system.Lh4_u * cos(S.dynamic.q(4)) / S.system.L4_u sin(S.dynamic.q(4))/S.system.L4_u; S.system.Lh4_u*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3_u      = [-S.system.Lh3_u * cos(S.dynamic.q(3)) / S.system.L3_u sin(S.dynamic.q(3))/S.system.L3_u; S.system.Lh3_u*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2_u      = [-S.system.Lh2_u * cos(S.dynamic.q(2)) / S.system.L2_u sin(S.dynamic.q(2))/S.system.L2_u; S.system.Lh2_u*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1_u      = [-S.system.Lh1_u * cos(S.dynamic.q(1)) / S.system.L1_u sin(S.dynamic.q(1))/S.system.L1_u; S.system.Lh1_u*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J_u       = [S.system.J4_u, S.system.J3_u, S.system.J2_u, S.system.J1_u];
            S.dynamic.f_rhs_u  = [ [1 0]*(eye(2)-S.system.J1_u) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2_u)*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3_u)*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4_u)*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.N+2))] * S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+2))] * S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+3))] * S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+3))] * S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+4))] * S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+4))] * S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(2*S.config.N+1))] * S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.N+1))] * S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u];
            S.system.Lh      = [    S.system.Lh4_u,S.system.L4_u;...
                                    S.system.Lh3_u,S.system.L3_u;...
                                    S.system.Lh2_u,S.system.L2_u;...
                                    S.system.Lh1_u,S.system.L1_u];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 5
            S.system.J5      = [-S.system.Lh5 * cos(S.dynamic.q(5)) / S.system.L5 sin(S.dynamic.q(5))/S.system.L5; S.system.Lh5*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J5, S.system.J4, S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.N+2))] *S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+2))] *S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+3))] *S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+3))] *S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+4))] *S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+4))] *S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+5))] *S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+5))] *S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(2*S.config.N+1))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.N+1))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u];
            S.system.J5_u      = [-S.system.Lh5_u * cos(S.dynamic.q(5)) / S.system.L5_u sin(S.dynamic.q(5))/S.system.L5_u; S.system.Lh5_u*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4_u      = [-S.system.Lh4_u * cos(S.dynamic.q(4)) / S.system.L4_u sin(S.dynamic.q(4))/S.system.L4_u; S.system.Lh4_u*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3_u     = [-S.system.Lh3_u * cos(S.dynamic.q(3)) / S.system.L3_u sin(S.dynamic.q(3))/S.system.L3_u; S.system.Lh3_u*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2_u      = [-S.system.Lh2_u * cos(S.dynamic.q(2)) / S.system.L2_u sin(S.dynamic.q(2))/S.system.L2_u; S.system.Lh2_u*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1_u      = [-S.system.Lh1_u * cos(S.dynamic.q(1)) / S.system.L1_u sin(S.dynamic.q(1))/S.system.L1_u; S.system.Lh1_u*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J_u       = [S.system.J5_u, S.system.J4_u, S.system.J3_u, S.system.J2_u, S.system.J1_u];
            %
            S.dynamic.f_rhs_u  = [ [1 0]*(eye(2)-S.system.J1_u) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2_u)*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3_u)*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4_u)*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5_u)*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.N+2))] *S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+2))] *S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+3))] *S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+3))] *S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+4))] *S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+4))] *S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+5))] *S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+5))] *S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(2*S.config.N+1))]*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.N+1))]*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u];

            S.system.Lh      = [    S.system.Lh5_u,S.system.L5_u;...
                                    S.system.Lh4_u,S.system.L4_u;...
                                    S.system.Lh3_u,S.system.L3_u;...
                                    S.system.Lh2_u,S.system.L2_u;...
                                    S.system.Lh1_u,S.system.L1_u];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 6
            S.system.J6      = [-S.system.Lh6 * cos(S.dynamic.q(6)) / S.system.L6 sin(S.dynamic.q(6))/S.system.L6; S.system.Lh6*sin(S.dynamic.q(6)) cos(S.dynamic.q(6))];
            S.system.J5      = [-S.system.Lh5 * cos(S.dynamic.q(5)) / S.system.L5 sin(S.dynamic.q(5))/S.system.L5; S.system.Lh5*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J6, S.system.J5, S.system.J4, S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.N+2))]*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+2))]*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+3))]*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+3))]*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+4))]*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+4))]*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(2*S.config.N+1))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.N+1))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u];

            S.system.J6_u      = [-S.system.Lh6_u * cos(S.dynamic.q(6)) / S.system.L6_u sin(S.dynamic.q(6))/S.system.L6_u; S.system.Lh6_u*sin(S.dynamic.q(6)) cos(S.dynamic.q(6))];
            S.system.J5_u      = [-S.system.Lh5_u * cos(S.dynamic.q(5)) / S.system.L5_u sin(S.dynamic.q(5))/S.system.L5_u; S.system.Lh5_u*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4_u      = [-S.system.Lh4_u * cos(S.dynamic.q(4)) / S.system.L4_u sin(S.dynamic.q(4))/S.system.L4_u; S.system.Lh4_u*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3_u      = [-S.system.Lh3_u * cos(S.dynamic.q(3)) / S.system.L3_u sin(S.dynamic.q(3))/S.system.L3_u; S.system.Lh3_u*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2_u      = [-S.system.Lh2_u * cos(S.dynamic.q(2)) / S.system.L2_u sin(S.dynamic.q(2))/S.system.L2_u; S.system.Lh2_u*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1_u      = [-S.system.Lh1_u * cos(S.dynamic.q(1)) / S.system.L1_u sin(S.dynamic.q(1))/S.system.L1_u; S.system.Lh1_u*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J_u       = [S.system.J6_u, S.system.J5_u, S.system.J4_u, S.system.J3_u, S.system.J2_u, S.system.J1_u];
            %
            S.dynamic.f_rhs_u  = [ [1 0]*(eye(2)-S.system.J1_u) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2_u)*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3_u)*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4_u)*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5_u)*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J6_u)*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J6_u*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.N+2))]*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+2))]*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+3))]*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+3))]*S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+4))]*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+4))]*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+5))]*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+5))]*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+6))]*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+6))]*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u                                   
                                   [0 cos(S.dynamic.q(2*S.config.N+1))]*S.system.J6_u*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.N+1))]*S.system.J6_u*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u];
            %
            S.system.Lh      = [    S.system.Lh6_u,S.system.L6_u;...
                                    S.system.Lh5_u,S.system.L5_u;...
                                    S.system.Lh4_u,S.system.L4_u;...
                                    S.system.Lh3_u,S.system.L3_u;...
                                    S.system.Lh2_u,S.system.L2_u;...
                                    S.system.Lh1_u,S.system.L1_u];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 7
            S.system.J7      = [-S.system.Lh7 * cos(S.dynamic.q(7)) / S.system.L7 sin(S.dynamic.q(7))/S.system.L7; S.system.Lh7*sin(S.dynamic.q(7)) cos(S.dynamic.q(7))];
            S.system.J6      = [-S.system.Lh6 * cos(S.dynamic.q(6)) / S.system.L6 sin(S.dynamic.q(6))/S.system.L6; S.system.Lh6*sin(S.dynamic.q(6)) cos(S.dynamic.q(6))];
            S.system.J5      = [-S.system.Lh5 * cos(S.dynamic.q(5)) / S.system.L5 sin(S.dynamic.q(5))/S.system.L5; S.system.Lh5*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J7, S.system.J6, S.system.J5, S.system.J4, S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J7)*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.N+2))]*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+2))]*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+3))]*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+3))]*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+4))]*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+4))]*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+5))]*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+6))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+7))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+7))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+8))]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+8))]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u];

            S.system.J7_u      = [-S.system.Lh7_u * cos(S.dynamic.q(7)) / S.system.L7_u sin(S.dynamic.q(7))/S.system.L7_u; S.system.Lh7_u*sin(S.dynamic.q(7)) cos(S.dynamic.q(7))];
            S.system.J6_u      = [-S.system.Lh6_u * cos(S.dynamic.q(6)) / S.system.L6_u sin(S.dynamic.q(6))/S.system.L6_u; S.system.Lh6_u*sin(S.dynamic.q(6)) cos(S.dynamic.q(6))];
            S.system.J5_u      = [-S.system.Lh5_u * cos(S.dynamic.q(5)) / S.system.L5_u sin(S.dynamic.q(5))/S.system.L5_u; S.system.Lh5_u*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4_u      = [-S.system.Lh4_u * cos(S.dynamic.q(4)) / S.system.L4_u sin(S.dynamic.q(4))/S.system.L4_u; S.system.Lh4_u*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3_u      = [-S.system.Lh3_u * cos(S.dynamic.q(3)) / S.system.L3_u sin(S.dynamic.q(3))/S.system.L3_u; S.system.Lh3_u*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2_u      = [-S.system.Lh2_u * cos(S.dynamic.q(2)) / S.system.L2_u sin(S.dynamic.q(2))/S.system.L2_u; S.system.Lh2_u*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1_u      = [-S.system.Lh1_u * cos(S.dynamic.q(1)) / S.system.L1_u sin(S.dynamic.q(1))/S.system.L1_u; S.system.Lh1_u*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J_u       = [S.system.J7_u, S.system.J6_u, S.system.J5_u, S.system.J4_u, S.system.J3_u, S.system.J2_u, S.system.J1_u];
            %
            S.dynamic.f_rhs_u  = [ [1 0]*(eye(2)-S.system.J1_u) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2_u)*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3_u)*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4_u)*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5_u)*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J6_u)*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J7_u)*S.system.J6_u*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J6_u*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [1 0] * S.system.J7_u*S.system.J6_u*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+1))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(S.config.N+2))]*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+2))]*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+3))]*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+3))]*S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+4))]*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+4))]*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+5))]*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+5))]*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+6))]*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+6))]*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u                                   
                                   [0 cos(S.dynamic.q(S.config.N+7))]*S.system.J6_u*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+7))]*S.system.J6_u*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u
                                   [0 cos(S.dynamic.q(S.config.N+8))]*S.system.J7_u*S.system.J6_u*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(S.config.N+8))]*S.system.J7_u*S.system.J6_u*S.system.J5_u*S.system.J4_u*S.system.J3_u*S.system.J2_u*S.system.J1_u * S.dynamic.u];
            %
            S.system.Lh      = [    S.system.Lh7_u,S.system.L7_u;...
                                    S.system.Lh6_u,S.system.L6_u;...
                                    S.system.Lh5_u,S.system.L5_u;...
                                    S.system.Lh4_u,S.system.L4_u;...
                                    S.system.Lh3_u,S.system.L3_u;...
                                    S.system.Lh2_u,S.system.L2_u;...
                                    S.system.Lh1_u,S.system.L1_u];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));            
        case 8
            S.system.J8      = [-S.system.Lh8 * cos(S.dynamic.q(8)) / S.system.L8 sin(S.dynamic.q(8))/S.system.L8; S.system.Lh8*sin(S.dynamic.q(8)) cos(S.dynamic.q(8))];
            S.system.J7      = [-S.system.Lh7 * cos(S.dynamic.q(7)) / S.system.L7 sin(S.dynamic.q(7))/S.system.L7; S.system.Lh7*sin(S.dynamic.q(7)) cos(S.dynamic.q(7))];
            S.system.J6      = [-S.system.Lh6 * cos(S.dynamic.q(6)) / S.system.L6 sin(S.dynamic.q(6))/S.system.L6; S.system.Lh6*sin(S.dynamic.q(6)) cos(S.dynamic.q(6))];
            S.system.J5      = [-S.system.Lh5 * cos(S.dynamic.q(5)) / S.system.L5 sin(S.dynamic.q(5))/S.system.L5; S.system.Lh5*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J8, S.system.J7, S.system.J6, S.system.J5, S.system.J4, S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J7)*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J8)*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(9))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(9))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(2*S.config.N+1))]*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.N+1))]*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u];%;...
                                   %
%                                    S.dynamic.u(1) / 1;
%                                    S.dynamic.u(2) / 1 ];
            %
            S.system.Lh      = [    S.system.Lh8_u,S.system.L8_u;...
                                    S.system.Lh7_u,S.system.L7_u;...
                                    S.system.Lh6_u,S.system.L6_u;...
                                    S.system.Lh5_u,S.system.L5_u;...
                                    S.system.Lh4_u,S.system.L4_u;...
                                    S.system.Lh3_u,S.system.L3_u;...
                                    S.system.Lh2_u,S.system.L2_u;...
                                    S.system.Lh1_u,S.system.L1_u];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));            
        case 9
            S.system.J9      = [-S.system.Lh9 * cos(S.dynamic.q(9)) / S.system.L9 sin(S.dynamic.q(9))/S.system.L9; S.system.Lh9*sin(S.dynamic.q(9)) cos(S.dynamic.q(9))];
            S.system.J8      = [-S.system.Lh8 * cos(S.dynamic.q(8)) / S.system.L8 sin(S.dynamic.q(8))/S.system.L8; S.system.Lh8*sin(S.dynamic.q(8)) cos(S.dynamic.q(8))];
            S.system.J7      = [-S.system.Lh7 * cos(S.dynamic.q(7)) / S.system.L7 sin(S.dynamic.q(7))/S.system.L7; S.system.Lh7*sin(S.dynamic.q(7)) cos(S.dynamic.q(7))];
            S.system.J6      = [-S.system.Lh6 * cos(S.dynamic.q(6)) / S.system.L6 sin(S.dynamic.q(6))/S.system.L6; S.system.Lh6*sin(S.dynamic.q(6)) cos(S.dynamic.q(6))];
            S.system.J5      = [-S.system.Lh5 * cos(S.dynamic.q(5)) / S.system.L5 sin(S.dynamic.q(5))/S.system.L5; S.system.Lh5*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J9, S.system.J8, S.system.J7, S.system.J6, S.system.J5, S.system.J4, S.system.J3, S.system.J2, S.system.J1];
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J7)*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J8)*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J9)*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(10))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(10))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(2*S.config.N+1))]*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.N+1))]*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u];%;...
                                   %
%                                    S.dynamic.u(1) / 1;
%                                    S.dynamic.u(2) / 1 ];
            %
            S.system.Lh      = [    S.system.Lh9_u,S.system.L9_u;...
                                    S.system.Lh8_u,S.system.L8_u;...
                                    S.system.Lh7_u,S.system.L7_u;...
                                    S.system.Lh6_u,S.system.L6_u;...
                                    S.system.Lh5_u,S.system.L5_u;...
                                    S.system.Lh4_u,S.system.L4_u;...
                                    S.system.Lh3_u,S.system.L3_u;...
                                    S.system.Lh2_u,S.system.L2_u;...
                                    S.system.Lh1_u,S.system.L1_u];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));            
        case 10
            S.system.J10     = [-S.system.Lh10 * cos(S.dynamic.q(10)) / S.system.L10 sin(S.dynamic.q(10))/S.system.L10; S.system.Lh10*sin(S.dynamic.q(10)) cos(S.dynamic.q(10))];
            S.system.J9      = [-S.system.Lh9 * cos(S.dynamic.q(9)) / S.system.L9 sin(S.dynamic.q(9))/S.system.L9; S.system.Lh9*sin(S.dynamic.q(9)) cos(S.dynamic.q(9))];
            S.system.J8      = [-S.system.Lh8 * cos(S.dynamic.q(8)) / S.system.L8 sin(S.dynamic.q(8))/S.system.L8; S.system.Lh8*sin(S.dynamic.q(8)) cos(S.dynamic.q(8))];
            S.system.J7      = [-S.system.Lh7 * cos(S.dynamic.q(7)) / S.system.L7 sin(S.dynamic.q(7))/S.system.L7; S.system.Lh7*sin(S.dynamic.q(7)) cos(S.dynamic.q(7))];
            S.system.J6      = [-S.system.Lh6 * cos(S.dynamic.q(6)) / S.system.L6 sin(S.dynamic.q(6))/S.system.L6; S.system.Lh6*sin(S.dynamic.q(6)) cos(S.dynamic.q(6))];
            S.system.J5      = [-S.system.Lh5 * cos(S.dynamic.q(5)) / S.system.L5 sin(S.dynamic.q(5))/S.system.L5; S.system.Lh5*sin(S.dynamic.q(5)) cos(S.dynamic.q(5))];
            S.system.J4      = [-S.system.Lh4 * cos(S.dynamic.q(4)) / S.system.L4 sin(S.dynamic.q(4))/S.system.L4; S.system.Lh4*sin(S.dynamic.q(4)) cos(S.dynamic.q(4))];
            S.system.J3      = [-S.system.Lh3 * cos(S.dynamic.q(3)) / S.system.L3 sin(S.dynamic.q(3))/S.system.L3; S.system.Lh3*sin(S.dynamic.q(3)) cos(S.dynamic.q(3))];
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2 sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)) cos(S.dynamic.q(2))];
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1 sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)) cos(S.dynamic.q(1))];
            %
            S.system.J       = [S.system.J10, S.system.J9, S.system.J8, S.system.J7, S.system.J6, S.system.J5, S.system.J4, S.system.J3, S.system.J2, S.system.J1];            
            %
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J7)*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J8)*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J9)*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0]*(eye(2)-S.system.J10)*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [1 0] * S.dynamic.u;...
                                   [1 0] * S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [1 0] * S.system.J10*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   %
                                   [0 cos(S.dynamic.q(11))] * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(11))] * S.dynamic.u;...
                                   [0 cos(S.dynamic.q(2*S.config.N+1))] * S.system.J10*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u;...
                                   [0 sin(S.dynamic.q(2*S.config.N+1))] * S.system.J10*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.u];%;...
                                   %
%                                    S.dynamic.u(1) / 1;
%                                    S.dynamic.u(2) / 1 ];
            %
            S.system.Lh      = [    S.system.Lh10_u,S.system.L10_u;...
                                    S.system.Lh9_u,S.system.L9_u;...
                                    S.system.Lh8_u,S.system.L8_u;...
                                    S.system.Lh7_u,S.system.L7_u;...
                                    S.system.Lh6_u,S.system.L6_u;...
                                    S.system.Lh5_u,S.system.L5_u;...
                                    S.system.Lh4_u,S.system.L4_u;...
                                    S.system.Lh3_u,S.system.L3_u;...
                                    S.system.Lh2_u,S.system.L2_u;...
                                    S.system.Lh1_u,S.system.L1_u];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));            
    end
    % Generate CASADI functions and integrators ***************************
    S.dynamic.f     = casadi.Function('f_rhs', {S.dynamic.q,S.dynamic.u}, {S.dynamic.f_rhs});    
    opts            = struct('main',true,'mex',true);
    S.dynamic.f.generate('f.c',opts)
    mex f.c -largeArrayDims;
    % RK4 -----------------------------------------------------------------
    k1              = S.dynamic.f(S.dynamic.q, S.dynamic.u);
    k2              = S.dynamic.f(S.dynamic.q + S.config.Ts / 2 * k1, S.dynamic.u);
    k3              = S.dynamic.f(S.dynamic.q + S.config.Ts / 2 * k2, S.dynamic.u);
    k4              = S.dynamic.f(S.dynamic.q + S.config.Ts * k3, S.dynamic.u);
    x_rk4           = S.dynamic.q + S.config.Ts / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    S.dynamic.FNt   = casadi.Function('FNt', {S.dynamic.q, S.dynamic.u}, {x_rk4});
    S.dynamic.FNt.generate('FNt.c',opts);
    mex FNt.c -largeArrayDims;
    % Output of the system ------------------------------------------------
    S.dynamic.h_rhs  = S.dynamic.q(S.config.outputs);
    S.dynamic.h      = casadi.Function('h', {S.dynamic.q}, {S.dynamic.h_rhs});
end

function S = init_ekfCasadi(S)

% S.Mtxs.Q    = diag([ones(1,S.config.N),1,0.1.*ones(1,S.config.N),1,1,0.1,0.1]);%,0.2,0.2]);
% S.Mtxs.R    = diag([ones(1,S.config.N),2,4,4]);%,1,1]);
% S.Mtxs.P    = 1e6.*eye(S.system.nq);

S.Mtxs.Q    = diag([ones(1,S.config.N),1,0.1.*ones(1,S.config.N-1),1,1,1,1,1]);%,0.2,0.2]);
S.Mtxs.R    = diag([ones(1,S.config.N),0.05,1,1]);%,1,1]);
S.Mtxs.P    = 1e6.*eye(S.system.nq);

S.algorithms.ekfCasadi = ekfCasadi(S.Mtxs, S.system.nq, S.system.ny, S.system.nu, S.dynamic.q, S.dynamic.u,[], S.dynamic.f_rhs_u, S.dynamic.h_rhs, S.config.Ts);

set_x0bar(S.algorithms.ekfCasadi, S.init_condition.x0bar);
%
updateMeasurement(S.algorithms.ekfCasadi, full(S.dynamic.h(S.init_condition.x0bar)));

end

function S = init_mheCasadi(S)
    S.box_constraints               = struct;
    S.box_constraints.QluBounds     = [];
    S.box_constraints.WluBounds     = [];
    S.box_constraints.VluBounds     = [];
    S.box_constraints.ZluBounds     = [repmat(-0.01,S.config.N,1), repmat(0.01,S.config.N,1)];
    S.box_constraints.UluBounds     = [];
    %
    S.Mtxs                          = struct;
    
    if S.config.SIM == false % G2T
        S.Mtxs.Q    = diag([0.5.*ones(1,S.config.N),1,0.25.*ones(1,S.config.N),ones(1,2),0.25.*ones(1,2*S.config.N)]);%,0.2,0.2]);
        S.Mtxs.R    = diag([0.25.*ones(1,S.config.N),0.65,1,1]);%,1,1]);
        S.Mtxs.P    = 1e6.*eye(S.system.nq);
        S.Mtxs.Z    = 1e0.*eye(S.config.N);
        S.Mtxs.dU   = eye(S.system.nu);
    else
        S.Mtxs.Q    = diag([0.5.*ones(1,S.config.N),1,0.25.*ones(1,S.config.N),ones(1,2),0.25.*ones(1,2*S.config.N)]);%,0.2,0.2]);
        S.Mtxs.R    = diag([0.25.*ones(1,S.config.N),0.65,1,1]);%,1,1]);
        S.Mtxs.P    = 1e6.*eye(S.system.nq);
        S.Mtxs.Z    = 1e0.*eye(S.config.N);
        S.Mtxs.dU   = eye(S.system.nu);
    end
    % nlmheCasadiNt(Ne,x,u,Nt,f_rhs,h_rhs,Mtxs,nq,nu,ny,boxConst,q0bar,Ts,dimensions)
    S.algorithms.mheCasadi = nlmheCasadiNt(S.config.Ne,S.dynamic.q,S.dynamic.u,S.config.N,S.dynamic.f_rhs_u,S.dynamic.h_rhs,...
        S.Mtxs,S.system.nq,S.system.nu,S.system.ny,S.box_constraints,S.init_condition.x0bar,S.config.Ts,S.system.Lh);
    %
    setSigma(S.algorithms.mheCasadi, 1);
    setC(S.algorithms.mheCasadi, 1e7);
end

function S = fill_mhe(S)
    for imhe=1:S.config.Ne
        updateMeasurement(S.algorithms.mheCasadi, full(S.dynamic.h(S.init_condition.x0bar)));
        updateInput(S.algorithms.mheCasadi,zeros(S.system.nu,1));       
    end
    updateMeasurement(S.algorithms.mheCasadi, full(S.dynamic.h(S.init_condition.x0bar)));
end

function S = init_mpc(S)
    % PARAMETERS MPC
    S.box_constraints               = struct;    
    S.box_constraints.QluBounds     = [repmat([-70*pi/180 70*pi/180],S.config.N,1); repmat([-inf inf],S.config.N+1,1); repmat([-inf inf],2*(S.config.N+1),1)];%; [-1 1; -5 5]];
    S.box_constraints.QNluBounds    = [repmat([-70*pi/180 70*pi/180],S.config.N,1); repmat([-inf inf],S.config.N+1,1); repmat([-inf inf],2*(S.config.N+1),1)];%; [-1 1; -5 5]];
    S.box_constraints.UluBounds     = [-2   2;...
                                       -0.5 0.5];
    S.box_constraints.dUluBounds    = [-1.5 1.5;...
                                       -0.5 0.5];
    %
    S.Mtxs                          = struct;
    
    if S.config.SIM == false % G2T
        if strcmp(S.path.path,'infinity')
            % Nc = 10 works fine
            S.Mtxs.indxMposesqToy = [2*S.config.N+2:S.system.nq];                                         % indices that point to those states that are being controlled or steered
            S.Mtxs.Q        = diag([0.5.*ones(1,S.config.N), 0.5.*ones(1,S.config.N+1), ones(1,2*(S.config.N+1)) ]); % output tracking
            S.Mtxs.R        = diag([0.75 1]);
            S.Mtxs.QrefOpt  = 0.75.*diag([0.5.*ones(1,S.config.N), 0.5.*ones(1,S.config.N+1), ones(1,2*(S.config.N+1))]);
%             S.Mtxs.QrefOpt  = diag([zeros(1,S.config.N), zeros(1,S.config.N+1), ones(1,2*(S.config.N+1))]);
            S.Mtxs.UrefOpt  = eye(S.system.nu);
            S.Mtxs.wq_q0    = 100.*blkdiag(eye(S.config.N), eye(S.config.N+1), eye(2*(S.config.N+1)));
            S.Mtxs.w_sigma  = 20;
            S.Mtxs.G        = 5.*eye(2);
            S.Mtxs.s_prm    = eye(S.config.N+1);
        elseif strcmp(S.path.path,'flat_infinity')
%             % Nc = 10 works fine
            S.Mtxs.indxMposesqToy = S.system.nq-1:S.system.nq;%[2*S.config.N+2:S.system.nq];                                         % indices that point to those states that are being controlled or steered
    %         S.Mtxs.qToqSteered = blkdiag([zeros(S.config.N), zeros(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
%             S.Mtxs.qToqSteered = blkdiag([eye(S.config.N), eye(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
            S.Mtxs.Q        = diag([0.*ones(1,S.config.N), 0*ones(1,S.config.N+1), ones(1,2), ones(1,2), ones(1,2) ]); % output tracking
            S.Mtxs.R        = diag([0.01 0.01]);
            S.Mtxs.QrefOpt  = diag([0.03.*eye(1,S.config.N), 0.1.*eye(1,S.config.N+1), 0.75.*ones(1,2), 0.75.*ones(1,2), 0.75.*ones(1,2)]);
            S.Mtxs.UrefOpt  = eye(S.system.nu);
            S.Mtxs.wq_q0    = 20.*blkdiag(eye(S.config.N), eye(S.config.N+1), eye(2*(S.config.N+1)));
            S.Mtxs.w_sigma  = 402;
            S.Mtxs.G        = 5.*eye(2);
            S.Mtxs.s_prm    = 300.*eye(S.config.N+1);
            %
            % S.Mtxs.indxMposesqToy = S.config.N+1:S.system.nq;%[2*S.config.N+2:S.system.nq];                                         % indices that point to those states that are being controlled or steered
            % S.Mtxs.qToqSteered = blkdiag([eye(S.config.N), eye(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
            % S.Mtxs.Q        = diag([0.5.*ones(1,S.config.N), 0.5.*ones(1,S.config.N+1), ones(1,2), ones(1,2), 5.*ones(1,2)]); % output tracking
            % S.Mtxs.R        = diag([0.1 0.01]);
            % S.Mtxs.QrefOpt  = diag([0.1.*ones(1,S.config.N), 0.1.*ones(1,S.config.N+1), 1.*ones(1,2), 1.*ones(1,2), 2.*ones(1,2)]);% 0.7.*eye(S.system.nq);
            % S.Mtxs.UrefOpt  = eye(S.system.nu);
            % S.Mtxs.wq_q0    = 20.*blkdiag(eye(S.config.N), eye(S.config.N+1), eye(2*(S.config.N+1)));%1e6.*eye(S.system.nq);
            % S.Mtxs.w_sigma  = 100;
            % S.Mtxs.G        = 5.*eye(2);
            % S.Mtxs.s_prm    = 500.*eye(S.config.N+1);
        elseif strcmp(S.path.path,'segment') || strcmp(S.path.path,'invert_segment')
            % Nc = 10 works fine
            S.Mtxs.indxMposesqToy = 1:S.system.nq;%[2*S.config.N+2:S.system.nq];                                         % indices that point to those states that are being controlled or steered
    %         S.Mtxs.qToqSteered = blkdiag([zeros(S.config.N), zeros(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
%             S.Mtxs.qToqSteered = blkdiag([eye(S.config.N), eye(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
            S.Mtxs.Q        = diag([0.5.*ones(1,S.config.N), 0.5.*ones(1,S.config.N+1), 5.*ones(1,2*(S.config.N+1)) ]); % output tracking
            S.Mtxs.R        = diag([0.1 0.05]);
            S.Mtxs.QrefOpt  = diag([0.25.*eye(1,S.config.N), 0.5.*eye(1,S.config.N+1), ones(1,2*(S.config.N+1))]);
            S.Mtxs.UrefOpt  = eye(S.system.nu);
            S.Mtxs.wq_q0    = 20.*blkdiag(eye(S.config.N), eye(S.config.N+1), eye(2*(S.config.N+1)));
            S.Mtxs.w_sigma  = 152;
            S.Mtxs.G        = 5.*eye(2);
            S.Mtxs.s_prm    = 20.*eye(S.config.N+1);
        elseif strcmp(S.path.path,'rectangular')
            % Nc = 6 works fine
            S.Mtxs.indxMposesqToy = S.system.nq-1:S.system.nq;%[2*S.config.N+2:S.system.nq];                                         % indices that point to those states that are being controlled or steered
            % S.Mtxs.qToqSteered = blkdiag([eye(S.config.N), eye(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
            S.Mtxs.Q        = diag([0.5.*ones(1,S.config.N), 0.5.*ones(1,S.config.N+1), ones(1,2), ones(1,2), 8.*ones(1,2)]); % output tracking
            S.Mtxs.R        = diag([0.075 0.01]);
            S.Mtxs.QrefOpt  = diag([0.1.*ones(1,S.config.N), 0.1.*ones(1,S.config.N+1), 1.*ones(1,2), 1.*ones(1,2), 2.*ones(1,2)]);% 0.7.*eye(S.system.nq);
            S.Mtxs.UrefOpt  = eye(S.system.nu);
            S.Mtxs.wq_q0    = 20.*blkdiag(eye(S.config.N), eye(S.config.N+1), eye(2*(S.config.N+1)));%1e6.*eye(S.system.nq);
            S.Mtxs.w_sigma  = 100;
            S.Mtxs.G        = 5.*eye(2);
            S.Mtxs.s_prm    = 500.*eye(S.config.N+1);
        end
    else
        if strcmp(S.path.path,'infinity')
            % Nc = 10 works fine
            S.Mtxs.indxMposesqToy = [2*S.config.N+2:S.system.nq];                                         % indices that point to those states that are being controlled or steered
    %         S.Mtxs.qToqSteered = blkdiag([zeros(S.config.N), zeros(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
%             S.Mtxs.qToqSteered = blkdiag([eye(S.config.N), eye(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
            S.Mtxs.Q        = diag([0.5.*ones(1,S.config.N), 0.5.*ones(1,S.config.N+1), ones(1,2*(S.config.N+1)) ]); % output tracking
            S.Mtxs.R        = diag([1 1]);
            S.Mtxs.QrefOpt  = diag([0.5.*ones(1,S.config.N), 0.5.*ones(1,S.config.N+1), ones(1,2*(S.config.N+1))]);% 0.7.*eye(S.system.nq);
            S.Mtxs.UrefOpt  = eye(S.system.nu);
            S.Mtxs.wq_q0    = 20.*blkdiag(eye(S.config.N), eye(S.config.N+1), eye(2*(S.config.N+1)));%1e6.*eye(S.system.nq);
            S.Mtxs.w_sigma  = 1000;
            S.Mtxs.G        = 5.*eye(2);
            S.Mtxs.s_prm    = eye(S.config.N+1);
        elseif strcmp(S.path.path,'flat_infinity')
            % Nc = 10 works fine
            S.Mtxs.indxMposesqToy = [2*S.config.N+2:S.system.nq];                                         % indices that point to those states that are being controlled or steered
    %         S.Mtxs.qToqSteered = blkdiag([zeros(S.config.N), zeros(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
%             S.Mtxs.qToqSteered = blkdiag([eye(S.config.N), eye(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
            S.Mtxs.Q        = diag([0.5.*ones(1,S.config.N), 0.5.*ones(1,S.config.N+1), ones(1,2*(S.config.N+1)) ]); % output tracking
            S.Mtxs.R        = diag([1 1]);
            S.Mtxs.QrefOpt  = diag([zeros(1,S.config.N), zeros(1,S.config.N+1), ones(1,2*(S.config.N+1))]);% 0.7.*eye(S.system.nq);
            S.Mtxs.UrefOpt  = eye(S.system.nu);
            S.Mtxs.wq_q0    = 20.*blkdiag(eye(S.config.N), eye(S.config.N+1), eye(2*(S.config.N+1)));%1e6.*eye(S.system.nq);
            S.Mtxs.w_sigma  = 10;
            S.Mtxs.G        = 5.*eye(2);
            S.Mtxs.s_prm    = eye(S.config.N+1);
        elseif strcmp(S.path.path,'segment') || strcmp(S.path.path,'invert_segment')
            % Nc = 10 works fine
            S.Mtxs.indxMposesqToy = [2*S.config.N+2:S.system.nq];                                         % indices that point to those states that are being controlled or steered
    %         S.Mtxs.qToqSteered = blkdiag([zeros(S.config.N), zeros(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
%             S.Mtxs.qToqSteered = blkdiag([eye(S.config.N), eye(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
            S.Mtxs.Q        = diag([0.5.*ones(1,S.config.N), 0.5.*ones(1,S.config.N+1), 10.*ones(1,2*(S.config.N+1)) ]); % output tracking
            S.Mtxs.R        = diag([0.1 0.1]);
            S.Mtxs.QrefOpt  = diag([zeros(1,S.config.N), zeros(1,S.config.N+1), 0.5.*ones(1,2*(S.config.N+1))]);
            S.Mtxs.UrefOpt  = eye(S.system.nu);
            S.Mtxs.wq_q0    = 20.*blkdiag(eye(S.config.N), eye(S.config.N+1), eye(2*(S.config.N+1)));
            S.Mtxs.w_sigma  = 12;
            S.Mtxs.G        = 5.*eye(2);
            S.Mtxs.s_prm    = 20.*eye(S.config.N+1);
        elseif strcmp(S.path.path,'rectangular')
            % Nc = 5 works fine
            S.Mtxs.indxMposesqToy = [2*S.config.N+2:S.system.nq];                                         % indices that point to those states that are being controlled or steered
            S.Mtxs.qToqSteered = blkdiag([eye(S.config.N), eye(S.config.N, S.config.N+1)], eye(2*(S.config.N+1)));
            S.Mtxs.Q        = diag([0.5.*ones(1,S.config.N), 0.5.*ones(1,S.config.N+1), ones(1,2), ones(1,2), 2.*ones(1,2)]); % output tracking
            S.Mtxs.R        = diag([1 1]);
            S.Mtxs.QrefOpt  = diag([0.1.*ones(1,S.config.N), ones(1,S.config.N+1), 0.5.*ones(1,2), 0.5.*ones(1,2), ones(1,2)]);% 0.7.*eye(S.system.nq);
            S.Mtxs.UrefOpt  = eye(S.system.nu);
            S.Mtxs.wq_q0    = 20.*blkdiag(eye(S.config.N), eye(S.config.N+1), eye(2*(S.config.N+1)));%1e6.*eye(S.system.nq);
            S.Mtxs.w_sigma  = 400;
            S.Mtxs.G        = 5.*eye(2);
            S.Mtxs.s_prm    = eye(S.config.N+1);
        end
    end 
    S.Mtxs.QN = S.Mtxs.Q;
    %
    S.algorithms.mpcCasadi = nlmpcCasadiNt(S.config.Nc,S.dynamic.q,S.dynamic.u,S.config.N,S.dynamic.f_rhs_u,S.Mtxs,...
        S.system.nq,S.system.nu,S.box_constraints,S.config.Ts,S.path.s,S.path.st,S.path.r,S.system.Lh,S.config.obs_strategy,S.config.Ne);
    %
    setVref(S.algorithms.mpcCasadi,S.config.vNr);
    %
    S.algorithms.mpc.last_tgt  = 0;
    S.algorithms.Controls      = zeros(S.system.nu,1);
    %
    S.algorithms.Controls_w0   = 0;
    S.algorithms.Controls_v0   = 0;
    %
    S.algorithms.mpc.Husky_w0  = 0;
    S.algorithms.mpc.Husky_v0  = 0;
end

function S = init_Michalek2017_robust_nSNT(S)
% robustnSNT_michalek(NumTT,Li,Lhi,b,r,integ)
    S.Michalek2017 = robustnSNTMichalek(S.config.N,S.system.Lh(:,2),S.system.Lh(:,1),S.system.b,S.system.r,S.config.Ts);
    set_wNr(S.Michalek2017, S.config.wNr);
    set_vNr(S.Michalek2017, S.config.vNr);
end

function S = init()
    % init fields of the structure
    S                 = struct;
    S.config          = struct;
    S.algorithms      = struct;
    S.noises          = struct;
    S.init_condition  = struct;
    S.path            = struct;
    S.obstacles       = struct;
    S.data            = struct;    
    S.plots           = struct;
    S.acado           = struct;
    S.acado.system    = struct;
    S.acado.states    = struct;
    S.acado.controls  = struct;
    S.exec_time       = struct;
    %
    S                 = reserve_notemp_memory(S);
    S                 = build_setup(S);
    S                 = gen_dynamic_and_solvers(S);    
    S.path.path       ='infinity';% 'rectangular';
    S.path.name       = [S.path.path,'-N=',num2str(S.config.N),'-Ne=',num2str(S.config.Ne),'-Nc=',num2str(S.config.Nc),'-Ts=',num2str(S.config.Ts),'.mat'];
    S                 = gen_path(S,S.path.path);
    S.config.reference  = false;        
    len                 = length(S.path.coordinates);
    S.path.ref          = [zeros(S.config.N,len); zeros(S.config.N+1,len); S.path.coordinates; S.path.coordinates];
    S.path.obstacles  = [];                                           % No obstacles               
    if strcmp(S.path.path,'infinity')
        if S.config.SIM
            S = gen_obstacle(S,1e6, 1e6, 0.35);
%             S = gen_obstacle(S,5.5,4.5, 0.5);
%             S = gen_obstacle(S,9,8.5, 0.45);
%             S = gen_obstacle(S,0,6, 0.55);
        else
            S = gen_obstacle(S,1e6, 1e6, 0.35);
        end
    elseif strcmp(S.path.path,'flat_infinity')        
        if S.config.SIM
            S = gen_obstacle(S,1e6, 1e6, 0.35);
        else
            S = gen_obstacle(S,1e6, 1e6, 0.35);
        end
    elseif strcmp(S.path.path,'segment') || strcmp(S.path.path,'invert_segment')
        if S.config.SIM
            S = gen_obstacle(S,3, 6, 0.35);
            % S = gen_obstacle(S,5, 6, 0.5);
            S = gen_obstacle(S,8, 6, 0.45);
        else
            S = gen_obstacle(S,1e6, 1e6, 0.35);
        end        
    elseif strcmp(S.path.path,'rectangular')
        if S.config.SIM

        else           
            S = gen_obstacle(S,1e6, 1e6, 0.35);
        end
%         S = gen_obstacle(S,10,3, 0.25);
    elseif strcmp(S.path.path,'test')

    elseif strcmp(S.path.path,'monaco')

    elseif strcmp(S.path.path,'test_betas')

    elseif strcmp(S.path.path,'rectangular_terrace')        
        if ~S.config.reference
            S = gen_obstacle(S,1e6, 1e6, 0.35);
        else
            S = gen_obstacle(S,1e6, 1e6, 0.25);
%             S = gen_obstacle(S,6, 31.75, 0.25);
%             S = gen_obstacle(S,0.8,30.75, 0.25);
%             S = gen_obstacle(S,0.8,25.75, 0.25);
        end
    elseif strcmp(S.path.path,'complex1_terrace')
        
    elseif strcmp(S.path.path,'circular_terrace')           
        % S = gen_obstacle(S,1e6, 1e6, 0.35);
        S = gen_obstacle(S,5.5, 8.5, 0.35);
        S = gen_obstacle(S,1.5, 4.5, 0.35);
    elseif strcmp(S.path.path,'ellipsoidal1_terrace')           

    elseif strcmp(S.path.path,'ellipsoidal2_terrace')           

    else  
    end
end

function S = gen_init_conditions(S)
    S   = gen_x0(S); 
    S   = gen_x0bar(S);
    S   = init_mheCasadi(S);
%     S   = init_ekfCasadi(S);
    S   = fill_mhe(S);
    S   = init_mpc(S);
%     S   = init_Michalek2017_robust_nSNT(S);
end

function S = init_flags_and_counters(S)
%
    S.path.reach_end_mhempc       = false;
    %
    S.config.time     = 0;
    S.config.iters    = 0;
end

function S = reserve_notemp_memory(S)
    S.data.mhempc.performance                       = struct;
    S.data.mhempc.performance.Psi_e      = [];
    S.data.mhempc.performance.Psi_e_vec  = [];
    S.data.mhempc.performance.maxPsi_e   = [];
    S.data.mhempc.performance.Psi_u      = [];    
    S.data.mhempc.performance.Psi_u_vec  = [];    
    S.data.mhempc.performance.maxPsi_u   = [];    
    S.data.mhempc.performance.time_tot   = [];
    S.data.mhempc.performance.est_err    = [];
end

function S = reserve_temp_memory(S)
    % Variables: states, outputs, etc. ____________________________________
    S.data.xsim                   = [];
    S.data.ysim                   = [];
    S.sensors.velocities          = [];
    S.data.references             = [];
    %
    S.data.slip                   = [];
    S.data.procDist               = [];
    %
    S.data.measNoise              = [];
    S.data.UmeasNoise             = [];
    %
    S.data.test                   = [];
    %                        
    S.data.xest                   = [];
    %                        
    S.data.time_mhempc            = [];
    %
    S.path.references_mhempc      = [];%S.path.tgt;
    %
    S.path.reach_end_mhempc       = false;
    %
    S.path.indices_tgt_mhempc     = [];
    %
    S.path.vel_tgt                = [];
    S.path.angle_tgt              = [];
    % auxiliar indice shared by all algorithms
    S.path.indx_alg               = [1;1];
    S.path.last_tgt               = 0;
    %   
    S.exec_time.t_mhe             = [];
    S.exec_time.t_pfa             = [];
    S.exec_time.t_mpc             = [];
    S.exec_time.t_sensors         = [];
    S.exec_time.t_obsdetector     = [];
    S.exec_time.t_ctrl            = 0;
    S.exec_time.t_tot             = 0;
    S.exec_time.t_acum            = 0;
    S.exec_time.t_mis             = 0;
    %
    S.sensors.theta0              = [];
    S.sensors.theta2              = [];
    S.sensors.theta2              = [];
    S.sensors.vlp16               = {};
    % 
    S.path.counter                = 1;
    % Sensor data from bolie phones
    S.mobile.data                 = {};
end

function S = gen_x0(S)
    % DO NOT FORGET CHECKING FEASIBILITY OF INITIAL CONDITION!!!
    Dx                  = S.path.coordinates(1,2)-S.path.coordinates(1,1);
    Dy                  = S.path.coordinates(2,2)-S.path.coordinates(2,1);
    theta0              = atan2(Dy,Dx);
    thetas              = repmat(theta0,S.config.N+1,1);
    betas               = -diff(thetas);
    xy_0                = S.path.coordinates(:,1);
    %
    if S.config.N == 10
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_8                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_9                = xy_8 - [S.system.Lh9*cos(thetas(9))+S.system.L9*cos(thetas(10)); S.system.Lh9*sin(thetas(9))+S.system.L9*sin(thetas(10))];
        xy_N                = xy_9 - [S.system.Lh10*cos(thetas(10))+S.system.L10*cos(thetas(11)); S.system.Lh10*sin(thetas(10))+S.system.L10*sin(thetas(11))];
    elseif S.config.N == 9
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_8                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_N                = xy_8 - [S.system.Lh9*cos(thetas(9))+S.system.L9*cos(thetas(10)); S.system.Lh9*sin(thetas(9))+S.system.L9*sin(thetas(10))];
    elseif S.config.N == 8
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_N                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_6;xy_7;xy_N];
    elseif S.config.N == 7
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_N                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_6;xy_N];
    elseif S.config.N == 6
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_N                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_N];
    elseif S.config.N == 5
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_N                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_N];
    elseif S.config.N == 4
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_N                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_N];
    elseif S.config.N == 3
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_N                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_N];
    elseif S.config.N == 2
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_N                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_0toN = [xy_0;xy_1;xy_N];
    else
        xy_N                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_0toN = [xy_0;xy_N];
    end
    %
    S.init_condition.x0       = [ betas; thetas; xy_0toN];%; zeros(S.system.nu,1) ];
    %
    S.data.xsim(:,1)   = S.init_condition.x0;%[S.init_condition.x0; [1;1]];
%     S.data.xsim(:,1)   = S.path.ref(:,1);%[S.init_condition.x0; [1;1]];
end

function S = gen_x0bar(S)
    % DO NOT FORGET CHECKING FEASIBILITY OF INITIAL CONDITION!!!
    thetas  = repmat(S.init_condition.x0(S.config.N+1),S.config.N+1,1);
    thetas(2:end) = thetas(2:end) + (rand(S.config.N,1)-0.5) .* (2*S.config.initUncertainty);
    betas   = -diff(thetas);
    xy_0    = S.init_condition.x0(2*S.config.N+2:2*S.config.N+3) + S.config.noise_lvl(S.config.N+2:S.config.N+3).*randn(2,1);
    xy_0Clean = S.init_condition.x0(2*S.config.N+2:2*S.config.N+3);
    % Initial uncertainty set
    r       = sum(S.system.Li(1:S.config.N))+sum(S.system.Lhi(1:S.config.N));
    alpha   = S.config.initUncertainty*S.config.N;
    t       = -thetas(1)-alpha:0.05:-thetas(1)+alpha;
    x       = xy_0Clean(1) + r*cos(t);
    x       = [xy_0Clean(1) x];
    y       = xy_0Clean(2) + r*sin(t);
    y       = [xy_0Clean(2) y];
    S.config.initUncertaintySet.x = x;
    S.config.initUncertaintySet.y = y;
    %
    if S.config.N == 10
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_8                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_9                = xy_8 - [S.system.Lh9*cos(thetas(9))+S.system.L9*cos(thetas(10)); S.system.Lh9*sin(thetas(9))+S.system.L9*sin(thetas(10))];
        xy_N                = xy_9 - [S.system.Lh10*cos(thetas(10))+S.system.L10*cos(thetas(11)); S.system.Lh10*sin(thetas(10))+S.system.L10*sin(thetas(11))];
    elseif S.config.N == 9
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_8                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_N                = xy_8 - [S.system.Lh9*cos(thetas(9))+S.system.L9*cos(thetas(10)); S.system.Lh9*sin(thetas(9))+S.system.L9*sin(thetas(10))];
    elseif S.config.N == 8
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_N                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_6;xy_7;xy_N];
    elseif S.config.N == 7
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_N                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_6;xy_N];
    elseif S.config.N == 6
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_N                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_5;xy_N];
    elseif S.config.N == 5
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_N                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_4;xy_N];
    elseif S.config.N == 4
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_N                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_3;xy_N];
    elseif S.config.N == 3
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_N                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_0toN = [xy_0;xy_1;xy_2;xy_N];
    elseif S.config.N == 2
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_N                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_0toN = [xy_0;xy_1;xy_N];
    else
        xy_N                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_0toN = [xy_0;xy_N];
    end   
    
    if S.config.SIM
        S.init_condition.x0bar = [ betas; thetas; xy_0toN];%; zeros(S.system.nu,1) ];
% S.init_condition.x0bar = S.init_condition.x0;
    else    
        S                       = read_sensors(S);
        S                       = measurements_vector(S);
        
        xy_0                    = [S.ROS.sensors.rtk.x0; S.ROS.sensors.rtk.y0];
        xy_1                    = xy_0 - [S.system.Lh1*cos(S.ROS.sensors.vectornav_theta0)+S.system.L1*cos(S.ROS.sensors.vectornav_theta0); S.system.Lh1*sin(S.ROS.sensors.vectornav_theta0)+S.system.L1*sin(S.ROS.sensors.vectornav_theta0)];
        xy_N                    = xy_1 - [S.system.Lh2*cos(S.ROS.sensors.vectornav_theta0)+S.system.L2*cos(S.ROS.sensors.vectornav_theta0); S.system.Lh2*sin(S.ROS.sensors.vectornav_theta0)+S.system.L2*sin(S.ROS.sensors.vectornav_theta0)];

        S.init_condition.x0bar = [ S.ROS.sensors.encoders.beta1;    % The vehicle must start aligned, i.e., all segment with the same attitude
                                   S.ROS.sensors.encoders.beta2;
                                   repmat(S.ROS.sensors.vectornav_theta0,S.config.N+1,1);...
                                   xy_0;...
                                   xy_1;...
                                   xy_N];
        
    end    
% ####    
S.data.xest = S.init_condition.x0bar;
% ####
end
