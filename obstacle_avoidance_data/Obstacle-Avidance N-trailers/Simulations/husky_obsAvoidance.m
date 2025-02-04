% #########################################################################
% 
% Estimation and Control for Generalised N-Trailer Systems
% Author: Nestor. N. Deniz - 2022
%
% #########################################################################
% TODO:
% 1) CHECK AC
%
% 
% #########################################################################


function S = husky_obsAvoidance()
    clear all; clc;
    import casadi.*
    % ---------------------------------------------------------------------
    S = init();
    % Init ROS's things ---------------------------------------------------
    if S.config.dataOnline == true; figure('units','normalized','outerposition',[0 0 1 1]); hold on; 
        dim = [.2 .5 .3 .3];
%         ann = annotation('textbox',dim,'String','','verticalalignment','top','horizontalalignment','left','FitBoxToText','on','Interpreter','Latex','FontSize',30);
    end
    if ~S.config.SIM == true; S = ROS(S); end
    for num_sims = 1:S.config.NUM_SIMS
        S = call_init_functions(S);
        %
        S.config.iters = 1;
        while (S.config.time(end) < S.config.tf) && ~(S.path.reach_end_mhempc || ~S.config.mhe_mpc)
            % *****************************************************************
            % #################################################################
            % Our approach NMHE+NMPC
            % #################################################################            
            %
            if S.config.mhe_mpc == true
            if check_end_condition(S); break; else
                % Solve estimation problem ________________________________
                S = call_MHE(S);
                % Path-tracking problem ***********************************        
                S = call_PFA(S);
                % Solve control problem ___________________________________
                S = call_MPC(S);
                % Apply controls to the Husky _____________________________
                if S.config.SIM == true
                    simulation_input.x = S.data.xsim_mhempc(:,end);
                    simulation_input.u = S.algorithms.mpcCasadi.u_k;
                    states             = S.dynamic.FNt(simulation_input.x,simulation_input.u);
                    S.data.xsim_mhempc = [S.data.xsim_mhempc, full(states)];
                else
                    S = write_ctrl(S, S.algorithms.mpc.Husky_w0, S.algorithms.mpc.Husky_v0);
                end
                % Update measurements _____________________________________
                S = update_measurements(S);
                % Perform real-time iteration _____________________________
                tic;
%                 S.exec_time.t_tot   = [S.exec_time.t_tot, S.exec_time.t_mhe(end)+S.exec_time.t_pfa(end)+S.exec_time.t_mpc(end)+S.exec_time.t_ctrl(end)+S.exec_time.t_mis(end)];
                S.exec_time.t_tot   = [S.exec_time.t_tot, S.exec_time.t_mhe(end)+S.exec_time.t_mpc(end)+S.exec_time.t_ctrl(end)+S.exec_time.t_mis(end)];
                S.exec_time.t_acum  = S.exec_time.t_tot(end);
                while ((S.exec_time.t_acum + toc) < S.config.Ts*0.8) && ~S.config.SIM; end
                % 
            end
            end            
            % #############################################################
            %
            % #############################################################        
            S.config.time   = [S.config.time, S.config.time(end)+S.config.Ts];
            S.config.iters  = S.config.iters+1;
            a               = [num_sims, S.config.iters];
            b = [norm(S.algorithms.mpcCasadi.qref-S.algorithms.mpcCasadi.qRefOpt1), norm(S.algorithms.mpcCasadi.qref-S.algorithms.mpcCasadi.qRefOpt2)]
            % -------------------------------------------------------------
            if S.config.dataOnline == true
                plot_mono2(S, S.data.xsim_mhempc(:,end));
                %
                plot_mono2(S, S.data.xest_mhempc(:,end));
%                 str = {['$\omega_0=\,$',num2str(S.data.xest_mhempc(2*S.config.num_trailers+6,end))],...
%                     ['$v_0=\,$',num2str(S.data.xest_mhempc(2*S.config.num_trailers+7,end))],...
%                     ['$\theta_0=\,$',num2str(S.data.xest_mhempc(S.config.num_trailers+1,end))],...
%                     ['$x_N=\,$',num2str(S.data.xest_mhempc(2*S.config.num_trailers+4,end))],...
%                     ['$x_{N_{ref}}=\,$',num2str(S.controller.ref(2*S.config.num_trailers+4,end))],...
%                     ['$y_N=\,$',num2str(S.data.xest_mhempc(2*S.config.num_trailers+5,end))],...      
%                     ['$y_{N_{ref}}=\,$',num2str(S.controller.ref(2*S.config.num_trailers+5,end))],...
%                     ['$x_0=\,$',num2str(S.data.xest_mhempc(2*S.config.num_trailers+2,end))],...
%                     ['$y_0=\,$',num2str(S.data.xest_mhempc(2*S.config.num_trailers+3,end))],...
%                     ['$rtf=\,$',num2str((S.exec_time.t_acum + toc)/S.config.Ts)] };
%                 set(ann,'String',str);
            end
            S.exec_time.t_mis = [S.exec_time.t_mis, toc];
            % -------------------------------------------------------------
            if S.config.mhe_mpc == true
                S.data.mhempc.performance.xfut{num_sims, S.config.iters-1}  = S.algorithms.mpcCasadi.Qtraj;
                S.data.mhempc.performance.J{num_sims, S.config.iters-1}     = S.algorithms.mpcCasadi.Jnum;
                S.data.mhempc.performance.qref{num_sims, S.config.iters-1}  = S.algorithms.mpcCasadi.qref;
                S.data.mhempc.performance.uref{num_sims, S.config.iters-1}  = S.algorithms.mpcCasadi.uref;
                S.data.mhempc.performance.XYR{num_sims, S.config.iters-1}   = S.algorithms.mpcCasadi.XYRobst;
            end
        end
        if ~S.config.SIM
            S = write_ctrl(S, 0, 0); % Stop the vehicle
        end
        % -----------------------------------------------------------------
        if S.config.mhe_mpc == true
%             S.data.mhempc.performance.time_tot      = [S.data.mhempc.performance.time_tot, mean(S.data.time_mhempc)/S.config.Ts];
%             [ctrl_eff, max_ctrl_eff, eff_vec]       = control_effort(S, S.data.ysim_mhempc(S.config.num_trailers+4:S.config.num_trailers+5,:));
%             S.data.mhempc.performance.Psi_u_vec     = [S.data.mhempc.performance.Psi_u_vec, eff_vec, -1];
%             S.data.mhempc.performance.Psi_u         = [S.data.mhempc.performance.Psi_u, ctrl_eff];
%             S.data.mhempc.performance.maxPsi_u      = [S.data.mhempc.performance.maxPsi_u, max_ctrl_eff];
%             [err, max_err, err_vec]                 = compute_tracking_on_Nt_error(S, S.path.indices_tgt_mhempc, S.data.ysim_mhempc, 'y');
%             S.data.mhempc.performance.Psi_e_vec     = [S.data.mhempc.performance.Psi_e_vec, err_vec, -1]; % the -1 is used as delimiter
%             S.data.mhempc.performance.Psi_e         = [S.data.mhempc.performance.Psi_e, err];
%             S.data.mhempc.performance.maxPsi_e      = [S.data.mhempc.performance.maxPsi_e, max_err];
%             S.data.mhempc.performance.est_err       = [S.data.mhempc.performance.est_err, compute_estimation_error(S, S.data.xsim_mhempc, S.data.xest_mhempc)];
            S.data.mhempc.performance.xsim{num_sims} = S.data.xsim_mhempc;            
            S.data.mhempc.performance.ysim{num_sims} = S.data.ysim_mhempc;
            S.data.mhempc.performance.xest{num_sims} = S.data.xest_mhempc;
        end
        if S.config.save_workspace == true; clk = clock; save([date,'-',num2str(clk(1:5))]); end
        if ~S.config.reference; ref = S.data.xsim_mhempc; save(S.path.name,'ref'); end
    end
end

function S = build_setup(S)
    % SIMULATION PARAMETERS ===============================================
    S.ROS.IMU_ZEROING_VEC1  = 1.16;
    S.ROS.IMU_ZEROING_VEC2  = 2.71;
    S.ROS.IMU_ZEROING_MIC1  = -0.21;
    S.ROS.ENCODERS_PULSE_TO_RAD = 0.6 * pi / 180;
    % Solver ______________________________________________________________
    S.config.solver         = 'casadi'; % options: 'casadi', 'acado'
    % Simulation or field experiment (Husky)_______________________________
    S.config.SIM            = true;
    %
    S.config.num_trailers   = 2;                       % number of trailers
    %
    S.config.verbose        = false;
    S.config.calcMtxConvCoord = true; if ~S.config.calcMtxConvCoord; warning('WARNING!!! The matrix for correcting x-y coordinates is not being computed...'); end;
    S.config.t0             = 0;                                        % initial time of sym. [Seg]
    S.config.tf             = 1000;                                     % final time of sym. [Seg]
    S.config.Ts             = 0.3;                                      % sampling period, 0=>discrete-time system
    S.config.N              = [];                                       % number of steps of sym.
    S.config.same_seed      = false;
    S.config.EXPORT         = false; if ~S.config.EXPORT; warning('WARNING!!! MHE and MPC were not compiled...'); end
    S.config.iters          = 0;
    S.config.time           = 0;
    S.config.NUM_SIMS       = 1;
    S.config.model.uncty    = false;
    S.config.model.dev      = 15;                                       % porcentual value of uncertainty
    S.config.updtAC         = false;
    S.config.outputs        = [1:S.config.num_trailers+1,2*S.config.num_trailers+2:2*S.config.num_trailers+3,2*S.config.num_trailers+6:2*S.config.num_trailers+7];
    S.config.Nc             = 5;                                       % Lenght of the control horizon
    S.config.Ne             = 5;                                        % Lenght of the estimation window
    % Algorithms to simulate ______________________________________________
    S.config.mhe_mpc        = true;
    % Noise amplitudes ____________________________________________________
    S.config.noise          = 'gaussian'; % 'gaussian' or 'uniform'
    S.config.noise_lvl      = 0.*[(1.2*pi/180).*ones(S.config.num_trailers, 1); 2*pi/180; [0.05; 0.05]; [0.1; 0.15]]; % measurement noise amplitude: [betas'; theta_0; x_0; y_0; w0; v0]
    % Reference velocities ________________________________________________
    S.config.vNr            = 0.05;%0.15;%/S.config.Ts;%0.1;
    S.config.wNr            = 0;
    S.config.uref           = [S.config.wNr; S.config.vNr];
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
    % Obstacles ___________________________________________________________
    S.config.distToObstacle = 10;       % Distance in meters to be aware of obstacles
    S.config.distFarAwayObs = 1.05;
    S.config.delta          = 2;
    S.config.r_safe         = 1;
    S.config.Qbar0          = 0.5;  %max([S.algorithms.mpcCasadi.Q(2*S.config.num_trailers+2,2*S.config.num_trailers+2),S.algorithms.mpcCasadi.Q(2*S.config.num_trailers+3,2*S.config.num_trailers+3)]);
    S.config.QbarN          = 2;    %max([S.algorithms.mpcCasadi.Q(2*S.config.num_trailers+4,2*S.config.num_trailers+4),S.algorithms.mpcCasadi.Q(2*S.config.num_trailers+5,2*S.config.num_trailers+5)]);
    % Distance for updating next coordinate _______________________________
    S.config.updt_coord     = 1.05;%0.25;
    % Distance to next coordinate _________________________________________
    S.config.coord_ahead    = 0;%1.25;
    % Distance to last target for finishing navigation and control ________
    S.config.dist_last_tgt  = 0.15;                                    % Distance in meters to next coordinate    
    % Plot data online during experiments. USeful for debbuging purposes __
    S.config.dataOnline     = true;
    % Save workspace flag _________________________________________________
    S.config.save_workspace = false;
end

function plot_mono(S, qk, scale, clr)
    if nargin ==4
        clrHb = clr;
        clrHr = clr;
        clrTb = clr;
        clrTr = clr;
        clrTN = clr;
        clrL  = clr;
        clrAw = clr;
        noRef = true;
    else
        clrHb = 'y';
        clrHr = 'k';
        clrTb = 'c';
        clrTr = 'k';
        clrTN = 'r';
        clrL  = 'k';
        clrAw = 'k';
        noRef = false;
    end
    %
    if ~noRef
        plot(S.path.coordinates(1,:),S.path.coordinates(2,:),'color',[0.6 0.6 0.6],'LineWidth',8); grid on; daspect([1 1 1]); hold on;
        plot(S.algorithms.mpcCasadi.Qtraj(2*S.config.num_trailers+2,:),S.algorithms.mpcCasadi.Qtraj(2*S.config.num_trailers+3,:),'y','LineWidth',5);
        plot(S.algorithms.mpcCasadi.Qtraj(2*S.config.num_trailers+4,:),S.algorithms.mpcCasadi.Qtraj(2*S.config.num_trailers+5,:),'r','LineWidth',5);
%         xlim([-1 10]); ylim([-15 35]); daspect([1 1 1]);
        xlim([min(S.path.coordinates(1,:))-1 max(S.path.coordinates(1,:))+1]); ylim([min(S.path.coordinates(2,:))-1 max(S.path.coordinates(2,:))+1]);
        if~isempty(S.path.obstacles)
            nroObs = size(S.path.obstacles,1);
            for i=1:nroObs
                circles(S.path.obstacles(i,1),S.path.obstacles(i,2),S.path.obstacles(i,3),'color','red')
            end
        end
    end    
    % Dimensions of the vehice
    widthH = 0.42;
    lengthH = 0.990;
    widthHr = 0.125;
    lengthHr = 0.33;
    lengthTr = 0.33;
    widthTr = 0.125;
    posHrCenY = 0.272;
    trailer_width = 0.32;
    % Husky's Wheel
    Praux  = [-widthHr/2, -widthHr/2, widthHr/2, widthHr/2; -lengthHr/2, lengthHr/2, lengthHr/2, -lengthHr/2].*scale;
    % Trailer's Wheel
    Prt  = [-widthTr/2, -widthTr/2, widthTr/2, widthTr/2; -lengthTr/2, lengthTr/2, lengthTr/2, -lengthTr/2].*scale;
    %
    Phaux  = [-widthH/2, -widthH/2, widthH/2, widthH/2; -lengthH/2, lengthH/2, lengthH/2, -lengthH/2].*scale;
    Prh    = [];
    Ph     = [];
    % Trailer's body
    % Rotate the Husky and wheels according to their attitude
    theta  = qk(S.config.num_trailers+1:2*S.config.num_trailers+1);
    R      = [cos(theta(1)-pi/2), -sin(theta(1)-pi/2); sin(theta(1)-pi/2), cos(theta(1)-pi/2)];
    for i=1:4
        Prh = [Prh, R*Praux(:,i)];
        Ph = [Ph, R*Phaux(:,i)];
    end
    % now rotate the trailers and their wheels
    Ptaux  = [];
    Ptaux2 = [];        
    Ptraux = [];
    for j=1:S.config.num_trailers
        % rotation matrix
        R      = [cos(theta(j+1)-pi/2), -sin(theta(j+1)-pi/2); sin(theta(j+1)-pi/2), cos(theta(j+1)-pi/2)];
        % Trailer's bodies
%         Ptaux2 = [Ptaux2; [-trailer_width/2, -trailer_width/2, trailer_width/2, trailer_width/2; -S.system.Lhi(j+1), S.system.Li(j), S.system.Li(j), -S.system.Lhi(j+1) ]];
        Ptaux2 = [Ptaux2; [-trailer_width/2, -trailer_width/2, trailer_width/2, trailer_width/2; 0, S.system.Li(j), S.system.Li(j), 0 ].*scale];
        Ptaux3 = [];
        Ptraux2 = [];
        for i=1:4
            Ptaux3 = [Ptaux3, R*Ptaux2((j-1)*2+1:j*2,i)];
            Ptraux2 = [Ptraux2, R*Prt(:,i)];
        end
        Ptaux = [Ptaux; Ptaux3];
        Ptraux = [Ptraux; Ptraux2];
    end
    % Husky's coordinates   
    P0c     = [qk(2*S.config.num_trailers+2); qk(2*S.config.num_trailers+3)];
    % Inter segments coordinates
    xy_i    = zeros(2,S.config.num_trailers+1);
    xy_i(:,1) = P0c;
    for i=1:S.config.num_trailers
        xy_i(:,i+1) = xy_i(:,i) - [S.system.Lhi(i)*cos(theta(i)) + S.system.Li(i)*cos(theta(i+1)); S.system.Lhi(i)*sin(theta(i)) + S.system.Li(i)*sin(theta(i+1))]; 
    end
    % Last trailer coordinates
    PNc     = [qk(2*S.config.num_trailers+4); qk(2*S.config.num_trailers+5)];    
    % angle that form each wheel respectthe centre of the Husky
    alpha   = atan2(posHrCenY, (widthH+widthHr)/2);
    % Husky's wheels coordinates
    Xrc     = -(widthH/2 + widthHr/2);
    Yrc     = (widthH/2 + widthHr/2);
    %
    Prh1    = Prh + repmat(P0c,1,4) + repmat([cos(theta(1)-pi/2), -sin(theta(1)-pi/2); sin(theta(1)-pi/2), cos(theta(1)-pi/2)] * [-Xrc; -Yrc],1,4);
    Prh1    = [Prh1, Prh1(:,1)];
    Prh2    = Prh + repmat(P0c,1,4) + repmat([cos(theta(1)-pi/2), -sin(theta(1)-pi/2); sin(theta(1)-pi/2), cos(theta(1)-pi/2)] * [-Xrc; Yrc],1,4);
    Prh2    = [Prh2, Prh2(:,1)];
    Prh3    = Prh + repmat(P0c,1,4) + repmat([cos(theta(1)-pi/2), -sin(theta(1)-pi/2); sin(theta(1)-pi/2), cos(theta(1)-pi/2)] * [Xrc; Yrc],1,4);
    Prh3    = [Prh3, Prh3(:,1)];
    Prh4    = Prh + repmat(P0c,1,4) + repmat([cos(theta(1)-pi/2), -sin(theta(1)-pi/2); sin(theta(1)-pi/2), cos(theta(1)-pi/2)] * [Xrc; -Yrc],1,4);
    Prh4    = [Prh4, Prh4(:,1)];
    %
    Phusky  = Ph + repmat(P0c,1,4);
    Phusky  = [Phusky, Phusky(:,1)];
    % Wheels of the trailers
    Pitrailer    = [];
    Pr1itrailer  = [];
    Pr2itrailer  = [];
    PrtrailerAxe = [];
    alpha        = atan2(S.system.Lhi(i+1)/2, trailer_width/1.75);
    for i=1:S.config.num_trailers
        %
        Pitrailer       = [Pitrailer; Ptaux((i-1)*2+1:i*2, :) + repmat(xy_i(:,i+1),1,4)];
        PrtrailerAxe    = [PrtrailerAxe, xy_i(:,i+1)+[-(trailer_width/2+S.system.Lhi(i+1)/2)*cos(theta(i+1)-pi/2+alpha); -(trailer_width/2+S.system.Lhi(i+1)/2)*sin(theta(i+1)-pi/2+alpha)],xy_i(:,i+1)+[(trailer_width/2+S.system.Lhi(i+1)/2)*cos(theta(i+1)-pi/2+alpha); (trailer_width/2+S.system.Lhi(i+1)/2)*sin(theta(i+1)-pi/2+alpha)] ];
        Pr1itrailer     = [Pr1itrailer; Ptraux((i-1)*2+1:i*2, :) + repmat(xy_i(:,i+1),1,4) + repmat([-trailer_width/1.75*cos(theta(i+1)-pi/2+alpha); -trailer_width/1.75*sin(theta(i+1)-pi/2+alpha)],1,4)];
        Pr2itrailer     = [Pr2itrailer; Ptraux((i-1)*2+1:i*2, :) + repmat(xy_i(:,i+1),1,4) + repmat([trailer_width/1.75*cos(theta(i+1)-pi/2+alpha); trailer_width/1.75*sin(theta(i+1)-pi/2+alpha)],1,4)];
    end
    PrtrailerAxe = [PrtrailerAxe, xy_i(:,S.config.num_trailers+1)+[-(trailer_width/2+S.system.Lhi(S.config.num_trailers+1)/2)*cos(theta(S.config.num_trailers+1)-pi/2+alpha); -(trailer_width/2+S.system.Lhi(S.config.num_trailers+1)/2)*sin(theta(S.config.num_trailers+1)-pi/2+alpha)],xy_i(:,S.config.num_trailers+1)+[(trailer_width/2+S.system.Lhi(S.config.num_trailers+1)/2)*cos(theta(S.config.num_trailers+1)-pi/2+alpha); (trailer_width/2+S.system.Lhi(S.config.num_trailers+1)/2)*sin(theta(S.config.num_trailers+1)-pi/2+alpha)] ];
    Pr1itrailer = [Pr1itrailer; Ptraux((S.config.num_trailers-1)*2+1:S.config.num_trailers*2, :) + repmat(xy_i(:,S.config.num_trailers+1),1,4) + repmat([-trailer_width/1.75*cos(theta(S.config.num_trailers+1)-pi/2+alpha); -trailer_width/1.75*sin(theta(S.config.num_trailers+1)-pi/2+alpha)],1,4)];
    Pr2itrailer = [Pr2itrailer; Ptraux((S.config.num_trailers-1)*2+1:S.config.num_trailers*2, :) + repmat(xy_i(:,S.config.num_trailers+1),1,4) + repmat([trailer_width/1.75*cos(theta(S.config.num_trailers+1)-pi/2+alpha); trailer_width/1.75*sin(theta(S.config.num_trailers+1)-pi/2+alpha)],1,4)];
    %
    Pitrailer = [Pitrailer, Pitrailer(:,1)];
    Pr1itrailer = [Pr1itrailer, Pr1itrailer(:,1)];
    Pr2itrailer = [Pr2itrailer, Pr2itrailer(:,1)];
    % Now the last trailer, fro whom the coordiantes are directly estimated
    % from the nmhe
    PNtrailer = Ptaux(end-1:end,:) + repmat(PNc,1,4);
    PNtrailer = [PNtrailer, PNtrailer(:,1)];
    % plt    
    if ~ noRef
        plot(S.path.references_mhempc(2*S.config.num_trailers+2,end),S.path.references_mhempc(2*S.config.num_trailers+3,end),'y+','linewidth',3,'markersize',20)
        plot(S.path.references_mhempc(2*S.config.num_trailers+2,end),S.path.references_mhempc(2*S.config.num_trailers+3,end),'yo','linewidth',3,'markersize',20)
        plot(S.path.references_mhempc(2*S.config.num_trailers+4,end),S.path.references_mhempc(2*S.config.num_trailers+5,end),'rx','linewidth',3,'markersize',20)
        plot(S.path.references_mhempc(2*S.config.num_trailers+4,end),S.path.references_mhempc(2*S.config.num_trailers+5,end),'ro','linewidth',3,'markersize',20)
        %
        plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+2,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+3,end),'c+','linewidth',3,'markersize',20)
        plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+2,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+3,end),'co','linewidth',3,'markersize',20)
        plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+4,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+5,end),'mx','linewidth',3,'markersize',20)
        plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+4,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+5,end),'mo','linewidth',3,'markersize',20)
    end
    %
    patch(Phusky(1,:),Phusky(2,:),clrHb);
    patch(Prh1(1,:),Prh1(2,:),clrHr);
    patch(Prh2(1,:),Prh2(2,:),clrHr);
    patch(Prh3(1,:),Prh3(2,:),clrHr);
    patch(Prh4(1,:),Prh4(2,:),clrHr);
    %
    for i=1:S.config.num_trailers-1        
        % Wheels of trailers 1 to N-1
        patch(Pr1itrailer((i-1)*2+1,:),Pr1itrailer(i*2,:),clrTr);
        patch(Pr2itrailer((i-1)*2+1,:),Pr2itrailer(i*2,:),clrTr);
        % Axe connecting adyacent trailers
        line(xy_i(1,i:i+1), xy_i(2,i:i+1),'color',clrL,'linewidth',2)
        % Body of trailer 1 to N-1
        patch(Pitrailer((i-1)*2+1,:),Pitrailer(i*2,:),clrTb);
        line(Pitrailer((i-1)*2+1,:),Pitrailer(i*2,:),'color','k','linewidth',1);
        % Axe between wheels of each trailer
        line(PrtrailerAxe(1,(i-1)*2+1:i*2), PrtrailerAxe(2,(i-1)*2+1:i*2),'color',clrAw,'linewidth',2)        
    end
    % Wheels of last trailer
    patch(Pr1itrailer((S.config.num_trailers-1)*2+1,:),Pr1itrailer(S.config.num_trailers*2,:),clrTr);
    patch(Pr2itrailer((S.config.num_trailers-1)*2+1,:),Pr2itrailer(S.config.num_trailers*2,:),clrTr);
    % Axe connecting adyacent trailers
    line(xy_i(1,S.config.num_trailers:S.config.num_trailers+1), xy_i(2,S.config.num_trailers:S.config.num_trailers+1),'color',clrL,'linewidth',2)
    % Body of the last trailer
    patch(PNtrailer(1,:),PNtrailer(2,:),clrTN);
    line(PNtrailer(1,:),PNtrailer(2,:),'color','k','linewidth',1);
    % Axe between wheels of the last trailer
    line(PrtrailerAxe(1,(S.config.num_trailers-1)*2+1:S.config.num_trailers*2), PrtrailerAxe(2,(S.config.num_trailers-1)*2+1:S.config.num_trailers*2),'color',clrAw,'linewidth',2)
    %
%     if ~ noRef
%         plot(S.path.references_mhempc(2*S.config.num_trailers+2,end),S.path.references_mhempc(2*S.config.num_trailers+3,end),'y+','linewidth',3,'markersize',20)
%         plot(S.path.references_mhempc(2*S.config.num_trailers+2,end),S.path.references_mhempc(2*S.config.num_trailers+3,end),'yo','linewidth',3,'markersize',20)
%         plot(S.path.references_mhempc(2*S.config.num_trailers+4,end),S.path.references_mhempc(2*S.config.num_trailers+5,end),'rx','linewidth',3,'markersize',20)
%         plot(S.path.references_mhempc(2*S.config.num_trailers+4,end),S.path.references_mhempc(2*S.config.num_trailers+5,end),'ro','linewidth',3,'markersize',20)
%         %
%         plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+2,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+3,end),'c+','linewidth',3,'markersize',20)
%         plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+2,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+3,end),'co','linewidth',3,'markersize',20)
%         plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+4,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+5,end),'mx','linewidth',3,'markersize',20)
%         plot(S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+4,end),S.algorithms.mpcCasadi.qRefOpt2(2*S.config.num_trailers+5,end),'mo','linewidth',3,'markersize',20)
%     end
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
    if ~noRef
        plot(S.path.coordinates(1,:),S.path.coordinates(2,:),'color',[0.6 0.6 0.6],'LineWidth',8); grid on; daspect([1 1 1]); hold on;
        xlim([min(S.path.coordinates(1,:))-1 max(S.path.coordinates(1,:))+1]); ylim([min(S.path.coordinates(2,:))-1 max(S.path.coordinates(2,:))+1]);
        if~isempty(S.path.obstacles)
            nroObs = size(S.path.obstacles,1);
            for i=1:nroObs
                circles(S.path.obstacles(i,1),S.path.obstacles(i,2),S.path.obstacles(i,3),'color','red')
            end
        end
%         patch(S.config.initUncertaintySet.x,S.config.initUncertaintySet.y,'r','edgecolor','none','facealpha',0.2)
    end    
    % Rotate the Husky and wheels according to their attitude
    thetas  = qk(S.config.num_trailers+1:2*S.config.num_trailers+1);
    betas   = qk(1:S.config.num_trailers);    
    xy0     = qk(2*S.config.num_trailers+2:2*S.config.num_trailers+3);
    xyi     = xy0;
    for i=1:S.config.num_trailers
        xyi(:,i+1) = xyi(:,i) - [S.system.Lhi(i)*cos(thetas(i)) + S.system.Li(i)*cos(thetas(i+1)); S.system.Lhi(i)*sin(thetas(i)) + S.system.Li(i)*sin(thetas(i+1))]; 
    end

xyi(:,end) = qk(2*S.config.num_trailers+4:2*S.config.num_trailers+5);

    R       = [cos(thetas(1)-pi/2), -sin(thetas(1)-pi/2); sin(thetas(1)-pi/2), cos(thetas(1)-pi/2)];
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
    for i=1:S.config.num_trailers
        R = [cos(thetas(i+1)-pi/2), -sin(thetas(i+1)-pi/2); sin(thetas(i+1)-pi/2), cos(thetas(i+1)-pi/2)];
        trailerAxePlt = R * S.system.XYtrailerAxe + xyi(:,i+1);
        trailerWheelLeftPlt = R * S.system.XYtrailerWheelLeft + xyi(:,i+1);
        trailerWheelRightPlt = R * S.system.XYtrailerWheelRight + xyi(:,i+1);
        trailerLongAxePlt = R * S.system.XYtrailerLongAxe((i-1)*2+1:i*2,:) + xyi(:,i+1);
        trailerLoadPlt = R * S.system.XYtrailerLoad((i-1)*2+1:i*2,:) + xyi(:,i+1);        
        if i~= S.config.num_trailers
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
    plot(S.path.references_mhempc(2*S.config.num_trailers+2,end),S.path.references_mhempc(2*S.config.num_trailers+3,end),'y+','linewidth',1.5,'markersize',20)
    plot(S.path.references_mhempc(2*S.config.num_trailers+2,end),S.path.references_mhempc(2*S.config.num_trailers+3,end),'yo','linewidth',1.5,'markersize',20)
    plot(S.path.references_mhempc(2*S.config.num_trailers+4,end),S.path.references_mhempc(2*S.config.num_trailers+5,end),'rx','linewidth',1.5,'markersize',20)
    plot(S.path.references_mhempc(2*S.config.num_trailers+4,end),S.path.references_mhempc(2*S.config.num_trailers+5,end),'ro','linewidth',1.5,'markersize',20)
    %
    if ~noRef
        hold off;
    end
    %
    drawnow limitrate
end

function S = plot_fo(S)
figure(2);
    xmin    = min(S.path.coordinates(1,:))-1;
    xmax    = max(S.path.coordinates(1,:))+1;
    ymin    = min(S.path.coordinates(2,:))-1;
    ymax    = max(S.path.coordinates(2,:))+3;
    xx      = xmin:0.25:xmax;
    yy      = ymin:0.25:ymax;
    [X,Y]   = meshgrid(xx,yy);
    for i=1:length(S.data.mhempc.performance.xsim{1})-1
        Z0posL = [];
        ZobsL  = [];
        ZNposL  = [];
        for j=1:S.config.Nc
            Z0pos = full(S.algorithms.mpcCasadi.Q(2*S.config.num_trailers+2,2*S.config.num_trailers+2)*(X-S.data.mhempc.performance.xfut{i}(2*S.config.num_trailers+2,j)).^2 + ...
                         S.algorithms.mpcCasadi.Q(2*S.config.num_trailers+3,2*S.config.num_trailers+3)*(Y-S.data.mhempc.performance.xfut{i}(2*S.config.num_trailers+3,j )).^2);
            Z0posL = [Z0posL, full(S.algorithms.mpcCasadi.Q(2*S.config.num_trailers+2,2*S.config.num_trailers+2)*(S.algorithms.mpcCasadi.qref(2*S.config.num_trailers+2)-S.data.mhempc.performance.xfut{i}(2*S.config.num_trailers+2,j)).^2 + ...
                         S.algorithms.mpcCasadi.Q(2*S.config.num_trailers+3,2*S.config.num_trailers+3)*(S.algorithms.mpcCasadi.qref(2*S.config.num_trailers+3)-S.data.mhempc.performance.xfut{i}(2*S.config.num_trailers+3,j )).^2);];
            Zobs = full(15*exp(-( (X - S.data.mhempc.performance.XYR{i}(1,1)).^2 /(2 * (S.data.mhempc.performance.XYR{i}(1,3))^2) + (Y - S.data.mhempc.performance.XYR{i}(1,2)).^2 /(2 * (S.data.mhempc.performance.XYR{i}(1,3))^2) )) +... % First obstacle
                        15*exp(-( (X - S.data.mhempc.performance.XYR{i}(2,1)).^2 /(2 * (S.data.mhempc.performance.XYR{i}(2,3))^2) + (Y - S.data.mhempc.performance.XYR{i}(2,2)).^2 /(2 * (S.data.mhempc.performance.XYR{i}(2,3))^2) )) +... % Second obstacle
                        15*exp(-( (X - S.data.mhempc.performance.XYR{i}(3,1)).^2 /(2 * (S.data.mhempc.performance.XYR{i}(3,3))^2) + (Y - S.data.mhempc.performance.XYR{i}(3,2)).^2 /(2 * (S.data.mhempc.performance.XYR{i}(3,3))^2) )));
            ZobsL = [ZobsL, full(15*exp(-( (S.algorithms.mpcCasadi.qref(2*S.config.num_trailers+2) - S.data.mhempc.performance.XYR{i}(1,1)).^2 /(2 * (S.data.mhempc.performance.XYR{i}(1,3))^2) + (S.algorithms.mpcCasadi.qref(2*S.config.num_trailers+3) - S.data.mhempc.performance.XYR{i}(1,2)).^2 /(2 * (S.data.mhempc.performance.XYR{i}(1,3))^2) )) +... % First obstacle
                        15*exp(-( (S.algorithms.mpcCasadi.qref(2*S.config.num_trailers+2) - S.data.mhempc.performance.XYR{i}(2,1)).^2 /(2 * (S.data.mhempc.performance.XYR{i}(2,3))^2) + (S.algorithms.mpcCasadi.qref(2*S.config.num_trailers+3) - S.data.mhempc.performance.XYR{i}(2,2)).^2 /(2 * (S.data.mhempc.performance.XYR{i}(2,3))^2) )) +... % Second obstacle
                        15*exp(-( (S.algorithms.mpcCasadi.qref(2*S.config.num_trailers+2) - S.data.mhempc.performance.XYR{i}(3,1)).^2 /(2 * (S.data.mhempc.performance.XYR{i}(3,3))^2) + (S.algorithms.mpcCasadi.qref(2*S.config.num_trailers+3) - S.data.mhempc.performance.XYR{i}(3,2)).^2 /(2 * (S.data.mhempc.performance.XYR{i}(3,3))^2) )));];
            
            ZNpos = full(S.algorithms.mpcCasadi.Q(2*S.config.num_trailers+4,2*S.config.num_trailers+4)*(X-S.data.mhempc.performance.xfut{i}(2*S.config.num_trailers+4,j )).^2 + ...
                         S.algorithms.mpcCasadi.Q(2*S.config.num_trailers+5,2*S.config.num_trailers+5)*(Y-S.data.mhempc.performance.xfut{i}(2*S.config.num_trailers+5,j )).^2);
            ZNposL = [ZNposL, full(S.algorithms.mpcCasadi.Q(2*S.config.num_trailers+4,2*S.config.num_trailers+4)*(S.algorithms.mpcCasadi.qref(2*S.config.num_trailers+4)-S.data.mhempc.performance.xfut{i}(2*S.config.num_trailers+4,j )).^2 + ...
                                   S.algorithms.mpcCasadi.Q(2*S.config.num_trailers+5,2*S.config.num_trailers+5)*(S.algorithms.mpcCasadi.qref(2*S.config.num_trailers+5)-S.data.mhempc.performance.xfut{i}(2*S.config.num_trailers+5,j )).^2);];
            %
            if j==1
                surf(X,Y,Z0pos+Zobs+ZNpos); zlim([0 300]);
            end
        end
        line(S.data.mhempc.performance.xfut{i}(2*S.config.num_trailers+2,:),S.data.mhempc.performance.xfut{i}(2*S.config.num_trailers+3,:),Z0posL,'color','y');
        line(S.data.mhempc.performance.xfut{i}(2*S.config.num_trailers+4,:),S.data.mhempc.performance.xfut{i}(2*S.config.num_trailers+5,:),Z0posL,'color','r');
        drawnow;
    end
end

function S = call_init_functions(S)
    S = reserve_temp_memory(S);
    S = gen_init_conditions(S);
    S = create_algorithms(S);
    S = init_flags_and_counters(S);
end

function S = compute_tras_rot_mtx(S)
    sdpvar a11 a12 a21 a22;
    sdpvar x1bar y1bar x2bar y2bar;
    % Besides the reference point, two more are needed to obtaint the local
    % reference frame
    % AZOTEA AC3E *********************************************************
    % Coordinates of point (x, 0)
    S.ROS.local_coord_1.lat     = -33.0343966667;       % hand coded value, measurement from rtk
    S.ROS.local_coord_1.long    = -71.5920216667;   % hand coded value, measurement from rtk
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
    S.ROS.local_coord_2.lat     = -33.03416;        % hand coded value, measurement from rtk
    S.ROS.local_coord_2.long    = -71.5921883333;         % hand coded value, measurement from rtk
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
        S.data.ysim_mhempc = [S.data.ysim_mhempc, S.data.xsim_mhempc(S.config.outputs,end) + noise];
    else
        S                  = read_sensors(S);
        S                  = measurements_vector(S);
        S.data.ysim_mhempc = [S.data.ysim_mhempc, S.ROS.sensors.measurement]; % <-- MEASUREMENTS are stored in this vector in this order: [beta_1, beta_2, theta_0, x_0, y_0]
%         S                  = update_sensorData_MHE(S);
    end
    S = update_sensorData_MHE(S);
    S.exec_time.t_sensors = [S.exec_time.t_sensors, toc];
end

function flag = check_end_condition(S)
    if (size(S.data.xest_mhempc,2) == 1) || isempty(S.controller.ref)
        flag = false;
    else       
        d0totgt = norm(S.data.xest_mhempc([2*S.config.num_trailers+2,2*S.config.num_trailers+3],end)-...
            S.controller.ref([2*S.config.num_trailers+2,2*S.config.num_trailers+3]));      
        flag = (S.algorithms.mpc.last_tgt && (d0totgt < S.config.dist_last_tgt)) || sum(sum(isnan(S.data.xest_mhempc)));
        
        v1 = [S.data.xest_mhempc(2*S.config.num_trailers+2)-S.data.xest_mhempc(2*S.config.num_trailers+4), S.data.xest_mhempc(2*S.config.num_trailers+3)-S.data.xest_mhempc(2*S.config.num_trailers+5)];
        v2 = [S.data.xest_mhempc(2*S.config.num_trailers+2)-S.controller.ref(2*S.config.num_trailers+2), S.data.xest_mhempc(2*S.config.num_trailers+3)-S.controller.ref(2*S.config.num_trailers+3)];
        v1dotv2 = v1*v2';

        if (v1dotv2 > 0) && S.algorithms.mpc.last_tgt
            flag = flag || true;
        end
    end
end

function S = ROS(S)
    S = init_ROS(S);
    S = create_obj_sens(S);
    S = create_obj_vel(S);    
end

function S = init_ROS(S)
    % Init coordinates of my local 2D plane
    S = get_reference(S);
    
    if S.config.calcMtxConvCoord == true;
        S = compute_tras_rot_mtx(S);
    else
        S.ROS.Mtx = eye(2);
    end
    
    S.ROS.sensors.measurements = [];
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
    fprintf('roslaunch RTK.launch\nroslaunch vectornav vectornav.launch\nroslaunch vectornav vectornav2.launch\nroslaunch imu_3dm_gx4 imu.launch\nroslaunch Encoders.launch\nroslaunch husky_base base.launch\npython speedholder.py\nroslaunch VP12 lidar.launch\n')
    aux1 = input('Presione enter...\n');
end

function S = create_obj_sens(S)
    sub_rtk         = rossubscriber('/fix');
    sub_vec1        = rossubscriber('/vectornav/IMU');
    sub_vec1_vel    = rossubscriber('/imu_INS');
    sub_vec2        = rossubscriber('/vectornav2/IMU');
    sub_encoders    = rossubscriber('/enc');
    sub_micro1      = rossubscriber('/imu/pose');
    sub_lidar       = rossubscriber('/angles');
    %
    S.ROS.rtk            = sub_rtk;
    S.ROS.vectornav1     = sub_vec1;
    S.ROS.vectornav2     = sub_vec2;
    S.ROS.vectornav1_vel = sub_vec1_vel;
    S.ROS.microstrain    = sub_micro1;
    S.ROS.IncEncoders    = sub_encoders;
    S.ROS.betasFromLidar = sub_lidar;

    % COrection factor for unwraping phase on real-time
    S.ROS.pose.correction_vec1              = 0;      % This value should be substracted to every pose measurement
    S.ROS.sensors.vectornav_euler_vec1      = [];
    S.ROS.sensors.vectornav_euler_vec1_Old  = 0;     % Value for computinng the difference. If it is bigger in absolute value than pi, a correction is needed
    %
    S.ROS.pose.correction_vec2              = 0;
    S.ROS.sensors.vectornav_euler_vec2      = [];
    S.ROS.sensors.vectornav_euler_vec2_Old  = 0;
    %
    S.ROS.pose.correction_micro1            = 0;      % This value should be substracted to every pose measurement
    S.ROS.sensors.microstrain_euler_micro1    = [];
    S.ROS.sensors.microstrain_euler_micro1_Old  = 0;
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
    S.ROS.LAT0   = -33.0344116667;
    S.ROS.LON0  = -71.59211133333;
    % CANCHA DE FUTBOL ****************************************************
%     S.ROS.LAT0   = -33.03453;
%     S.ROS.LON0  = -71.594498;
    %
end

function S = read_rtk(S)
    S.ROS.sensors.gpsrtk = receive(S.ROS.rtk);
    
    if isnan(S.ROS.sensors.gpsrtk.Latitude) || isnan(S.ROS.sensors.gpsrtk.Longitude)
        S = write_ctrl(S, 0, 0);
        while isnan(S.ROS.sensors.gpsrtk.Latitude) || isnan(S.ROS.sensors.gpsrtk.Longitude)
            fprintf('NO GPS SIGNAL...')
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
    S.ROS.sensors.rtk.xN = xy_cor(1);
    S.ROS.sensors.rtk.yN = xy_cor(2);
end

function S = write_ctrl(S, w0, v0)
    tic;
    S.ROS.CTRL.msg.Linear.X    = v0;
    S.ROS.CTRL.msg.Angular.Z   = w0;
    send(S.ROS.CTRL.pub,S.ROS.CTRL.msg);
    S.exec_time.t_ctrl = [S.exec_time.t_ctrl, toc];
end

function S = read_vectornav(S) % measure theta0 and the speeds
    % pose
    S.ROS.sensors.vectornav             = receive(S.ROS.vectornav1);
    quat                                = S.ROS.sensors.vectornav.Orientation;
    S.ROS.sensors.vectornav_euler_vec1  = quat2eul([quat.X,quat.Y,quat.Z,quat.W]);    
    % Unwrap pahse from online data _______________________________________
    if (S.ROS.sensors.vectornav_euler_vec1(3)-S.ROS.sensors.vectornav_euler_vec1_Old) >= pi
        S.ROS.pose.correction_vec1 = S.ROS.pose.correction_vec1 + 2*pi;
    elseif (S.ROS.sensors.vectornav_euler_vec1(3)-S.ROS.sensors.vectornav_euler_vec1_Old) <= -pi
        S.ROS.pose.correction_vec1 = S.ROS.pose.correction_vec1 - 2*pi;         
    end         
    %
    S.ROS.sensors.vectornav_euler_vec1_Old = S.ROS.sensors.vectornav_euler_vec1(3);    
    % No compute the attitude angle in my reference frama -----------------
    S.ROS.sensors.vectornav_theta0 = -S.ROS.sensors.vectornav_euler_vec1(3) + S.ROS.IMU_ZEROING_VEC1 + S.ROS.pose.correction_vec1;
    % Measure the speed ---------------------------------------------------
    S.ROS.sensors.vectornav_vel     = receive(S.ROS.vectornav1_vel);
    S.ROS.sensors.vectornav_NedVelX = S.ROS.sensors.vectornav_vel.Data(1);
    S.ROS.sensors.vectornav_NedVelY = S.ROS.sensors.vectornav_vel.Data(2);
    S.ROS.sensors.vectornav_u0      = sqrt(S.ROS.sensors.vectornav_NedVelX^2 + S.ROS.sensors.vectornav_NedVelY^2);
    S.ROS.sensors.vectornav_w0      = -1 * S.ROS.sensors.vectornav.AngularVelocity.Z;
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
    S = read_vectornav2(S);
    S = read_microstrain(S);
    S = read_encoders(S);
    S = read_rtk(S);
end

function S = measurements_vector(S)
    S.ROS.sensors.measurement   = [S.ROS.sensors.encoders.beta1; S.ROS.sensors.encoders.beta2; S.ROS.sensors.vectornav_theta0; S.ROS.sensors.rtk.xN; S.ROS.sensors.rtk.yN; S.ROS.sensors.vectornav_w0; S.ROS.sensors.vectornav_u0];
    %
    S.ROS.sensors.measurements  = [S.ROS.sensors.measurements, S.ROS.sensors.measurement];   
    % Store the attitude of each trailer
    S.data.theta0      = [S.data.theta0, S.ROS.sensors.vectornav_theta0];
    S.data.theta1      = [S.data.theta1, S.ROS.sensors.microstrain_theta1];
    S.data.theta2      = [S.data.theta2, S.ROS.sensors.vectornav_theta2];
end

function S = call_MHE(S)
    tic;        
    solve(S.algorithms.mheCasadi);
    S.data.xest_mhempc = [S.data.xest_mhempc, S.algorithms.mheCasadi.q_k];
% S.data.xest_mhempc = [S.data.xest_mhempc, S.data.xsim_mhempc(:,end) ];
    %
    S.exec_time.t_mhe               = [S.exec_time.t_mhe, toc];
end

function S = call_PFA(S)
    tic;
%     S                              = PF_observer(S,S.data.xest_mhempc(:,end),S.path.ind_tgt_mhempc,S.algorithms.mpc.last_tgt);
    q_k             = S.data.xest_mhempc(:,end);
    pos0_k          = q_k(2*S.config.num_trailers+2:2*S.config.num_trailers+3);
    % =====================================================================
    % Search for obstacles in the neighbourhood and pass coordinates in
    % case of finding any
    % NEED TO BE REIMPLEMENTED FOR ACADO SOLVERS!!!! **********************
    if isempty(S.path.obstacles)
        setXYRobstacle(S.algorithms.mpcCasadi, [1e6, 1e6, 1]); % no obstacle found
    else
        dis             = pdist2(S.path.obstacles(:,1:2),pos0_k');
        [disort, indx]  = sort(dis);
        nroObs          = numel(indx);
        if nroObs<3
            setXYRobstacle(S.algorithms.mpcCasadi, [S.path.obstacles(indx,:);repmat([1e6, 1e6, 1],3-nroObs,1)])
        else
            setXYRobstacle(S.algorithms.mpcCasadi, S.path.obstacles(indx(1:3),:));
        end
    end
    %
    indx                = min(floor(S.path.counter),length(S.path.ref));
    S.path.counter      = S.path.counter+1;
    S.controller.ref    = S.path.ref(:,indx);
    %
    if indx >= length(S.path.ref)
        S.path.last_tgt             = true;
        S.algorithms.mpc.last_tgt   = 1;
    end
    %
    S.exec_time.t_pfa = [S.exec_time.t_pfa, toc];
end

function S = call_MPC(S)
    tic;
    % Solve control problem _______________________________________________
    setq0(S.algorithms.mpcCasadi,S.data.xest_mhempc(:,end));

%     setQref(S.algorithms.mpcCasadi,[S.controller.ref; S.config.wNr; S.config.vNr]);
%     S.controller.ref = S.path.ref(:,min(floor(S.config.iters+5),length(S.path.ref)));
    if length(S.controller.ref) < 2*S.config.num_trailers+7
        S.path.references_mhempc    = [S.path.references_mhempc, S.controller.ref];
        S.controller.ref            = [S.controller.ref; S.config.wNr; S.config.vNr];
    else
        S.path.references_mhempc    = [S.path.references_mhempc, S.controller.ref(1:end-2)];
    end
    setQref(S.algorithms.mpcCasadi, S.controller.ref);

    S.path.ind_tgt_mhempc          = S.path.indx_alg;
    S.path.indices_tgt_mhempc      = [S.path.indices_tgt_mhempc, reshape(S.path.ind_tgt_mhempc,2,1)];
    
%     S.path.references_mhempc       = [S.path.references_mhempc, S.controller.ref];    

    solve(S.algorithms.mpcCasadi);

    S.algorithms.mpc.controls_MPC  = [S.algorithms.mpc.controls_MPC, S.algorithms.mpcCasadi.u_k];
    % To the last estimated velocties, add the last acceleration computed _
    S.algorithms.mpc.controls_MPC_w0 = [S.algorithms.mpc.controls_MPC_w0, S.data.xest_mhempc(2*S.config.num_trailers+6,end)+S.algorithms.mpc.controls_MPC(1,end)];
    S.algorithms.mpc.controls_MPC_v0 = [S.algorithms.mpc.controls_MPC_v0, S.data.xest_mhempc(2*S.config.num_trailers+7,end)+S.algorithms.mpc.controls_MPC(2,end)];
    % Velocities to apply to the Husky ____________________________________
    S.algorithms.mpc.Husky_w0       = S.algorithms.mpc.controls_MPC_w0(end);
    S.algorithms.mpc.Husky_v0       = S.algorithms.mpc.controls_MPC_v0(end);
    %
    S.exec_time.t_mpc              = [S.exec_time.t_mpc, toc];
end

function S = update_sensorData_MHE(S)
    updateMeasurement(S.algorithms.mheCasadi, double(S.data.ysim_mhempc(:,end)));
    updateInput(S.algorithms.mheCasadi, double(S.algorithms.mpc.controls_MPC(:,end)));
end

function S = gen_path(S,type)
    if strcmp(type,'infinity')
        a = 5;
        t = 0:0.01:2*pi;
        x = (a*sqrt(2).*cos(t))./(sin(t).^2+1);
        y = (a*sqrt(2).*cos(t).*sin(t))./(sin(t).^2+1);
        
        d = -max(S.system.Lh(:,2))*S.config.num_trailers*log(tan(0.1/2));
        y1 = y(1)-d:0.1:y(1);
        x1 = x(1).*ones(size(y1));
        
        S.path.coordinates = [x1,x;y1,y];        
    elseif strcmp(type,'test')
        % TEST TRAJECTORY #################################################
%         t = 0:3500;
%         x = 3*sin(2*pi.*t/t(end)) + 10;
%         y = 6*cos(2*pi.*t/t(end)) + 10;
%         S.path.coordinates = [x;y];
%         %
        S.config.N  = (S.config.tf-S.config.t0)/S.config.Ts;
        t           = linspace(0,2*pi,S.config.N);
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

        d = -max(S.system.Lh(:,2))*S.config.num_trailers*log(tan(0.1/2));
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
           
            d = -max(S.system.Lh(:,2))*S.config.num_trailers*log(tan(0.1/2));
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
           y = 14:0.1:23.5;
           x = repmat(8,1,length(y));
           
           t = 0:0.01:1;
           x = [x, 4.5 + 3.5.*cos(2*pi*t)];
           y = [y, 23.5 + 4.5.*sin(2*pi*t)];
           
           x = x-0.5;
           y = y+2;

           S.path.coordinates = [x; y];
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
    else    
        % REAL TRAJECTORY #####################################################
        % First segment
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
%         trajX = x9-S.path.deltaL:-S.path.deltaL:S.init_condition.x0(2*S.config.num_trailers+2)-2*(1+S.config.num_trailers)*(S.system.Lh1+S.system.L1);
%         trajY = y9.*ones(size(trajX));
%         S.path.coordinates = [S.path.coordinates, [trajX; trajY]];
%         %
%         x10 = S.path.coordinates(1,end);
%         y10 = S.path.coordinates(2,end);    
%         trajY = y10-S.path.deltaL:-S.path.deltaL:S.init_condition.x0(2*S.config.num_trailers+3);
%         trajX = x10.*ones(size(trajY));
%         S.path.coordinates = [S.path.coordinates, [trajX; trajY]];
%         %
%         x11 = S.path.coordinates(1,end);
%         y12 = S.path.coordinates(2,end);
%         trajX = x11+S.path.deltaL:S.path.deltaL:S.init_condition.x0(2*S.config.num_trailers+2);
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

function [effort, max_effort, effort_vec] = control_effort(S, U)
    l           = length(U);
    effort      = 0;
    max_effort  = 0;
    effort_vec  = zeros(1,l);    
    for i=1:l
        effort_step     = U(1,i)^2 + U(2,i)^2;
        effort_vec(i)   = effort_step;
        if effort_step  > max_effort
            max_effort  = effort_step;
        end
        effort      = effort + effort_step;
    end
    effort = effort/l;
end

function est_err = compute_estimation_error(S, x, xhat)
    err = [];
    l = min(length(x), length(xhat));
    rows = min(numel(x(:,1)), numel(xhat(:,1)));
    for i=1:l
        err = [err, norm(x(1:rows,i)-xhat(1:rows,i))];
    end
    %
    est_err = sum(err)/l;
end

function S = gen_dynamic_and_solvers(S)
    % System's dimension __________________________________________________        
    S.system.nu      = 2;                                                              % input dimension: [w0,v0]
    S.system.ny      = S.config.num_trailers + 3 + S.system.nu;      % output dimension: [beta_i,theta_0,xy_0]
    S.system.nv      = S.system.ny;  
    S.system.nq      = 2*S.config.num_trailers + 7;
    % System's parameters -------------------------------------------------
    S.system.Lh1     = -12;%0.3;%0.342;
    S.system.L1      = 0.75;%0.38;
    S.system.Lh2     = -12;%0.3;%0;
    S.system.L2      = 0.85;%1.08;
    % ---------------------------------------------------------------------
    S.system.Lh3     = -12;%0.3;%0;
    S.system.L3      = 1.05;%0.78;%0.229;
    S.system.Lh4     = 0;%0.3;%0.342;%0;
    S.system.L4      = 0.25;%0.78;%0.229;
    S.system.Lh5     = 0;%0.3;%-0.342;%0.048;
    S.system.L5      = 0.25;%0.78;%0.229;
    S.system.Lh6     = 0;%0.3;%-0.342;%0.048;
    S.system.L6      = 0.25;%0.78;%0.229;
    S.system.Lh7     = 0;%0.3;%-0.342;%0.048;
    S.system.L7      = 0.25;%0.78;%0.229;
    S.system.Lh8     = 0;%0.3;%0.342;%0.048;
    S.system.L8      = 0.25;%0.78;%0.229;
    S.system.Lh9     = 0;%0.3;%0.342;%0.048;
    S.system.L9      = 0.25;%0.78;%0.229;
    S.system.Lh10    = 0;%0.3;%0.048;
    S.system.L10     = 0.25;%0.78;%0.229;
    %
    S.system.Lhi     = [S.system.Lh1;S.system.Lh2;S.system.Lh3;S.system.Lh4;S.system.Lh5;S.system.Lh6;S.system.Lh7;S.system.Lh8;S.system.Lh9;S.system.Lh10]; 
    S.system.Li      = [S.system.L1;S.system.L2;S.system.L3;S.system.L4;S.system.L5;S.system.L6;S.system.L7;S.system.L8;S.system.L9;S.system.L10]; 
    % Casadi variabes ----------------------------------------------------
    S.dynamic.q      = casadi.MX.sym('q',S.system.nq);
    S.dynamic.u      = casadi.MX.sym('u',S.system.nu);                  % accelerations (angular and linear)
%     S.dynamic.p      = casadi.MX.sym('beta_0',S.config.num_trailers);   % initial joint-angle values
    
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
    %
    S.system.r       = 0.165;                                                           % wheel radius
    S.system.b       = 0.55;%0.229;                                                           % vehicles's width    
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
    for i=1:S.config.num_trailers
        S.system.XYtrailerLongAxe = [S.system.XYtrailerLongAxe; [0 0; 0 S.system.Li(i)]];
        S.system.XYtrailerLoad = [S.system.XYtrailerLoad; [-S.system.b/2*0.9 -S.system.b/2*0.9 S.system.b/2*0.9 S.system.b/2*0.9 -S.system.b/2*0.9; -S.system.r/2*1.2 S.system.Li(i)*0.65 S.system.Li(i)*0.65 -S.system.r/2*1.2 -S.system.r/2*1.2]];
    end
    % System's variables  -------------------------------------------------
    % Since ACADO does not allows to incorporate directly rate of change
    % constraints, I've add the controls w0 and v0 as additional states and
    % define as controls the variaton of w0 and v0 and add constraints on
    % these new inputs. Proper differential equations were added too.
    
    Control delta_w0 delta_v0;
    %
    switch(S.config.num_trailers)
        case 1            
            % CASADI FORMULATION ******************************************
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1, sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)), cos(S.dynamic.q(1))];
            S.system.J       = [S.system.J1];
            %
            S.dynamic.f_rhs  = [ [1 0]*(eye(2)-S.system.J1) * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                 %
                                 [1 0] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                 [1 0] * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                 %
                                 [0 cos(S.dynamic.q(2))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                 [0 sin(S.dynamic.q(2))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                 [0 cos(S.dynamic.q(3))] * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                 [0 sin(S.dynamic.q(3))] * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                 %
                                 S.dynamic.u(1);...
                                 S.dynamic.u(2) ];
            S.system.Lh      = [S.system.Lh1_u,S.system.L1_u];
            S.system.LLh     = sum(S.system.Lh(:,1)+S.system.Lh(:,2));
        case 2
            % CASADI FORMULATION ******************************************
            S.system.J2      = [-S.system.Lh2 * cos(S.dynamic.q(2)) / S.system.L2, sin(S.dynamic.q(2))/S.system.L2; S.system.Lh2*sin(S.dynamic.q(2)), cos(S.dynamic.q(2))];            
            S.system.J1      = [-S.system.Lh1 * cos(S.dynamic.q(1)) / S.system.L1, sin(S.dynamic.q(1))/S.system.L1; S.system.Lh1*sin(S.dynamic.q(1)), cos(S.dynamic.q(1))];            
            S.system.J       = [S.system.J2, S.system.J1];
            f1               = simplify([1 0] * (eye(2)-S.system.J1) * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7));
            f2               = simplify([1 0] * (eye(2)-S.system.J2) * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7));
            f3               = simplify([1 0] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7));
            f4               = simplify([1 0] * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7));
            f5               = simplify([1 0] * S.system.J2 * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7));
            f6               = simplify([0 cos(S.dynamic.q(3))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7));
            f7               = simplify([0 sin(S.dynamic.q(3))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7));
            f8               = simplify([0 cos(S.dynamic.q(2*S.config.num_trailers+1))] * S.system.J2 * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7));
            f9               = simplify([0 sin(S.dynamic.q(2*S.config.num_trailers+1))] * S.system.J2 * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7));

            S.dynamic.f_rhs  = [ f1;...
                                 f2;...
                                 %
                                 f3;...
                                 f4;...
                                 f5;...
                                 %
                                 f6;...
                                 f7;...
                                 f8;...
                                 f9;...
                                 %
                                 S.dynamic.u(1);...
                                 S.dynamic.u(2) ];
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
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [1 0] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [0 cos(S.dynamic.q(4))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(4))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...S.dynamic.q
                                   [0 cos(S.dynamic.q(2*S.config.num_trailers+1))]*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(2*S.config.num_trailers+1))]*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   S.dynamic.u(1) / 1;
                                   S.dynamic.u(2) / 1 ];
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
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [1 0] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [0 cos(S.dynamic.q(5))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(5))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 cos(S.dynamic.q(2*S.config.num_trailers+1))] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(2*S.config.num_trailers+1))] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   S.dynamic.u(1) / 1;
                                   S.dynamic.u(2) / 1 ];
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
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [1 0] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [0 cos(S.dynamic.q(6))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(6))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 cos(S.dynamic.q(2*S.config.num_trailers+1))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(2*S.config.num_trailers+1))]*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   S.dynamic.u(1) / 1;
                                   S.dynamic.u(2) / 1 ];
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
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [1 0] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [0 cos(S.dynamic.q(7))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(7))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 cos(S.dynamic.q(2*S.config.num_trailers+1))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(2*S.config.num_trailers+1))]*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   S.dynamic.u(1) / 1;
                                   S.dynamic.u(2) / 1 ];
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
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J7)*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [1 0] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [0 cos(S.dynamic.q(8))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(8))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 cos(S.dynamic.q(2*S.config.num_trailers+1))]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(2*S.config.num_trailers+1))]*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   S.dynamic.u(1) / 1;
                                   S.dynamic.u(2) / 1 ];
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
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J7)*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J8)*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [1 0] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [0 cos(S.dynamic.q(9))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(9))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 cos(S.dynamic.q(2*S.config.num_trailers+1))]*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(2*S.config.num_trailers+1))]*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   S.dynamic.u(1) / 1;
                                   S.dynamic.u(2) / 1 ];
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
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J7)*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J8)*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J9)*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [1 0] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [0 cos(S.dynamic.q(10))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(10))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 cos(S.dynamic.q(2*S.config.num_trailers+1))]*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(2*S.config.num_trailers+1))]*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   S.dynamic.u(1) / 1;
                                   S.dynamic.u(2) / 1 ];
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
            S.dynamic.f_rhs  = [   [1 0]*(eye(2)-S.system.J1) * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J2)*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J3)*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J4)*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J5)*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J6)*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J7)*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J8)*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J9)*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0]*(eye(2)-S.system.J10)*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [1 0] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [1 0] * S.system.J10*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   [0 cos(S.dynamic.q(11))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(11))] * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 cos(S.dynamic.q(2*S.config.num_trailers+1))] * S.system.J10*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   [0 sin(S.dynamic.q(2*S.config.num_trailers+1))] * S.system.J10*S.system.J9*S.system.J8*S.system.J7*S.system.J6*S.system.J5*S.system.J4*S.system.J3*S.system.J2*S.system.J1 * S.dynamic.q(2*S.config.num_trailers+6:2*S.config.num_trailers+7);...
                                   %
                                   S.dynamic.u(1) / 1;
                                   S.dynamic.u(2) / 1 ];
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
    %
    % RK4 -----------------------------------------------------------------
    k1              = S.dynamic.f(S.dynamic.q, S.dynamic.u);
    k2              = S.dynamic.f(S.dynamic.q + S.config.Ts / 2 * k1, S.dynamic.u);
    k3              = S.dynamic.f(S.dynamic.q + S.config.Ts / 2 * k2, S.dynamic.u);
    k4              = S.dynamic.f(S.dynamic.q + S.config.Ts * k3, S.dynamic.u);
    x_rk4           = S.dynamic.q + S.config.Ts / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    S.dynamic.FNt   = casadi.Function('FNt', {S.dynamic.q, S.dynamic.u}, {x_rk4});
    % Output of the system ------------------------------------------------
    S.dynamic.h_rhs  = S.dynamic.q(S.config.outputs);
    S.dynamic.h      = casadi.Function('h', {S.dynamic.q}, {S.dynamic.h_rhs});
    %
    x0 = S.dynamic.q(2*S.config.num_trailers+2);
    y0 = S.dynamic.q(2*S.config.num_trailers+3);
    xN = S.dynamic.q(2*S.config.num_trailers+4);
    yN = S.dynamic.q(2*S.config.num_trailers+5);
end

function S = updateAC(S,q,v)
    S.algorithms.mhe.Mk     = exp( -((q' * S.algorithms.mhe.P * q) * ( (v * v' ) / S.algorithms.mhe.B )) -1e-3);
    S.algorithms.mhe.alphak = 1 - S.algorithms.mhe.Mk;
    l                       = real(eig(S.algorithms.mhe.P));
    L                       = max(l);
    S.algorithms.mhe.Wk     = S.algorithms.mhe.P * ( eye(S.system.nq+2*S.config.num_trailers) - ( (q * q') * S.algorithms.mhe.P) ./ (1 + q' * q * L) );
    traceWoverAlphak        = trace(S.algorithms.mhe.Wk) / S.algorithms.mhe.alphak;

    if traceWoverAlphak >= S.algorithms.mhe.A
        S.algorithms.mhe.P = S.algorithms.mhe.Wk .* (S.algorithms.mhe.alphak)^(1/S.system.nq);
    else
        S.algorithms.mhe.P = S.algorithms.mhe.Wk;
    end
    
end

function S = init_mheCasadi(S)
    S.box_constraints               = struct;
    S.box_constraints.QluBounds     = [];
    S.box_constraints.WluBounds     = [];
    S.box_constraints.VluBounds     = [];
    S.box_constraints.ZluBounds     = [];
    S.box_constraints.UluBounds     = [];
    %
    S.Mtxs                          = struct;
    
    if S.config.SIM == false % G2T
        S.Mtxs.Q    = eye(2*S.config.num_trailers+7);
        S.Mtxs.R    = diag([0.6.*ones(1,S.config.num_trailers) 1.5 5 5 0.5 0.5]);
        S.Mtxs.P    = 1e7.*eye(S.system.nq);
        S.Mtxs.Z    = 100.*eye(S.config.num_trailers);
        S.Mtxs.dU   = 100.*eye(S.system.nu);
    else
        S.Mtxs.Q    = diag([ones(1,S.config.num_trailers),1,ones(1,S.config.num_trailers),1,1,1,1,1,1]);
        S.Mtxs.R    = diag([ones(1,S.config.num_trailers),2,4,4,1,1]);
        S.Mtxs.P    = 1e6.*eye(S.system.nq);
        S.Mtxs.Z    = 1e3.*eye(S.config.num_trailers);
        S.Mtxs.dU   = eye(S.system.nu);
    end
    % nlmheCasadiNt(Ne,x,u,Nt,f_rhs,h_rhs,Mtxs,nq,nu,ny,boxConst,q0bar,Ts,dimensions)
    S.algorithms.mheCasadi = nlmheCasadiNt(S.config.Ne,S.dynamic.q,S.dynamic.u,S.config.num_trailers,S.dynamic.f_rhs,S.dynamic.h_rhs,...
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

function S = compute_Lambda(S)
    S.algorithms.mhe.RAS.lrho = min(eig(S.algorithms.mhe.input_MHE.SAC));
    S.algorithms.mhe.RAS.Lrho = max(eig(S.algorithms.mhe.input_MHE.SAC));
    S.algorithms.mhe.RAS.lvarphi = min(eig(S.algorithms.mhe.input_MHE.W));
    S.algorithms.mhe.RAS.Lvarphi = max(eig(S.algorithms.mhe.input_MHE.W));
    %
    S.algorithms.mhe.RAS.Lambda = (sqrt(S.algorithms.mhe.RAS.lrho*S.algorithms.mhe.RAS.lvarphi)+sqrt(S.algorithms.mhe.RAS.Lrho*S.algorithms.mhe.RAS.lvarphi)+sqrt(S.algorithms.mhe.RAS.Lrho*S.algorithms.mhe.RAS.lrho))/sqrt(S.algorithms.mhe.RAS.lrho*S.algorithms.mhe.RAS.lvarphi);
    %
    sdpvar lrho Lrho lvarphi Lvarphi;
    obj = (lrho*lvarphi + Lrho*lvarphi + Lrho*lrho)/(lrho*lvarphi);
    con = [lrho >= eps; Lrho >= eps; lvarphi >= eps; Lvarphi >= eps; Lvarphi>=lvarphi; Lrho>=lrho];
    optimize(con,obj);
    %
    S.algorithms.mhe.RAS.lrho_sym = value(lrho);
    S.algorithms.mhe.RAS.Lrho_sym = value(Lrho);
    S.algorithms.mhe.RAS.lvarphi_sym = value(lvarphi);
    S.algorithms.mhe.RAS.Lvarphi_sym = value(Lvarphi);
    S.algorithms.mhe.RAS.Lambda_sym = (sqrt(S.algorithms.mhe.RAS.lrho_sym*S.algorithms.mhe.RAS.lvarphi_sym)+sqrt(S.algorithms.mhe.RAS.Lrho_sym*S.algorithms.mhe.RAS.lvarphi_sym)+sqrt(S.algorithms.mhe.RAS.Lrho_sym*S.algorithms.mhe.RAS.lrho_sym))/sqrt(S.algorithms.mhe.RAS.lrho_sym*S.algorithms.mhe.RAS.lvarphi_sym);
end

function S = init_MHE_AC(S)
    S.algorithms.mhe.B                = 1000.5;
    S.algorithms.mhe.C                = 1e6;
    S.algorithms.mhe.A                = S.algorithms.mhe.C;
    S.algorithms.mhe.P                = S.algorithms.mhe.C.*eye(4*S.config.num_trailers+7);
    S.algorithms.mhe.P(1:S.config.num_trailers,1:S.config.num_trailers) = S.algorithms.mhe.P(1:S.config.num_trailers,1:S.config.num_trailers) .* 10;
    %
    S.algorithms.mhe.input_MHE.xAC    = S.init_condition.x0bar;
    S.algorithms.mhe.input_MHE.SAC    = inv(S.algorithms.mhe.P);
    S.algorithms.mhe.input_MHE.WL     = chol(S.algorithms.mhe.input_MHE.SAC);
    S.algorithms.mhe.countAC          = 1;
    S.algorithms.mhe.MAXcountAC       = 10;
end

function S = init_mpc(S)
    % PARAMETERS MPC
    Vref                               = [0; 0];
    
    S.box_constraints               = struct;    
    S.box_constraints.QluBounds     = [repmat([-pi/2 pi/2],S.config.num_trailers,1); repmat([-inf inf],S.config.num_trailers+1,1); [-inf inf;-inf inf; -inf inf; -inf inf]; [-1 1; -0.75 0.75]];
    S.box_constraints.QNluBounds    = [repmat([-pi/2 pi/2],S.config.num_trailers,1); repmat([-inf inf],S.config.num_trailers+1,1); [-inf inf;-inf inf; -inf inf; -inf inf]; [-1 1; -0.75 0.75]];
    S.box_constraints.UluBounds     = [-0.25 0.25;...
                                       -0.15 0.15];

    %
    S.Mtxs                          = struct;
    
    if S.config.SIM == false % G2T
        S.Mtxs.Q    = eye(2*S.config.num_trailers+7);
        S.Mtxs.QN   = eye(2*S.config.num_trailers+7);
        S.Mtxs.R    = eye(2);
    else
        if ~S.config.reference
            S.Mtxs.Q    = diag([0.01.*ones(1,S.config.num_trailers),0.0001.*ones(1,S.config.num_trailers+1),[0.5 0.5 0.00002 0.00002], 0.01 10]);
            S.Mtxs.R    = 0.1.*eye(2);
        else
            S.Mtxs.Q    = diag([0.1.*ones(1,S.config.num_trailers),0.1.*ones(1,S.config.num_trailers+1),[0.5 0.5 2 2], 0.5 10]);
            S.Mtxs.R    = 0.1.*eye(2);
        end
        S.Mtxs.QN   = S.Mtxs.Q;        
        S.Mtxs.QrefOpt  = 1e6.*eye(S.system.nq);
        S.Mtxs.UrefOpt  = 1e3.*eye(S.system.nu);
    end  

    % nlmpcCasadi(Nc,x,u,Nt,f_rhs,Ji,Mtxs,nq,nu,boxConst,Ts)
    S.algorithms.mpcCasadi = nlmpcCasadiNt(S.config.Nc,S.dynamic.q,S.dynamic.u,S.config.num_trailers,S.dynamic.f_rhs,S.Mtxs,S.system.nq,S.system.nu,S.box_constraints,S.config.Ts);
    %
    setUref(S.algorithms.mpcCasadi,Vref);
    setQref(S.algorithms.mpcCasadi,[[]; S.config.wNr; S.config.vNr]);
    %
    S.algorithms.mpc.last_tgt          = 0;
    S.algorithms.mpc.controls_MPC      = zeros(S.system.nu,1);
    %
    S.algorithms.mpc.controls_MPC_w0   = 0;
    S.algorithms.mpc.controls_MPC_v0   = 0;
    %
    S.algorithms.mpc.Husky_w0          = 0;
    S.algorithms.mpc.Husky_v0          = 0;
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
    S.path.path       = 'infinity';
    S.path.name       = [S.path.path,'-N=',num2str(S.config.num_trailers),'-Ne=',num2str(S.config.Ne),'-Nc=',num2str(S.config.Nc),'-Ts=',num2str(S.config.Ts),'.mat'];
    S                 = gen_path(S,S.path.path);
    if exist(S.path.name,'file')
        S.config.reference  = true;
        load(S.path.name);
        ref(end,:)          = S.config.vNr;
        S.path.ref          = ref;
    else
        S.config.reference  = false;        
        len                 = length(S.path.coordinates);
        S.path.ref          = [zeros(S.config.num_trailers,len); zeros(S.config.num_trailers+1,len); S.path.coordinates; S.path.coordinates; repmat(S.config.wNr,1,len); repmat(S.config.vNr,1,len)];
    end
    S.path.obstacles  = [];                                           % No obstacles               
    S.path.counter    = 1;
    %
    if strcmp(S.path.path,'infinity')

    elseif strcmp(S.path.path,'test')

    elseif strcmp(S.path.path,'monaco')

    elseif strcmp(S.path.path,'test_betas')

    elseif strcmp(S.path.path,'rectangular_terrace')        
        if ~S.config.reference
            S = gen_obstacle(S,1e6, 1e6, 0.25);
        else
%             S = gen_obstacle(S,1e6, 1e6, 0.25);
            S = gen_obstacle(S,6, 31.75, 0.25);
            S = gen_obstacle(S,0.8,30.75, 0.25);
            S = gen_obstacle(S,0.8,25.75, 0.25);
        end
    elseif strcmp(S.path.path,'complex1_terrace')
        
    elseif strcmp(S.path.path,'circular_terrace')           

    elseif strcmp(S.path.path,'ellipsoidal1_terrace')           

    elseif strcmp(S.path.path,'ellipsoidal2_terrace')           

    else  
    end   
%     %        
%     S                 = gen_xend(S);                                % come back to starting point
%     S                 = gen_constraints(S);

end

function S = gen_init_conditions(S)
    S                 = gen_x0(S); 
    S                 = gen_x0bar(S);
    S                 = init_mheCasadi(S);
    S                 = fill_mhe(S);
    S                 = init_mpc(S);
end

function S = create_algorithms(S)
    % Init MHE and MPC ****************************************************
%     S = init_mhe(S);
%     S = init_mpc(S);
    
    %
    % Create an instance of the Michalek's 2017 controller ****************
%     S                                      = gen_Michalek2017_robust_nSNT(S);
%     S                                      = gen_EKF_michalekS(S);
%     set_x0bar(S.algorithms.michalek.ekf, S.init_condition.x0bar);
    % Create an instance of the Keymasi's 2019 controller *****************
%     S                                      = gen_Keymasi2019_nSNT(S);
%     S                                      = gen_EKF_keymasi(S);
%     set_x0bar(S.algorithms.keymasi.ekf, S.init_condition.x0bar);
    % Create an instance of the Morales's 2013 controller *****************
%     S                                      = gen_Morales2013_nSNT(S);
%     S                                      = gen_EKF_morales(S);
%     set_x0bar(S.algorithms.morales.ekf, S.init_condition.x0bar);
    % Create an instance of the de Wit's 1996 controller ******************
%     S                                      = gen_deWit1996_nSNT(S);
%     S                                      = gen_EKF_dewit(S);
%     set_x0bar(S.algorithms.dewit.ekf, S.init_condition.x0bar);
    % Create an instance of the Morin's 2013 controller *******************
%     S                                      = gen_Morin2008_nSNT(S);
%     S                                      = gen_EKF_morin(S);
%     set_x0bar(S.algorithms.morin.ekf, S.init_condition.x0bar);
end

function S = compute_init_reference(S)
    new_pos_tgtXYN      = S.path.coordinates(:,2);
    dis                 = 0;
    for i=1:length(S.path.coordinates)-2
        dis = dis + norm(S.path.coordinates(:,2+i)-S.path.coordinates(:,2+i-1));
        if dis >= S.system.LLh
            break;
        end
    end
    dxN                 = S.path.coordinates(1,2)-S.path.coordinates(1,1);
    dyN                 = S.path.coordinates(2,2)-S.path.coordinates(2,1);        
    new_angle_tgtXYN    = atan2c(dxN, dyN, 0);
    %
    new_pos_tgtXY0      = S.path.coordinates(:,2+i);
    dx0                 = S.path.coordinates(1,2+i)-S.path.coordinates(1,2+i-1);
    dy0                 = S.path.coordinates(2,2+i)-S.path.coordinates(2,2+i-1);        
    new_angle_tgtXY0    = atan2c(dx0, dy0, 0);    
    %
    beta_i              = (new_angle_tgtXY0-new_angle_tgtXYN)/S.config.num_trailers;
    
    S.path.tgt          = [repmat(beta_i,S.config.num_trailers,1); repmat(new_angle_tgtXY0,S.config.num_trailers+1,1)-(0:S.config.num_trailers)'.*beta_i;new_pos_tgtXY0; new_pos_tgtXYN];   
    S.path.ind_tgt_mhempc = [2+i; 2];
    % Compute intial coeffcients of the rect ______________________________
    ux0                  = S.path.coordinates(1,2+i) - S.path.coordinates(1,1+i);
    uy0                  = S.path.coordinates(2,2+i) - S.path.coordinates(2,1+i);
    if uy0~=0
        vx0              = 1;
        vy0              = -ux0/uy0;
    else
        vy0              = 1;
        vx0              = 0;
    end                        
    v_norm0              = norm([vx0, vy0]);
    vx0                  = vx0/v_norm0;
    vy0                  = vy0/v_norm0;
    A0                   = vy0;
    B0                   = -1*vx0;
    C0                   = S.path.coordinates(2,2+i)*vx0 - S.path.coordinates(1,2+i)*vy0;
    norm_AB0             = sqrt(A0^2 + B0^2);
    %
    uxN                  = S.path.coordinates(1,2) - S.path.coordinates(1,1);
    uyN                  = S.path.coordinates(2,2) - S.path.coordinates(2,1);    
    if uyN~=0
        vxN              = 1;
        vyN              = -uxN/uyN;
    else
        vyN              = 1;
        vxN              = 0;
    end                        
    v_normN              = norm([vxN, vyN]);
    vxN                  = vxN/v_normN;
    vyN                  = vyN/v_normN;
    AN                   = vyN;
    BN                   = -1*vxN;
    CN                   = S.path.coordinates(2,2)*vx0 - S.path.coordinates(1,2)*vy0;
    norm_ABN             = sqrt(AN^2 + BN^2);
    %
    S.path.rect.mhe_mpc.A0 = A0;
    S.path.rect.mhe_mpc.B0 = B0;
    S.path.rect.mhe_mpc.C0 = C0;
    S.path.rect.mhe_mpc.norm_AB0 = norm_AB0;
    S.path.rect.mhe_mpc.ux0 = ux0;
    S.path.rect.mhe_mpc.uy0 = uy0;
    %
    S.path.rect.mhe_mpc.AN = AN;
    S.path.rect.mhe_mpc.BN = BN;
    S.path.rect.mhe_mpc.CN = CN;
    S.path.rect.mhe_mpc.norm_ABN = norm_ABN;
    S.path.rect.mhe_mpc.uxN = uxN;
    S.path.rect.mhe_mpc.uyN = uyN;

    S.path.tg_rect.A0        = A0;
    S.path.tg_rect.B0        = B0;
    S.path.tg_rect.C0        = C0;
    S.path.tg_rect.norm_AB0  = norm_AB0;
    S.path.tg_rect.ux0       = ux0;
    S.path.tg_rect.uy0       = uy0;
    %
    S.path.tg_rect.AN        = AN;
    S.path.tg_rect.BN        = BN;
    S.path.tg_rect.CN        = CN;
    S.path.tg_rect.norm_ABN  = norm_ABN;
    S.path.tg_rect.uxN       = uxN;
    S.path.tg_rect.uyN       = uyN;
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
    S.data.xsim_mhempc            = [];
    S.data.ysim_mhempc            = [];
    S.data.Uk_mhempc              = [];    
    %
    S.data.slip                   = [];
    S.data.slipest_mhempc         = [];
    %
    S.data.meas_noise             = [];
    %
    S.data.test                   = [];
    %                        
    S.data.xest_mhempc            = [];
    %                        
    S.data.time_mhempc            = [];
    %
%     S                             = compute_init_reference(S);
    S.path.references_mhempc      = [];%S.path.tgt;
    %
    S.path.reach_end_mhempc       = false;
    %
    S.path.slipp_mhempc_correctCtrl = [0;0];
    %
    S.path.indices_tgt_mhempc     = [];
    %
    S.path.vel_tgt                = [];
    S.path.angle_tgt              = [];
    %
%     S.path.inds_tgt_mhempc        = 1;
%     S.path.inds_tgt_michalek      = 1;
%     S.path.inds_tgt_keymasi       = 1;
%     S.path.inds_tgt_morales       = 1;
%     S.path.inds_tgt_dewit         = 1;
%     S.path.inds_tgt_morin         = 1;
    % auxiliar indice shared by all algorithms
    S.path.indx_alg               = [1;1];
    S.path.uref                   = [];
    S.path.last_tgt               = 0;
    %   
    S.exec_time.t_mhe             = [];
    S.exec_time.t_pfa             = [];
    S.exec_time.t_mpc             = [];
    S.exec_time.t_sensors         = [];
    S.exec_time.t_ctrl            = 0;
    S.exec_time.t_tot             = 0;
    S.exec_time.t_acum            = 0;
    S.exec_time.t_mis             = 0;
    %
    S.debug.flags_mhe             = [];
    S.debug.flags_mpc             = [];
    S.debug.flags_index           = [];
    %
    S.init_flag.mhe               = false;
    S.init_flag.mpc               = false;
    S.init_flag.pfa               = false;
    %
    S.data.lidar.beta1            = [];
    S.data.lidar.beta2            = [];
    %
    S.data.theta0                 = [];
    S.data.theta1                 = [];
    S.data.theta2                 = [];
    % Flag to indicate which MHE should use
    S.config.switchMHE.mhe1       = true; % Start with mhe1(the one that estimates the join angles, then swith to mhe2)
    S.data.initValBetas           = zeros(S.config.num_trailers,1);
end

function S = gen_x0(S)
    % DO NOT FORGET CHECKING FEASIBILITY OF INITIAL CONDITION!!!
    Dx                  = S.path.coordinates(1,2)-S.path.coordinates(1,1);
    Dy                  = S.path.coordinates(2,2)-S.path.coordinates(2,1);
    theta0              = atan2(Dy,Dx);
    thetas              = [theta0; repmat(theta0,S.config.num_trailers,1)+(rand(S.config.num_trailers,1)-0.5)*pi/2 ];
    betas               = -diff(thetas);
%     if S.config.SIM == true;
%         S.config.beta1_0        = betas(1);
%         S.config.beta2_0        = betas(2);
%     end
    %
    xy_0                = S.path.coordinates(:,1);
    xy_0(2) = xy_0(2)-1; % start 
    %
    if S.config.num_trailers == 10
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
    elseif S.config.num_trailers == 9
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_8                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_N                = xy_8 - [S.system.Lh9*cos(thetas(9))+S.system.L9*cos(thetas(10)); S.system.Lh9*sin(thetas(9))+S.system.L9*sin(thetas(10))];
    elseif S.config.num_trailers == 8
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_N                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
    elseif S.config.num_trailers == 7
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_N                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
    elseif S.config.num_trailers == 6
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_N                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
    elseif S.config.num_trailers == 5
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_N                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
    elseif S.config.num_trailers == 4
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_N                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
    elseif S.config.num_trailers == 3
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_N                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
    elseif S.config.num_trailers == 2
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_N                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
    else
        xy_N                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
    end
    %
    S.init_condition.x0       = [ betas; thetas; xy_0; xy_N; zeros(S.system.nu,1) ];
    %
    S.data.xsim_mhempc(:,1)   = S.init_condition.x0;%[S.init_condition.x0; [1;1]];
end

function S = gen_x0bar(S)
    % DO NOT FORGET CHECKING FEASIBILITY OF INITIAL CONDITION!!!
    max_b   = inf;
    i       = 1;
%     while max_b >= pi/2
%         thetas  = (5./(5+i)).*[0.1*randn; (rand(S.config.num_trailers,1).*2-1).*(89*pi/180)+0.15.*randn(S.config.num_trailers,1)];
%         thetas  = S.init_condition.x0(S.config.num_trailers+1:2*S.config.num_trailers+1);
%         thetas  = S.init_condition.x0(S.config.num_trailers+1:2*S.config.num_trailers+1);%+(rand-0.5)*pi/4.*ones(S.config.num_trailers+1,1);
        thetas  = repmat(S.init_condition.x0(S.config.num_trailers+1),S.config.num_trailers+1,1);
%         betas   = zeros(S.config.num_trailers,1);
        betas   = -diff(thetas);
%         max_b   = max(abs(betas));
%     end    
    %
    xy_0                = S.init_condition.x0(2*S.config.num_trailers+2:2*S.config.num_trailers+3)+0.2.*randn(2,1);
    %
    if S.config.num_trailers == 10
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
    elseif S.config.num_trailers == 9
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_8                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
        xy_N                = xy_8 - [S.system.Lh9*cos(thetas(9))+S.system.L9*cos(thetas(10)); S.system.Lh9*sin(thetas(9))+S.system.L9*sin(thetas(10))];
    elseif S.config.num_trailers == 8
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_7                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
        xy_N                = xy_7 - [S.system.Lh8*cos(thetas(8))+S.system.L8*cos(thetas(9)); S.system.Lh8*sin(thetas(8))+S.system.L8*sin(thetas(9))];
    elseif S.config.num_trailers == 7
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_6                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
        xy_N                = xy_6 - [S.system.Lh7*cos(thetas(7))+S.system.L7*cos(thetas(8)); S.system.Lh7*sin(thetas(7))+S.system.L7*sin(thetas(8))];
    elseif S.config.num_trailers == 6
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_5                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
        xy_N                = xy_5 - [S.system.Lh6*cos(thetas(6))+S.system.L6*cos(thetas(7)); S.system.Lh6*sin(thetas(6))+S.system.L6*sin(thetas(7))];
    elseif S.config.num_trailers == 5
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_4                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
        xy_N                = xy_4 - [S.system.Lh5*cos(thetas(5))+S.system.L5*cos(thetas(6)); S.system.Lh5*sin(thetas(5))+S.system.L5*sin(thetas(6))];
    elseif S.config.num_trailers == 4
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_3                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
        xy_N                = xy_3 - [S.system.Lh4*cos(thetas(4))+S.system.L4*cos(thetas(5)); S.system.Lh4*sin(thetas(4))+S.system.L4*sin(thetas(5))];
    elseif S.config.num_trailers == 3
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_2                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
        xy_N                = xy_2 - [S.system.Lh3*cos(thetas(3))+S.system.L3*cos(thetas(4)); S.system.Lh3*sin(thetas(3))+S.system.L3*sin(thetas(4))];
    elseif S.config.num_trailers == 2
        xy_1                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
        xy_N                = xy_1 - [S.system.Lh2*cos(thetas(2))+S.system.L2*cos(thetas(3)); S.system.Lh2*sin(thetas(2))+S.system.L2*sin(thetas(3))];
    else
        xy_N                = xy_0 - [S.system.Lh1*cos(thetas(1))+S.system.L1*cos(thetas(2)); S.system.Lh1*sin(thetas(1))+S.system.L1*sin(thetas(2))];
    end   
    
    if S.config.SIM
%         S.init_condition.x0bar = [ betas; thetas; xy_0; xy_N; zeros(S.system.nu,1) ];
S.init_condition.x0bar = S.init_condition.x0;
    else    
        S                      = read_sensors(S);
        S                      = measurements_vector(S);
        if S.config.sensors.incremental_encoder == true
            S.init_condition.x0bar = [ S.ROS.sensors.encoders.beta1; S.ROS.sensors.encoders.beta2; repmat(S.ROS.sensors.vectornav_theta0,S.config.num_trailers+1,1);...
                                   S.ROS.sensors.rtk.xN+S.system.LLh*cos(S.ROS.sensors.vectornav_theta0); S.ROS.sensors.rtk.yN+S.system.LLh*sin(S.ROS.sensors.vectornav_theta0);...
                                   S.ROS.sensors.rtk.xN; S.ROS.sensors.rtk.yN; zeros(S.system.nu,1); zeros(2*S.config.num_trailers,1) ];
            
        else
            S.init_condition.x0bar = [ S.ROS.sensors.lidar.beta1; S.ROS.sensors.lidar.beta2; repmat(S.ROS.sensors.vectornav_theta0,S.config.num_trailers+1,1);...
                                   S.ROS.sensors.rtk.xN+S.system.LLh*cos(S.ROS.sensors.vectornav_theta0); S.ROS.sensors.rtk.yN+S.system.LLh*sin(S.ROS.sensors.vectornav_theta0);...
                                   S.ROS.sensors.rtk.xN; S.ROS.sensors.rtk.yN; zeros(S.system.nu,1); zeros(2*S.config.num_trailers,1) ];
        end
        
    end    
% ####    
S.data.xest_mhempc = S.init_condition.x0bar;
% ####
end
