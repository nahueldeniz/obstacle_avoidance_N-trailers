clear all;
%
%
traj_field_exp = 'rect'; % 'flat', 'rect', 'flinfit'

if strcmp(traj_field_exp, 'flinfit')
    num_exps = [2:3,6:8,10:17,21:23];
elseif strcmp(traj_field_exp, 'flat')
    num_exps = [1:3,5:8];
elseif strcmp(traj_field_exp, 'rect')
    num_exps = [1:21];%[2,5,7,8:12,15:19,21];%[1:21]; % exps_15_nov
%     num_exps = [1:3,6,8:9,11,13:20]; % exps_13_nov
%     num_exps = [2:8,11:12,14:18];
end

% traj_field_exp  = 'flat'; % 'flat'
% num_exps        = [1 6 7 8];

pltNomPath  = true;
pltEst      = true; % plt estimated position or measured

LAT0    = -33.034115;
LON0    = -71.592205;

devRefGlobal_0      = [];
devRefGlobal_N      = [];
devPathGlobal_0     = [];
devPathGlobal_N     = [];
devRefPathGlobal_0  = [];
devRefPathGlobal_N  = [];

figure(2); grid on; hold on;

for ii=num_exps    
    % N-trailer GPS data
    lats   = [];
    longs  = [];
    hdop        = [];
    h_accuracy  = [];
    x           = [];
    y           = [];
    gpsFileName = strcat([traj_field_exp, num2str(ii),'.sbp.json']);
    fileName    = strcat([traj_field_exp, num2str(ii),'.mat']);    
    load(fileName);
    try
        fileID  = fopen(gpsFileName, 'r');
        gpsData = fread(fileID, '*char')';
        fclose(fileID);
        lines   = strsplit(gpsData, '\n');
        computeCorrFac = true;
        for j = 1:length(lines)
            try
                jsonData        = jsondecode(lines{j});    
                a=1;
                if isfield(jsonData, 'lat') && isfield(jsonData, 'lon')
                    lats   = [lats, jsonData.lat];
                    longs  = [longs, jsonData.lon];
                    %
                    [xx, yy]    = latlon2xy(lats(end), longs(end), LAT0, LON0);
                    xy          = S.ROS.Mtx * [xx*1000; yy*1000];
                    x           = [x, xy(1)];
                    y           = [y, xy(2)];
                elseif isfield(jsonData, 'hdop')
                    hdop = [hdop, jsonData.hdop];
                elseif isfield(jsonData, 'h_accuracy')
                    if jsonData.h_accuracy ~= 0
                        h_accuracy = [h_accuracy, jsonData.h_accuracy];
                    end
                end
            catch
                % Ignorar líneas que no se pueden decodificar como JSON
%                 fprintf('Something went wrong...\n')
            end
        end
        if computeCorrFac
            xCorr = S.data.mhempc.performance.xest{1}(end-1,end) - mean(x(end-15:end));
            yCorr = S.data.mhempc.performance.xest{1}(end,end) - mean(y(end-15:end));
            %
            computeCorrFac = false;
        end
        x = x + xCorr;
        y = y + yCorr + 0.5;
        %
        if ii~=120
            figure(1)
%             patchline(x,y,'linestyle','-','edgecolor','r','linewidth',0.5,'edgealpha',0.3);
        end
    catch
        fprintf('No gps file...\n')
    end
    %
    startP  = 200;
    endP    = 200;
    x       = x(startP: end-endP);
    y       = y(startP: end-endP);
    xaux    = [];
    yaux    = [];
    dres    = 0.1;
    indxToDel = [];
    i = 1;
    while i<length(x)-50
        v1 = [x(i); y(i)];        
        for j=min(i+1,length(x)-50):length(x)-1
            v2 = [x(j); y(j)];
            if norm(v2-v1) <= dres
                indxToDel = [indxToDel, j];
            else
                i = j;
                break,
            end
        end
    end
    x(indxToDel) = [];
    y(indxToDel) = [];
    %
%     if ~isempty(x)        
%         errUp   = [];
%         errDown = [];
%         for i=2:length(x)
%             v       = [x(i)-x(i-1); y(i)-y(i-1)];
%             errUp   = [errUp, [-v(2); v(1)].*((1/100)*h_accuracy(i)/norm(v))];
%             errDown = [errDown, [v(2); -v(1)].*((1/100)*h_accuracy(i)/norm(v))];
%         end
%         xUp     = smooth(x(1:end-1)+errUp(1,:),1);
%         yUp     = smooth(y(1:end-1)+errUp(2,:),1);
%         xDown   = smooth(x(1:end-1)+errDown(1,:),1);
%         yDown   = smooth(y(1:end-1)+errDown(2,:),1);
% 
%         figure(1); hold on;
% %         plot(xUp,yUp,'r','LineWidth',2);
% %         plot(xDown,yDown,'r','LineWidth',2);
%         
%         figure(1)
%         for i=1:length(errUp)-1
% %             patch([xDown(i:i+1) xUp(i:i+1)],[yDown(i:i+1) yUp(i:i+1)],'r')
%             patch([min([xDown(i) xDown(i+1) xUp(i) xUp(i+1)]), max([xDown(i) xDown(i+1) xUp(i) xUp(i+1)]),...
%                 max([xDown(i) xDown(i+1) xUp(i) xUp(i+1)]), min([xDown(i) xDown(i+1) xUp(i) xUp(i+1)])],...
%                 [min([yDown(i) yDown(i+1) yUp(i) yUp(i+1)]), min([yDown(i) yDown(i+1) yUp(i) yUp(i+1)]),...
%                 max([yDown(i) yDown(i+1) yUp(i) yUp(i+1)]), max([yDown(i) yDown(i+1) yUp(i) yUp(i+1)])],'r','edgecolor','none')
%             clear alpha;
%             alpha(0.05)
%         end
%     end
    %
    if ii>0
        xyref_0     = [];
        xyref_1     = [];
        xyref_N     = [];
        thetaref_0  = [];
        thetaref_1  = [];
        thetaref_N  = [];
        indx_ref    = 2;
        for jj=2:size(S.data.mhempc.performance.qref,2)
            xyref_0     = [xyref_0, [S.data.mhempc.performance.qref{jj}(2*S.config.N+2,indx_ref); S.data.mhempc.performance.qref{jj}(2*S.config.N+3,indx_ref)]];
            xyref_1     = [xyref_1, [S.data.mhempc.performance.qref{jj}(2*S.config.N+4,indx_ref); S.data.mhempc.performance.qref{jj}(2*S.config.N+5,indx_ref)]];
            xyref_N     = [xyref_N, [S.data.mhempc.performance.qref{jj}(end-1,indx_ref); S.data.mhempc.performance.qref{jj}(end,indx_ref)]];
            %
            thetaref_0  = [thetaref_0, S.data.mhempc.performance.qref{jj}(S.config.N+1,indx_ref)];
            thetaref_1  = [thetaref_1, S.data.mhempc.performance.qref{jj}(S.config.N+2,indx_ref)];
            thetaref_N  = [thetaref_N, S.data.mhempc.performance.qref{jj}(S.config.N+3,indx_ref)];
        end    
        % #####################################################################
        % ERROR RESPECTO DE LA REFERENCIA
        %
        for j=1:S.config.N+1
            devRef_0       = [];
            devRef_N       = [];
            for k=1:length(xyref_0)
                if j==1
                    dist_aux    = sqrt((xyref_0(1,k)-S.data.ysim(4,:)).^2 +...
                                       (xyref_0(2,k)-S.data.ysim(5,:)).^2);
                    [val, ~]    = min(dist_aux);
                    devRef_0    = [devRef_0, val];                
                elseif j==S.config.N+1
                    dist_aux    = sqrt((xyref_N(1,k)-x).^2 +...
                                       (xyref_N(2,k)-y).^2);
                    [val, ~]    = min(dist_aux);
                    devRef_N    = [devRef_N, val]; 
                end            
            end
            if j==1
                devRefGlobal_0 = [devRefGlobal_0, mean(devRef_0)];
            elseif j==S.config.N+1
                devRefGlobal_N = [devRefGlobal_N, mean(devRef_N)];
            end
            % #####################################################################
            % ERROR RESPECTO DEL CAMNINO
            %
            devPath_0       = [];
            devPath_N       = [];
            for k=1:length(S.path.coordinates)
                if j==1
                    dist_aux    = sqrt((S.path.coordinates(1,k)-S.data.ysim(4,:)).^2 +...
                                       (S.path.coordinates(2,k)-S.data.ysim(5,:)).^2);
                    [val, ~]    = min(dist_aux);
                    devPath_0   = [devPath_0, val];                
                elseif j==S.config.N+1
                    dist_aux    = sqrt((S.path.coordinates(1,k)-x).^2 +...
                                       (S.path.coordinates(2,k)-y).^2);
                    [val, ~]    = min(dist_aux);
                    devPath_N   = [devPath_N, val]; 
                end            
            end
            if j==1
                devPathGlobal_0 = [devPathGlobal_0, mean(devPath_0)];
            elseif j==S.config.N+1
                devPathGlobal_N = [devPathGlobal_N, mean(devPath_N)];
            end
            % #####################################################################
            % DESVIACION DE LA REFERENCIA RESPECTO DEL CAMINO
            %
            devRefPath_0       = [];
            devRefPath_N       = [];
            for k=1:length(xyref_0)
                if j==1
                    dist_aux    = sqrt((xyref_0(1,k)-S.path.coordinates(1,:)).^2 +...
                                       (xyref_0(2,k)-S.path.coordinates(2,:)).^2);
                    [val, ~]    = min(dist_aux);
                    devRefPath_0 = [devRefPath_0, val];                
                elseif j==S.config.N+1
                    dist_aux    = sqrt((xyref_N(1,k)-S.path.coordinates(1,:)).^2 +...
                                       (xyref_N(2,k)-S.path.coordinates(2,:)).^2);
                    [val, ~]    = min(dist_aux);
                    devRefPath_N = [devRefPath_N, val]; 
                end            
            end
            if j==1
                devRefPathGlobal_0 = [devRefPathGlobal_0, mean(devRefPath_0)];
            elseif j==S.config.N+1
                devRefPathGlobal_N = [devRefPathGlobal_N, mean(devRefPath_N)];
            end
            %
        end  
    end
    drawnow;
end

%
boxplot([devRefGlobal_0; devRefGlobal_N; devPathGlobal_0; devPathGlobal_N; devRefPathGlobal_0; devRefPathGlobal_N]','Labels',{'Ref0','RefN','Path0','PathN','RefPath0','RefPathN'},'PlotStyle','traditional','BoxStyle','outline')

grid on;
ax              = gca;
ax.GridAlpha    = 0.4;
box on
ax.FontSize     = 12;
ax.GridAlpha    = 0.3;
ylabel({'$|err.|\,(m)$'},'interpreter','latex','fontsize',15);


% NO OLVIDAR AGREGAR EL TEXTO DE GPS ERROR

% exportgraphics(figure(1),'exp3_boxplot.eps','contenttype','vector')


%%
for ii=num_exps
    fileName    = strcat([traj_field_exp, num2str(ii),'.mat']);    
    load(fileName);
    for j=1:S.config.N+1
        figure(1)
        if j==1
            if pltEst
                patchline(S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j)*2,:),'linestyle','-','edgecolor','y','linewidth',1.5,'edgealpha',0.3);
            else
                patchline(S.data.ysim(end-1,:), S.data.ysim(end,:),'linestyle','-','edgecolor','y','linewidth',1.5,'edgealpha',0.3);
            end            
        elseif j~=S.config.N+1
            patchline(S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j)*2,:),'linestyle','-','edgecolor','b','linewidth',1.5,'edgealpha',0.3);
        else
            if pltEst
                patchline(S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j)*2,:),'linestyle','-','edgecolor','r','linewidth',1.5,'edgealpha',0.45);
            else
                patchline(x,y,'linestyle','-','edgecolor','r','linewidth',1.5,'edgealpha',0.45);
%                 patchline(S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j-1)*2+1,:), S.data.mhempc.performance.xest{1}(2*S.config.N+1+(j)*2,:),'linestyle','-.','edgecolor','r','linewidth',0.5,'edgealpha',1);
            end            
        end        
    end  
end
%%
grid on;
ax              = gca;
ax.GridAlpha    = 0.4;
box on
ax.FontSize     = 12;
ax.GridAlpha    = 0.3;
daspect([1 1 1])
xlim([2 12])
ylim([-0.25 8.5])
xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',15);
ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',15);

plot(S.path.coordinates(1,:),S.path.coordinates(2,:),'k-.','LineWidth',2);

plot_mono2(S, S.data.mhempc.performance.xest{1}(:,13))
plot_mono2(S, S.data.mhempc.performance.xest{1}(:,31))
plot_mono2(S, S.data.mhempc.performance.xest{1}(:,79))
plot_mono2(S, S.data.mhempc.performance.xest{1}(:,125))

% NO OLVIDAR AGREGAR EL TEXTO DE GPS ERROR

exportgraphics(figure(1),'exp3.eps','contenttype','vector')
%%
xyref_0 = [];
xyref_1 = [];
xyref_N = [];
thetaref_0 = [];
thetaref_1 = [];
thetaref_N = [];
for ii=2:size(S.data.mhempc.performance.qref,2)
    xyref_0 = [xyref_0, [S.data.mhempc.performance.qref{ii}(2*S.config.N+2,2); S.data.mhempc.performance.qref{ii}(2*S.config.N+3,2)]];
    xyref_1 = [xyref_1, [S.data.mhempc.performance.qref{ii}(2*S.config.N+4,2); S.data.mhempc.performance.qref{ii}(2*S.config.N+5,2)]];
    xyref_N = [xyref_N, [S.data.mhempc.performance.qref{ii}(end-1,2); S.data.mhempc.performance.qref{ii}(end,2)]];
    %
    thetaref_0 = [thetaref_0, S.data.mhempc.performance.qref{ii}(S.config.N+1,2)];
    thetaref_1 = [thetaref_1, S.data.mhempc.performance.qref{ii}(S.config.N+2,2)];
    thetaref_N = [thetaref_N, S.data.mhempc.performance.qref{ii}(S.config.N+3,2)];
end
for ii=1:size(S.data.mhempc.performance.qref,2)    
    patchline(xyref_0(1,:), xyref_0(2,:),'linestyle','-.','edgecolor','y','linewidth',2,'edgealpha',1);
    patchline(xyref_1(1,:), xyref_1(2,:),'linestyle','-.','edgecolor','b','linewidth',2,'edgealpha',1);
    patchline(xyref_N(1,:), xyref_N(2,:),'linestyle','-.','edgecolor','r','linewidth',2,'edgealpha',1);   
end


len = length(S.data.mhempc.performance.xest{1});        
% infty traj
indx1 = 1;
indx2 = ceil(0.45*(len/4));
indx3 = ceil(2.3*(len/4));
indx4 = ceil(1.8*(len/4));

plot_mono2(S, S.data.mhempc.performance.xest{iplt}(:,indx1));hold on;
plot_mono2(S, S.data.mhempc.performance.xest{iplt}(:,indx2));hold on;
plot_mono2(S, S.data.mhempc.performance.xest{iplt}(:,indx3));hold on;
plot_mono2(S, S.data.mhempc.performance.xest{iplt}(:,indx4));hold on;

plot_mono2(S, S.data.mhempc.performance.xfut{iplt,indx1}(:,1),'m');hold on;
plot_mono2(S, S.data.mhempc.performance.xfut{iplt,indx2}(:,1),'m');hold on;
plot_mono2(S, S.data.mhempc.performance.xfut{iplt,indx3}(:,1),'m');hold on;
plot_mono2(S, S.data.mhempc.performance.xfut{iplt,indx4}(:,1),'m');hold on;

plot_mono2(S, S.data.mhempc.performance.xfut{iplt,indx1}(:,end),'c');hold on;
plot_mono2(S, S.data.mhempc.performance.xfut{iplt,indx2}(:,end),'c');hold on;
plot_mono2(S, S.data.mhempc.performance.xfut{iplt,indx3}(:,end),'c');hold on;
plot_mono2(S, S.data.mhempc.performance.xfut{iplt,indx4}(:,end),'c');hold on;
%


% circ traj
indx = 1;
% for i=1:S.config.N+1
%     if i==1
%         plot(S.data.mhempc.performance.xfut{iplt,indx}(2*S.config.N+1+(i-1)*2+1,:),S.data.mhempc.performance.xfut{iplt,indx}(2*S.config.N+1+i*2,:),'y','LineWidth',5); hold on;
%     elseif i==S.config.N+1
%         plot(S.data.mhempc.performance.xfut{iplt,indx}(2*S.config.N+1+(i-1)*2+1,:),S.data.mhempc.performance.xfut{iplt,indx}(2*S.config.N+1+i*2,:),'r','LineWidth',5); hold on;
%     else    
%         plot(S.data.mhempc.performance.xfut{iplt,indx}(2*S.config.N+1+(i-1)*2+1,:),S.data.mhempc.performance.xfut{iplt,indx}(2*S.config.N+1+i*2,:),'b','LineWidth',5); hold on;
%     end
% end
% indx_aux = ceil(2.3*(len/4));
% 
% % plot_mono2(S, S.data.mhempc.performance.xsim{1}(:,indx));hold on;
% % plot_mono2(S, S.data.mhempc.performance.xsim{1}(:,indx_aux));hold on;
% 
% plot_mono2(S, S.data.mhempc.performance.xfut{1,indx}(:,1),'m');hold on;
% plot_mono2(S, S.data.mhempc.performance.xfut{1,indx_aux}(:,1),'m');hold on;
% 
% plot_mono2(S, S.data.mhempc.performance.xfut{1,indx}(:,end),'c');hold on;
% plot_mono2(S, S.data.mhempc.performance.xfut{1,indx_aux}(:,end),'c');hold on;




% exportgraphics(axes,'inftyNoAttN7.eps','contenttype','vector')
% exportgraphics(axes,'circSiAttN7-refs.eps','contenttype','vector')

% InSet = get(axes, 'TightInset');
% set(axes, 'Position', [InSet(1:2), InSet(1)-InSet(3), InSet(2)-InSet(4)])

figure; hold on; grid on;
patchline(1:length(thetaref_0),thetaref_0.*180/pi,'linestyle','-.','edgecolor','y','linewidth',2,'edgealpha',0.75);
patchline(1:length(thetaref_1),thetaref_1.*180/pi,'linestyle','-.','edgecolor','b','linewidth',2,'edgealpha',0.75);
patchline(1:length(thetaref_N),thetaref_N.*180/pi,'linestyle','-.','edgecolor','r','linewidth',2,'edgealpha',0.75);
%%

% traj = 'circ';
traj    = 'infty';
figure; hold on; grid on;
devs    = {};
dk      = 1;
mean_devs = zeros(S.config.N+1,1);

if strcmp(traj, 'circ')
    start_val           = 1;%2700;
    end_val             = 2700;%1800;
else
    start_val           = 700;%9000;%2500
    end_val             = 200;
end

avoid_val1 = 70;
avoid_val2 = 70;
[curvSorted, indx]  = sort(S.path.curvature(1:end));
for ii=1:S.config.NUM_SIMS
    if ii~=avoid_val1 && ii~=avoid_val2
        for j=[2:S.config.N,S.config.N+1,1]%1:S.config.N+1
            devRef_0       = [];
            indx_min    = [];
            curvature   = [];
            xval        = [];
            yval        = [];

            for k=1:length(S.path.coordinates)
                dist_aux    = sqrt((S.path.coordinates(1,k)-S.data.mhempc.performance.xsim{ii}(2*S.config.N+1+(j-1)*2+1,:)).^2 +...
                                   (S.path.coordinates(2,k)-S.data.mhempc.performance.xsim{ii}(2*S.config.N+1+j*2,:)).^2);
                [val, ~]    = min(dist_aux);
                devRef_0       = [devRef_0, val];                
            end
            devRef_0               = devRef_0(1:length(S.path.curvature));
            t1                  = linspace(start_val*S.config.Ts, S.config.Ts*(start_val+length(devRef_0)),length(devRef_0));
            devs{ii,j}           = devRef_0;          

            if j==1                
                if strcmp(traj, 'circ')
                    plot(t1(start_val:end-end_val), devRef_0(start_val:end-end_val),'y-'); hold on;
                else
                    plot(curvSorted(start_val:end-end_val), devRef_0(start_val:end-end_val),'y-'); hold on;
%                     plot(t1(start_val:end-end_val), dev_i(start_val:end-end_val),'y-'); hold on;
                end                
            elseif j==S.config.N+1
                if strcmp(traj, 'circ')
                    plot(t1(start_val:end-end_val), devRef_0(start_val:end-end_val),'r-'); hold on; grid on;
                else
                    plot(curvSorted(start_val:end-end_val), devRef_0(start_val:end-end_val),'r-'); hold on;
%                     plot(t1(start_val:end-end_val), dev_i(start_val:end-end_val),'r-'); hold on;
                end                
            else                   
                if strcmp(traj, 'circ')
                    plot(t1(start_val:end-end_val), devRef_0(start_val:end-end_val),'b-'); hold on; grid on;
                else
                    plot(curvSorted(start_val:end-end_val), devRef_0(start_val:end-end_val),'b-'); hold on;
%                     plot(t1(start_val:end-end_val), dev_i(start_val:end-end_val),'b-'); hold on;
                end            
            end            
        end
    end
end

for j=1:S.config.N+1
    sum_aux = 0;
    for ii=1:S.config.NUM_SIMS
        if ii~=avoid_val1 && ii~=avoid_val2
            sum_aux = sum_aux + mean(devs{ii,j}(start_val:end-end_val));
        end
    end
    mean_devs(j) = sum_aux/(S.config.NUM_SIMS);
end
mean_devs
%%
% if strcmp(traj, 'circ')
    % exp1                  Lh2    Lh2    Lh2
    data1_stack(:,:,1) = [1720 1608 1449
                          1606 1464 1279
                          1314 1194 1386].*1e-4;  % Lh1
    data1_stack(:,:,2) = [1055 1176 1344
                          1289 1440 1644
                          1569 1615 1452].*1e-4;  % Lh1
    data1_stack(:,:,3) = [3891 4141 4204
                          4150 4437 4549
                          4455 4635 4291].*1e-4;  % Lh1
    % exp2                  Lh2    Lh2    Lh2
    data2_stack(:,:,1) = [4718 5092 5220
                          5313 5559 5755
                          5969 6196 5677].*1e-4;  % Lh1
    data2_stack(:,:,2) = [2229 2611 2749
                          2733 3000 3194
                          3525 3738 3255].*1e-4;  % Lh1
    data2_stack(:,:,3) = [0430 0237 0156
                          0141 0292 0594
                          0945 1292 1456].*1e-4;  % Lh1
% else
    % exp1                  Lh2    Lh2    Lh2
    data3_stack(:,:,1) = [1619 1606 1572
                          1621 1596 1543
                          1613 1585 1530].*1e-4;  % Lh1
    data3_stack(:,:,2) = [1708 1755 1857
                          1827 1904 2009
                          1721 1802 1905].*1e-4;  % Lh1
    data3_stack(:,:,3) = [4866 5058 5022
                          5000 5236 5214
                          4869 5112 5114].*1e-4;  % Lh1
    % exp2                  Lh2    Lh2    Lh2
    data4_stack(:,:,1) = [2790 2903 2893
                          2856 2985 2930
                          2860 2998 2921].*1e-4;  % Lh1
    data4_stack(:,:,2) = [1375 1467 1565
                          1414 1588 1841
                          1533 1806 2168].*1e-4;  % Lh1
    data4_stack(:,:,3) = [2271 2218 2065
                          2295 2176 2100
                          2118 2014 2333].*1e-4;  % Lh1
% end

figure;
stacked_bar3(data1_stack);
xticks([1 2 3]);
set(gca, 'XTickLabel',{-0.3 0 0.3},'fontsize',26)
yticks([1 2 3]);
set(gca, 'YTickLabel',{-0.3 0 0.3},'fontsize',26)
ax1 = gca;
ax1.Projection = 'perspective';
ax1.ZGrid = 'on';
ax1.GridAlpha  = 0.4;
zlim([0 1.2])
xlabel({'$L_{h_2}$'},'interpreter','latex')
ylabel({'$L_{h_1}$'},'interpreter','latex')
box on
ax1.FontSize   = 20;
ax1.GridAlpha  = 0.3;
daspect([1 1 1])
exportgraphics(ax1,'bars3dcirc-1.eps','contenttype','vector')

figure
stacked_bar3(data2_stack);
xticks([1 2 3]);
set(gca, 'XTickLabel',{-0.3 0 0.3},'fontsize',26)
yticks([1 2 3]);
set(gca, 'YTickLabel',{-0.3 0 0.3},'fontsize',26)
ax2 = gca;
ax2.Projection = 'perspective';
ax2.ZGrid = 'on';
ax2.GridAlpha  = 0.4;
zlim([0 1.2])
xlabel({'$L_{h_2}$'},'interpreter','latex')
ylabel({'$L_{h_1}$'},'interpreter','latex')
box on
ax2.FontSize   = 20;
ax2.GridAlpha  = 0.3;
daspect([1 1 1])
exportgraphics(ax2,'bars3dcirc-2.eps','contenttype','vector')

figure;
stacked_bar3(data3_stack);
xticks([1 2 3]);
set(gca, 'XTickLabel',{-0.3 0 0.3},'fontsize',26)
yticks([1 2 3]);
set(gca, 'YTickLabel',{-0.3 0 0.3},'fontsize',26)
ax3 = gca;
ax3.Projection = 'perspective';
ax3.ZGrid = 'on';
ax3.GridAlpha  = 0.4;
zlim([0 1.2])
xlabel({'$L_{h_2}$'},'interpreter','latex')
ylabel({'$L_{h_1}$'},'interpreter','latex')
box on
ax3.FontSize   = 20;
ax3.GridAlpha  = 0.3;
daspect([1 1 1])
exportgraphics(ax3,'bars3dinfty-1.eps','contenttype','vector')

figure
stacked_bar3(data4_stack);
xticks([1 2 3]);
set(gca, 'XTickLabel',{-0.3 0 0.3},'fontsize',26)
yticks([1 2 3]);
set(gca, 'YTickLabel',{-0.3 0 0.3},'fontsize',26)
ax4 = gca;
ax4.Projection = 'perspective';
ax4.ZGrid = 'on';
ax4.GridAlpha  = 0.4;
zlim([0 1.2])
xlabel({'$L_{h_2}$'},'interpreter','latex')
ylabel({'$L_{h_1}$'},'interpreter','latex')
box on
ax4.FontSize   = 20;
ax4.GridAlpha  = 0.3;
daspect([1 1 1])
exportgraphics(ax4,'bars3dinfty-2.eps','contenttype','vector')



%%

axis('tight')
xlabel({'$curv.\,(m^{-1})$'},'interpreter','latex','fontsize',34)
ylabel({'$dev.\,(m)$'},'interpreter','latex','fontsize',34)
box on
% ylim([0 2.5])
axes            = gca;
axes.FontSize   = 35;
axes.GridAlpha  = 0.45;




if strcmp(traj,'circ')
    clear axes;
    axes('position',[.4 .4 .45 .25])
    box on % put box around new pair of axes
    ylim([0 0.6])
    t2 = (t1 < 60) & (t1 > 30); % range of t near perturbation
    for ii=1:S.config.NUM_SIMS
        if ii~=7
            for j=1:S.config.N+1
                if j==1
                    plot(t1(t2),devs{ii,j}(t2),'y'); hold on;
                elseif j==S.config.N+1
                    plot(t1(t2),devs{ii,j}(t2),'r'); hold on;
                else
                    plot(t1(t2),devs{ii,j}(t2),'b'); hold on;
                end            
            end
        end
    end
    box on; grid on;
    axes            = gca;
    axes.FontSize   = 20;
    axes.GridAlpha  = 0.35;
    xlim([30 60])
    ylim([0 0.6])
else
%     clear axes;
%     axes('position',[.57 .15 .45 .25])
%     box on % put box around new pair of axes
% 
%     plot(S.path.coordinates(1,:),S.path.coordinates(2,:),'color',[0.6 0.6 0.6],'LineWidth',8); grid on; daspect([1 1 1]); hold on;
%         xlim([min(S.path.coordinates(1,:))-1 max(S.path.coordinates(1,:))+1]); ylim([min(S.path.coordinates(2,:))-1 max(S.path.coordinates(2,:))+1]);
% 
%     box on; grid on;
%     axes            = gca;
%     axes.FontSize   = 20;
%     axes.GridAlpha  = 0.35;
end


%%

dt  = 0.1;
t   = 1:dt:2;
f   = 4.*sin(1.*t);
R   = [0 -1;1 0];

figure; hold on; daspect([1 1 1]); grid on;
plot(t(2:end-1),f(2:end-1),'ko','MarkerSize',20)

for ii=2:length(t)-2
    dx1     = dt;
    dy1     = f(ii)-f(ii-1);
    alpha   = f(ii+1)-f(ii);
    vec_aux = R * [dx1; dy1];
    dx2     = vec_aux(1);
    dy2     = vec_aux(2);
    quiver(t(ii),f(ii),dt,alpha,'Color','b','LineWidth',2)
    if ii==2
        quiver(t(ii),f(ii),dt,0,'Color','k','LineWidth',1)
        quiver(t(ii),f(ii),0,dt,'Color','k','LineWidth',1)
        text(t(ii)+dt,f(ii),'$x$','Interpreter','latex','fontsize',20);
        text(t(ii),f(ii)+dt,'$y$','Interpreter','latex','fontsize',20);
        text(t(ii),f(ii)-0.02,['$p_',num2str(ii-2),'$'],'Interpreter','latex','fontsize',20);
        theta0 = atan(dy1/dx1);
        w = 0:0.001:theta0-0.06;
        xarx = t(ii)+0.05*cos(w);
        yarx = f(ii)+0.05*sin(w);
        plot(xarx,yarx,'k')
        text(t(ii)+0.0157,f(ii)+0.015,'$\theta_0$','Interpreter','latex','fontsize',20);
    else
        quiver(t(ii),f(ii),dt,dy1,'Color',[0.07,0.62,1.00],'LineWidth',1)
        quiver(t(ii),f(ii),dx2,dy2,'Color',[0.07,0.62,1.00],'LineWidth',1)
        text(t(ii)+dt,f(ii)+dy1,['$x^',num2str(ii-2),'$'],'Interpreter','latex','fontsize',20);
        text(t(ii)+dx2,f(ii)+dy2,['$y^',num2str(ii-2),'$'],'Interpreter','latex','fontsize',20);
        text(t(ii),f(ii)-0.02,['$p_',num2str(ii-2),'$'],'Interpreter','latex','fontsize',20);
    end    
end
ax = gca;
ax.GridAlpha  = 0.4;
ylim([3.5 4.15])
xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',20);
% xticks(1:0.1:1.9);
% set(gca, 'XTickLabel',{0:10},'fontsize',26)
ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',20);
box on
ax.FontSize   = 20;
ax.GridAlpha  = 0.3;
daspect([1 1 1])
% exportgraphics(ax,'scheme_change_basis.eps','contenttype','vector')


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
% plot the vehicle reference generatio system
figure1=figure(1); grid on; hold on;
patchline(S.path.coordinates(1,1:end/4),S.path.coordinates(2,1:end/4),'linestyle','-.','edgecolor','k','linewidth',2,'edgealpha',0.5);
xlim([4 12]); ylim([4 14]);
daspect([1 1 1]);

N   = 7;
N1  = [1:N,S.config.N+1:S.config.N+1+N,2*S.config.N+2:2*S.config.N+1+2*(N+1)];


len = length(S.data.mhempc.performance.xest{1});        
% infty traj
indx1 = 1;
indx2 = ceil(0.25*(len/4));

plot_mono3(S, S.data.mhempc.performance.xest{1}(N1,indx2),[], N,[1 1 1 1 1 1 1 0.3]);hold on;

ax              = gca;
ax.GridAlpha    = 0.4;
box on
ax.FontSize     = 20;
ax.GridAlpha    = 0.3;
xlim([6 12])
ylim([9 14])
xlabel({'$x\,(m)$'},'interpreter','latex','fontsize',25);
ylabel({'$y\,(m)$'},'interpreter','latex','fontsize',25);
daspect([1 1 1])

a1 = annotation(figure1,'doublearrow',[0.37142857142857 0.392410714285711],...
    [0.573809523809525 0.35978725519435],'LineWidth',1,'LineStyle','--');

a2 = annotation(figure1,'doublearrow',[0.3125 0.326339285714282],...
    [0.671428571428572 0.521692017099113],'LineWidth',1,'LineStyle','--');

ax              = gca;
%  exportgraphics(figure1,'trac.eps','contenttype','vector')


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