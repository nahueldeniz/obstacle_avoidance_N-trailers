function S = handleStaticAndDynamicObs(S)
%     S.path.staticObs = S.path.listOfObsStr{end};

if S.config.iters==104
    a=1
end
    %
    if numel(S.path.listOfObsStr)<S.config.nroFrDetMObs        
        return;
    end
    % *********************************************************************
    % Here I compute how much the obstacles have moved in succesive
    % frames.If any reachs some treshould, it could be a moving obstacle
    diagDists = [];
    for i=S.config.nroFrDetMObs-1:-1:1
        dists       = pdist2(S.path.listOfObsStr{end-i}(:,1:2), S.path.listOfObsStr{end-i+1}(:,1:2));
        diagAux     = diag(dists);
%         diagDists   = [diagDists, diagAux(1:S.config.maxStaticObs)];
        diagDists   = [diagDists, diagAux(1:S.config.totNumObs)];
    end
    [f1,~] = find(diagDists(:,end)>S.config.thrMovingObs); % From all detected obstacles, this variable contains the indices of those that have moved.
%     nroFrmMov = [];
%     if ~isempty(f1)
%         for i=1:f1
%             nroFrmMov = [nroFrmMov; sum(diagDists(f1,:)>S.config.thrMovingObs,2)];
%         end
%     else
%         return;
%     end    
%     [~,indx] = sort(nroFrmMov);
indx=f1;    
    if numel(indx)>S.config.maxMovingObs
        indx = indx(1:S.config.maxMovingObs);
    end
    S.path.dynamicObs = S.path.listOfObsStr{end}(indx,:);
    %
    indxAux           = 1:numel(S.path.listOfObsStr{end}(:,1));
    indxStatic        = setdiff(indxAux', indx');
    S.path.staticObs  = S.path.listOfObsStr{end}(indxStatic,:); % after isolating the dynamic obstacles, the list of static obstacles is updated
    %
    % Find the ellipse's equations that model each moving obstacle
    %
    if ~isempty(S.path.dynamicObs)
        if isempty(S.path.dynObsEllipsePrms.a)  % If moving obstacles were found but I do not have yet modelled the ellipses
            for i=1:size(S.path.dynamicObs,1)
                indxMovObs  = pdist2(S.path.listOfObsStr{end}(:,1:2),S.path.dynamicObs(i,1:2));
                indx        = find(indxMovObs==0);
                xy          = [];
                for j=S.config.numMeasEll-1:-1:0
                    xy = [xy; S.path.listOfObsStr{end-j}(indx,1:2)];
                end
                S                               = findCoeffEllipse(S,xy(:,1),xy(:,2));
                S.path.dynObsEllipsePrms.a      = [S.path.dynObsEllipsePrms.a; S.path.dynObsEllipsePrmsAux.a];
                S.path.dynObsEllipsePrms.b      = [S.path.dynObsEllipsePrms.b; S.path.dynObsEllipsePrmsAux.b];
                S.path.dynObsEllipsePrms.xc     = [S.path.dynObsEllipsePrms.xc; S.path.dynObsEllipsePrmsAux.xc];
                S.path.dynObsEllipsePrms.yc     = [S.path.dynObsEllipsePrms.yc; S.path.dynObsEllipsePrmsAux.yc];
                S.path.dynObsEllipsePrms.t0     = [S.path.dynObsEllipsePrms.t0; S.path.dynObsEllipsePrmsAux.t0];
                S.path.dynObsEllipsePrms.kt     = [S.path.dynObsEllipsePrms.kt; S.path.dynObsEllipsePrmsAux.kt];
                S.path.dynObsEllipsePrms.theta  = [S.path.dynObsEllipsePrms.theta; S.path.dynObsEllipsePrmsAux.theta];
                S.path.dynObsEllipsePrms.t      = [S.path.dynObsEllipsePrms.t; S.path.dynObsEllipsePrmsAux.t];
                S.path.dynObsEllipsePrms.tOld   = [S.path.dynObsEllipsePrms.tOld; S.path.dynObsEllipsePrmsAux.tOld];
            end            
            return;
        elseif size(S.path.dynObsEllipsePrms.a,1)>size(S.path.dynamicObs,1) % If I have more ellipses than moving obstacles
            while size(S.path.dynObsEllipsePrms.a,1)>size(S.path.dynamicObs,1)
                D = [];
                for i=1:size(S.path.dynamicObs,1)
                    d = [];
                    for j=1:size(S.path.dynObsEllipsePrms.a,1)
                        xObs = S.path.dynamicObs(i,1);
                        yObs = S.path.dynamicObs(i,2);
                        d    = [d; sqrt( (xObs-S.path.dynObsEllipsePrms.xc(j)-S.path.dynObsEllipsePrms.a(j)*cos(S.path.dynObsEllipsePrms.t(i)))^2+...
                                         (yObs-S.path.dynObsEllipsePrms.yc(j)-S.path.dynObsEllipsePrms.b(j)*sin(S.path.dynObsEllipsePrms.t(i)))^2  )];
                    end
                    D = [D, d];
                end
                if size(D,2)==1
                    [~,indx] = max(D);                    
                else
                    diffD    = abs(diff(D,1,2));
                    [~,indx] = min(diffD);
                end
                S.path.dynObsEllipsePrms.a(indx)    = [];
                S.path.dynObsEllipsePrms.b(indx)    = [];
                S.path.dynObsEllipsePrms.xc(indx)   = [];
                S.path.dynObsEllipsePrms.yc(indx)   = [];
                S.path.dynObsEllipsePrms.t0(indx)   = [];
                S.path.dynObsEllipsePrms.kt(indx)   = [];
                S.path.dynObsEllipsePrms.theta(indx)= [];
                S.path.dynObsEllipsePrms.t(indx)    = [];
                S.path.dynObsEllipsePrms.tOld(indx) = [];
            end
        else % If I have more dynamic obstacles than ellipses
            while size(S.path.dynObsEllipsePrms.a,1)<size(S.path.dynamicObs,1)
                for i=1:size(S.path.dynamicObs,1)
                    d = [];
                    for j=1:size(S.path.dynObsEllipsePrms.a,1)
                        xObs = S.path.dynamicObs(i,1);
                        yObs = S.path.dynamicObs(i,2);
                        d = sqrt( (xObs-S.path.dynObsEllipsePrms.xc(j)-S.path.dynObsEllipsePrms.a(j)*cos(S.path.dynObsEllipsePrms.t(j)))^2+...
                                  (yObs-S.path.dynObsEllipsePrms.yc(j)-S.path.dynObsEllipsePrms.b(j)*sin(S.path.dynObsEllipsePrms.t(j)))^2  );
                        if any(d<S.config.distToUpdEll)
                            break;
                        end
                    end
                    [v,~] = min(d);
                    if v > S.config.distToUpdEll
                        xy   = [];
                        for k=S.config.numMeasEll-1:-1:0
                            xy = [xy; S.path.listOfObsStr{end-k}(i,1:2)];
                        end
                        S = findCoeffEllipse(S,xy(:,1),xy(:,2));
                        S.path.dynObsEllipsePrms.a      = [S.path.dynObsEllipsePrms.a; S.path.dynObsEllipsePrmsAux.a];
                        S.path.dynObsEllipsePrms.b      = [S.path.dynObsEllipsePrms.b; S.path.dynObsEllipsePrmsAux.b];
                        S.path.dynObsEllipsePrms.xc     = [S.path.dynObsEllipsePrms.xc; S.path.dynObsEllipsePrmsAux.xc];
                        S.path.dynObsEllipsePrms.yc     = [S.path.dynObsEllipsePrms.yc; S.path.dynObsEllipsePrmsAux.yc];
                        S.path.dynObsEllipsePrms.t0     = [S.path.dynObsEllipsePrms.t0; S.path.dynObsEllipsePrmsAux.t0];
                        S.path.dynObsEllipsePrms.kt     = [S.path.dynObsEllipsePrms.kt; S.path.dynObsEllipsePrmsAux.kt];
                        S.path.dynObsEllipsePrms.theta  = [S.path.dynObsEllipsePrms.theta; S.path.dynObsEllipsePrmsAux.theta];
                        S.path.dynObsEllipsePrms.t      = [S.path.dynObsEllipsePrms.t; S.path.dynObsEllipsePrmsAux.t];
                        S.path.dynObsEllipsePrms.tOld   = [S.path.dynObsEllipsePrms.tOld; S.path.dynObsEllipsePrmsAux.tOld];
                    end
                end
            end
        end
        % Now I have the same number of dynamic obstacles and ellipses
        indxObs = 0;   % I have to check if the current ellipses are good approximation to the moving obstacles, otherwise, I have to update the ellipse's equations
        D       = [];
        t0Temp  = [];
        for j=1:size(S.path.dynObsEllipsePrms.a,1)
            d = [];
            for i=1:size(S.path.dynamicObs,1)
                xObs = S.path.dynamicObs(i,1);
                yObs = S.path.dynamicObs(i,2);
                dis = sqrt( (xObs-S.path.dynObsEllipsePrms.xc(j)-S.path.dynObsEllipsePrms.a(j)*cos(S.path.dynObsEllipsePrms.t(j)))^2+...
                            (yObs-S.path.dynObsEllipsePrms.yc(j)-S.path.dynObsEllipsePrms.b(j)*sin(S.path.dynObsEllipsePrms.t(j)))^2  );
                if any(dis<S.config.distToUpdEll)
                    d = [d;0];
                    t0Temp = [];
                else
                    d = [d;dis];
                end                
            end
            D = [D, d];
        end
        [f,~] = find(D==0);
        f     = unique(f);
        if isempty(f)
            S.path.dynObsEllipsePrms.a      = [];
            S.path.dynObsEllipsePrms.b      = [];
            S.path.dynObsEllipsePrms.xc     = [];
            S.path.dynObsEllipsePrms.yc     = [];
            S.path.dynObsEllipsePrms.t0     = [];
            S.path.dynObsEllipsePrms.kt     = [];
            S.path.dynObsEllipsePrms.theta  = [];
            S.path.dynObsEllipsePrms.t      = [];
            S.path.dynObsEllipsePrms.tOld   = [];
            for i=1:size(S.path.dynamicObs,1)
                indxMovObs  = pdist2(S.path.listOfObsStr{end}(:,1:2),S.path.dynamicObs(i,1:2));
                indx        = find(indxMovObs==0);
                xy          = [];
%                 indxsEllip  = randperm(length(S.path.listOfObsStr),S.config.numMeasEll);
                for j=S.config.numMeasEll-1:-1:0
                    xy = [xy; S.path.listOfObsStr{end-(j)}(indx,1:2)];
                end
                S                               = findCoeffEllipse(S,xy(:,1),xy(:,2));
                S.path.dynObsEllipsePrms.a      = [S.path.dynObsEllipsePrms.a; S.path.dynObsEllipsePrmsAux.a];
                S.path.dynObsEllipsePrms.b      = [S.path.dynObsEllipsePrms.b; S.path.dynObsEllipsePrmsAux.b];
                S.path.dynObsEllipsePrms.xc     = [S.path.dynObsEllipsePrms.xc; S.path.dynObsEllipsePrmsAux.xc];
                S.path.dynObsEllipsePrms.yc     = [S.path.dynObsEllipsePrms.yc; S.path.dynObsEllipsePrmsAux.yc];
                S.path.dynObsEllipsePrms.t0     = [S.path.dynObsEllipsePrms.t0; S.path.dynObsEllipsePrmsAux.t0];
                S.path.dynObsEllipsePrms.kt     = [S.path.dynObsEllipsePrms.kt; S.path.dynObsEllipsePrmsAux.kt];
                S.path.dynObsEllipsePrms.theta  = [S.path.dynObsEllipsePrms.theta; S.path.dynObsEllipsePrmsAux.theta];
                S.path.dynObsEllipsePrms.t      = [S.path.dynObsEllipsePrms.t; S.path.dynObsEllipsePrmsAux.t];
                S.path.dynObsEllipsePrms.tOld   = [S.path.dynObsEllipsePrms.tOld; S.path.dynObsEllipsePrmsAux.tOld];
            end
        else
            indxTot     = [1:size(S.path.dynObsEllipsePrms.a,1)]';
            indxToLoop  = setdiff(indxTot, f');
% indxToLoop(find(indxToLoop>size(S.path.dynamicObs,1)))=[];            
            if ~isempty(indxToLoop)
                for i=indxToLoop
                    indxMovObs  = pdist2(S.path.listOfObsStr{end}(:,1:2),S.path.dynamicObs(i,1:2));
                    indx        = find(indxMovObs(:,1)==0);
                    xy          = [];
                    for j=S.config.numMeasEll-1:-1:0
                        xy = [xy; S.path.listOfObsStr{end-j}(indx,1:2)];
                    end
                    S                                  = findCoeffEllipse(S,xy(:,1),xy(:,2));
                    S.path.dynObsEllipsePrms.a(i)      = S.path.dynObsEllipsePrmsAux.a;
                    S.path.dynObsEllipsePrms.b(i)      = S.path.dynObsEllipsePrmsAux.b;
                    S.path.dynObsEllipsePrms.xc(i)     = S.path.dynObsEllipsePrmsAux.xc;
                    S.path.dynObsEllipsePrms.yc(i)     = S.path.dynObsEllipsePrmsAux.yc;
    %                     S.path.dynObsEllipsePrms.t0(i)     = S.path.dynObsEllipsePrmsAux.t0;
                    S.path.dynObsEllipsePrms.kt(i)     = S.path.dynObsEllipsePrmsAux.kt;
                    S.path.dynObsEllipsePrms.theta(i)  = S.path.dynObsEllipsePrmsAux.theta;
                    S.path.dynObsEllipsePrms.t(i)      = S.path.dynObsEllipsePrmsAux.t;
                    S.path.dynObsEllipsePrms.tOld(i)   = S.path.dynObsEllipsePrmsAux.tOld;
                end
            end
        end
        % Update the phsae for eahc moving obstacle as it moves along the
        % ellipse
        for j=1:length(S.path.dynObsEllipsePrms.a)
            xObs = S.path.dynamicObs(j,1);
            yObs = S.path.dynamicObs(j,2);
            S.path.dynObsEllipsePrms.t(j)     = atan2c((xObs-S.path.dynObsEllipsePrms.xc(j))/S.path.dynObsEllipsePrms.a(j), (yObs-S.path.dynObsEllipsePrms.yc(j))/S.path.dynObsEllipsePrms.b(j), S.path.dynObsEllipsePrms.t(j));
            S.path.dynObsEllipsePrms.tOld(j)  = S.path.dynObsEllipsePrms.t(j);                    
        end
        setT0MovingObs(S.mpc.mpcCasadi,S.path.dynObsEllipsePrms.t);
        






    elseif ~isempty(S.path.dynObsEllipsePrms.a) % If I do not have any moving obstacle
        for i=1:size(S.path.dynObsEllipsePrms.a,1)
            S.path.dynObsEllipsePrms.a      = [];
            S.path.dynObsEllipsePrms.b      = [];
            S.path.dynObsEllipsePrms.xc     = [];
            S.path.dynObsEllipsePrms.yc     = [];
            S.path.dynObsEllipsePrms.t0     = [];
            S.path.dynObsEllipsePrms.kt     = [];
            S.path.dynObsEllipsePrms.theta  = [];
            S.path.dynObsEllipsePrms.t      = [];
            S.path.dynObsEllipsePrms.tOld   = 0;
        end
    end



%     if ~isempty(S.path.dynObsAux)
% %         [~,indxs] = setdiff(S.path.listOfObs, S.path.dynObsAux,'rows');
% %         if ~isempty(indxs)
% %             S.path.dynObsXYcoords = [S.path.dynObsXYcoords; S.path.dynamicObs(indxs,1:2)];
% S.path.dynObsXYcoords = [S.path.dynObsXYcoords; S.path.dynamicObs(1,1:2)];
% %         end
%     end
%     if size(S.path.dynObsXYcoords,1) >= S.config.numMeasEll
%         d = sqrt((S.path.dynObsXYcoords(end,1)-S.path.dynObsEllipsePrms.xc-S.path.dynObsEllipsePrms.a*cos(S.path.dynObsEllipsePrms.t))^2+(S.path.dynObsXYcoords(end,2)-S.path.dynObsEllipsePrms.yc-S.path.dynObsEllipsePrms.b*sin(S.path.dynObsEllipsePrms.t))^2);
%         if any(d>0.2) || ~any(d) % When the measurement is beyond that 0.2 (m) with respect to the equation, it is computed again
% %             indxsEllip  = randperm(length(S.path.dynObsXYcoords),S.config.numMeasEll);
%             S               = findCoeffEllipse(S,S.path.dynObsXYcoords(end-S.config.numMeasEll+1:end,1),S.path.dynObsXYcoords(end-S.config.numMeasEll+1:end,2));
            setDynObsPrms(S.mpc.mpcCasadi,100*ones(length(S.path.dynObsEllipsePrms.xc),1),S.path.dynObsEllipsePrms.xc,S.path.dynObsEllipsePrms.yc,S.path.dynObsEllipsePrms.kt,S.path.dynObsEllipsePrms.a,S.path.dynObsEllipsePrms.b,S.path.dynObsEllipsePrms.theta);
            %
            if isempty(S.path.dynamicObs)
                S.fnmppc.optiMPC.set_value(S.fnmppc.ObsD, repmat([1e6, 1e6, 1],S.config.maxMovingObs,1));
            else
                if size(S.path.dynamicObs,1) <= S.config.maxMovingObs
                    S.fnmppc.optiMPC.set_value(S.fnmppc.ObsD, [S.path.dynamicObs;repmat([1e6 1e6 1],S.config.maxMovingObs-size(S.path.dynamicObs,1),1)] );
                else
                    S.fnmppc.optiMPC.set_value(S.fnmppc.ObsD, S.path.dynamicObs(1:S.config.maxMovingObs,:) );
                end
            end
%         end
                
%     end
end