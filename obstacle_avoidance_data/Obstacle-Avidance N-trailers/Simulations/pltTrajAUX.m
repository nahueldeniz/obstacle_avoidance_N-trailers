hold on;
for i=1:S.config.Nt+1 
    if i==1 
        clr = 'y'; 
    elseif i==S.config.Nt+1
        clr = 'r'; 
    else 
        clr = 'b'; 
    end 
    plot(S.data.xsim(2*S.config.Nt+1+(i-1)*2+1,:),S.data.xsim(2*S.config.Nt+1+i*2,:),clr) 
end

%%
disSegToMovingObs       = [];
for i=1:min(length(S.path.dynObsXYcoords),length(S.data.xsim))
    disSegToMovingObsAux = [];
    for j=1:S.config.Nt+1
        disSegToMovingObsAux = [disSegToMovingObsAux; norm(S.data.xsim(2*S.config.Nt+1+(j-1)*2+1:2*S.config.Nt+1+j*2,i)-S.path.dynObsXYcoords(i,:)')];
    end
    disSegToMovingObs = [disSegToMovingObs, disSegToMovingObsAux];
end

figure; hold on;
for j=1:S.config.Nt+1
    if j==1
        clr = 'y';
    elseif j==S.config.Nt+1
        clr = 'r';
    else
        clr = 'b';
    end
    plot(disSegToMovingObs(j,:),clr)
end