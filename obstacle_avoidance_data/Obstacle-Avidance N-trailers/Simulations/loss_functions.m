clc;

e = -30:0.1:30;
a = -2;
c = 1;

figure; grid on; hold on,

for i=1:length(a)
    for j=1:length(c)
        f = (abs(a(i)-2)/a(i)) .* ( ( ((e/c(j)).^2)./abs((a(i)-2)) +1).^(a(i)/2) - 1);
        f2 = (e.^2)./(10 + e.^2);
        g = exp(-e.^2);
        subplot(1,2,1); hold on; grid on; plot(e, f,'y'); plot(e,2.*f2,'b'); plot(e,g,'r'); 
        subplot(1,2,2); hold on; grid on; plot(e,f2+g,'b')
        drawnow;
    end
end