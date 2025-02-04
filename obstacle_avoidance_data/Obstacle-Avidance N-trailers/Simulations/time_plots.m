

figure; hold on; grid on;
max_de_max      = 0;
mean_hard       = 0;
mean_soft       = 0;
mean_barrier    = 0;

for i=1:10
    m = max(time_hard{i}.t_tot);
    if m > max_de_max
        max_de_max = m;
    end    
    %
    m = max(time_soft{i}.t_tot);
    if m > max_de_max
        max_de_max = m;
    end    
    %
    m = max(time_barrier{i}.t_tot);
    if m > max_de_max
        max_de_max = m;
    end    
end

for i=1:10
    mean_hard = mean_hard + mean(time_hard{i}.t_tot./max_de_max);
    mean_soft = mean_soft + mean(time_soft{i}.t_tot./max_de_max);
    mean_barrier = mean_barrier + mean(time_barrier{i}.t_tot./max_de_max);

    plot(time_hard{i}.t_tot./max_de_max,'r')
    plot(time_soft{i}.t_tot./max_de_max,'b')
    plot(time_barrier{i}.t_tot./max_de_max,'g')
end
mean_hard = mean_hard/10;
mean_soft = mean_soft/10;
mean_barrier = mean_barrier/10;



