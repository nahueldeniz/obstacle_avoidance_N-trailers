

%Firefly

NUM_OF_ITERATIONS = 50;
NUM_OF_FIREFLIES = 20;
ALPHA = 0.3;
BETA = 0.5;
GAMMA = 0.5;
K1 = 1;
K2 = 1;

% Inicialización de las posiciones de las luciérnagas
fireflies = rand(num_fireflies, 2) * 20 - 10;
initial_tractor_position = [0 0];
fireflies(1,:) = initial_tractor_position;
goal = [9, 9];
obstacles = [4, 4; 5, 8; 8, 5; 7, 7; 8, 2];

% Evaluar la función objetivo para cada luciérnaga
brightness = zeros(num_fireflies, 1);
for i = 1:NUM_OF_FIREFLIES
    brightness(i) = objective_function(fireflies(i,:), goal, obstacles, K1, K2);
end

% Para visualización
figure;
xlim([0, 10]);
ylim([0, 10]);
hold on;
goal_plot = plot(goal(1), goal(2), 'go', 'MarkerSize', 10, 'DisplayName', 'Meta');
obstacles_plot = plot(obstacles(:,1), obstacles(:,2), 'ro', 'MarkerSize', 10, 'DisplayName', 'Obstáculo');
fireflies_plot = plot(fireflies(:,1), fireflies(:,2), 'bo');
grid on

% Algoritmo Firefly adaptado según paper
robot_path = [];
for t = 1:NUM_OF_ITERATIONS
    for i = 1:NUM_OF_FIREFLIES
        for j = 1:NUM_OF_FIREFLIES
            if brightness(j) < brightness(i) % Luciérnaga j es más brillante que luciérnaga i
                r = norm(fireflies(i,:) - fireflies(j,:));
                attractiveness = BETA * exp(-gamma * r^2);
                fireflies(i,:) = fireflies(i,:) + attractiveness * (fireflies(j,:) - fireflies(i,:)) + ALPHA * (rand(1,2) - 0.5);
                fireflies(i,:) = max(min(fireflies(i,:), 10), 0); % Limitar dentro del espacio de búsqueda
                brightness(i) = objective_function(fireflies(i,:), goal, obstacles, K1, K2);
            end
        end
    end
    robot_path = [robot_path; fireflies(1,:)]; % Guardar la trayectoria del primer robot
    update_plot(fireflies, robot_path);
end


function fi = objective_function(x,goal,obstacles,K1,K2)
    
    distance_firefly_to_goal = norm(x - goal);
    min_On = inf;
    for i = 1:length(obstacles)
        distance_firefly_to_obstacle = norm(x - obstacles(i,:));
        if distance_firefly_to_obstacle < min_On
           min_On = distance_firefly_to_obstacle;
        end
    end
    fi = K1*(1/min_On)+K2*distance_firefly_to_goal;
end


% Función de actualización de la visualización
function update_plot(fireflies, robot_path)
    plot(fireflies(:,1), fireflies(:,2), 'bo');
    path_x = robot_path(:,1);
    path_x = robot_path(:,1);
    path_x = robot_path(:,1);
    path_y = robot_path(:,2);
    plot(path_x, path_y, 'b--');
    drawnow;
    pause(0.1);
end