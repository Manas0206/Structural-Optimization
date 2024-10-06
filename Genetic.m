% Define x and y as symbolic variables
syms x y;

% Define expressions for P1, P2, u1, and u2
P1 = -32 + 859.7 * y + 26.2 * x - 1010.7 * y^2 + 77.3 * x * y - 9.7 * x^2;
P2 = -43.6459 + 145.0015 * y + 131.0431 * x - 218.7569 * y^2 + 22.4834 * x * y - 15.5354 * x^2;
u1 = 18.6153 - 13.1276 * y + 8.8413 * x + 17.5066 * y^2 + 2.1654 * x * y - 0.6867 * x^2;
u2 = 125.5877 - 24.6597 * y - 51.1271 * x - 36.5219 * y^2 + 27.0401 * x * y + 6.6905 * x^2;

% Define the toughness expression
Toughness = 0.5 * P1 * u1 + 0.5 * (P1 + P2) * (u2 - u1);

% Convert the symbolic expression to a MATLAB function
Toughness_fun = matlabFunction(Toughness);

% Define the objective function
objective_function = @(x) -Toughness_fun(x(1), x(2));

% Define initial guess
x0 = [2 0.7];

% Define bounds for X
lb = [1, 0.13]; % Lower bounds
ub = [3, 1]; % Upper bounds

% Call the GA solver
options = optimoptions('ga', 'Display', 'iter');
[X_opt, fval, exitflag, output] = ga(objective_function, 2, [], [], [], [], lb, ub, [], options);

% Display results
fprintf('Optimal X: %.4f, %.4f\n', X_opt(1), X_opt(2));
fprintf('Maximum toughness: %.4f\n', -fval);
disp(output);

% Plot the surface of the objective function
figure;
[X_mesh, Y_mesh] = meshgrid(linspace(1, 4, 100), linspace(0, 1, 100));
surf(X_mesh, Y_mesh, Toughness_fun(X_mesh, Y_mesh));
xlabel('X');
ylabel('Y');
zlabel('Toughness');
title('Surface Plot of Toughness Objective Function');

figure;
contour(X_mesh, Y_mesh, Toughness_fun(X_mesh, Y_mesh), 'LevelStep', 200, 'LineWidth', 1.5);
hold on;
plot(x0(1), x0(2), 'ro', 'MarkerSize', 10, 'LineWidth', 1.5);  % Initial guess
plot(X_opt(1), X_opt(2), 'go', 'MarkerSize', 10, 'LineWidth', 1.5);  % Optimal point
xlabel('X', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y', 'FontSize', 12, 'FontWeight', 'bold');
title('Contours of Toughness Objective Function', 'FontSize', 14, 'FontWeight', 'bold');
legend('Contours', 'Initial Guess', 'Optimal Point', 'Location', 'northwest', 'FontSize', 10);
xlim([1, 4]);
ylim([0, 0.6]);
grid on;
box on;
set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'FontWeight', 'bold');
