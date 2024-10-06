% Define objective function
fun = @(x1, x2) (x1 - 2)^2 + (x2 - 1)^2;

% Define constraints
x1 = linspace(-2, 4, 100);
x2 = linspace(-2, 4, 100);
[X1, X2] = meshgrid(x1, x2);

% Constraint: x1 + x2 <= 2
C1 = X1 + X2 - 2;

% Constraint: x1^2 <= x2
C2 = X1.^2 - X2;

% Plot contours of the objective function
figure;
contour(X1, X2, fun(X1, X2), 30);
hold on;

% Plot the feasible region
contour(X1, X2, C1, [0 0], 'LineWidth', 2, 'Color', 'blue');
contour(X1, X2, C2, [0 0], 'LineWidth', 2, 'Color', 'red');

% Add color bar to represent objective function
colorbar;

% Find and plot KKT points
x0 = [0; 0];
options = optimoptions('fmincon', 'Display', 'off');
[x, fval] = fmincon(@(x) fun(x(1), x(2)), x0, [], [], [], [], [], [], @nonlcon, options);
plot(x(1), x(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Mark local minima
[x_min, y_min] = fminsearch(@(x) fun(x(1), x(2)), x0);
plot(x_min(1), y_min, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');

% Labels and legend
xlabel('x_1');
ylabel('x_2');
title('Contour Plot with Feasible Region, KKT Points, and Local Minima');
legend('Objective Function Contours', 'Feasible Region (x1 + x2 <= 2)', 'Feasible Region (x1^2 <= x2)', 'KKT Points', 'Local Minima', 'Location', 'Best');
grid on;
axis equal;
hold off;

function [c, ceq] = nonlcon(x)
    % Nonlinear inequality constraint: x1^2 - x2 <= 0
    c = x(1)^2 - x(2);
    % No nonlinear equality constraints
    ceq = [];
end
