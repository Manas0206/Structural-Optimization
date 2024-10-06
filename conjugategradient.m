 clc
clear
% format long

% Define symbolic variables for the objective function:
syms x1 x2;

% Objective function:
f = 100*(x2 - (x1)^2)^2 + (1 - x1)^2;

% Initial guess:
x_1(1) = 0;
x_2(1) = 1;

Epsilon = 10^(-4); % Convergence criteria
i = 1; % Iteration counter

% Contour plot settings:
x1Label = linspace(-4,4,100);
x2Label = linspace(-3,6,100);
[x,y] = meshgrid(x1Label,x2Label);
f2 = 100.*(y - x.^2).^2 + (1 - x).^2; % Objective function for contour plot

figure(1)
contour(x,y,f2,'Fill','On');
hold on
plot(x_1(1), x_2(1), '*-k');
text(-3,5,['Initial Point (x1, x2) = (' num2str(x_1(1)), ', ', num2str(x_2(1)), ')'], 'Color', 'k')
xlabel('x1')
ylabel('x2')
title('Conjugate Gradient Method')
grid on
hold on

% Gradient computation:
df_dx1 = diff(f, x1);
df_dx2 = diff(f, x2);

% Evaluate the gradient at the initial guess:
J = [subs(df_dx1, [x1, x2], [x_1(1), x_2(1)]), subs(df_dx2, [x1, x2], [x_1(1), x_2(1)])];
S = -J; % Search direction

% Minimization condition:
while norm(double(S)) > Epsilon
    I = [x_1(i), x_2(i)]';
    
    % Step size determination:
    syms lambda; % Step size
    g = subs(f, [x1, x2], [x_1(i) + S(1)*lambda, x_2(i) + lambda*S(2)]);
    
    % Optimize the step length:
    dg_dlambda = diff(g, lambda) == 0;
    lambda_vals = double(solve(dg_dlambda, lambda));
    lambda_vals = lambda_vals(imag(lambda_vals) == 0); % Only real values
    
    % Evaluate the objective function for all real lambdas:
    for k = 1:length(lambda_vals)
        fun_value(k) = double(subs(f, [x1, x2], [x_1(i) + lambda_vals(k)*S(1), x_2(i) + lambda_vals(k)*S(2)]));
    end
    
    % Find the optimal step length:
    [~, index] = min(fun_value);
    lambda_opt = lambda_vals(index);
    
    % Update the plot:
    plot(x_1(i), x_2(i), '*-r');
    hold on
    
    % Update the points:
    x_1(i+1) = double(I(1) + lambda_opt*S(1)); % New x value
    x_2(i+1) = double(I(2) + lambda_opt*S(2)); % New y value
    
    % Update the gradient:
    J_old = [subs(df_dx1, [x1, x2], [x_1(i), x_2(i)]), subs(df_dx2, [x1, x2], [x_1(i), x_2(i)])];
    i = i + 1;
    J_new = [subs(df_dx1, [x1, x2], [x_1(i), x_2(i)]), subs(df_dx2, [x1, x2], [x_1(i), x_2(i)])];
    
    % Update the search direction:
    S = -double(J_new) + ((norm(double(J_new))^2) / (norm(double(J_old))^2)) * double(S);
end

% Plot the final point:
plot(x_1, x_2, '*-r');
hold on
plot(x_1(i), x_2(i), '*-k');
plot(x_1(1), x_2(1), '*-k');
text(0,-2,['Optimum point [x1*, x2*] = [', num2str(x_1(i)), ', ', num2str(x_2(i)), ']'], 'Color', 'k')

% Output the result:
fprintf('Initial Objective Function Value: %.4f\n\n', double(subs(f, [x1, x2], [x_1(1), x_2(1)])));

if (norm(S) < Epsilon)
    fprintf('Minimum successfully obtained...\n\n');
end

fprintf('Number of Iterations for Convergence: %d\n\n', i);
fprintf('Point of Minima: [%.4f, %.4f]\n\n', x_1(i), x_2(i));
fprintf('Objective Function Minimum Value: %.4f\n\n', double(subs(f, [x1, x2], [x_1(i), x_2(i)])));

% Result table:
Iter = (1:i)';
X_coordinate = x_1';
Y_coordinate = x_2';
T = table(Iter, X_coordinate, Y_coordinate);
disp(T)
