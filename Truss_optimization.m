%% 7 bars, 5 nodes, 1-point load
clear all;
clc;
format long;
%% Input values
rho = 7800; %kg/m^3
E = 210 * 1e9; % Pa
len = [0.075 0.1581 0.090 0.090 0.15 0.1581 0.075]; % m
W = 0.025; % kg
Amin = 1e-6; % m^2
Amax = 16 * 1e-6; % m^2
tol = 1e-10; % convergence tolerance
p = zeros(10, 1); % load vector
p(4) = -100; % Newtons
p = p(3:9); % remove rows 1,2,10 corresponding to support reactions
% initial guess
% Ai = 0.5 * (Amin + Amax) * ones(7, 1);
Ai = 0.5 * (Amin + Amax) * rand(7, 1);
% Set A3 and A4 to Amin
Ai(3) = 1e-6; % m^2
Ai(4) = 1e-6; % m^2
%% OC iteration
 Ak = Ai; % Reset Ak to initial guess
 res = 1;
 count = 1;
 res_vec = zeros(1, 10); % Initialize residual vector
 while res > tol && count <= 15 % Limit maximum iterations for the outerloop to 10
 [K, K1, K2, K3, K4, K5, K6, K7] = Calc_Kmatrix(E, len, Ak);
 u = inv(K) * p;
 vec = [u' * (K1 / Ak(1)) * u, u' * (K2 / Ak(2)) * u, u' * (K3 /Ak(3)) * u, u' * (K4 / Ak(4)) * u, u' * (K5 / Ak(5)) * u, u' * (K6 / Ak(6)) * u, u' * (K7 / Ak(7))* u]';
 mu = sum(vec .* Ak) / W;
 Ak1 = (Ak .* vec) ./ (mu * rho * len');

 % Set A3 and A4 to Amin
 Ak1(3) = 1e-6; % m^2
 Ak1(4) = 1e-6; % m^2
 count_in = 1; % Initialize inner loop counter
 while any(Ak1 < Amin) || any(Ak1 > Amax) % need to do binning
 if count_in > 15
 disp('Maximum iterations reached in inner loop.');
 break; % Exit the loop if maximum iterations reached
 end
 disp('Warning: Areas are out of bounds, running the inner loop');
 ind_min = find(Ak1 < Amin);
 ind_max = find(Ak1 > Amax);
 ind_mid = setdiff(1:7, [ind_min; ind_max]);
 sum1=0;
 sum2=0;
        for i=1:length(ind_min)
            k=ind_min(i);
            sum1=sum1+len(k);
        end
        for i=1:length(ind_max)
            k=ind_max(i);
            sum2=sum2+len(k);
        end
 mu = sum(vec(ind_mid).*Ak(ind_mid))/(W-(rho*(sum1*Amin + sum2*Amax )));
 Ak1 = (Ak .* vec) ./ (mu * rho * len');
 % Set A3 and A4 to Amin
 Ak1(3) = 1e-6; % m^2
 Ak1(4) = 1e-6; % m^2
  alpha = 0.1; % mixing parameter
  Ak1 = alpha * Ak1 + (1 - alpha) * Ak; % mixing

 fprintf('Inner loop no. %d \n', count_in);
 count_in = count_in + 1;
 end
 res = norm(Ak - Ak1) / norm(Ak);
 res_vec(count) = res;
 fprintf('Outer loop no. %d, res=%.8f \n', count, res);
 alpha = 0.65; % mixing parameter
 Ak = alpha * Ak1 + (1 - alpha) * Ak; % mixing
 count = count + 1;
 end
 Ak
%% convergence plot
semilogy(1:length(res_vec),res_vec,'-*')
xlabel('iterations'); ylabel('residual')
objective_function=p'*u;
%% Functions required for the above problem
function [Kr, K1, K2, K3, K4, K5, K6, K7] = Calc_Kmatrix(E, len, A)
K = zeros(10, 10);

% bar 1 (DOFs:1,2,5,6)
member_ind = 1;
dof_set = [1, 2, 5, 6];
theta = 0 * pi / 180;
K1 = Assemble_Ke_into_K(E, len, A, member_ind, dof_set, theta, K);
% bar 2 (DOFs:1,2,3,4)
member_ind = 2;
dof_set = [1, 2, 3, 4];
theta = 18.435 * pi / 180;
K2 = Assemble_Ke_into_K(E, len, A, member_ind, dof_set, theta, K);
% bar 3 (DOFs:3,4,5,6)
member_ind = 3;
dof_set = [3, 4, 5, 6];
theta = 33.7 * pi / 180;
K3 = Assemble_Ke_into_K(E, len, A, member_ind, dof_set, theta, K);
% bar 4 (DOFs:3,4,7,8)
member_ind = 4;
dof_set = [3, 4, 7, 8];
theta = -33.7 * pi / 180;
K4 = Assemble_Ke_into_K(E, len, A, member_ind, dof_set, theta, K);
% bar 5 (DOFs:5,6,7,8)
member_ind = 5;
dof_set = [5, 6, 7, 8];
theta = 0 * pi / 180;
K5 = Assemble_Ke_into_K(E, len, A, member_ind, dof_set, theta, K);
% bar 6 (DOFs:5,6,9,10)
member_ind = 6;
dof_set = [3, 4, 9, 10];
theta = -18.435 * pi / 180;
K6 = Assemble_Ke_into_K(E, len, A, member_ind, dof_set, theta, K);
% bar 7 (DOFs:7,8,9,10)
member_ind = 7;
dof_set = [7, 8, 9, 10];
theta = 0 * pi / 180;
K7 = Assemble_Ke_into_K(E, len, A, member_ind, dof_set, theta, K);
K = K1 + K2 + K3 + K4 + K5 + K6 + K7;
% Remove the rows and columns of K matrices corresponding to DOFs 1,2,10
Kr = K(3:9, 3:9);
K1 = K1(3:9, 3:9);
K2 = K2(3:9, 3:9);

K3 = K3(3:9, 3:9);
K4 = K4(3:9, 3:9);
K5 = K5(3:9, 3:9);
K6 = K6(3:9, 3:9);
K7 = K7(3:9, 3:9);
end
function K = Assemble_Ke_into_K(E, len, A, member_ind, dof_set, theta, K)
ke = (E * A(member_ind)/ len(member_ind)) * [1, -1; -1, 1];
l = cos(theta);
m = sin(theta);
L = [l, m, 0, 0; 0, 0, l, m];
Ke = L' * ke * L; % Element stiffness matrix
K(dof_set, dof_set) = K(dof_set, dof_set) + Ke;
end
