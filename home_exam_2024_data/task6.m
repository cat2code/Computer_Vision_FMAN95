%% Task 6
clearvars;close all;

% Given data
P1 = [1 0 0 0; 0 1 0 0; 0 0 1 0];
P2 = [1 0 4 1; 0 1 2 0; 0 0 1 1];
P3 = [-1 1 -2 1; -1 0 -1 0; 0 -2 0 0];
H23 = [0 0 2; -2 1 1; -2 0 0]; % homography between x3 ~ H23*2
x1 = [1; 1; 1]; % homogeneous coordinates for the point in the first image

% Find the 3D point X.

% Scale factors for the null space
lambda = sym('lambda');
mu = sym('mu');

% Equations from P1 and x1
eq1 = P1(1,:) * [lambda; lambda; lambda; 1] == x1(1) * P1(3,:) * [lambda; lambda; lambda; 1];
eq2 = P1(2,:) * [lambda; lambda; lambda; 1] == x1(2) * P1(3,:) * [lambda; lambda; lambda; 1];

% Equations from P2, P3 and the homography
eq3 = P2(1,:) * [mu; mu; mu; 1] == P3(1,:) * H23 * [mu; mu; 1];
eq4 = P2(2,:) * [mu; mu; mu; 1] == P3(2,:) * H23 * [mu; mu; 1];
eq5 = P2(3,:) * [mu; mu; mu; 1] == P3(3,:) * H23 * [mu; mu; 1];

% Solve the system
[sol_lambda, sol_mu] = solve([eq1, eq2, eq3, eq4, eq5], [lambda, mu]);

% Get the numeric solutions
num_lambda = double(sol_lambda);
num_mu = double(sol_mu);

% Use the first valid solution (if multiple)
valid_idx = find(~isnan(num_lambda) & ~isnan(num_mu), 1);
lambda = num_lambda(valid_idx);
mu = num_mu(valid_idx);

% Now, compute the 3D point using the found scale factors
X = [lambda; lambda; lambda; 1];
x2 = [mu; mu; 1];
x3 = H23 * x2;

% Display the results
disp('The 3D point X is:');
disp(vpa(X, 6)); % Display with 6 decimal places
disp('The projection x2 is:');
disp(vpa(x2, 6)); % Normalize and display
disp('The projection x3 is:');
disp(vpa(x3, 6)); % Normalize and display

%% Task 6
clearvars;close all;

% Given data
P1 = [1 0 0 0; 0 1 0 0; 0 0 1 0];
P2 = [1 0 4 1; 0 1 2 0; 0 0 1 1];
P3 = [-1 1 -2 1; -1 0 -1 0; 0 -2 0 0];
H23 = [0 0 2; -2 1 1; -2 0 0];
x1 = [1; 1; 1]; % homogeneous coordinates for the point in the first image

% Assume Z=1 for simplicity, making X = [x, y, z, 1]
X_assumed = [1; 1; 1; 1]; % Assumed for demonstration, directly using x1 with an assumed Z.

% Project X to find x2 and x3
x2_projected = P2 * X_assumed;
x3_projected = P3 * X_assumed;

% Normalize to homogeneous coordinates
x2 = x2_projected / x2_projected(3);
x3 = x3_projected / x3_projected(3);

disp('x2 = ');
disp(x2);
disp('x3 = ');
disp(x3);





























