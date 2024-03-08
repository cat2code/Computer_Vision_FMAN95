clearvars;close all;clc;
%%% Exercises 1-6 %%%
%% E1 new

% Define A and t from P2
A = [1 1 0; 0 2 0; 0 0 1];
t = [0; 2; 0];

% Compute the cross-product matrix [t]_x
t_x = [0, -t(3), t(2); t(3), 0, -t(1); -t(2), t(1), 0];

% Compute the fundamental matrix F
F = t_x * A;

% Given point x in homogeneous coordinates
x_homogeneous = [1; 1; 1];

% Compute the epipolar line in the second image for the point x
l2 = F * x_homogeneous;

% Given points to check
points = [2 0 1; 2 1 1; 4 2 1];

% Display the epipolar line
disp('Epipolar line in the second image:');
disp(l2);

% Check if the given points lie on the epipolar line l2
for i = 1:size(points, 1)
    point = points(i, :)';
    % The point lies on the epipolar line if the dot product is approximately zero
    if abs(l2' * point) < 1e-6  % Using a small threshold to account for floating-point errors
        disp(['Point (' num2str(point(1)) ',' num2str(point(2)) ') could be a projection of the same 3D point X into P2.']);
    else
        disp(['Point (' num2str(point(1)) ',' num2str(point(2)) ') could NOT be a projection of the same 3D point X into P2.']);
    end
end


%% E2

t1 = zeros(1,3)';
P1 = [eye(3) t1];
t2 = [2 2 0]';
A = [1 1 1; 0 2 0; 0 0 1];
P2 = [A t2];


% Compute the camera centers
C1 = -pinv(P1(:, 1:3)) * P1(:, 4);
C2 = -pinv(P2(:, 1:3)) * P2(:, 4);

% Compute the epipoles
e1 = P1(:, 1:3) * C2;
e2 = P2(:, 1:3) * C1;

% Display the results
disp('Epipole in the first image:');
disp(e1);
disp('Epipole in the second image:');
disp(e2);


% Compute the fundamental matrix
F = skew(t2) * A * pinv(A);

% Compute the determinant of the fundamental matrix
detF = det(F);

% Verify eT2 F = 0
eT2F = e2' * F;

% Verify F e1 = 0
Fe1 = F * e1;

% Display the results
disp('Fundamental matrix:');
disp(F);
disp('Determinant of the fundamental matrix:');
disp(detF);
disp('eT2 F:');
disp(eT2F);
disp('F e1:');
disp(Fe1);

%% E2.2
% TODO: Do this on paper

%For a general camera pair P1 = [I 0] and P2 = [A t]. Compute the epipoles, by projecting the camera centers. (You may assume that A is invertible.)

% Compute the camera centers
C1 = -pinv(P1(:, 1:3)) * P1(:, 4);
C2 = -pinv(P2(:, 1:3)) * P2(:, 4);

% Compute the epipoles
e1 = P1(:, 1:3) * C2;
e2 = P2(:, 1:3) * C1;

% Display the results
disp('Epipole in the first image:');
disp(e1);
disp('Epipole in the second image:');
disp(e2);

% Verify that for the fundamental matrix F = [t]×A the epipoles will always fulfill eT2 F = 0 and Fe1 = 0.

% Given the above result explain why the fundamental matrix has to have determinant 0.


%% E3

% TODO: Do this on paper

% Compute the fundamental matrix F using the 8-point algorithm
F_tilde = skew(t2) * A * pinv(A);

% Compute the normalization matrices N1 and N2
N1 = [1/size(x1, 2) 0 -mean(x1(1, :)); 0 1/size(x1, 2) -mean(x1(2, :)); 0 0 1];
N2 = [1/size(x2, 2) 0 -mean(x2(1, :)); 0 1/size(x2, 2) -mean(x2(2, :)); 0 0 1];

% Compute the fundamental matrix F for the original (un-normalized) points
F = N2' * F_tilde * N1;

%% E4

F = [0 1 1; 1 0 0; 0 1 1];

% Define the 3D scene points
X1 = [1; 2; 3];
X2 = [3; 2; 1];
X3 = [1; 0; 1];

% Define the camera matrices
P1 = [eye(3) zeros(3, 1)];
P2 = [skew(e2) * F e2];

% Compute the image projections
x1 = P1 * [X1; 1];
x2 = P2 * [X2; 1];
x3 = P2 * [X3; 1];

% Verify the epipolar constraint
epipolar_constraint = x2' * F * x1;

% Compute the camera center of P2
C2 = -pinv(P2(:, 1:3)) * P2(:, 4);

% Display the results
disp('Epipolar constraint (xT2 Fx1):');
disp(epipolar_constraint);
disp('Camera center of P2:');
disp(C2);

%% E4

% Define the fundamental matrix F
F = [0 1 1; 1 0 0; 0 1 1];

% Compute the epipole e2 as the null space of the transpose of F
e2 = null(F');

% Normalize e2 to ensure it is a homogeneous coordinate with the last entry as 1
e2 = e2 ./ e2(end);

% Display the result
disp('Epipole e2 in the second camera:');
disp(e2);

% Define the 3D scene points as homogeneous coordinates
X1 = [1, 2, 3, 1]';
X2 = [3, 2, 1, 1]';
X3 = [1, 0, 1, 1]';

% Define the camera matrices
P1 = [eye(3), zeros(3,1)]; % Camera P1 is the identity matrix and zero vector

E2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0]; % Skew-symmetric matrix for e2

% Define camera P2
P2 = [E2x * F, e2];

% Project points into each camera
x1_P1 = P1 * X1;
x2_P1 = P1 * X2;
x3_P1 = P1 * X3;

x1_P2 = P2 * X1;
x2_P2 = P2 * X2;
x3_P2 = P2 * X3;

% Normalize the projected points to convert from homogeneous to Cartesian coordinates
x1_P1 = x1_P1 / x1_P1(3);
x2_P1 = x2_P1 / x2_P1(3);
x3_P1 = x3_P1 / x3_P1(3);

%x1_P2 = x1_P2 / x1_P2(3);
x2_P2 = x2_P2 / x2_P2(3);
%x3_P2 = x3_P2 / x3_P2(3);

% Verify the epipolar constraint for each point
epipolar1 = x1_P2' * F * x1_P1;
epipolar2 = x2_P2' * F * x2_P1;
epipolar3 = x3_P2' * F * x3_P1;

% Output the results
fprintf('Epipolar constraint for X1: %f\n', epipolar1);
fprintf('Epipolar constraint for X2: %f\n', epipolar2);
fprintf('Epipolar constraint for X3: %f\n', epipolar3);

% The camera center for P2 can be found by null space of P2
camera_center_P2 = null(P2);

% Output the camera center of P2
fprintf('Camera center of P2: [%f, %f, %f]\n', camera_center_P2(1), camera_center_P2(2), camera_center_P2(3));





%% E6

% Define the matrices U and V
U = [1/sqrt(2), -1/sqrt(2), 0;
     1/sqrt(2), 1/sqrt(2), 0;
     0, 0, 1];

V = [1, 0, 0;
     0, 0, -1;
     0, 1, 0];

% Calculate the essential matrix E using singular value decomposition
E = U * diag([1, 1, 0]) * V';

% Verify that det(UV^T) = 1
det_UVT = det(U * V');

% Apply the essential matrix E to x2
x2 = [1; 1; 1]; % Homogeneous coordinates for (1,1)
Ex2 = E * x2;

% Display the results
disp('Essential Matrix E:');
disp(E);

disp('Determinant of UV^T:');
disp(det_UVT);

disp('Result of E * x2:');
disp(Ex2);

% Convert points to homogeneous coordinates
x1 = [0; 0; 1];
x2 = [1; 1; 1];

% Check the epipolar constraint
plausible_correspondence = x2' * E * x1;

% Display the result
disp('Plausibility of correspondence:');
disp(plausible_correspondence == 0);


P1 = [eye(3) zeros(3,1)];

pinv(P1)*x1 


% Define matrix W and vector u3
W = [0 -1 0;
     1 0 0;
     0 0 1];
 
% Extract u3 from previously defined U matrix
u3 = U(:, 3);

% Define possible camera matrix P2 configurations
P2_options = {
    U*W*V',            % [UWV']
    U*W*V' * -u3,      % [UWV' -u3]
    U*W'*V' * u3,      % [UW'V' u3]
    U*W'*V' * -u3,     % [UW'V' -u3]
};

% Point x2 in homogeneous coordinates
x2 = [1; 1; 1];

% Loop through all P2 options to find the correct s
for i = 1:length(P2_options)
    % Form the equation Ax = 0 where A is the matrix formed by stacking the skew-symmetric
    % form of x2 multiplied by the first two rows of P2_options{i} and x is the vector [s; 1]
    skew_x2 = [0 -x2(3) x2(2); x2(3) 0 -x2(1); -x2(2) x2(1) 0];
    A = skew_x2 * P2_options{i}(1:2, :);
    
    % Solve for s using singular value decomposition
    [~, ~, V_svd] = svd(A);
    s_candidate = V_svd(:, end);
    s_candidate = s_candidate / s_candidate(end); % Normalize to make the last entry 1
    
    % Calculate the 3D point X(s)
    X = [0; 0; s_candidate(1)];
    
    % Check if X(s) is in front of both cameras
    if P2_options{i}(3,:) * [X; 1] > 0 && U(3,:) * [X; 1] > 0
        s = s_candidate(1);
        P2 = P2_options{i};
        break;
    end
end

% Display the value of s and the chosen P2 configuration
disp('Value of s:');
disp(s);

disp('Chosen Camera Matrix P2:');
disp(P2);



%%% Done

