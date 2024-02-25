clearvars;close all;clc;
%%% Exercises 1-4 %%%

%% E1

t1 = zeros(1,3)';
P1 = [eye(3) t1];
t2 = [0 2 0]';
A = [1 1 0; 0 2 0; 0 0 1];
P2 = [A t2];

F = skew(t2) * A * pinv(A);

% Given point x
x = [1; 1];

% Compute the epipolar line in the second image
l2 = F * x;

% Display the result
disp('Epipolar line in the second image:');
disp(l2);

% Given points
points = [2 0; 2 1; 4 2];

% Check if the points could be projections of the same point X into P2
for i = 1:size(points, 1)
    x = [points(i, 1); points(i, 2); 1];
    l2 = F * x;
    disp(['Point (' num2str(points(i, 1)) ',' num2str(points(i, 2)) ') is on the epipolar line in the second image:']);
    disp(l2);
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


%%% Done

