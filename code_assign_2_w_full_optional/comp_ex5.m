%% 7 Triangulation using DLT ((OPTIONAL))
%%% Computer Exercise 5 %%%
clearvars;close all;clc;

% Load necessary data
load for_comp_e5.mat          % Loads x1 (2x1610), x2 - should be correct
load for_comp_e5_2.mat        % From computer exercise 3
load compEx3data.mat          % Contains stuff needed

% Read images
im1 = imread("cube1.jpg");
im2 = imread("cube2.jpg");

% Attempt without normalization
K1 = eye(3);
K2 = eye(3);

% Normalize SIFT points using inverse camera calibration matrices
x1 = pflat(inv(K1)*[x1; ones(1, size(x1, 2))], 0);
x2 = pflat(inv(K2)*[x2; ones(1, size(x1, 2))], 0);

% K-Normalize camera matrices
P1 = inv(K1) * P1;
P2 = inv(K2) * P2;

% Initialize matrix for 3D points
X = [];

% DLT loop for triangulation
for i = 1:size(x1, 2)
    % Form matrix M for SVD
    M = [P1 -x1(:,i) zeros(3, 1); P2 zeros(3, 1) -x2(:,i)];
    
    % Singular Value Decomposition
    [U, S, V] = svd(M);
    v = V(:, end); % Extract last column V
    X = [X v(1:4, 1)]; % Append solution
end
X = pflat(X, 0); % Project back to homogeneous coordinates

% K-Un-normalize camera matrices
P1 = K1 * P1;
P2 = K2 * P2;

% Project calculated 3D points onto image planes
x_proj_1 = pflat(P1 * X, 0);
x_proj_2 = pflat(P2 * X, 0);

% Re-project measured points onto image planes
x1 = pflat(K1 * x1, 0);
x2 = pflat(K2 * x2, 0);

% Plot for first image
figure; 
imagesc(im1); 
hold on; 
axis equal; 
plot(x_proj_1(1,:), x_proj_1(2,:), 'g*');  % Projected points
plot(x1(1,:), x1(2,:), 'ro');              % Measured points
%title('Image 1 projected points (with normalization)')
title('Image 1 projected points (no normalization)')

% Plot for second image
figure; 
imagesc(im2); 
hold on; 
axis equal; 
plot(x_proj_2(1,:), x_proj_2(2,:), 'g*');  % Projected points
plot(x2(1,:), x2(2,:), 'ro');              % Measured points
%title('Image 2 projected points (with normalization)')
title('Image 1 projected points (no normalization)')

% Find points with error less than 3 pixels
good_fit = (sqrt(sum((x1(1:2,:) - x_proj_1(1:2,:)).^2)) < 3 & ...
            sqrt(sum((x2(1:2,:) - x_proj_2(1:2,:)).^2)) < 3);
X = X(:, good_fit); % Filter points based on good fit

% Final 3D plotting of good points and Rubik's cube edges
figure;
axis equal;
plot3([Xmodel(1, startind); Xmodel(1, endind)], ...
      [Xmodel(2, startind); Xmodel(2, endind)], ...
      [Xmodel(3, startind); Xmodel(3, endind)], 'b-');
hold on;
plot3(X(1,:), X(2,:), X(3,:), 'r.');


% Plot camera positions and orientations (assuming plotcams is a predefined function)
plotcams({P1, P2});

title('3D plot of points with small error')



