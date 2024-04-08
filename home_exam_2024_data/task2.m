% Written by Eliot Montesino Petrén 
% el6183mo-s
clearvars;close all;clc;
%% Task 2
close all;clearvars;
load compEx3.mat

figure()
plot([startpoints(1,:); endpoints(1,:)], ... 
    [startpoints(2,:); endpoints(2,:)],'b-');
title("Original grid")
axis equal

%H1 = [sqrt(3) -1 1;1 sqrt(3) 1; 0 0 2];
%H2 = [1 -1 1 ;1 1 0; 0 0 1];

% Define the matrices H1 and H2
H1 = [1 0 0; 0 2 0; 0 0 1];
H2 = [0 -1 -1; 0 -1 0; -1 2 0];

% Define points in P^2
startpoints_hat = [startpoints; ones(1,length(startpoints))];
endpoints_hat = [endpoints; ones(1,length(endpoints))];


%%% H1 %%%

% Let H1 transform points. Keeps it a 3x42 matrix
startpoints_hat = H1*startpoints_hat;
endpoints_hat = H1*endpoints_hat;
% Calculate cartesian coordinates
startpoints_hat = pflat(startpoints_hat,0);
endpoints_hat = pflat(endpoints_hat,0);

figure()
plot([startpoints(1,:); endpoints(1,:)], ... 
    [startpoints(2,:); endpoints(2,:)],'b-');
hold on;
plot([startpoints_hat(1,:); endpoints_hat(1,:)], ... 
    [startpoints_hat(2,:); endpoints_hat(2,:)],'b-');
title("H1 figure")
axis equal

% Reset starting and ending points
startpoints_hat = [startpoints; ones(1,length(startpoints))];
endpoints_hat = [endpoints; ones(1,length(endpoints))];


%%% H2 %%%

% Let H2 transform points. Keeps it a 3x42 matrix
startpoints_hat = H2*startpoints_hat;
endpoints_hat = H2*endpoints_hat;
% Calculate cartesian coordinates
startpoints_hat = pflat(startpoints_hat,0);
endpoints_hat = pflat(endpoints_hat,0);

figure()
plot([startpoints(1,:); endpoints(1,:)], ... 
    [startpoints(2,:); endpoints(2,:)],'b-');
hold on;
plot([startpoints_hat(1,:); endpoints_hat(1,:)], ... 
    [startpoints_hat(2,:); endpoints_hat(2,:)],'b-');
title("H2 figure")
axis equal

% Reset starting and ending points
startpoints_hat = [startpoints; ones(1,length(startpoints))];
endpoints_hat = [endpoints; ones(1,length(endpoints))];


% b) Transform the points (1, 0), (1, 1) and (2, 1) using H1 and H2. Which 
% of the transformed points can be interpreted as regular points in R2? 
% What is the interpretation of the other transformed points?

% Define the points in homogeneous coordinates
points = [1 0 1; 1 1 1; 2 1 1]'; % Transpose to make columns

points_H1 = H1 * points;
points_H2 = H2 * points;

% Normalize against homogenous coord
points_H1 = points_H1 ./points_H1(3, :);
points_H2 = points_H2 ./ points_H2(3, :);

% Extract Cartesian coord 
points_H1_cart = points_H1(1:2, :);
points_H2_cart= points_H2(1:2, :);


disp('Transformed points using H1:');
disp(points_H1);
disp('Transformed points using H2:');
disp(points_H2);
disp('Cartesian Transformed points using H1:');
disp(points_H1_cart);
disp('Cartesian Transformed points using H2:');
disp(points_H2_cart);




%%% c)

% Compute H = H2_inv * H1
H2_inv = inv(H2);
H = H2_inv * H1;

% Find the eigenvectors and eigenvalues of H
[eigenvectors, eigenvalues] = eig(H);

% Find where eigenvalue is 1 
idx = abs(diag(eigenvalues) - 1) <= 0;

% The eigenvectors corresponding to the eigenvalue 1 the invariant points
invariant_points = eigenvectors(:, idx);
invariant_points = invariant_points./invariant_points(end,end);

disp('Invariant points for which H1x = H2x are:');
disp(invariant_points);

% Quick verification.
points = [2 0 -2; 0.7071 0 -0.7071; 1 0 -1]';
points_H1 = H1 * points;
points_H2 = H2 * points;














