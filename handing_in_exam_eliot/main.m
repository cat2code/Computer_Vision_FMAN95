% Written by Eliot Montesino Petrén 
% el6183mo-s
clearvars;close all;clc;
%% Task 1
clearvars;close all;

%%% a)

% Define the projection matrix P
P = [2 0 1 0; 1 -1 1 1; 2 1 1 0];

% Define the 3D points as homogeneous coordinates
points3D = [1 0 0 1; 0 1 0 1; 0 0 1 1];

% Compute the projections of the 3D points
proj2D = P * points3D';

% Normalize the points to convert from homogeneous coordinates to Cartesian coordinates
proj2D_cartesian = proj2D(1:2, :) ./ proj2D(3, :);

% Compute the projection of a point at infinity in the direction (1, 0, -3)
point_at_infinity = [1 0 -3 0]; % Homogeneous coordinates for a point at infinity
proj_infinity = P * point_at_infinity';

% Since the point is at infinity, we do not need to normalize
proj_infinity_cartesian = proj_infinity(1:2);

% Display the results
disp('2D projections of the 3D points:');
disp(proj2D_cartesian);

disp('2D projection of the point at infinity:');
disp(proj_infinity_cartesian);

%%% b) 

A = P(:, 1:3);
A3 = A(:, 3);

% Perform Singular Value Decomposition on P
[U, S, V] = svd(P);

% The camera center is given by the last column of V (in homogeneous coord)
C = V(:, end);

% Normalize against homogenous coordinate and obtain cartesian coord.
C = C(1:3) / C(4);


disp('Camera Center (Cartesian coordinates):');
disp(C);

% Compute the principal axis R3, normalize A3 by its norm (p.32 lec notes)
R3 = sign(det(A)) * A3 / norm(A3);

% Display the principal axis
disp('Principal Axis (normalized):');
disp(R3);

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


%% Task 3
clearvars;close all;

% Essential Matrix
E = [0 1 0; -1 0 0; 0 0 0];

% a) Determine which of the three following 2D point pairs could be 
% the projections of the same 3D point?

% Answer:
% By the definition of the Essential Matrix, two projections corresponding
% to the same point should result in zero (in the ideal case where there 
% is no noise and margin of error). Check this for all 3 pairs.

% Def 2D points in homogenous coordinates for all pairs
x1 = [1; 1; 1];
x1_prime = [-1; -1; 1];
result1 = x1_prime' * E * x1;

x2 = [0; 2; 1];
x2_prime = [2; 1; 1];
result2 = x2_prime' * E * x2;

x3 = [2; -2; 1];
x3_prime = [-1; 1; 1];
result3 = x3_prime' * E * x3;

disp(result1); % Zero!
disp(result2); % Non-zero!
disp(result3); % Zero!


% b) Compute the two epipoles.

% Essential Matrix
E = [0 1 0; -1 0 0; 0 0 0];

% Compute the epipoles
% For the left epipole (null space of E)
[U,S,V] = svd(E);
e_left = V(:,end); % The last column of V from SVD of E

% For the right epipole (null space of E transpose)
[U_t,S_t,V_t] = svd(E');
e_right = V_t(:,end);

disp('Left Epipole:');
disp(e_left);
disp('Right Epipole:');
disp(e_right);


% c) The matrix E has the singular value decomposition ****
% Use this to compute the four possible camera pairs from E.

% Define the matrices U, S, and V_trans (V') from the provided SVD of E
U = eye(3); 
S = [1 0 0; 0 1 0; 0 0 0]; 
V_trans = [0 1 0; -1 0 0; 0 0 1];
V = V_trans'; 

% Define the W matrix
W = [0 -1 0; 1 0 0; 0 0 1];

% Third column of U (which is just the last column of the identity matrix)
U3 = U(:,3);

% Compute the four possible configurations
R1 = U * W * V';
t1 = U3;

R2 = U * W * V';
t2 = -U3;

R3 = U * W' * V';
t3 = U3;

R4 = U * W' * V';
t4 = -U3;

% Output the results
disp('Rotation R1 and Translation t1:');
disp(R1);
disp(t1);

disp('Rotation R2 and Translation t2:');
disp(R2);
disp(t2);

disp('Rotation R3 and Translation t3:');
disp(R3);
disp(t3);

disp('Rotation R4 and Translation t4:');
disp(R4);
disp(t4);



% d) The point pair x4 = (1, 0) and x′4 = (−2, 0) is a real match. Use this 
% to select the real solution from d).

x4 = [1; 0; 1; 1]; 
x4_prime = [-2; 0; 1; 1]; 

P1 = [R1, t1]; % Plausible
P2 = [R2, t2]; 
P3 = [R3, t3]; % Plausible
P4 = [R4, t4]; 

% Check the depth for each configuration
validity = zeros(1,4);
validity(1) = checkDepth(R1, t1, x4, x4_prime);
validity(2) = checkDepth(R2, t2, x4, x4_prime);
validity(3) = checkDepth(R3, t3, x4, x4_prime);
validity(4) = checkDepth(R4, t4, x4, x4_prime);

% Output the valid configurations
disp('Validity of configurations:');
disp(validity);


% e) What kind of motion is the camera undergoing?




%% Task 4
clearvars;close all;

%%% a)

% The camera matrix has 11 variables, but since scale invariance is
% apparent that leaves us with 10 variables. Thus 5 point correspondences
% are required instead of what would have been 6.


%%% b) 

% Paramz
P = 0.9999; % Desired prob
epsilon = 0.25; % Prop of outliers
s = 6; % Minimal set of correspondences needed

% Calculate num RANSAC iterations req
N = log(1 - P) / log(1 - (1 - epsilon)^s);
N_iterations = ceil(N); % Ceiling to round up

fprintf('The number of RANSAC iterations req: %d\n', N_iterations);

save("eliot_task4.mat","N_iterations","s")

%%% c) 




% Solution 1
clearvars;close all;
load('ex4.mat');
load('eliot_task4.mat')

%%% Define "hyperparams"

% A point is considered to be an inlier if its projection 
% is less than 5 pixels from the corresponding image point.
threshold = 5;  
% From a), b)
iterations = N_iterations; 
point_correspondences = s;

% Step 2 & 3: Est camera matrix using RANSAC w. DLT
[P_best, inliers] = RANSAC_CameraMatrix(x, X, iterations, threshold, point_correspondences);

% Display the results
disp('Best Camera Matrix P:');
disp(P_best);

disp('Number of Inliers:');
disp(length(inliers));



% Optionally, visualize the inliers and outliers
figure;
plot(x(1, inliers), x(2, inliers), 'go', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
outliers = setdiff(1:size(x, 2), inliers);
plot(x(1, outliers), x(2, outliers), 'ro', 'LineWidth', 2, 'MarkerSize', 10);
title('Inliers and Outliers in the 2D points');
legend('Inliers', 'Outliers');
xlabel('x coordinate');
ylabel('y coordinate');
axis equal;
grid on;
hold off;



%% Task 5
clearvars;close all;

%%% a)

% Define the camera matrices
P1 = [1 0 0 0; 0 1 0 0; 0 0 1 0];
P2 = [1 0 4 1; 0 1 2 0; 0 0 1 1];
P3 = [-1 1 -2 1; -1 0 -1 0; 0 -2 0 0];

A1 = P1(1:3,1:3);
A2 = P2(1:3,1:3);
A3 = P3(1:3,1:3);

% Compute the epipoles for each camera matrix by solving P * e = 0
C1 = null(P1);
C2 = null(P2);
C3 = null(P3);

C1 = C1./C1(end,end);
C2 = C2./C2(end,end);
C3 = C3./C3(end,end);

e1 = P1*C1;
e2 = P2*C2;
e3 = P3*C3;

% Normalize epipoles
e1 = e1 / e1(end);
%e2 = e2 / e2(end); % Do not normalize
e3 = e3 / e3(end); 

% Convert the epipoles into skew symmetric matrices
e1x = [0 -e1(3) e1(2); e1(3) 0 -e1(1); -e1(2) e1(1) 0];
e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0];
e3x = [0 -e3(3) e3(2); e3(3) 0 -e3(1); -e3(2) e3(1) 0];

F12 = e2x * P2 * pinv(P1);
F13 = e3x * P3 * pinv(P1);
F23 = e3x * P3 * pinv(P2);

disp('Fundamental matrix F12:');
disp(F12);

disp('Fundamental matrix F13:');
disp(F13);

disp('Fundamental matrix F23:');
disp(F23);


%%% Compute the epipoles and verify equations

% Compute epipoles for each fundamental matrix
% The epipole e_ij is the null space of F_ij (right epipole)
% The epipole e_ji is the null space of F_ij' (left epipole)

e12 = null(F12); % Right 
e21 = null(F12'); % Left 

e13 = null(F13); % Right 
e31 = null(F13'); % Left 

e23 = null(F23); % Right 
e32 = null(F23'); % Left

% Normalize 
%e12 = e12 / e12(end);
%e21 = e21 / e21(end);
e13 = e13 / e13(end);
e31 = e31 / e31(end);
e23 = e23 / e23(end);
e32 = e32 / e32(end);

% Verify the equations
eq1 = e23' * F12 * e13;
eq2 = e32' * F13 * e12;
eq3 = e31' * F23 * e21;


disp('Verification of equations:');
disp(['e23'' * F12 * e13 = ', num2str(eq1)]);
disp(['e32'' * F13 * e12 = ', num2str(eq2)]);
disp(['e31'' * F23 * e21 = ', num2str(eq3)]);

tolerance = 1e-10; 
disp('Equations verification results:');
disp(['eq1 ≈ 0: ', num2str(abs(eq1) < tolerance)]);
disp(['eq2 ≈ 0: ', num2str(abs(eq2) < tolerance)]);
disp(['eq3 ≈ 0: ', num2str(abs(eq3) < tolerance)]);



%% Functions 

% Task 3
function isValid = checkDepth(R, t, x, x_prime)
    P_prime = [R, t]; % Projection matrix for the second camera
    
    % Ensure that x is homogeneous
    if size(x, 1) == 3
        x = [x; 1];
    end
   
    % Transform the point x into the second camera's frame
    x_transformed = P_prime * x;
    
    % Check if depth is positive in both cameras
    isValid = (x_transformed(3) > 0) && (x_prime(3) > 0);
end




function P = DLT(x2D, x3D)
    % Normalize the points.
    [x2D_normalized, T] = normalize2Dpoints(x2D);
    [x3D_normalized, U] = normalize3Dpoints(x3D);

    % Construct matrix A for DLT.
    A = [];
    for i = 1:size(x2D_normalized, 2)
        X = x3D_normalized(:, i);
        x = x2D_normalized(1, i);
        y = x2D_normalized(2, i);
        A = [A;
             -X', zeros(1, 4), x*X';
             zeros(1, 4), -X', y*X'];
    end

    % Solve for the homography matrix using SVD.
    [~, ~, V] = svd(A);
    P_normalized = reshape(V(:, end), 4, 3)';

    % Denormalize the projection matrix.
    P = T \ (P_normalized * U);
end

function [P_best, inliers] = RANSAC_CameraMatrix(x2D, x3D, iter, threshold, point_correspondences)
    % Initialize variables
    numPts = size(x2D, 2);
    maxInliers = 0;
    P_best = [];
    
    for i = 1:iter
        % Randomly select a subset of points for estimating P
        subset = randperm(numPts, point_correspondences); % Minimum number of points needed for P estimation
        P_temp = DLT(x2D(:, subset), x3D(:, subset));
        
        % Project 3D points to 2D using the temporary camera matrix
        x2D_proj = P_temp * [x3D; ones(1, numPts)];
        x2D_proj = x2D_proj ./ x2D_proj(3, :); % Convert to inhomogeneous coordinates
        
        % Calculate the reprojection error
        reprojectionError = sqrt(sum((x2D(1:2, :) - x2D_proj(1:2, :)).^2, 1));
        
        % Determine inliers based on the threshold
        inliers_temp = find(reprojectionError < threshold);
        
        % Update the best model if the current one is better
        if length(inliers_temp) > maxInliers
            maxInliers = length(inliers_temp);
            P_best = P_temp;
            inliers = inliers_temp;
        end
    end
end



function [normalized_points, T] = normalize2Dpoints(x2D)
    % Calculate the mean of the 2D points
    mean_points = mean(x2D(1:2,:), 2);
    
    % Calculate the standard deviation of the 2D points
    std_dev = std(x2D(1:2,:), 0, 2);
    
    % Construct the normalization matrix T
    T = [1/std_dev(1)   0           -mean_points(1)/std_dev(1)
         0              1/std_dev(2) -mean_points(2)/std_dev(2)
         0              0             1                     ];
    
    % Apply normalization matrix to the original points
    normalized_points = T * [x2D; ones(1, size(x2D,2))];
end


function [normalized_points, U] = normalize3Dpoints(x3D)
    % Calculate the mean of the 3D points
    mean_points = mean(x3D(1:3,:), 2);
    
    % Calculate the standard deviation of the 3D points
    std_dev = std(x3D(1:3,:), 0, 2);
    
    % Construct the normalization matrix U
    U = [1/std_dev(1)   0             0             -mean_points(1)/std_dev(1)
         0              1/std_dev(2)   0             -mean_points(2)/std_dev(2)
         0              0              1/std_dev(3)  -mean_points(3)/std_dev(3)
         0              0              0              1                       ];
    
    % Apply normalization matrix to the original points
    normalized_points = U * [x3D; ones(1, size(x3D,2))];
end






