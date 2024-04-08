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



%% Outline

% Plotting
figure;
plot3(points3D(:,1), points3D(:,2), points3D(:,3), '.');
hold on;
camera_center = -inv(P(:,1:3)) * P(:,4);
principal_axis = P(3, 1:3);
quiver3(camera_center(1), camera_center(2), camera_center(3), principal_axis(1), principal_axis(2), principal_axis(3), 50);
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Model with Camera Center and Principal Axis');
grid on;
hold off;



%% Functions
%%% Solution 1

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

function [P_best, inliers] = RANSAC_CameraMatrix2(x2D, x3D, iterations, threshold)
    numPoints = size(x2D, 2);
    P_best = [];
    inliers_best = [];

    for i = 1:iterations
        % Randomly select points for DLT.
        indices = randperm(numPoints, 6); % Minimum number for P matrix estimation.
        P = DLT(x2D(:, indices), x3D(:, indices));

        % Compute inliers where reprojection error is below threshold.
        x2D_projected = P * [x3D; ones(1, numPoints)];
        x2D_projected = x2D_projected ./ x2D_projected(3, :);
        
        
        errors = sqrt(sum((x2D(1:2, :) - x2D_projected(1:2, :)).^2, 1));
        inliers = find(errors < threshold);

        if length(inliers) > length(inliers_best)
            P_best = P;
            inliers_best = inliers;
        end
    end
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




%%% Old stuff
function P = DLT2(points2D, points3D)
    % Assuming points2D and points3D are Nx2 and Nx3 matrices respectively.
    
    N = size(points2D, 1);
    A = [];
    
    for i = 1:N
        X = points3D(i, :);
        x = points2D(i, :);
        A = [A; 
             x(1) * [0, 0, 0, -1, -X, X(2)];
             x(2) * [0, 0, 0, -1, -X, X(1)]];
    end
    
    [~, ~, V] = svd(A);
    P = reshape(V(:, end), 4, 3)';
end

function [P, inliers] = ransacDLT(points2D, points3D, iterations, threshold)
    numPoints = size(points2D, 1);
    maxInliers = [];
    bestP = [];
    
    for i = 1:iterations
        subset = randperm(numPoints, 5);
        P_est = DLT(points2D(subset, :), points3D(subset, :));
        
        inliers = [];
        for j = 1:numPoints
            projected = P_est * [points3D(j, :)'; 1];
            projected = projected ./ projected(3);
            if norm(projected(1:2) - points2D(j, :)') < threshold
                inliers = [inliers, j];
            end
        end
        
        if length(inliers) > length(maxInliers)
            maxInliers = inliers;
            bestP = P_est;
        end
    end
    
    % Recompute P using all inliers
    P = DLT(points2D(maxInliers, :), points3D(maxInliers, :));
    inliers = maxInliers;
end