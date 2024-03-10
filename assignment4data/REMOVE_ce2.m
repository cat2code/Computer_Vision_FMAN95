% Clear the workspace and close all figures to ensure a clean environment.
clear
close all

% Add paths to external datasets and the VLFeat toolbox for feature extraction.
addpath('../assignment4data')
load('compEx2data.mat') % Load variables for the computation, including 'K' and point matches.
addpath("../../vlfeat-0.9.21/toolbox");
vl_setup; % Setup VLFeat library.

% Read the two images for stereo matching.
im1 = imread('im1.jpg');
im2 = imread('im2.jpg');

% Append a row of ones to the homogeneous coordinates of matched points in both images.
x{1} = [x{1}; ones(1, size(x{1},2))];
x{2} = [x{2}; ones(1, size(x{2},2))];

% Assign the matches to more convenient variables.
x1 = x{1};
x2 = x{2};

% Normalize the points using the inverse of the camera calibration matrix.
x1_norm = K^-1 * x1;
x2_norm = K^-1 * x2;

% Define RANSAC parameters: minimum iterations and correspondences.
min_iter = 5;
min_corr = 5;

% Initialize variables to track the best inliers and Essential matrix.
nbest_inls = 0;
best_inls = [];
best_E = [];

% RANSAC loop to robustly estimate the Essential matrix.
for i=1:min_iter*20
    % Randomly select correspondences.
    r = randperm(size(x1_norm, 2), min_corr);
    rx1_norm = x1_norm(:,r);
    rx2_norm = x2_norm(:,r);
    
    % Solve for Essential matrices using the five-point solver.
    E = fivepoint_solver(rx1_norm, rx2_norm);
    
    % Test each Essential matrix solution.
    for j=1:size(E,2)
        % Convert Essential matrix to Fundamental matrix using calibration matrix.
        F = (K^-1)'*E{j}*K^-1;
        
        % Compute epipolar lines in both images.
        l1 = pflat(F'*x2); % For image 1
        l1 = l1 ./ sqrt(repmat(l1(1, :).^2 + l1(2, :).^2 ,[3 1]));
        l2 = pflat(F*x1); % For image 2
        l2 = l2 ./ sqrt(repmat(l2(1, :).^2 + l2(2, :).^2 ,[3 1]));
        
        % Compute distances to epipolar lines in both images.
        d1 = abs(sum(l1.*x1)); % Distance for image 1 points
        d2 = abs(sum(l2.*x2)); % Distance for image 2 points
        
        % Determine inliers based on the distance threshold.
        inls = (d1 < 5) & (d2 < 5);
        
        % Update the best Essential matrix and inliers if this model is better.
        nbr_inls = sum(inls(:));
        if nbr_inls > nbest_inls
            best_E = E{j};
            nbest_inls = nbr_inls;
            best_inls = inls;
        end 
    end 
end
disp("number of inliers = " + nbest_inls);

% From the best Essential matrix, extract camera matrices and triangulate points.
P1 = [eye(3) zeros(3, 1)]; % First camera matrix (canonical form).
P2 = create_cm(best_E); % Function to extract second camera matrix from the Essential matrix.

% Refine camera matrix and 3D points to ensure points are in front of both cameras.
[P2_best, X_best] = dlt_infront(P1, P2, x1_norm, x2_norm); % Triangulate and select points in front.

% Convert camera matrices to unnormalized form using calibration matrix.
P1 = K*P1;
P2_best = K*P2_best;

% Flatten the homogeneous coordinates of 3D points.
X_best = pflat(X_best);

% Filter matches and 3D points to only include inliers.
x ={x{1}(:,best_inls==1),x{2}(:,best_inls==1)};
X_best=X_best(:,best_inls==1);

% Compute and plot reprojection errors and RMS error.
[err, res] = ComputeReprojectionError({P1, P2_best}, X_best, x);
RMS = sqrt(err/size(res,2));
figure();
hist(res,100);

% Plot 3D points and camera positions.
figure();
plot3(X_best(1,:), X_best(2,:), X_best(3,:), '.', 'Markersize', 2.5);
axis equal;
hold on;
plotcams({P1, P2_best});

% Save variables for future use.
P2 = P2_best;
X = X_best;
save('ce2_variables', 'P1', 'P2', 'X', 'x');
