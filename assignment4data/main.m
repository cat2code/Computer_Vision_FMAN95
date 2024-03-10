% Assignment 4 in Computer Vision FMAN95
% Written by Eliot Petrén
clearvars;close all;clc;

%% Computer Exercise 1 (CE1) - Robust Homography Estimation and Stitching
clearvars; close all;

% Setup paths and VLFeat
addpath('../assignment4data', 'vlfeat-0.9.21/toolbox');
vl_setup;

% Load and display images
A = imread('a.jpg');
B = imread('b.jpg');
figure, imshow(A, 'InitialMagnification', 'fit'), title('Image A');
figure, imshow(B, 'InitialMagnification', 'fit'), title('Image B');

% Compute SIFT features and matches
imgA_gray = single(rgb2gray(A));
imgB_gray = single(rgb2gray(B));
[fA, dA] = vl_sift(imgA_gray);
[fB, dB] = vl_sift(imgB_gray);
[matches, ~] = vl_ubcmatch(dA, dB);

% Setup matched feature points
xA = [fA(1:2, matches(1,:)); ones(1, size(matches,2))];
xB = [fB(1:2, matches(2,:)); ones(1, size(matches,2))];

% RANSAC for homography estimation
threshold = 5;
H_best = [];
inls_max = 0;
corr_min = 4;
for i=1:30
    r = randperm(size(xA, 2), corr_min);
    H = estimateHomographyRANSAC(xA(:,r), xB(:,r), threshold);
    txA = homographyTransform(xA, H);
    inls = sum(sqrt(sum((txA(1:2,:) - xB(1:2,:)).^2)) < threshold);
    if inls > inls_max
        H_best = H;
        inls_max = inls;
    end 
end

% Prepare for stitching
Htform = projective2d(H_best');
Rout = imref2d(size(A), [-200 800], [-400 600]);
Atransf = imwarp(A, Htform, 'OutputView', Rout);
Btransf = imwarp(B, projective2d(eye(3)), 'OutputView', Rout);

% Stitch images
AB = Btransf;
AB(Btransf < Atransf) = Atransf(Btransf < Atransf);

% Display stitched image
figure, imshow(AB, 'XData', Rout.XWorldLimits, 'YData', Rout.YWorldLimits);
title('Stitched Image');


%% CE2 - Robust Essential Matrix Estimation
clearvars;close all;
% Setup
addpath('../assignment4data')
addpath("vlfeat-0.9.21/toolbox");
load compEx2data.mat
vl_setup;

% (Stereo matching these pictures later)
im1 = imread('im1.jpg');
im2 = imread('im2.jpg');

% Homogenous coord.
x{1} = [x{1}; ones(1, size(x{1},2))];
x{2} = [x{2}; ones(1, size(x{2},2))];

% Define
x1 = x{1};
x2 = x{2};

% K-normalize image points (K is camera calibration matrix)
x1_norm = K^-1 * x1;
x2_norm = K^-1 * x2;


%%% Estimate Essential matrix %%%

% Define RANSAC params: min iterations and min correspondences.
min_iter = 5;
corr_min = 5;

% Init
nbest_inls = 0;
best_inls = [];
best_E = [];

% RANSAC loop finds best estimate of the Essential Matrix (E)
for i=1:min_iter*20
    % Randomly select correspondences
    r = randperm(size(x1_norm, 2), corr_min);
    rx1_norm = x1_norm(:,r);
    rx2_norm = x2_norm(:,r);
    
    % Build E with five-point solver
    E = fivepoint_solver(rx1_norm, rx2_norm);
    
    % Test E solution
    for j=1:size(E,2)
        % Build Fundamental Matrix (F)
        F = (K^-1)'*E{j}*K^-1;
        
        % Epipolar lines in the images
        l1 = pflat(F'*x2);
        l1 = l1 ./ sqrt(repmat(l1(1, :).^2 + l1(2, :).^2 ,[3 1]));
        l2 = pflat(F*x1);
        l2 = l2 ./ sqrt(repmat(l2(1, :).^2 + l2(2, :).^2 ,[3 1]));
        
        % Epipolar distances in the images
        d1 = abs(sum(l1.*x1)); 
        d2 = abs(sum(l2.*x2)); 
        
        % Find inliers and count
        inls = (d1 < 5) & (d2 < 5);
        nbr_inls = sum(inls(:));

        % Update with the best E and number of inliers
        if nbr_inls > nbest_inls
            best_E = E{j};
            nbest_inls = nbr_inls;
            best_inls = inls;
        end 
    end 
end
disp("number of inliers = " + nbest_inls);

% Extract camera matrices
P1 = [eye(3) zeros(3, 1)];
P2 = E_to_P(best_E);

% Ensure points are in front of the cameras.
[P2_best, X_best] = selectVisible3DPointsDLT(P1, P2, x1_norm, x2_norm);

% K-Unnormalize
P1 = K * P1;
P2_best = K * P2_best;

% Normalize against homogenous so that they are 1
X_best = pflat(X_best);

% Filter matches and 3D points to only include inliers.
x = {x{1}(:,best_inls == 1), x{2}(:, best_inls == 1)};
X_best = X_best(:, best_inls == 1);

% Compute reprojection errors and RMS
[err, res] = ComputeReprojectionError({P1, P2_best}, X_best, x);
RMS = sqrt(err / size(res, 2));

% Plot error in hist
figure()
hist(res,100)

% Plot 3D points and camera centers
figure()
plot3(X_best(1,:), X_best(2,:), X_best(3,:), '.', 'Markersize', 2.5)
axis equal
hold on
plotcams({P1, P2_best})

% Save results
P2 = P2_best;
X = X_best;
save('ce2_results', 'P1', 'P2', 'X', 'x');



%% CE3 - Calibrated Structure from Motion and Local Optimization
clearvars;close all;
% Setup
addpath('../assignment4data')
load ce2_results.mat

% Step size in gradient descent
gamma_k = 10^-10;

% Takes two camera matrices and puts them in a cell.
P = {P1,P2};

% New variable names
U = X;
u = x;

% Computes the reprejection error and the values of all the residuals
% for the current solution P,U,u.
[err,res] = ComputeReprojectionError(P,U,u);

% New variable names
Pnew = P;
Unew = U;

% Plot the initial objective func value (sum of residuals).
figure();
plot(0, sum(res), 'r*')
hold on

% Perform gradient descent 10 times
for i=1:10
    % Linearize the reprojection error around the current solution to 
    % obtain residuals (r) and Jacobian (J).
    [r, J] = LinearizeReprojErr(Pnew, Unew, u);
    
    % Computes the LM update.
    deltav = -gamma_k*J'*r;
    
    % Updates the variabels
    [Pnew, Unew] = update_solution(deltav, Pnew, Unew);
    
    % Recompute the reprojection error and residuals
    [err,res] = ComputeReprojectionError(Pnew,Unew,u);
    
    % Plot the updated objective func value
    plot(i, sum(res), 'r*')
end 

% Compute the final and (hopefully) minimized RMS
RMS = sqrt(err / size(res, 2));
disp(['Calculated RMS CE3: ', RMS])

% Done wahoo
%% CE4 - Calibrated Structure from Motion and Local Optimization
clearvars;close all;
% Setup
addpath('../assignment4data')
%load compEx2data.mat
%load compEx4data.mat
load ce2_results.mat

%gamma_k = 10^-9; % never used

% damping factor in LM algorithm Gauss-Newton something something
lambda = 0.1;

% Takes two camera matrices and puts them in a cell.
P = {P1, P2};

% New variable names
U = X;
u = x;

% Computes the reprejection error and the values of all the residuals
% for the current solution P,U,u.
[err, res] = ComputeReprojectionError(P, U, u);

% New variable names
Pnew = P;
Unew = U;

% Plot the initial objective func value (sum of residuals).
figure();
plot(0, sum(res), 'r*')
hold on

% Levenberg-Marquardt method and perform 10 updates
for i=1:10
    % Linearize the reprojection error around the current solution to 
    % obtain residuals (r) and Jacobian (J).
    [r,J] = LinearizeReprojErr(Pnew,Unew,u);
    
    % Calculate LM update in Gauss-Newton with damping factor 
    C = J'*J + lambda*speye(size(J,2));
    c = J'*r;
    deltav = -C\c;
    
    % Updates the variables
    [Pnew, Unew] = update_solution(deltav,Pnew,Unew);
    
    % Recompute the reprojection error and residuals
    [err,res] = ComputeReprojectionError(Pnew,Unew,u);
    
    % Plot the updated objective func value
    plot(i, sum(res), 'r*')
end 

% Compute the final and (hopefully) minimized RMS
RMS = sqrt(err/size(res,2));
disp(['Calculated RMS CE4: ', RMS])

% DONE again wahoo
%% CE5 - Dynamic Objects and Factorization
clearvars;close all;
% Setup
%load compEx2data.mat
load compEx4data.mat

% Best approx with k = 7
k = 7;

% Draw the first face in A
drawface(A(1,:), A(2,:)); 

% Compute the mean over the rows of A.
m = mean(A,2);

% Subtract row means from A and compute the SVD.
[U,S,V] = svd(A - m*ones(1,size(A,2)),'econ');

% Use only the first k singular values to approximate A.
Uk = U(:,1:k);
Sk = S(1:k,1:k);
Vk = V(:,1:k);
X = Uk * Sk * Vk' + m*ones(1, size(A,2));

% Plot the original and the estimated faces from image 1
figure;
subplot(1,2,1);
drawface(A(1,:), A(2,:));
title('Original Face');
subplot(1,2,2);
drawface(X(1,:), X(2,:));
title('Approximated Face');

% 2. Determine the shape basis Bk
Bk = Uk * sqrt(Sk);

% 3. Coefficient Determination
% If we had a new image (new_points), we would compute its coefficients as follows:
% new_points_normalized = new_points - m;
% cik = Bk \ new_points_normalized;

% 4. Handle Noisy Data
figure;
for i = 1:size(Anoise,2)/2
    subplot(1,3,i);
    drawface(Anoise(2*i-1,:), Anoise(2*i,:));
    title(['Noisy Face ', num2str(i)]);
end

% 5. Handle Missing Data
Am = Amissing;
% Remove rows with NaN and solve for coefficients
for i = 1:size(Am,2)/2
    current_points = [Am(2*i-1,:); Am(2*i,:)];
    missing_indices = any(isnan(current_points));
    Bk_reduced = Bk(~missing_indices,:);
    current_points_reduced = current_points(~missing_indices);
    cik = Bk_reduced \ current_points_reduced;
    % Compute full points
    reconstructed_points = Bk * cik;
    figure;
    drawface(reconstructed_points(1:2:end), reconstructed_points(2:2:end));
    title(['Reconstructed Face ', num2str(i)]);
end

%% ver 2.0 CE5 - Dynamic Objects and Factorization
clearvars;close all;
% Setup
%load compEx2data.mat
load compEx4data.mat

% 1. Best Approximation with k = 7
k = 7;

% Draw the first face in A
figure;
subplot(1,2,1);
drawface(A(1,:), A(2,:)); 
title('Original Face');

% Compute the mean over the rows of A.
m = mean(A,2);

% Subtract row means from A and compute the SVD.
[U,S,V] = svd(A - m*ones(1,size(A,2)),'econ');

% Use only the first k singular values to approximate A.
Uk = U(:,1:k);
Sk = S(1:k,1:k);
Vk = V(:,1:k);
X = Uk * Sk * Vk' + m*ones(1, size(A,2));

% Plot the estimated faces from image 1
subplot(1,2,2);
drawface(X(1,:), X(2,:));
title('Approximated Face');

% 2. Determine the shape basis Bk
Bk = Uk * sqrt(Sk);

% Since this code does not include the unseen new points, the coefficients cik cannot be computed.
% We assume that this part will be run with the actual new points when available.

% 4. Handle Noisy Data
% Plotting faces with noise. We assume Anoise contains pairs of (x,y) coordinates for each face.
for i = 1:3
    figure
    %subplot(1,3,i);
    drawface(Anoise(2*i-1,:), Anoise(2*i,:));
    title(['Noisy Face ', num2str(i)]);
end

% 5. Handle Missing Data
% Plotting faces with missing data, assuming Amissing is structured similar to Anoise.
Am = Amissing;
for i = 1:3
    disp(i)
    current_points = [Am(2*i-1,:); Am(2*i,:)]; 
    missing_indices = any(isnan(current_points), 2); 
    current_points_reduced = current_points(~missing_indices);
    current_points_reduced = current_points_reduced(:); 
    Bk_reduced = Bk(~missing_indices(:),:); % Ensure Bk_reduced matches non-missing x and y indices
    
    cik = Bk_reduced \ current_points_reduced; % Solve for coefficients
    
    % Reshape the result of Bk * cik to 2 rows for x and y coordinates, respectively
    reconstructed_points = reshape(Bk * cik, 2, []) + reshape(m,2,114);
    figure;
    % Use the first row for x and the second row for y when calling drawface
    drawface(reconstructed_points(1,:), reconstructed_points(2,:));
    title(['Reconstructed Face ', num2str(i)]);
end




%% Functions:


