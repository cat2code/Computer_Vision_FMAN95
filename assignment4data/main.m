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

disp(size(fA)) % 947 matched features
disp(size(fB)) % 865 matched features
disp(size(matches)) % 204 matches

% Setup matched feature points
xA = [fA(1:2, matches(1,:)); ones(1, size(matches,2))];
xB = [fB(1:2, matches(2,:)); ones(1, size(matches,2))];

% RANSAC for homography estimation
threshold = 5;
H_best = [];
inls_max = 0;
corr_min = 4;
best_inlier_mask = [];
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
disp(['Found inliers: ', num2str(inls_max)])

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


%% CE2 - Robust Essential Matrix Estimation (STABLE NOW)
clearvars;close all;
% Setup
addpath('../assignment4data', 'vlfeat-0.9.21/toolbox');
load compEx2data.mat
vl_setup;

% (Stereo matching these pictures later)
im1 = imread('im1.jpg');
im2 = imread('im2.jpg');
figure, imshow(im1, 'InitialMagnification', 'fit'), title('Image 1');
figure, imshow(im2, 'InitialMagnification', 'fit'), title('Image 2');

% Homogenous coord.
x{1} = [x{1}; ones(1, size(x{1},2))];
x{2} = [x{2}; ones(1, size(x{2},2))];

% Define
x1 = x{1};
x2 = x{2};

% K-normalize image points (K is camera calibration matrix)
x1_n = K^-1 * x1;
x2_n = K^-1 * x2;

%%% Estimate Essential matrix %%%

% Init
corr_min = 5;
nbest_inls = 0;
best_inls = [];
best_E = [];
threshold = 5;

% RANSAC loop finds best estimate of the Essential Matrix (E)
for i=1:100
    % Randomly select correspondences
    r = randperm(size(x1_n, 2), corr_min);
    rx1_n = x1_n(:,r);
    rx2_n = x2_n(:,r);
    
    % Build E with five-point solver
    E = fivepoint_solver(rx1_n, rx2_n);
    
    % Test E solution
    for j=1:size(E,2)
        % Build Fundamental Matrix (F)
        F = (K^-1)' * E{j} * K^-1;
        
        % Epipolar lines in the images
        l1 = pflat(F' * x2, 0);
        l1 = l1 ./ sqrt(repmat(l1(1, :).^2 + l1(2, :).^2 ,[3 1]));
        l2 = pflat(F * x1, 0);
        l2 = l2 ./ sqrt(repmat(l2(1, :).^2 + l2(2, :).^2 ,[3 1]));
        
        % Epipolar distances in the images
        d1 = abs(sum(l1.*x1)); 
        d2 = abs(sum(l2.*x2)); 
        
        % Find inliers and count
        inls = (d1 < threshold) & (d2 < threshold);
        nbr_inls = sum(inls(:));

        % Update with the best E and number of inliers
        if nbr_inls > nbest_inls
            best_E = E{j};
            nbest_inls = nbr_inls;
            best_inls = inls;
        end 
    end 
end
disp("Found total of inliers = " + nbest_inls);

%%% Best inliers and their indices found, continue

P1 = K*[eye(3) zeros(3,1)];

% Extracting possible cameras from the essential matrix E
P2 = E_to_P(best_E);

% Choose solution (camera) with the most points in front of them
[P2_best, X] = selectVisible3DPointsDLT(P2, x, x1_n, x2_n);

% K-transform P_best
P2_best = K * P2_best;

% Prepare for plots
x = {x{1}(:,best_inls==1),x{2}(:,best_inls==1)};
X = pflat(X,0);
X = X(:,best_inls==1);

% Compute reprojection errors and RMS
[err, res] = ComputeReprojectionError({P1, P2_best}, X, x);
RMS = sqrt(err / size(res, 2));
disp(['RMS: ', num2str(RMS)])

figure()
hist(res,100)
title('Error histogram')

figure()
plot3(X(1,:),X(2,:),X(3,:),'.','Markersize',3)
axis equal
hold on
plotcams({P2_best,P1})

P2 = P2_best;

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
title('Objective value vs iteration')
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
disp(['Calculated RMS CE3: ', num2str(RMS)])

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
title('Objective value vs iteration')
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
disp(['Calculated RMS CE4: ', num2str(RMS)])

% DONE again wahoo

%% ver 2.0 CE5 - Dynamic Objects and Factorization (first draft not done)
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

Uk = U(:,1:k);
Sk = S(1:k,1:k);
Vk = V(:,1:k);
X = Uk * Sk * Vk' + m*ones(1, size(A,2));

% Plot face from image 1
subplot(1,2,2);
drawface(X(1,:), X(2,:));
title('Approximated Face');

Bk = Uk;  % The shape basis is just Uk

for i = 1:3
    figure
    %subplot(1,3,i);
    drawface(Anoise(2*i-1,:), Anoise(2*i,:));
    title(['Noisy Face ', num2str(i)]);
end

Am = Amissing;
for i = 1:3
    current_points = [Am(2*i-1,:); Am(2*i,:)]; 
    missing_indices = isnan(current_points); 
    current_points_reduced = current_points(~missing_indices);
    
    % Adjust Bk and m to match the non-missing data
    Bk_reduced = Bk(~missing_indices,:);
    m_reduced = m(~missing_indices);
    
    
    cik = Bk_reduced \ current_points_reduced;
    
    % Reconstruct the full set of points
    reconstructed_points = Bk * cik;
    reconstructed_points(missing_indices) = NaN;  % Mark the missing points as NaN
    
    % Adjust the mean back to the reconstructed points
    reconstructed_points = reconstructed_points + m;
    
    % Reshape for plotting
    reconstructed_points = reshape(reconstructed_points, 2, []);
    
    figure;
    drawface(reconstructed_points(1,:), reconstructed_points(2,:));
    title(['Reconstructed Face ', num2str(i)]);
end


%% Functions:

function H = estimateHomographyRANSAC(xA, xB, threshold)
    % DLT Homography estimation simplified for 4 points
    A = [];
    for i = 1:size(xA,2)
        xs = xA(1,i); ys = xA(2,i);
        xd = xB(1,i); yd = xB(2,i);
        A = [A; -xs, -ys, -1, 0, 0, 0, xs*xd, ys*xd, xd;
                0, 0, 0, -xs, -ys, -1, xs*yd, ys*yd, yd];
    end
    [~, ~, V] = svd(A);
    H = reshape(V(:,end), 3, 3)';
end

function txA = homographyTransform(xA, H)
    % Apply homography transform to points xA
    txA = H * xA;
    txA = txA ./ txA(3,:);
end
