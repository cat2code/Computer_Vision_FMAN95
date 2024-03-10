% Clear workspace and close all figures
clear
close all

% Add paths to external datasets and VLFeat toolbox
addpath('../assignment4data')
addpath("../../vlfeat-0.9.21/toolbox");
vl_setup;

% Read two images from files
A = imread('a.jpg');
B = imread('b.jpg');

% Compute SIFT features for both images. 'fA' and 'fB' contain the feature
% frames and 'dA' and 'dB' contain the descriptors.
[fA, dA] = vl_sift(single(rgb2gray(A)));
[fB, dB] = vl_sift(single(rgb2gray(B)));

% Find matching features between the two sets of descriptors.
[matches, scores] = vl_ubcmatch(dA, dB);

% Prepare matched feature points for estimating homography, adding a row
% of ones for homogeneous coordinates.
xA = [fA(1:2, matches(1,:)); ones(1, size(matches,2))];
xB = [fB(1:2, matches(2,:)); ones(1, size(matches,2))];

% Initialize variables for RANSAC.
best_H = [];
most_inls = 0;
min_corr = 4; % Minimum number of correspondences

% RANSAC loop to estimate the best homography
for i=1:30
    % Randomly select 'min_corr' correspondences
    r = randperm(size(xA, 2), min_corr);
    rxA = xA(:,r);
    rxB = xB(:,r);
    
    % Initialize the matrix for DLT
    M = zeros(3*size(rxA,2), 3*size(rxA,1)+min_corr);
    
    % Build the matrix M for DLT
    row = 1;
    col = 3*size(rxA,1)+1;
    for j=1:size(rxA, 2)
        m = [rxA(:,j)'    zeros(1,6); 
             zeros(1,3)  rxA(:,j)' zeros(1,3); 
             zeros(1,6)  rxA(:,j)'];
        M(row:row+2, 1:9) = m;
        M(row:row+2, col) = -rxB(:,j);
        
        row = row + 3;
        col = col + 1;
    end
    
    % Solve for homography using SVD
    [U, S, V] = svd(M);
    v = V(:, end);
    H = reshape(v(1:9), 3, 3)';
    
    % Apply the estimated homography and normalize
    txA = zeros(size(xA));
    for j=1:size(xA,2)
        txA(:,j) = H*xA(:,j);
        txA(:,j) = txA(:,j) ./ txA(3,j); % Normalize
    end 
    
    % Count inliers based on a distance threshold
    inls = sqrt(sum((txA(1:2,:) - xB(1:2,:)).^2))<5;
    nbr_inls = sum(inls);
    if nbr_inls > most_inls
        best_H = H;
        most_inls = nbr_inls;
    end 
end

% Convert the best homography matrix to a format compatible with MATLAB's imwarp function
Htform = projective2d(best_H');

% Define the output view for the stitched images
Rout = imref2d(size(A), [-200 800], [-400 600]);

% Warp image A using the estimated homography
[Atransf] = imwarp(A, Htform, 'OutputView', Rout);

% Warp image B using an identity transformation (i.e., no change)
Idtform = projective2d (eye(3));
[Btransf] = imwarp(B, Idtform, 'OutputView', Rout);

% Combine the two transformed images into a single panorama
AB = Btransf;
AB(Btransf < Atransf) = Atransf(Btransf < Atransf);

% Display the resulting stitched image
imagesc(Rout.XWorldLimits, Rout.YWorldLimits ,AB);
