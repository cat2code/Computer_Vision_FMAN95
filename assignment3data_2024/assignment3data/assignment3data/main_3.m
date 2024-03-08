% Written by Eliot M.P el6183mo-s

%% 2. 2 The Fundamental Matrix
clearvars;close all;

%%% Computer Exercise 1 %%%
load compEx1data.mat

%%%% SET UP %%%%

% plot images 
pic1 = "kronan1.jpg";
figure('Name',pic1)
im1 = imread(pic1); %Reads the imagefile with name in imfiles{i}
imagesc(im1)

pic2 = "kronan2.jpg";
figure('Name',pic2)
im2 = imread(pic2); %Reads the imagefile with name in imfiles{i}
imagesc(im2)
%%% --- %%%

% Declare the matching points
x1 = x{1};
x2Knorm = x{2};

% Original x1 and x2
x1_og = x1;     
x2_og = x2Knorm;
%%%% END SET UP %%%%


%%% NORMALIZE %%%

for i = 1:2
    sigma(:,i) = std(x{i}(1:2,:),0,2); %Standard deviation of the 1st and 2nd rows ow x{i}
    tilde(:,i) = mean(x{i}(1:2,:),2);  %Computes the mean value of the 1st and 2nd rows ow x{i}
end

% Declare N
N1 = [1/sigma(1,1) 0 -tilde(1,1)/sigma(1,1); 0 1/sigma(2,1) -tilde(2,1)/sigma(2,1); 0 0 1];
N2 = [1/sigma(1,2) 0 -tilde(1,2)/sigma(1,2); 0 1/sigma(2,2) -tilde(2,2)/sigma(2,2); 0 0 1];

% Normalize x1 and x2
x1 = N1*x1;
x2Knorm = N2*x2Knorm;

% Save normalized x1, x2
x1Knorm = x1;
x2Knorm = x2Knorm;
%%% END NORMALIZE %%% 


%%% SET UP M-matrix %%%
M = []; % Initialize M matrix
for i = 1:length(x1Knorm)
    row = [x2Knorm(1,i)*x1Knorm(:,i)', x2Knorm(2,i)*x1Knorm(:,i)', x2Knorm(3,i)*x1Knorm(:,i)'];
    M = [M; row];
end


% Solve the homogeneous least squares system using SVD
[U,S,V] = svd(M); 

S_squared = S'*S;

% Check that minimum singular value is small
smallest_egeinval_S = S_squared(length(S_squared),length(S_squared));

% Solution to SVD problem
V_star = V(:,length(V));    

% Calculate ||Mv|| and check that it is small
Mv = norm(M*V_star);

% Display check if they are small
T1 = table(smallest_egeinval_S);
T2 = table(Mv);
disp(T1); 
disp(T2);

%%% Answer
% The min S is indeed small and so is ||Mv|| as well. This indicates that
% the method is working as intended and I can proceed.


%%% FORM THE NORMALIZED F-MATRIX %%%

% Construct approx F hat 
Fhat = reshape(V_star(1:9), [3,3])';

% Exact Fhat
[U,S,V] = svd(Fhat);

% Set minimum S to zero
S(3,3) = 0;

% Construct the normalized fundamental matrix from the solution v
Fn = U*S*V';

% Check if Fn is approximately zero
disp("Is Fn approximately zero?")
tolerance = 1e-6; % Define the tolerance for approximation
check = det(Fn) < tolerance;
if check
    disp("Fn is approximately zero (below 1e-6)")
else
    disp("Fn is not approximately zero: " + det(Fn))
end

epipolar_constraints=[];
for k=1:1:length(x1Knorm)
    epipolar_constraints=[epipolar_constraints;x2Knorm(:,k)'*Fn*x1Knorm(:,k)];
end
disp('max epipolar constr: ')
max(epipolar_constraints) 
disp('min epipolar constr: ')
mean(epipolar_constraints)

% Answer: The maximum and minimum epipolar constraints are small and we
% find the epipolar constraint sufficiently satisfied (they 
% are sufficiently close to zero)
%%% END FORM THE NORMALIZED F-MATRIX %%%



%%% Transform back to un-normalized F %%%
F = N2'*Fn*N1;

%Computes and plots all the epipolar constraints (should be roughly 0)
l = F*x{1}; %Computes the epipolar lines
l = l./sqrt(repmat(l(1,:).^2 + l(2,:).^2,[3 1]));

lnn = Fn*x{1}; %Computes the epipolar lines
lnn = lnn./sqrt(repmat(lnn(1,:).^2 + lnn(2,:).^2,[3 1]));

% Pick 20 random points in the second image
numPoints = 20;
randIndices = randperm(size(x{2}, 2), numPoints);
randPoints = x{2}(:, randIndices);

% Plot the image and randomly selected points
figure;
imagesc(im1); % Assuming pic2 is the second image
hold on;
plot(randPoints(1, :), randPoints(2, :), 'r*', MarkerSize = 7);

for i = 1:numPoints
    rital(l(:, randIndices(i)));
end

%
%Makes sure that the line has a unit normal
%(makes the distance formula easier) 

% Compute epipolar distances
l_dist = abs(sum(l.*x{2}));
lnn_dist = abs(sum(lnn.*x{2}));

figure()
hist(l_dist,100);
title('Histogram 100 bins of epipolar distances')
xlabel('Distance')
ylabel('Number of occurances')
%Computes all the the distances between the
%and there corresponding lines, and plots in a histogram

%%% Answer
F = F./F(3,3);
disp('Here come F: ')
disp(F)
save('F_matrix.mat', 'F');

mean_epipolar_dist_F = mean(l_dist);
disp('mean_epipolar_dist_F: ')
disp(mean_epipolar_dist_F)
save('mean_epipolar_dist_F_not_normalized.mat','mean_epipolar_dist_F')

mean_epipolar_dist_Fn = mean(lnn_dist);
disp('mean_epipolar_dist_F_normalized')
disp(mean_epipolar_dist_Fn )
save('mean_epipolar_dist_F_normalized.mat','mean_epipolar_dist_Fn')

% Necessary values for the next task
save('for_next_ex.mat','N1','N2',"x1_og","x2_og","x1Knorm","x2Knorm","x1","x2Knorm")




%% 2. 2 The Fundamental Matrix
clearvars;close all;
%%% Computer Exercise 2 %%%

load for_next_ex.mat
load F_matrix.mat

% Compute camera matrices P1 and P2
P1 = [eye(3,3), zeros(3,1)];
% Compute epipole
e2 = null(F');
% Cross prod matrix
e2x = [0, -e2(3), e2(2); e2(3), 0, -e2(1); -e2(2), e2(1), 0];
P2 = [e2x*F, e2];
P2 = P2./P2(3,4);

% Normalization is regained from last exercise
% x1n, with N1
% x2n, with N2

%%% Triangulation with DLT %%%

% Normalized stereo pair for triangular
normTriangMats = [N1*P1; N2*P2];

% Initialize variables
X = []; % 3D points
%d = [];
absMV = [];

for i = 1:length(x1Knorm)
    xh = [x1Knorm(:,i) zeros(3,1); zeros(3,1) x2Knorm(:,i)];
    M = [normTriangMats -xh];
    [~, ~, V] = svd(M);
    absMV = [absMV; norm(M*V(1:6,end))];
    X(:,i) = V(1:4,end);
    %d(:,i) = V(5:6,end);
end

% Check if Mv is approximately zero
disp("Is Mv approx zero?")
tolerance = 1e-5; % Define the tolerance for approximation
check = max(absMV) < tolerance;
if check
    disp("Mv is approximately zero (below 1e-5)")
else
    disp("Mv is not approximately zero: " + check + max(absMV))
end
%%% Traingulation with DLT END %%%

% Normalize on homogenous
%X = X ./ X(4, :);


%%% Plot: Image, Image Points, Projected 3D points in same fig %%% 

% Project 3D world points onto camera
P1img = pflat((P1 * X));
P2img = pflat((P2 * X));

%%% Plot: Image, Image Points, Projected 3D points in same fig END%%% 

% Prep images
pic1 = "kronan1.jpg";
im1 = imread(pic1); %Reads the image file with name in imfiles{i}
pic2 = "kronan2.jpg";
im2 = imread(pic2); %Reads the image file with name in imfiles{i}

% Project 3D world points onto camera
P1img = pflat(P1 * X);
P2img = pflat(P2 * X);

%%% Plot: Image 1, Image Points, and Projected 3D points %%%
figure()
imagesc(im1);
hold on
plot(x1_og(1,:), x1_og(2,:), 'bo', 'DisplayName', 'Original Points')
plot(P1img(1,:), P1img(2,:), 'r.', 'DisplayName', 'Projected Points')
hold off
title('Image 1: Original and Projected Points')
legend('show')
legend('Location', 'best') % Adjusts the location of the legend to the best location

%%% Plot: Image 2, Image Points, and Projected 3D points %%%
figure()
imagesc(im2);
hold on
plot(x2_og(1,:), x2_og(2,:), 'bo', 'DisplayName', 'Original Points')
plot(P2img(1,:), P2img(2,:), 'r.', 'DisplayName', 'Projected Points')
hold off
title('Image 2: Original and Projected Points')
legend('show')
legend('Location', 'best') % Adjusts the location of the legend to the best location

%%% 3D Reconstruction Visualization %%%
XX = pflat(X);

figure()
plot3(XX(1,:), XX(2,:), XX(3,:), '.');
title('3D Reconstruction')
%%% END %%% 

%% 3 The Essential Matrix
clearvars;close all;clc;

%%% Computer Exercise 3 %%%
load compEx3data.mat;
load for_next_ex.mat

% K-Normalize the image points
x1Knorm = K\x1_og; % Inverse of K times x1_og
x2Knorm = K\x2_og;


%%% SET UP M MATRIX USING SVD 8 POINT ALGO %%%

M=[];
%for k=1:1:length(x1n)
 %   row=[(x1Knorm(:,k)')*x2Knorm(1,k),(x1Knorm(:,k)')*x2Knorm(2,k),(x1Knorm(:,k)')*x2Knorm(3,k)]; 
 %   M=[M;row];
%end

% Append rows to M
M=[];
for k = 1:length(x1Knorm)
    x1k = x1Knorm(:,k)';
    x2k = x2Knorm(:,k);
    row = [x1k*x2k(1), x1k*x2k(2), x1k*x2k(3)]; 
    M = [M; row];
end

% Perform Singular Value Decomposition on M
[~, S, V] = svd(M);

%%% Check that the minimum singular value and Mv are both small. %%%
% Extract the minimum eigenvalue and the corresponding eigenvector
smin = S(end, end)^2;
Vmin = V(:, end)
disp(min(abs(Vmin)))
% Compute the norm of M times the eigenvector corresponding to the minimum eigenvalue
Mv = norm(M * Vmin)

% Both Vmin and Mv are small
%%% END CHECK SMALL %%% 


%%% Construct the Essential matrix from the solution v %%%

% make sure that E has two equal singular values and the third one zero

% Check that the epipolar constraints x ̃T2 Ex ̃1 = 0 are roughly fulfilled.

% init
E = zeros(3,3); 

% Reshape the smallest eigenvector to a 3x3 matrix and transpose
Etilde = reshape(Vmin(1:9), [3,3])';

% Perform SVD on Etilde to extract E
[U, ~, V] = svd(Etilde); 

% Construct E and choose positive
if det(U * V') > 0
    E = U * diag([1, 1, 0]) * V';
else
    V = -V;
    E = U * diag([1, 1, 0]) * V';
end

% Check the epipolar constraint for each point
epipolarConstraints = arrayfun(@(k) x2Knorm(1:3,k)' * E * x1Knorm(1:3,k), 1:length(x2Knorm));

% Display maximum and mean of the epipolar constraints
maxConstraint = max(epipolarConstraints);
meanConstraint = mean(epipolarConstraints);
fprintf('Maximum epipolar constraint: %e\n', maxConstraint);
fprintf('Mean epipolar constraint: %e\n', meanConstraint);

% Comment: The epipolar constraint is approx 0 and sufficiently satisfied 

% Computing the Fundamental Matrix for unnormalized points
F = K' * E * K;



% Compute the epipolar lines for the first set of points
l = F * x1_og;
l = l ./ sqrt(repmat(l(1,:).^2 + l(2,:).^2, [3, 1])); % Normalize lines

% Randomly select 20 indices
randomIndices = randi(length(l), 20, 1);
randPoint = x2_og(:, randomIndices);
randEpi = l(:, randomIndices);


%%% PLOT: Epipolar Lines on Image 2 %%%

% Prep images
pic1 = "kronan1.jpg";
im1 = imread(pic1);
pic2 = "kronan2.jpg";
im2 = imread(pic2);

figure()
imagesc(im2);
hold on
plot(randPoint(1,:), randPoint(2,:), 'r.', 'MarkerSize', 7, 'DisplayName', 'Random Points')
for k = 1:20
    % Plot epipolar lines for the random points
    rital(randEpi(:,k)); % Assuming 'rital' plots the epipolar lines correctly
end
title('Epipolar Lines with Random Points on Image 2');
legend('show');
hold off

%%% Histogram of Total Epipolar Line Distances %%%

figure()
totales = abs(sum(l .* x2_og));
hist(totales, 100);
title('Histogram of Total Epipolar Line Distances');
xlabel('Distance');
ylabel('Frequency');
% Additional line to display the mean of totales
meanValue = mean(totales);
hold on;
plot([meanValue, meanValue], ylim, 'r--', 'LineWidth', 2, 'DisplayName', 'Mean Distance');
legend('show');
hold off

% Normalize the Essential matrix E for further use or saving
E = E./E(3,3)

% Save the Essential matrix to a file
save("Essential_matrix.mat", "E");


save("Essential_matrix.mat","E")
%%% END ESSENTIAL MATRIX %%%

%% 3. The Essential Matrix
%%% Computer Exercise 4. %%%
close all;

% From E6
% rotation ambiguity resolution in Essential Matrix decomposition
W = [0 -1 0;
     1 0 0;
     0 0 1];

% Camera configs
P2 = {};

% Extract the third column of U from SVD(E) as camera center (prev task)
U3 = U(:,3);

% Generate four possible configurations for the second camera matrix P2
P2{1} = [U * W * V', U3];
P2{2} = [U * W * V', -U3];
P2{3} = [U * W' * V', U3];
P2{4} = [U * W' * V', -U3];

% Define the camera projection matrix for the first camera (reference camera)
P1 = [eye(3,3), zeros(3,1)];

% Initialize cell arrays 
xproj = {};
lamdas = {};
X3d = {};

% Initialize an array to keep count of inliers or valid points during triangulation
count = [];

% K-Normalize the image points
x1Knorm = K\x1_og; % Inverse of K times x1_og
x2Knorm = K\x2_og;

% Iterate over each configuration of the second camera matrix
for k = 1:4
    P = P2{k}; % Select the k-th configuration of P2
    
    % Combine the projection matrices of the first and k-th second camera
    Pcombo = [P1; P]; % Combined projection matrix
    
    % Initialize matrices for storing the homogeneous coordinates of 3D points
    X = [];
    
    % Iterating over each pair of corresponding normalized image points
    for i = 1:length(x1Knorm)
        % Construct the matrix for the homogeneous system of equations
        xh = [x1Knorm(:,i), zeros(3,1); zeros(3,1), x2Knorm(:,i)];
        M = [Pcombo, -xh];
        
        % Solve the system using SVD
        [Unew, Snew, Vnew] = svd(M);
        
        % Extract the homogeneous coordinates of the 3D point
        X(:,i) = Vnew(1:4, end);
    end
    
    % Store 3D points
    X3d{k} = X;
    
    % Project the 3D points back 
    xPic = P * X;
    
    % Store the projections for the second camera
    xPics{k} = xPic;
    
    % Determine the number of points with positive depth relative to the second camera
    pos = xPics{k}(3,:) >= 0;
    count = [count; sum(pos)];
end


% Plot the original image points on the second image
figure();
plot(pflat(xPics{1}(1,:)), pflat(xPics{1}(2,:)), 'o'); % Plot original points
hold on;
imagesc(im2); % Display the second image
hold off;
axis equal; % Equal aspect ratio
title('Original Image Points on Second Image'); % Add figure title

% Process and plot 3D points for all camera configurations
X1 = pflat(X3d{1});
X2 = pflat(X3d{2});
X3 = pflat(X3d{3});
X4 = pflat(X3d{4});

% Process and plot 3D points for all camera configurations
figure();
plot3(X1(1,:), X1(2,:), X1(3,:), '.'); % Configuration 1
hold on;
plotcams(P2); % Plot camera positions for all configurations
plot3(X2(1,:), X2(2,:), X2(3,:), 'o'); % Configuration 2
plot3(X3(1,:), X3(2,:), X3(3,:), '+'); % Configuration 3
plot3(X4(1,:), X4(2,:), X4(3,:), '*'); % Configuration 4
hold off;
title('3D Points for All Camera Configurations'); % Add figure title

% Selecting the correct camera configuration after visualization
Cam = {};
Cam{1} = P2{2}; % The second configuration is selected as correct
Cam{2} = P1;

% Plot 3D points from the correct configuration
figure();
plot3(X2(1,:), X2(2,:), X2(3,:), '.');
hold on;
plotcams(Cam); % Plot the selected and reference cameras
hold off;
title('3D Points with Selected Camera Configuration'); % Add figure title

% Project the 3D points using the correct camera matrix
xPics = pflat(K * Cam{1} * X2); % Corrected projection code

% Overlay the projected points on the original second image
figure();
imagesc(im2);
hold on;
plot(xPics(1,:), xPics(2,:), '.'); % Projected points
plot(x2_og(1,:), x2_og(2,:), 'o'); % Original points
axis equal;
hold off;
title('Projected 3D Points on Second Image'); % Add figure title

% Compute and display the discrepancy between original and projected points
dist = sqrt((x2_og(1,:) - xPics(1,:)).^2 + (x2_og(2,:) - xPics(2,:)).^2);
disp(['Max distance: ', num2str(max(dist))]);
disp(['Mean distance: ', num2str(mean(dist))]);

























