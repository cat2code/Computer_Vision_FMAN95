% Assignment 2 
% Written by Eliot M.P el6183mo-s
% I hereby state that this is my own and original work


%% 2. Calibrated vs. Uncalibrated Reconstruction - part 1
clearvars;close all;clc;

%%% Computer Exercise 1 %%%
load compEx1data.mat

% Description 1
%the 3D points of the reconstruction X, 
%the camera matrices P, 
%the image points x 
%and the filenames imfiles of the images

% Description 2
% X is a 4 × 9471 matrix containing the homogeneous coordinates for all 3D points, 
% x{i} is a 3 × 9471 matrix containing the homogeneous coordinates of the image points seen in image i (NaN means that the point has not been detected in this image). 
% P{i} contains the camera matrix of image i 
% and imfiles{i} contains the name of that image.


%%% 1. Plot the 3D points of the reconstruction. %%%

X_cartesian = pflat(X);
hold on
axis equal
plotcams(P);
title("3D Reconstruction")
legend("3D reconstruction","9 Cameras" )

%%% Answer 1 %%%
% No, the reconstruction looks distorted. The corners of the wall should be
% parallel. 
% ---


%%% 2. Project the 3D points into one of the cameras (only those that have
%%% been detected in that camera) %%%

%%% 3. Plot the image, the projected points, and the image points in the
%%% same figure. %%%

%%% 4. Do the projections appear to be close to the corresponding image points? 
%%% (If not: Did you forget to divide by the third coordinate?) %%%


% For loop to do it with i number of cameras. Here only camera 1 is chosen
for i = 1:1
    figure('Name',"Camera P"+i)
    pic = imfiles{i};
    im1 = imread(pic); %Reads the imagefile with name in imfiles{i}
    imagesc(im1)
    axis equal
    hold on
    
    % Project 3D points onto camera i
    P_actual = cell2mat(P(i));
    xtilde = P_actual*X_cartesian;
    xtilde = pflat(xtilde,0);
    
    visible = isfinite(x{i}(1,:)); % Determines which of the points are visible in image i
    
    % Plot projected points and image points
    plot(x{i}(1,visible), x{i}(2,visible),'*'); %Plots a  *’ at each point coordinate
    plot(xtilde(1,visible), xtilde(2,visible),'ro'); %Plots a red  ot at each visible point in xproj
    
    title("2D projection, camera = "+i, pic)
    legend("Projected points","Image points")
end

%%% Answer 2-4 %%%
% Yes, the image points seem to be corresponding to the projections
% exactly. This is to be expected since the projection and image points are
% suppossed to represent the same 2D image of the same 3D world. 




%% 2. Calibrated vs. Uncalibrated Reconstruction - part 2 
clearvars;close all;clc;
%%% Computer Exercise 1 %%%
load compEx1data.mat

%%% 5. modify the 3D points and cameras with T1 and T2 %%%

% T1 transform - 3D %
figure('Name', "3D Reconstruction T1 transformed")

T1 = [1 0 0 0; 0 4 0 0; 0 0 1 0; 0.1 0.1 0 1];

X_tilde1 = pflat(T1*X);   % calc Xtilde, divide with 4 coord and 3D plot

for i = 1:9
    P{i} = P{i}*inv(T1);    % transforming camera positions
end

hold on
axis equal
plotcams(P);
title("3D Reconstruction T1 transformed")
legend("3D reconstruction","9 Cameras" )


% 2D projections and so on...
for i = 1:1
    figure('Name',"2D projection T1 transformed")
    pic = imfiles{i};
    im1 = imread(pic); %Reads the imagefile with name in imfiles{i}
    imagesc(im1)
    axis equal
    hold on
    
    % Project 3D points onto camera i
    P_actual = cell2mat(P(i));
    x_proj = P_actual*X_tilde1;
    x_proj = pflat(x_proj,0);
    
    visible = isfinite(x{i}(1,:)); % Determines which of the points are visible in image i
    
    % Plot projected points and image points
    plot(x{i}(1,visible), x{i}(2,visible),'*'); %Plots a  *’ at each point coordinate
    plot(x_proj(1,visible), x_proj(2,visible),'ro'); %Plots a red  ot at each visible point in xproj
    
    title(" 2D projection T1 transformed, camera = "+i, pic)
    legend("Projected points","Image points")
end

 
% T2 transform - 3D %
clearvars;
load compEx1data.mat % Reset P

figure('Name', "3D Reconstruction T2 transformed")

T1 = [1 0 0 0; 0 4 0 0; 0 0 1 0; 0.1 0.1 0 1];
T2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 1/16 1/16 0 1];

X_tilde2 = pflat(T2*X);   % calc Xtilde, divide with 4 coord and 3D plot
for i = 1:9
    P{i} = P{i}*inv(T2);    % transforming camera positions
end

hold on
axis equal
plotcams(P);
title("3D Reconstruction T2 transformed")
legend("3D reconstruction","9 Cameras" )


% 2D projections and so on...
for i = 1:1
    figure('Name',"2D projection T2 transformed")
    pic = imfiles{i};
    im1 = imread(pic); %Reads the imagefile with name in imfiles{i}
    imagesc(im1)
    axis equal
    hold on
    
    % Project 3D points onto camera i
    P_actual = cell2mat(P(i));
    x_proj = P_actual*X_tilde2;
    x_proj = pflat(x_proj,0);
    
    visible = isfinite(x{i}(1,:)); % Determines which of the points are visible in image i
    
    % Plot projected points and image points
    plot(x{i}(1,visible), x{i}(2,visible),'*'); %Plots a  *’ at each point coordinate
    plot(x_proj(1,visible), x_proj(2,visible),'ro'); %Plots a red  ot at each visible point in xproj
    
    title(" 2D projection T2 transformed, camera = "+i, pic)
    legend("Projected points","Image points")
end


%%% Answer 5 %%%
% It appears that the T2 transformation creates a 3D representation that is
% more accurate to the real world (no / small distortion and corners are 
% parallel). The T1 transformation does not look good and the camera
% positions does not seem to make sense. 
% However, when it comes to the projections and the 2D image the projection
% points and the image points are on top of each other and the result is
% this preserved compared to the last one.

save('for_comp_e2.mat', 'T1', 'T2');


%% 4 RQ Factorization and Computation of K

%%% Computer Exercise 2 %%%
clc;clearvars;close all;

load for_comp_e2.mat
load compEx1data.mat

% Task:
% Compute K for one camera in each of the solutions you obtained in computer exercise 1.
% Choosing camera P(1)
% Submit the K matrices. (Make sure that element K(3,3) = 1 by division; K = K./K(3, 3).)


% Re-generate camera matrices from Ex1
P1 = P{2}; % The assignment in 2024 wants you to choose P{2}
PT1 = P{2}*inv(T1);
PT2 = P{2}*inv(T2);   

%P1 = P{1};
%PT1 = P{1}*inv(T1);
%PT2 = P{1}*inv(T2);   

% Calculate K
% RQ factorization
% r upper tri, q unit matrix
[K1,q1]     = rq(P1);
[K1T1,qT1]   = rq(PT1);
[K1T2,qT2]   = rq(PT2);

% Make sure K(3,3) = 1
K1  = K1./K1(3,3);
K1T1 = K1T1./K1T1(3,3);
K1T2 = K1T2./K1T2(3,3);


% Do they represent the same transformation?
test1 = K1 == K1T1  % not the same
test2 = K1T1 == K1T2 % not the same
test3 = K1 == K1T2  % the same!


save('K_matrices.mat', 'K1', 'K1T1', 'K1T2');


%%% Answer %%%
% In the case of camera matrix P(1) and P(1) with T2 transformation, yes
% they do represent the same transformation. This is not true for all
% combinations of P(1) with T1 transformation. This is not due to a scale
% factor not being correct. 





%% 5. Direct Linear Transformation DLT
clc;clearvars;close all;

%%% Computer Exercise 3 %%%
im1 = imread('cube1.JPG');
figure(1)
imagesc(im1)
title("Picture: cube1.JPG")
im2 = imread('cube2.JPG');
figure(2)
imagesc(im2)
title("Pciture: cube2.JPG'")
load compEx3data.mat;



% Plot 3D model (Xmodel)
figure()
%Plots the lines of the cube model (works only if all points are included)
plot3([Xmodel(1,startind); Xmodel(1,endind)],... 
    [Xmodel(2,startind); Xmodel(2,endind)],... 
    [Xmodel(3,startind); Xmodel(3,endind)],'b-'); 
title('3D model plot')
legend('Xmodel')


x1 = x{1};
x2 = x{2};

%%% Optional exercise

% Define indices of points to keep
points_to_keep = [1, 4, 13, 16, 25, 28, 31];

% Filter points
%Xmodel = Xmodel(:, points_to_keep);
%x1 = x1(:, points_to_keep);
%x2 = x2(:, points_to_keep);

%%% %%% %%% 

for i = 1:2
    sigma(:,i) = std(x{i}(1:2,:),0,2); %Standard deviation of the 1st and 2nd rows ow x{i}
    tilde(:,i) = mean(x{i}(1:2,:),2);  %Computes the mean value of the 1st and 2nd rows ow x{i}
end


N1 = [1/sigma(1,1) 0 -tilde(1,1)/sigma(1,1); 0 1/sigma(2,1) -tilde(2,1)/sigma(2,1); 0 0 1];
N2 = [1/sigma(1,2) 0 -tilde(1,2)/sigma(1,2); 0 1/sigma(2,2) -tilde(2,2)/sigma(2,2); 0 0 1];

%%% Optional exercise 

% do not normalize:
%N1 = eye(3);
%N2 = eye(3);

%%% %%% %%% %%%


% Normalization to original points
x1tilde = N1*x1;
x2tilde = N2*x2;


figure()
pflat(x1tilde);
hold on
pflat(x2tilde);
title("Projection of the normalized points ")
legend('Normalized x1','Normalized x1')

%%% Answer %%%
% Plot the normalized points in a new figure. Does it look like the points 
% are centered around (0, 0) with mean distance 1 to (0, 0)?
% 
% -> Yes when just observing the points seem to be centered 
%       around origo (0,0), though I have not checked exactly. 
% 


% normalize original image points and save the orginal in another
x1_og = x1;
x2_og = x2;
x1 = N1*x1;
x2 = N2*x2;


%%% SET UP DLT USING SVD %%%
Xmodel(4,:) = ones(1,length(Xmodel));


% General method:
%• Set up the linear homogeneous system
%• Compute the singular value decomposition.
%• Extract the solution v∗ from the last column of V .


% GENERATE M1-MATRIX
% Generate left side X^T
Mleft = [];
for i = 1:length(Xmodel)
    number_of_blocks = 3;
    left = zeros(3,4*number_of_blocks);
    for j = 0:2
        left(j+1,j*4+1:j*4+4) = Xmodel(:,i)';
    end
    Mleft = [Mleft;left];
end

% Generate right side
Mright = [];
xlen = length(x1(1,:));
for i = 1:xlen
    v = [zeros(1,3*(i-1)) -x1(:,i)' zeros(1,(xlen-i)*3)];
    Mright = [Mright v'];
end


M1 = [Mleft Mright];


% GENERATE M2-MATRIX
% (left side is the same)

% Generate right side
Mright = [];
xlen = length(x2(1,:));
for i = 1:xlen
    v = [zeros(1,3*(i-1)) -x2(:,i)' zeros(1,(xlen-i)*3)];
    Mright = [Mright v'];
end

M2 = [Mleft Mright];


[U1,S1,V1] = svd(M1); %Computes the singular value decomposition of M
[U2,S2,V2] = svd(M2); 

S1_squared = S1'*S1;
S2_squared = S2'*S2;

smallest_egeinval_S1 = S1_squared(length(S1_squared),length(S1_squared));
smallest_egeinval_S2 = S2_squared(length(S2_squared),length(S2_squared));

V1_star = V1(:,length(V1));
V2_star = V2(:,length(V2));


% Calculate ||Mv||
M1v1 = norm(M1*V1_star);
M2v2 = norm(M2*V2_star);


% Create a table
T2 = table(M1v1,M2v2);
T1 = table(smallest_egeinval_S1,smallest_egeinval_S2);


% Print the table
disp(T1); 
disp(T2); 


%%% Answer %%% 
% * The smallest egeinvalue is found from the last element in the diagonal of 
% S'*S according to theory, yes they are very small for both S1 and S2.
% They are though not zero due to unavoidable noise from the picture. 
%
% * ||Mv|| =/= 0. This is also expected due to the same reason as mentioned
% above. We see that M1v1 = 0.0092149 and M2v2 = 0.0085228. This should
% though according to theory be the minimal solution to the problem of
% finding a solution to the DLT problem "Mv = 0".



%%% SET UP THE CAMERA MATRICES %%%

% THE CHECK LATER GIVES THAT THE V2 SOLUTION HAS TO BE TURNED AROUND
V1_star = -V1_star;
V2_star = -V2_star;
%%%%%


P1 = reshape(V1_star(1:12),[4 3])'; %Takes the first 12 entries of sol and row-stacks them in a 3x4 matrix
P2 = reshape(V2_star(1:12),[4 3])';



% Check that they are in front of camera
P1X = P1*Xmodel;
P2X = P2*Xmodel;

negative_elements_P1X = P1X(3,:) < 0;
negative_elements_P2X = P2X(3,:) < 0;

disp(negative_elements_P1X); 
disp(negative_elements_P2X); 
% -> V2 has to be turned around, all elements are negative in the third row
% when it is not turned around by negative multiplication.




% Un-normalize camera matrices 
P1 = inv(N1)*P1;
P2 = inv(N2)*P2;



x1tilde = P1*Xmodel;
x2tilde = P2*Xmodel;

%%% 2D projections start %%%
i = 1;
im = im1;
P = P1;
X = Xmodel;
x = x1_og;

figure('Name',"Camera P"+i)
imagesc(im)
axis equal
hold on

% Project 3D points onto camera calculation
xtilde = P*X;
xtilde = pflat(xtilde,0);

x = pflat(x,0);

% Plot projected points and image points 
plot(xtilde(1,:), xtilde(2,:),'*');  % Plots a  *
plot(x(1,:), x(2,:),'ro');           % Plots a red dot

title("2D projection, camera = "+i)
legend("Projected points","Image points")


i = 2;
im = im2;
P = P2;
X = Xmodel;
x = x2_og;

figure('Name',"Camera P"+i)
imagesc(im)
axis equal
hold on

% Project 3D points onto camera calculation
xtilde = P*X;
xtilde = pflat(xtilde,0);

x = pflat(x,0);

% Plot projected points and image points 
plot(xtilde(1,:), xtilde(2,:),'*');  % Plots a  *
plot(x(1,:), x(2,:),'ro');           % Plots a red dot

title("2D projection, camera = "+i)
legend("Projected points","Image points")


%%% 2D projections end %%%




% Make third row P3 larger (!!!)
P1 = P1./P1(3,3);
P2 = P2./P2(3,3);
% This fixes the problem of the principal axis being to small to plot the
% viewing angle...



%%% Plot camera centers in 3D %%%

center1 = null(P1);
center1 = pflat(center1,0);     % Center in (0,0,0)
center2 = null(P2);
center2 = pflat(center2,0);     % Center in ≈(6,6;14,8;-15,0)

principal1 = P1(3,1:3);
principal2 = P2(3,1:3);

figure()
pflat(Xmodel);
axis equal;
hold on;
title("3D plot") 

% Plot camera centers
plot3(center1(1),center1(2),center1(3),'*','Color','r');   %Same as plot but with smaller points
plot3(center2(1),center2(2),center2(3),'*','Color','r');   %Same as plot but with smaller points

% Plot viewing angles
s = 10;
quiver3(center1(1),center1(2),center1(3),principal1(1),principal1(2),principal1(3),s) 
quiver3(center2(1),center2(2),center2(3),principal2(1),principal2(2),principal2(3),s) 

legend('3d model','P1 camera center','P2 camera center','P1 viewing angle','P2 viewing angle')

%%% End of plot camera centers in 3D %%%



%%% ANSWERS AND DISCUSSION %%%
% * As for the 2D projection. Yes, the results are satisfactory I
% believe? The projection in its form seems to be correct and they
% correspond exactly to the image points. The normalization seems to be
% working as intended as well.
%
% * The 3D projections. Yes, aside from the rather wierd rotation the form
% seems to be correct, nice! There was also another fix I implemented that 
% I am unsure of why it works. The third row of P, P3, was very small so
% that the principal axis could not be interpreted by Matlab. To fix this I
% normalized the size of the elements by calculating P = P./P(3,3). This
% fixed the plot for the viewing angle even though I am unsure what is in
% detail going on here. Reasonable results.


%%% CALCULATE INNER PARAMETERS %%%

% Make P smaller again to original. Unnecessary?
%P1 = P1./P1(3,3);
%P2 = P2./P2(3,3);


[K1, q1] = rq(P1);
K1 = K1./K1(3,3);

[K2, q2] = rq(P2);
K2 = K2./K2(3,3);

T3 = table(K1,q1);
disp(T3);

T4 = table(K2,q2);
disp(T4);

%%% ANSWER
% We know that these are "true", or at least approximatively true, due to
% that the camera matrix has another origin. I have based it on the camera 
% data and made a good approximation of P1 with the DLT method. 
% It is affected of noise though, which was the reason of using SVD earlier
% as well. 


save('for_comp_e5_2.mat', 'P1', 'P2', 'K1', 'K2');


%%% OPTIONAL PART OF EXERCISE 3 %%%

% Given that Xmodel is already [4 N], no need to add a row of ones
Xmodel_homogeneous = Xmodel; % This is just to clarify; you can use Xmodel directly

% Now call calculateRMS with Xmodel without adding an extra row
RMS_error_1 = calculateRMS(Xmodel, x1_og, P1); % Use Xmodel directly if it's already [4 N]
RMS_error_2 = calculateRMS(Xmodel, x2_og, P2);

disp(['RMS error for Camera 1: ', num2str(RMS_error_1)]);
disp(['RMS error for Camera 2: ', num2str(RMS_error_2)]);

save('RMS_errors.mat', 'RMS_error_1', 'RMS_error_2');
%save('RMS_errors_subset.mat', 'RMS_error_1', 'RMS_error_2');
%save('RMS_errors_not_normalized.mat', 'RMS_error_1', 'RMS_error_2');
%save('RMS_errors_subset_not_normalized.mat', 'RMS_error_1', 'RMS_error_2');

%% Analyze difference normalization and not
clearvars;

% Load RMS errors for non-normalized points
load('RMS_errors_not_normalized.mat', 'RMS_error_1', 'RMS_error_2');
RMS_error_1_nonorm = RMS_error_1;
RMS_error_2_nonorm = RMS_error_2;

% Load RMS errors for normalized points
load('RMS_errors.mat', 'RMS_error_1', 'RMS_error_2');


% Display RMS errors for comparison
disp(['RMS error for Camera 1 with normalization: ', num2str(RMS_error_1)]);
disp(['RMS error for Camera 2 with normalization: ', num2str(RMS_error_2)]);
disp('---'); % Separator
disp(['RMS error for Camera 1 without normalization: ', num2str(RMS_error_1_nonorm)]);
disp(['RMS error for Camera 2 without normalization: ', num2str(RMS_error_2_nonorm)]);

% Analyze the difference
difference_1 = abs(RMS_error_1 - RMS_error_1_nonorm);
difference_2 = abs(RMS_error_2 - RMS_error_2_nonorm);

disp(['Difference in RMS error for Camera 1: ', num2str(difference_1)]);
disp(['Difference in RMS error for Camera 2: ', num2str(difference_2)]);


%% 6. Feature Extraction and Matching using SIFT
clear;clc;close all;

%%% Computer Exercise 4 %%%

im1 = imread('cube1.JPG');
figure(1)
imagesc(im1)
im2 = imread('cube2.JPG');
figure(2)
imagesc(im2)


[f1 d1] = vl_sift( single(rgb2gray(im1)), 'PeakThresh', 1);
[f2 d2] = vl_sift( single(rgb2gray(im2)), 'PeakThresh', 1);


figure(1)
vl_plotframe(f1);

figure(2)
vl_plotframe(f2);


[matches ,scores] = vl_ubcmatch(d1,d2);


x1 = [f1(1,matches(1,:));f1(2,matches(1,:))]; 
x2 = [f2(1,matches(2,:));f2(2,matches(2,:))];


perm = randperm(size(matches ,2)); figure;
imagesc([im1 im2]);
hold on;
plot([x1(1,perm(1:10)); x2(1,perm(1:10))+size(im1,2)], ... 
    [x1(2,perm(1:10)); x2(2,perm(1:10))],'-');
hold off;

save('for_comp_e5.mat', 'x1', 'x2');

%%% Answer
% How many of the matches appear to be correct?
% -> All of them actually.


%% 7 Triangulation using%% 7 Triangulation using DLT ((OPTIONAL))
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

