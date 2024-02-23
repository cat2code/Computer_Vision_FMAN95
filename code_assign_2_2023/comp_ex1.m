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


