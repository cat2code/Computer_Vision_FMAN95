% Assignment 1 
% Written by Eliot Montesino Petrén el6183mo-s
% I hereby state that this is my own and original work
close all;clearvars;clc;

%% Points in Homogeneous Coordinates.

load compEx1.mat
A2D = pflat(x2D);
title("2D cartesian coordinates transformation")
A3D = pflat(x3D);
title("3D cartesian coordinates transformation")

%% Lines 
close all;clearvars;clc;

figure(1)
set(1,'DefaultLineLineWidth',1)
im = imread('compEx2.JPG');
imagesc(im)
colormap gray
hold on

load compEx2.mat

%Plotting lines traditionally to verify
plot(p1(1,1:2),p1(2,1:2), 'LineWidth', 3);
plot(p2(1,1:2),p2(2,1:2), 'LineWidth', 3);
plot(p3(1,1:2),p3(2,1:2), 'LineWidth', 3);

% Solving nullspace of p1' finds parameters for the line going through the
% pair of points in the plane in p1, p2 and p3.
% The null space of two homogenous coordinates finds the lines that intersects
% the corresponding ones in cartesian coordinates.
all_lines=[null(p1') null(p2') null(p3')];
rital(all_lines)
legend('p3','p2','p1','p3','p2','p1', 'intersection')



% Find intersection with nullspace again
intersection = pflat(null([all_lines(:,2) all_lines(:,3)]'),0);
plot(intersection(1,1),intersection(2,1),'y*')

%Answer: 
% The lines are not parallell, they all intersect. They intersect close to
% each other but still distinctly not all intersecting in the same point.
% This is interesting as the lines should be almost parallell to each other
% in the 3D world that they were taken photo of. There are many instances
% that introduce room for error that the lines are not exactly parallell 
% so that they do not intersect in the same point in our plot.
% What is interesting is that they are still enough parallel-alike to 
% almost intersect in the same point.


% Calculate distance
distance = abs(all_lines(:,1)'*intersection)/hypot(all_lines(1,1),all_lines(2,1))
%distance = abs(a*x1 + b*x2 + c)/sqrt(a^2 + b^2); % d = 8.1950
save('distance','distance')

% Answer: Distinctly not zero. The distance is relatively small though.
% This ties into the discussion above. 


%% Projective Transformations
close all;clearvars;clc;

H = [1 1 0; 0 1 0; -1 0 1];
x1 = [1 0 1]';
x2 = [0 1 1]';

% Compute transformations
y1 = H*x1;
y2 = H*x2;

% Compute line intersections
X = [x1' ; x2'];
l1 = null(X);

Y = [y1' ; y2'];
l2 = null(Y);

% Compute inverse and compare
compare1 = inv(H)'*l1;


%%%% Comp. Exercise 3 %%%%
%close all;clearvars;clc;

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

H3 = [1 1 0; 0 2 0; 0 0 1];
H4 = [sqrt(3) -1 1; 1 sqrt(3) 1; 1/4 1/2 2];


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


%%% H3 %%%

% Let H3 transform points. Keeps it a 3x42 matrix
startpoints_hat = H3*startpoints_hat;
endpoints_hat = H3*endpoints_hat;
% Calculate cartesian coordinates
startpoints_hat = pflat(startpoints_hat,0);
endpoints_hat = pflat(endpoints_hat,0);

figure()
plot([startpoints(1,:); endpoints(1,:)], ... 
    [startpoints(2,:); endpoints(2,:)],'b-');
hold on;
plot([startpoints_hat(1,:); endpoints_hat(1,:)], ... 
    [startpoints_hat(2,:); endpoints_hat(2,:)],'b-');
title("H3 figure")
axis equal

% Reset starting and ending points
startpoints_hat = [startpoints; ones(1,length(startpoints))];
endpoints_hat = [endpoints; ones(1,length(endpoints))];


%%% H4 %%%

% Let H4 transform points. Keeps it a 3x42 matrix
startpoints_hat = H4*startpoints_hat;
endpoints_hat = H4*endpoints_hat;
% Calculate cartesian coordinates
startpoints_hat = pflat(startpoints_hat,0);
endpoints_hat = pflat(endpoints_hat,0);

figure()
plot([startpoints(1,:); endpoints(1,:)], ... 
    [startpoints(2,:); endpoints(2,:)],'b-');
hold on;
plot([startpoints_hat(1,:); endpoints_hat(1,:)], ... 
    [startpoints_hat(2,:); endpoints_hat(2,:)],'b-');
title("H4 figure")
axis equal



%%% Answer %%%
% Classification of transformations:
% H1 - Euclidian transformation
%       Preserves the size and shape of the original grid. Only rotattion.
% H2 - Similarity Transform
%       Rotation but does not preserve size of the shape. 
% H3 - Affine transformation
%       Typical shape for affine transformation and preserving paralel
%       lines. It also has an appropiate transformation matrix with [A t;0 1]
% H4 - Projective Transformation
%       All other were chosen so this is the only one left. We can see 
%       though that all straight lines are preserved as it should but the 
%       angles between them are not preserved which is not required. 
%       The transformation is otherwise arbitrary.



%% The Pinhole Camera
close all;clearvars;clc;

%%% Exercise 5 - computing camera projections %%%
P = [eye(3) [0;0;1]];

x1 = [1 2 3 1]';
x2 = [1 1 1 1]';
x3 = [1 1 -1 1]';

proj1 = P*x1;
proj2 = P*x2;
proj3 = P*x3;

center = null(P);


%%% Comp exercise 4 - 3D plotting %%% 
load compEx4.mat

im1 = imread('compEx4im1.JPG');
figure()
imagesc(im1)
colormap gray
title("im1")

im2 = imread('compEx4im2.JPG');
figure()
imagesc(im2)
colormap gray
title("im2")


% Define camera matrices 
P1 = K*[R1 t1];
P2 = K*[R2 t2];

% Camera centers and principal axis
center1 = null(P1);
center1 = pflat(center1,0);     % Center in (0,0,0)
center2 = null(P2);
center2 = pflat(center2,0);     % Center in ≈(6,6;14,8;-15,0)

principal1 = P1(3,1:3);
principal2 = P2(3,1:3);

U_cartesian = pflat(U);
axis equal;
hold on;
title("Original point model") 

% Plot camera centers
plot3(center1(1),center1(2),center1(3),'*','Color','r');   %Same as plot but with smaller points
plot3(center2(1),center2(2),center2(3),'*','Color','r');   %Same as plot but with smaller points


% Plot viewing angles
s = 10;
quiver3(center1(1),center1(2),center1(3),principal1(1),principal1(2),principal1(3),s) 
quiver3(center2(1),center2(2),center2(3),principal2(1),principal2(2),principal2(3),s) 

% ANSWER:
% The results are not entirely reasonable, yet! The point figure looks good
% and the center points for the cameras look very reasonable comparing with
% the images. The same goes for the angle. Very cool!
% The plot could be improved so that the statue is vertical, that would
% have looked more nice...


% Correction: Adding the 3D points projected in camera 1 and camera 2
U1 = pflat(P1*U, 0);
U2 = pflat(P2*U, 0);

%plot of the transformations together with the images
figure(1)
imagesc(im1)
colormap gray
hold on
plot(U1(1,:),U1(2,:),'.','MarkerSize',2)


figure(2)
imagesc(im2)
colormap gray
hold on
plot(U2(1,:),U2(2,:),'.','MarkerSize',2)



%% Optional stuff
close all;clearvars;clc;

im = imread('compEx5.JPG');
load CompEx5.mat


P1 = [eye(3) zeros(3,1)];
%P2 = [R t]

%% Plottinggg
figure(1)
imagesc(im)
colormap gray
hold on

%Plots the cornerpoints and connects them with lines.
plot(corners(1,[1:end 1]), corners(2,[1:end 1]),'*-'); 

axis ij %Makes the y-axis point down (as in an image)



%%

tform = maketform('projective',Htot );
%Creates a projective transformation that can be used in imtransform 
%NOTE: Matlab uses the transposed version of the homografi.

[new_im,xdata,ydata] = imtransform(im,tform,'size',size(im)); 
%Creastes a transformed image (using tform)
%of the same size as the original one.

imagesc(xdata,ydata,new_im);
%plots the new image with xdata and ydata on the axes



















