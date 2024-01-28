% Assignment 1 
% Written by Eliot M.P el6183mo-s
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

im = imread('compEx2.JPG');

figure(1)
imagesc(im)
colormap gray
hold on

load compEx2.mat
p1_no_Z = p1(1:2,:);
p2_no_Z = p2(1:2,:);
p3_no_Z = p3(1:2,:);

%Plotting lines traditionally to verify
plot(p1_no_Z(1,:),p1_no_Z(2,:),'-')
plot(p2_no_Z(1,:),p2_no_Z(2,:),'-')
plot(p3_no_Z(1,:),p3_no_Z(2,:),'-')

legend('p3','p2','p1')

% Solving nullspace of p1' finds parameters for the line going through the
% pair of points in the plane in p1, p2 and p3.
p1_intersect_line = null(p1');
p2_intersect_line = null(p2');
p3_intersect_line = null(p3');

p1_intersect_line(2) = -p1_intersect_line(2);
p2_intersect_line(2) = -p2_intersect_line(2);
p3_intersect_line(2) = -p3_intersect_line(2);


% Store in one matrix
P = [p1_intersect_line p2_intersect_line p3_intersect_line];


% PLOT WITH RITAL
% Rital was unclear how to use I found :( had to add "hold on" in the 
% function's code
rital(P)        % They all intersect
hold on;


%Answer: 
% The lines are not parallell, they all intersect. They intersect close to
% each other but still distinctly not all intersecting in the same point.
% This is interesting as the lines should be almost parallell to each other
% in the 3D world that they were taken photo of. There are many instances
% that introduce room for error that the lines are not exactly parallell 
% so that they do not intersect in the same point in our plot.
% What is interesting is that they are still enough parallel-alike to 
% almost intersect in the same point.

% Calculate intersection point bewteen p2 and p3 lines
p2_p3_point = null([p2_intersect_line p3_intersect_line]');
p2_p3_point = p2_p3_point./p2_p3_point(3);

p2_p3_point_c = p2_p3_point(1:2);        % Intersection point in cartesian

x1 = p2_p3_point_c(1);
x2 = p2_p3_point_c(2);

plot(x1,x2,'*')
legend('p1','p2','p3','intersection p2 p3');

% Calculate distance between intersection point and p1 line
a = p1_intersect_line(1);
b = p1_intersect_line(2);
c = p1_intersect_line(3);

d = abs(a*x1 + b*x2 + c)/sqrt(a^2 + b^2); % d = 8.1950
save('distance','d')

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

H1 = [sqrt(3) -1 1;1 sqrt(3) 1; 0 0 2];
H2 = [1 -1 1 ;1 1 0; 0 0 1];
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





















