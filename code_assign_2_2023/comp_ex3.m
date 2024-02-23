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

for i = 1:2
    sigma(:,i) = std(x{i}(1:2,:),0,2); %Standard deviation of the 1st and 2nd rows ow x{i}
    tilde(:,i) = mean(x{i}(1:2,:),2);  %Computes the mean value of the 1st and 2nd rows ow x{i}
end


N1 = [1/sigma(1,1) 0 -tilde(1,1)/sigma(1,1); 0 1/sigma(2,1) -tilde(2,1)/sigma(2,1); 0 0 1];
N2 = [1/sigma(1,2) 0 -tilde(1,2)/sigma(1,2); 0 1/sigma(2,2) -tilde(2,2)/sigma(2,2); 0 0 1];

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
% above. We see that M1v1 = 0.015061 and M2v2 = 0.012153. This should
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

T3 = table(K1,q1);
disp(T3);

%%% ANSWER
% We know that these are "true", or at least approximatively true, due to
% that the camera matrix has another origin. I have based it on the camera 
% data and made a good approximation of P1 with the DLT method. 
% It is affected of noise though, which was the reason of using SVD earlier
% as well. 




save('for_comp_e5_2.mat', 'P1', 'P2');
