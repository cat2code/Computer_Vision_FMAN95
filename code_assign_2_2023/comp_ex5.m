%% 7 Triangulation using DLT ((OPTIONAL))
%%% Computer Exercise 5 %%%
clearvars;close all;clc;
load for_comp_e5.mat
load for_comp_e5_2.mat

x1(3,:) = ones(1,length(x1));
x2(3,:) = ones(1,length(x2));

% Set up the DLT equations for triangulation, and solve the homogeneous 
% least squares system.

% lambda * x = P * X


% GENERATE M1-MATRIX
% Generate left side X^T
Mleft = [P1;P2];


% Generate right side
Mright = [];
xlen = length(x1(1,:));
for i = 1:2
    v = [zeros(3,xlen*(abs(1-i))) -x1 zeros(3,xlen*(abs(2-i)))];
    Mright = [Mright; v];
end

M = [Mleft Mright];

[U,S,V] = svd(M); %Computes the singular value decomposition of M

S_squared = S'*S;

smallest_egeinval_S = S_squared(length(S_squared),length(S_squared));

V_star = V(:,length(V));


% Calculate ||Mv||
Mv = norm(M*V_star);


T1 = table(Mv,smallest_egeinval_S);
disp(T1); 


P = reshape(V_star(1:12),[4 3])'; %Takes the first 12 entries of sol and row-stacks them in a 3x4 matrix
% This does not look right...?


% Check that they are in front of camera
%PX = P*Xmodel;
%negative_elements_PX = PX(3,:) < 0;
%disp(negative_elements_PX); 
%xtilde = P*Xmodel;





% Project the computed points into the two images and compare with the 
% corresponding SIFT-points.







%Useful matlab commands:

%xproj1 = pflat(Pa*X); 
%xproj2 = pflat(Pb*X); %Computes the projections

%good_points = (sqrt(sum((x1-xproj1(1:2,:)).^2)) < 3 & ... 
%    sqrt(sum((x2-xproj2(1:2,:)).^2)) < 3);
%Finds the points with reprojection error less than 3 pixels in both images

%X = X(:,good_points);
%Removes points that are not good enough.

