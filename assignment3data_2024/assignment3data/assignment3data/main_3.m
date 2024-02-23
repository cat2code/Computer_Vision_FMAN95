% Assignment 3
% Written by Eliot M.P el6183mo-s
% I hereby state that this is my own and original work

%% 2. 2 The Fundamental Matrix
clearvars;close all;clc;

%%% Computer Exercise 1 %%%
load compEx1data.mat

%%% plot images %%% 
pic = "kronan1.jpg";
figure('Name',pic)
im1 = imread(pic); %Reads the imagefile with name in imfiles{i}
imagesc(im1)

pic = "kronan2.jpg";
figure('Name',pic)
im2 = imread(pic); %Reads the imagefile with name in imfiles{i}
imagesc(im2)
%%% --- %%%

% Declare the matching points
x1 = x{1};
x2 = x{2};

% Original x1 and x2
x1_og = x1;     
x2_og = x2;

%%% Normalization %%%
for i = 1:2
    sigma(:,i) = std(x{i}(1:2,:),0,2); %Standard deviation of the 1st and 2nd rows ow x{i}
    tilde(:,i) = mean(x{i}(1:2,:),2);  %Computes the mean value of the 1st and 2nd rows ow x{i}
end

% Declare N
N1 = [1/sigma(1,1) 0 -tilde(1,1)/sigma(1,1); 0 1/sigma(2,1) -tilde(2,1)/sigma(2,1); 0 0 1];
N2 = [1/sigma(1,2) 0 -tilde(1,2)/sigma(1,2); 0 1/sigma(2,2) -tilde(2,2)/sigma(2,2); 0 0 1];

% Normalize x1 and x2
x1 = N1*x1;
x2 = N2*x2;

%%% Form M-matrix %%%
for i = 1:3
    xx = x2(:,i)*x1(:,i)';      %Computes a 3x3 matrix containing all multiplications 
                                %of coordinates from x1n(:,i) and x2n(:,i).
    M(i,:) = xx(:)';            %Reshapes the matrix above and adds to the M matrix
end


%%% Solve the homogeneous least squares system using SVD %%%

[U,S,V] = svd(M); %Computes the singular value decomposition of M

S_squared = S'*S;

smallest_egeinval_S = S_squared(length(S_squared),length(S_squared));

V_star = V(:,length(V));    % Solution to SVD problem

% Calculate ||Mv||
Mv = norm(M*V_star);

% Create a table
T2 = table(Mv);
T1 = table(smallest_egeinval_S);

% Print the table
disp(T1); 
disp(T2); 


%%% Answer
% The min S is indeed small and so is ||Mv|| as well. This indicates that
% the method is working as intended and I can proceed.


%%% Form the normalized F-matrix %%%

%Forms an F-matrix from the solution v of the least squares problem
Fn = reshape(V_star,[3 3]);  

% Check that epipolar contraint is fullfilled
disp("Is det(F) = 0?")
check = det(Fn) == 0;
if check
    disp("true")
else
    disp("det(F) = " + det(Fn))
end

% Check that det(F) = 0

%Computes and plots all the epipolar constraints (should be roughly 0)
figure()
plot(diag(x2'*Fn*x1))


%%% Answer
% det(F) is zero and the epipolar contraints are roughly filled since they
% are all about 0.


%%% Transform back to un-normalized F %%%
F = N2'*Fn*N1;

%Computes and plots all the epipolar constraints (should be roughly 0)
l = F*x{1}; %Computes the epipolar lines
l = l./sqrt(repmat(l(1,:).^2 + l(2,:).^2,[3 1])); 

%Makes sure that the line has a unit normal
%(makes the distance formula easier) 


hist(abs(sum(l.*x{2})),100);
%Computes all the the distances between the
%and there corresponding lines, and plots in a histogram


%%% Answer

F = F./F(3,3);
save('F_matrix.mat', 'F');
%save('mean_epipolar_dist_F_non_normalized.mat','mean_epipolar_dist_F')
%save('for_other_ex.mat','N1','N2',"x1_og","x2_og","x1tilde","x2tilde")

%% No normalization attempt

N1 = eye(3);
N2 = N1;


%save('mean_epipolar_dist_F_normalized.mat','mean_epipolar_dist_Fn')


%% 2. 2 The Fundamental Matrix
clearvars;close all;clc;
%%% Computer Exercise 2 %%%


load F_matrix.mat
load for_other_ex.mat














%% 3 The Essential Matrix
clearvars;close all;clc;

%%% Computer Exercise 3 %%%
load compEx3data.mat;
load for_other_ex.mat

% K-Normalize the image points
x1 = inv(K)*x1tilde;
x2 = inv(K)*x2tilde;


% Set up the matrix M in the eight point algorithm, and solve the homogeneous least squares system
% using SVD.

M = 0;


E = 0;
%E = E./E(3,3);
save("Essential_matrix.mat","E")


%%
%Useful matlab commands:

[U,S,V] = svd(Eapprox); 
if det(U*V')>0
    E = U*diag([1 1 0])*V'; 
else
    V = -V;
    E = U*diag([1 1 0])*V'; 
end

% Creates a valid essential matrix from an approximate solution.
% Note: Computing svd on E may still give U and V that does not fulfill % det(U*V') = 1 since the svd is not unique.
% So don¨t recompute the svd after this step.¨


%% 3. The Essential Matrix
%%% For E6 %%%
clearvars;close all;clc;


U = [1/sqrt(2), -1/sqrt(2), 0; 1/sqrt(2), 1/sqrt(2), 0; 0, 0, 1];
V = [1 0 0; 0 0 -1; 0 1 0];

check = det(U*V');
table(check)        % It is 1!


% Compute the essential matrix and verify that x1 = (0, 0) (in camera 1) and x2 = (1, 1) (in camera 2)
% is a plausible correspondence.





%% 3. The Essential Matrix
%%% Computer Exercise 4. %%%
clearvars;close all;clc;








































