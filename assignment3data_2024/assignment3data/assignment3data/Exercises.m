% By Eliot

%% 2. 2 The Fundamental Matrix
clearvars;close all;clc;
%%% Exercises 1-4 %%%

t1 = zeros(1,3)';
P1 = [eye(3) t1];
t2 = [2 2 c0]';
A = [1 1 1; 0 2 0; 0 0 1];
P2 = [A t2];

% Calculate camera centers
C1 = null(P1)
C2 = null(P2);
C2(3) = 0;
C2 = C2./C2(4)  % normalizing C2

% Project camera centers onto camera 1 and camera 2
e1 = P1*C2
e2 = P2*C1

% Compute fundamental matrix 
y = [0 -t2(3) t2(2); t2(3) 0 -t2(1); -t2(2) t2(1) 0]
F = y*A

% ...its determinant
determinant_of_F = det(F)

% ...and verify that e2'F = 0 and Fe1 = 0;

e2F = e2'*F
Fe1 = F*e1

% Define F
F = [0, 1, 1; 1, 0, 0; 0, 1, 1];

% F transpose for finding e2
F';

% Find e_{2}
e2 = null(F');
e2 = e2./e2(3)

% Define vector product operator [e_{2}]x
y = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0]

% Calculate A for P2
A = y*F

% Define P2
P2 = [A e2]

% Define 3D scene points
scene1 = [1 2 3 1]';
scene2 = [3 2 1 1]';

% Camera center och P2
CP2 = null(P2);
CP2 = CP2./CP2(2)

% Project scene1
x1a = P1*scene1;
x2a = P2*scene1;

% Project scene2
x1b = P1*scene2;
x2b = P2*scene2;

% epipolar constraint check
scene1_check = x2a'*F*x1a == 0
scene2_check = x2b'*F*x1b == 0

% They fulfill the constraint, as expected.
