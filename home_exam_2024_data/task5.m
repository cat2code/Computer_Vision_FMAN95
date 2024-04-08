%% Task 5
clearvars;close all;

%%% a)

% Define the camera matrices
P1 = [1 0 0 0; 0 1 0 0; 0 0 1 0];
P2 = [1 0 4 1; 0 1 2 0; 0 0 1 1];
P3 = [-1 1 -2 1; -1 0 -1 0; 0 -2 0 0];

A1 = P1(1:3,1:3);
A2 = P2(1:3,1:3);
A3 = P3(1:3,1:3);

% Compute the epipoles for each camera matrix by solving P * e = 0
C1 = null(P1);
C2 = null(P2);
C3 = null(P3);

C1 = C1./C1(end,end);
C2 = C2./C2(end,end);
C3 = C3./C3(end,end);

e1 = P1*C1;
e2 = P2*C2;
e3 = P3*C3;

% Normalize epipoles
e1 = e1 / e1(end);
%e2 = e2 / e2(end); % Do not normalize
e3 = e3 / e3(end); 

% Convert the epipoles into skew symmetric matrices
e1x = [0 -e1(3) e1(2); e1(3) 0 -e1(1); -e1(2) e1(1) 0];
e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0];
e3x = [0 -e3(3) e3(2); e3(3) 0 -e3(1); -e3(2) e3(1) 0];

F12 = e2x * P2 * pinv(P1);
F13 = e3x * P3 * pinv(P1);
F23 = e3x * P3 * pinv(P2);

disp('Fundamental matrix F12:');
disp(F12);

disp('Fundamental matrix F13:');
disp(F13);

disp('Fundamental matrix F23:');
disp(F23);


%% Task 6: Compute the epipoles and verify equations

% Compute epipoles for each fundamental matrix
% The epipole e_ij is the null space of F_ij (right epipole)
% The epipole e_ji is the null space of F_ij' (left epipole)

e12 = null(F12); % Right 
e21 = null(F12'); % Left 

e13 = null(F13); % Right 
e31 = null(F13'); % Left 

e23 = null(F23); % Right 
e32 = null(F23'); % Left

% Normalize 
%e12 = e12 / e12(end);
%e21 = e21 / e21(end);
e13 = e13 / e13(end);
e31 = e31 / e31(end);
e23 = e23 / e23(end);
e32 = e32 / e32(end);

% Verify the equations
eq1 = e23' * F12 * e13;
eq2 = e32' * F13 * e12;
eq3 = e31' * F23 * e21;


disp('Verification of equations:');
disp(['e23'' * F12 * e13 = ', num2str(eq1)]);
disp(['e32'' * F13 * e12 = ', num2str(eq2)]);
disp(['e31'' * F23 * e21 = ', num2str(eq3)]);

tolerance = 1e-10; 
disp('Equations verification results:');
disp(['eq1 ≈ 0: ', num2str(abs(eq1) < tolerance)]);
disp(['eq2 ≈ 0: ', num2str(abs(eq2) < tolerance)]);
disp(['eq3 ≈ 0: ', num2str(abs(eq3) < tolerance)]);



%% b) 

% Compute epipoles
e12 = null(F12); % Right null space of F12, gives e12
e21 = null(F12'); % Left null space of F12, gives e21

% Compute the other epipoles
e23 = null(F23); 
e32 = null(F23');
e13 = null(F13); 
e31 = null(F13');

% Normalize 
e12 = e12 / e12(end);
e21 = e21 / e21(end);
e23 = e23 / e23(end);
e32 = e32 / e32(end);
e13 = e13 / e13(end);
e31 = e31 / e31(end);

% Verify the equations
verification1 = e23' * F12 * e13;
verification2 = e32' * F13 * e12;
verification3 = e31' * F23 * e21;

% Display the results
fprintf('Verification of e23''*F12*e13 = %f\n', verification1);
fprintf('Verification of e32''*F13*e12 = %f\n', verification2);
fprintf('Verification of e31''*F23*e21 = %f\n', verification3);



%% Alternate solution
%%% a) 
clearvars; close all;

% Define the camera matrices
P1 = [1 0 0 0; 0 1 0 0; 0 0 1 0];
P2 = [1 0 4 1; 0 1 2 0; 0 0 1 1];
P3 = [-1 1 -2 1; -1 0 -1 0; 0 -2 0 0];

% Compute the epipoles for each camera matrix by solving P * e = 0
[~, ~, V1] = svd(P1);
[~, ~, V2] = svd(P2);
[~, ~, V3] = svd(P3);

C1 = V1(:,end) / V1(end,end);
C2 = V2(:,end) / V2(end,end);
C3 = V3(:,end) / V3(end,end);

e1 = P1 * C1;
e2 = P2 * C2;
e3 = P3 * C3;

% Normalize epipoles
e1 = e1 / e1(end);
%e2 = e2 / e2(end); % skip normalizing
e3 = e3 / e3(end); 

% Convert the epipoles into skew symmetric matrices
e1x = [0 -e1(3) e1(2); e1(3) 0 -e1(1); -e1(2) e1(1) 0];
e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0];
e3x = [0 -e3(3) e3(2); e3(3) 0 -e3(1); -e3(2) e3(1) 0];

% Compute the fundamental matrices using corrected formulas
F12 = e2x * P2 * pinv(P1);
F13 = e3x * P3 * pinv(P1);
F23 = e3x * P3 * pinv(P2);

disp('Fundamental matrix F12:');
disp(F12);

disp('Fundamental matrix F13:');
disp(F13);

disp('Fundamental matrix F23:');
disp(F23);

