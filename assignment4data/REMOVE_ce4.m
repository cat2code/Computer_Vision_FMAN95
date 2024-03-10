% Clear the workspace and close all figures to start fresh.
clear
close all

% Add the path to the dataset and load the variables obtained from a previous computer exercise. These variables include camera matrices (P1, P2), 3D points (X), and image points (x).
addpath('../assignment4data')
load('ce2_variables.mat')

% Set a small step size for the update and initialize the lambda parameter for the LM algorithm. Lambda controls the balance between the gradient descent and Gauss-Newton method behaviors.
gammak = 10^-10;
lambda = 0.1;

% Organize the camera matrices into a cell array and assign the 3D and image points to variables for easier manipulation.
P = {P1,P2};
U = X;
u = x;

% Compute the initial reprojection error and residuals for the current solution (cameras and 3D points).
[err,res] = ComputeReprojectionError(P,U,u);

% Initialize variables to hold the updated solution.
Pnew = P;
Unew = U;

% Plot the initial objective value, which is the sum of residuals, indicating the total reprojection error.
figure();
plot(0, sum(res), 'r*')
hold on

% Iterate 10 times to perform updates using the Levenberg-Marquardt method.
for i=1:10
    % Linearize the reprojection error around the current solution to obtain the residuals (r) and the Jacobian (J).
    [r,J] = LinearizeReprojErr(Pnew,Unew,u);
    
    % Calculate the LM update. Unlike in Steepest Descent, here we adjust the update formula to include a damping factor (lambda) which helps in navigating both steep and flat regions of the error surface.
    C = J'*J + lambda*speye(size(J,2));
    c = J'*r;
    deltav = -C\c;
    
    % Apply the computed update to refine the solution (camera matrices and 3D points).
    [Pnew, Unew] = update_solution(deltav,Pnew,Unew);
    
    % Recompute the reprojection error and residuals with the updated solution.
    [err,res] = ComputeReprojectionError(Pnew,Unew,u);
    
    % Plot the updated objective value (sum of residuals) for each iteration to visualize the optimization progress.
    plot(i, sum(res), 'r*')
end

% Calculate the final Root Mean Square (RMS) error using the total error and the number of residuals.
RMS = sqrt(err/size(res,2))
