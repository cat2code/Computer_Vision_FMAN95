% Clear the workspace and close all figures to start with a clean state.
clear
close all

% Add the path to the dataset and load variables from a previous exercise, including camera matrices (P1, P2), 3D points (X), and image points (x).
addpath('../assignment4data')
load('ce2_variables.mat')

% Set a very small step size for the update (gamma_k), crucial for convergence in Steepest Descent methods.
gammak = 10^-10;

% Organize the camera matrices into a cell array for easier manipulation.
P = {P1,P2};

% Assign the 3D points and image points to more conveniently named variables.
U = X;
u = x;

% Compute the initial reprojection error and residuals for the current solution.
[err,res] = ComputeReprojectionError(P,U,u);

% Initialize new variables for the updated camera matrices and 3D points.
Pnew = P;
Unew = U;

% Plot the initial objective value (sum of residuals).
figure();
plot(0, sum(res), 'r*')
hold on

% Iterate 10 times to perform Steepest Descent updates.
for i=1:10
    % Linearize the reprojection error at the current solution to obtain the residuals (r) and Jacobian (J).
    [r,J] = LinearizeReprojErr(Pnew,Unew,u);
    
    % Calculate the update delta_v using the Steepest Descent formula: delta_v = -gamma_k * J' * r.
    deltav = -gammak*J'*r;
    
    % Update the solution (camera matrices and 3D points) using the calculated delta_v.
    [Pnew, Unew] = update_solution(deltav,Pnew,Unew);
    
    % Recompute the reprojection error and residuals with the updated solution.
    [err,res] = ComputeReprojectionError(Pnew,Unew,u);
    
    % Plot the updated objective value (sum of residuals) for each iteration.
    plot(i, sum(res), 'r*')
end 

% Compute the final RMS error using the total error and number of residuals.
RMS = sqrt(err/size(res,2))

% The script concludes by plotting the objective value vs. iteration number to visualize the optimization process. The use of Steepest Descent method aims to minimize the sum of residuals, which represents the reprojection error, by iteratively updating the solution in the direction that decreases this error the most at each step. The choice of γ_k is critical for the success of the convergence, and in this script, a very small value is used to ensure a careful approach towards the minimum.
