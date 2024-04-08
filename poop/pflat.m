function A = pflat(X, makeplot)
%PFLAT Projects homogenous coordinates to cartesian coordinates.
%
% This function computes and optionally plots the cartesian coordinates 
% of the input homogenous coordinates. It handles both 2D and 3D 
% homogenous coordinates.
%
% Inputs:
%   X        - A matrix representing points in homogenous coordinates. 
%              The function supports both 2D and 3D homogenous coordinates.
%   makeplot - A binary flag to determine if the function should plot the 
%              output. 1 (or true) to plot, 0 (or false) to not plot. 
%              If not specified, the function defaults to plotting.
%
% Output:
%   A - The cartesian representation of the input points. This is the 
%       "flattened" version of the input homogenous coordinates.
%
% Example:
%   A = pflat(X); % Defaults to plotting
%   A = pflat(X, 0); % No plot will be produced

% Check if 'makeplot' argument is given, default to 1 (true) if not
if nargin < 2
    makeplot = 0;
end

% Extract the scale factor (last row of X)
Z = X(end, :);  

% Element-wise division to convert to cartesian coordinates
A = X ./ Z;

% Check if plotting is requested
if makeplot == 1
    disp('new fig')
    % Create a new figure
    figure();
    
    % Determine whether to plot in 2D or 3D based on the dimensions of A
    if size(A, 1) == 3
        % 2D plot for 3 rows in A (2D homogenous coordinates)
        plot(A(1, :), A(2, :), '.'); % Plots a point for each column of A
    elseif size(A, 1) == 4
        % 3D plot for 4 rows in A (3D homogenous coordinates)
        plot3(A(1, :), A(2, :), A(3, :), '.'); % Plots a 3D point for each column of A
    end
    
    % Ensure equal scaling of axes
    axis equal;
end

end
