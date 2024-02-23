function A = pflat(X,makeplot)
%PFLAT Summary of this function goes here
%   Computes and plots the cartesian coordinates of the homogenous
%   coordinates input. Handles 2D and 3D homogenous coordinates. 
%
%   @ A - the cartesian representation. The "flattened" object on homogenous
%   coordinates
%   @ X - input matrix in homogenous coordinates. Allowed to be 2D and 3D. 
%   @ makeplot - Choose if it should plot the output or not, 0 if no plot
%                should be produced, no input assumes making a plot. 

if nargin < 2
    makeplot = 1;
end

% Divide by last homogenous coordinate
Z = X(end,:);                  
A = X./Z;                       

if makeplot == 1
    %figure()
    if length(A(:,1)) == 3                      
        plot(A(1,:),A(2,:),'.')         %Plots a point at (a(1,i),a(2,i)) for each i. 
    elseif length(A(:,1)) == 4  
        plot3(A(1,:),A(2,:),A(3,:),'.') %Same as above but 3D.
    end
    
    axis equal                      %Makes sure that all axes have the same scale.
end

end

