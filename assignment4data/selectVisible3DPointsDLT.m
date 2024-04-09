function [a,b] = selectVisible3DPointsDLT(P, data, x1_n, x2_n)
    % automated comments with GPT, hope that is OK :)
    % Initialize a counter for points in front of cameras for each configuration
    count_front = [0,0,0,0];
    
    % Define the projection matrix for the first camera, assuming it's located at the origin
    P1 = [eye(3) zeros(3,1)];
    
    % Initialize an empty array to store 3D points for each camera configuration
    X =[];
    
    % Loop over the 4 possible camera configurations derived from the essential matrix decomposition
    for i = 1:4
        % Loop over each pair of corresponding points in the two images
        for j = 1:size(data{1},2)
            % Setup the matrix M for triangulation using DLT. This involves creating a system of equations that combines the projections from both cameras.
            M = [P1 -x1_n(:,j) zeros(3,1);
                 P{i} zeros(3,1) -x2_n(:,j)];
            % Perform Singular Value Decomposition (SVD) on M to solve for the 3D point
            [~,~,V] = svd(M);
            % The 3D point is given by the last column of V (homogeneous coordinates)
            v = V(:,end);
            % Store the 3D point in an array
            X{i}(1:4,j) = v(1:4,1);
            
            % Calculate the depth of the point relative to both camera centers. This involves checking the point's position along the camera's viewing direction.
            d0 = depth(P1,v(1:4,1)); % Depth in the first camera
            d1 = depth(P{i},v(1:4,1)); % Depth in the ith possible second camera
            
            % If the point is in front of both cameras (positive depth), increment the counter for visible points
            if sign(d0)>0 && sign(d1)>0
                count_front(i) = count_front(i)+1;
            end 
        end
    end
    
    % Find the camera configuration with the maximum number of visible 3D points
    [~,index] = max(count_front);
    
    % Return the best camera matrix (a) and its corresponding 3D points (b)
    a = P{index};
    b = X{index};
end
