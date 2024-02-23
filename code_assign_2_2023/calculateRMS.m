function eRMS = calculateRMS(Xmodel, x_measured, P)
    %%% Function to calculate RMS error %%%
    x_projected = P * Xmodel; % Project the model points
    x_projected = pflat(x_projected, 0); % Normalize the projected points, assuming pflat can handle this without plotting

    % Ensure x_measured is compatible for difference calculation
    % Assuming x_measured is already 2xN, matching the projected points' dimensions
    differences = x_measured(1:2,:) - x_projected(1:2,:); % Compute differences
    eRMS = sqrt(mean(sum(differences.^2, 1))); % Calculate RMS error
end
