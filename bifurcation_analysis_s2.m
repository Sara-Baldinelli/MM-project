% Parameters
k1 = 720; 
k2 = 720; 
k3 = 0.007; 
k4 = 0.007;
d_alpha_syn = 15;
K_alpha_syn = 8.5;

S3 = 1;   % Genetic predisposition factor
S4 = 1;   % Protein clearance factor
S1 = 0;   % Toxins exposure

% Function for steady-state ROS calculation
SteadyStateROS = @(ROS, S2) ...
    k1 * (1 + S1 + d_alpha_syn * ((k3 * ROS * S3 / (k4 * S4 * K_alpha_syn))^4) / ...
    (1 + (k3 * ROS * S3 / (k4 * S4 * K_alpha_syn))^4)) - k2 * ROS * S2;

S2_vals = linspace(0, 3, 1000); 
ROS_grid = linspace(0, 20, 5000); 

ROS_vals = nan(length(S2_vals), 3); % For up to three steady states
for i = 1:length(S2_vals)
    S2 = S2_vals(i);
    
    % Evaluate the function on the grid
    f_vals = arrayfun(@(ROS) SteadyStateROS(ROS, S2), ROS_grid);
    
    % Find where the function crosses zero
    crossings = find(diff(sign(f_vals)) ~= 0);
    
    % Refine roots using fzero and sort them
    refined_roots = nan(size(crossings));
    for j = 1:length(crossings)
        refined_roots(j) = fzero(@(ROS) SteadyStateROS(ROS, S2), ROS_grid(crossings(j)));
    end
    roots = sort(refined_roots); % Sort roots
    ROS_vals(i, 1:length(roots)) = roots; % Store up to three roots
end

% Remove discontinuities in branches
threshold = 0.5; % Threshold to detect jumps
for j = 1:3
    for i = 2:length(S2_vals)
        if abs(ROS_vals(i, j) - ROS_vals(i-1, j)) > threshold
            ROS_vals(i, j) = NaN; % Remove discontinuous points
        end
    end
end

% Plot bifurcation diagram
figure;
hold on;
plot(S2_vals, ROS_vals(:, 1), 'k-', 'LineWidth', 1.5); % Stable lower branch
plot(S2_vals, ROS_vals(:, 2), 'r--', 'LineWidth', 1.5); % Unstable branch
plot(S2_vals, ROS_vals(:, 3), 'k-', 'LineWidth', 1.5); % Stable upper branch
xlabel('S2 (Anti-oxidant mechanisms)');
ylabel('Steady-State ROS');
title('Bifurcation Diagram');
legend('Stable State', 'Unstable State');
grid on;
hold off;