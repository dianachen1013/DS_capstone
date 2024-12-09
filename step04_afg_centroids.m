% Load the second_level_centers and VV
load('D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing\valid_files\result_matlab_storage\second_level_centers.mat');
load('D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing\valid_files\result_matlab_storage\VV.mat');

% Extract the 7th column from VV, representing the temporal clustering result
temporal_clustering_result = VV(:, 7);

% Number of states (in your case, 3 states marked by 1, 2, and 3)
num_states = 3;

% Initialize the centroids matrix
[num_groups, num_ROIs] = size(second_level_centers);
centroids = zeros(num_states, num_ROIs);

% Calculate the centroid for each state
for state = 1:num_states
    % Find the indices of the groups belonging to the current state
    state_indices = find(temporal_clustering_result == state);
    
    % Compute the average pattern for the current state
    centroids(state, :) = mean(second_level_centers(state_indices, :), 1);
end

% Display the centroids
disp('Centroids for each state:');
disp(centroids);

% Define the file path where you want to save second_level_centers.mat
save_path = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files/result_matlab_storage/centroids.mat';

% Save the variable second_level_centers to the specified file
save(save_path, 'centroids');