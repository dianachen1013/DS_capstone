% Start timing the script
total_tic = tic;

% Define the directory containing the data files
data_dir = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files';
pattern = 'sub-*_ses-*_rsfmri_BP_space-fsnative_atlas-schaefer-400_desc-timeseries.txt';
files = dir(fullfile(data_dir, pattern));

% Adjacency matrix construction for GCAFG
adj_matrix_tic = tic;
adj_matrix = corr(second_level_centers');
adj_matrix_time = toc(adj_matrix_tic);
fprintf('Adjacency matrix construction completed in %.2f seconds\n', adj_matrix_time);

% Scale parameter
scale = 0.1:0.1:1.5; % Example scale parameter, with smaller steps

% Apply GCAFG
GCAFG_tic = tic;
VV = GCAFG(adj_matrix, scale);
GCAFG_time = toc(GCAFG_tic);
fprintf('GCAFG application completed in %.2f seconds\n', GCAFG_time);

% Save the GCAFG results
save('GCAFG_results.mat', 'VV', 'scale');

% End timing the script
total_time = toc(total_tic);

% Print timing information for each step
fprintf('Adjacency matrix construction time: %.2f seconds\n', adj_matrix_time);
fprintf('GCAFG application time: %.2f seconds\n', GCAFG_time);
fprintf('Total script time: %.2f seconds\n', total_time);

% Define the number of clusters
num_clusters_second_level = size(VV, 2);

% Initialize variables to store ROIs for each cluster
cluster_rois = cell(1, num_clusters_second_level);

% Assign ROIs to clusters based on the VV matrix
for comm_idx = 1:num_clusters_second_level
    cluster_rois{comm_idx} = find(VV(:, comm_idx) == 1);
end

% Display the results for each community
for comm_idx = 1:num_clusters_second_level
    community_states = cluster_rois{comm_idx};
    fprintf('Community %d: %d states\n', comm_idx, length(community_states));
    
    % Initialize mean DLI values for each cluster
    mean_dli_cluster = [];

    % Collect DLI values for the states in this community
    for state_idx = community_states'
        if state_idx <= num_clusters_second_level % Ensure state_idx is within valid bounds
            mean_dli_cluster = [mean_dli_cluster; group_dli_values{state_idx}(:, community_states)];
        end
    end
    
    % Calculate and display the mean DLI values for this community
    fprintf('Community %d Mean DLI value: %.4f\n', comm_idx, mean(mean_dli_cluster(:)));
end

% Define the file path where you want to save second_level_centers.mat
save_path = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing//valid_files/result_matlab_storage/VV.mat';

% Save the variable second_level_centers to the specified file
save(save_path, 'VV');

% Define the file path where you want to save second_level_centers.mat
save_path = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files/result_matlab_storage/second_level_centers.mat';

% Save the variable second_level_centers to the specified file
save(save_path, 'second_level_centers');

disp(VV)
histcounts(VV)/15