% Define the directory containing the data files
data_dir = 'D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing\valid_files';
pattern = 'sub-*_ses-*_rsfmri_BP_space-fsnative_atlas-schaefer-400_desc-timeseries.txt';
files = dir(fullfile(data_dir, pattern));

% Check if files are found
if isempty(files)
    error('No files found matching the pattern.');
end

% Display all file paths
fprintf('List of files to be processed:\n');
for file_idx = 1:length(files)
    fprintf('%s\n', fullfile(files(file_idx).folder, files(file_idx).name));
end

num_samples = length(files);
num_rois = 400;
num_sessions = 11; 
num_participants = ceil(num_samples / num_sessions);

% Initialize the matrix to store the spatial clusters
spatial_clusters = zeros(num_participants, num_sessions, num_rois);

% Initialize DLI storage matrix (454*22*10)*400 = 99880*400
total_windows = 454; % Total number of windows per session
DLI_storage = zeros(total_windows * num_sessions * num_participants, num_rois);

% Start overall timer
overall_tic = tic;

% Process each file
for file_idx = 1:num_samples
    % Start file processing timer
    file_tic = tic;

    % Load the data from the current file
    load_tic = tic;
    file_path = fullfile(files(file_idx).folder, files(file_idx).name);
    data = load(file_path);
    load_time = toc(load_tic);

    % Ensure the file has 450 columns
    if size(data, 2) ~= 450
        error('Expected 450 columns for ROIs, but found %d', size(data, 2));
    end

    % Remove the first 50 columns
    remove_tic = tic;
    tc = data(:, 51:end);
    remove_time = toc(remove_tic);

    % Split data into left and right hemispheres
    split_tic = tic;
    left_brain = tc(:, 1:200);
    right_brain = tc(:, 201:400);
    split_time = toc(split_tic);

    % Calculate the global signal for each hemisphere
    global_signal_tic = tic;
    global_signal_left = mean(left_brain, 2);
    global_signal_right = mean(right_brain, 2);
    global_signal_time = toc(global_signal_tic);

    % Define the window length and step size
    window_length = 22; % Adjust based on TR and analysis requirement
    step_size = 1; % Typically 1 TR step

    % Initialize variables to store DLI results
    num_windows = floor((length(global_signal_left) - window_length) / step_size) + 1;
    DLI_matrix = zeros(num_windows, num_rois); % 400 ROIs

    % Calculate DLI
    DLI_tic = tic;
    for win = 1:num_windows
        win_start = (win - 1) * step_size + 1;
        win_end = win_start + window_length - 1;

        % Extract the window data
        window_roi_signal = tc(win_start:win_end,:);
        window_global_signal_left = global_signal_left(win_start:win_end);
        window_global_signal_right = global_signal_right(win_start:win_end);

        % Calculate correlation coefficients
        corr_left = corr(window_roi_signal, window_global_signal_left);
        corr_right = corr(window_roi_signal, window_global_signal_right);

        % Apply Fisher's z-transformation
        z_corr_left = 0.5 * log((1 + corr_left) ./ (1 - corr_left));
        z_corr_right = 0.5 * log((1 + corr_right) ./ (1 - corr_right));

        % Calculate DLI for the ROI within the window
        DLI_matrix(win, :) = z_corr_left - z_corr_right;
    end
    DLI_time = toc(DLI_tic);

    % Store the DLI data in DLI_storage
    participant_idx = ceil(file_idx / num_sessions);
    session_idx = mod(file_idx - 1, num_sessions) + 1;
    start_idx = (participant_idx - 1) * num_sessions * total_windows + (session_idx - 1) * total_windows + 1;
    end_idx = start_idx + num_windows - 1;
    DLI_storage(start_idx:end_idx, :) = DLI_matrix;

    % Create the correlation matrix from the DLI time series
    corr_matrix_tic = tic;
    correlation_matrix = corr(DLI_matrix);
    corr_matrix_time = toc(corr_matrix_tic);

    % Run Louvain algorithm 100 times to find the most balanced result
    best_std = Inf;
    most_balanced_community_labels = [];
    most_balanced_Q = -Inf;

    for i = 1:100
        [community_labels, Q] = community_louvain(correlation_matrix, 1, [], 'negative_asym');
        unique_labels = unique(community_labels);
        cluster_sizes = arrayfun(@(x) sum(community_labels == x), unique_labels);
        current_std = std(cluster_sizes);

        if current_std < best_std
            best_std = current_std;
            most_balanced_community_labels = community_labels;
            most_balanced_Q = Q;
        end
    end

    % Store the community labels for this sample
    spatial_clusters(participant_idx, session_idx, :) = most_balanced_community_labels';

    % Display the time taken for each step and for this file
    fprintf('Time taken to process %s:\n', file_path);
    fprintf('  Load data: %.2f seconds\n', load_time);
    fprintf('  Remove columns: %.2f seconds\n', remove_time);
    fprintf('  Split hemispheres: %.2f seconds\n', split_time);
    fprintf('  Calculate global signals: %.2f seconds\n', global_signal_time);
    fprintf('  Calculate DLI: %.2f seconds\n', DLI_time);
    fprintf('  Create correlation matrix: %.2f seconds\n', corr_matrix_time);
    fprintf('  Total time for file: %.2f seconds\n\n', toc(file_tic));
end

% Save the spatial clusters and DLI storage
save('spatial_clusters.mat', 'spatial_clusters');
save('DLI_storage.mat', 'DLI_storage');

% Display the total time taken
fprintf('Total time taken: %.2f seconds\n', toc(overall_tic));

% Display the results
disp('Most balanced community partitioning:');
disp(most_balanced_community_labels');
disp(['Modularity (Q) of most balanced result: ', num2str(most_balanced_Q)]);

% Display the number of clusters
num_clusters = numel(unique(most_balanced_community_labels));
disp(['Number of clusters: ', num_clusters]);

% Adjust the number of colors in the colormap to match the number of clusters
cluster_colors = distinguishable_colors(num_clusters);

% Initialize cell arrays to store the ROIs for each cluster
cluster_rois = cell(num_clusters, 1);

% Save the ROIs for each cluster
unique_labels = unique(most_balanced_community_labels);
for i = 1:numel(unique_labels)
    % Find the ROIs that belong to this cluster
    cluster_rois{i} = find(most_balanced_community_labels == unique_labels(i));
    
    % Dynamically create variable names for clusters
    cluster_name = sprintf('cluster%d', i);
    
    % Assign the ROIs to the variable
    assignin('base', cluster_name, cluster_rois{i});
    
    % Save each cluster's ROIs into a separate file
    save(fullfile('D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/result_matlab_storage', ...
        sprintf('%s.mat', cluster_name)), cluster_name);
    
    % Display the ROIs in the command window
    disp(['Cluster ', num2str(i), ': ', num2str(numel(cluster_rois{i})), ' ROIs']);
    fprintf('ROIs in Cluster %d: ', i);
    fprintf('%d ', cluster_rois{i});
    fprintf('\n');
end

% Optionally, save all clusters together in one file
save('D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files/result_matlab_storage/all_clusters.mat', 'cluster_rois');

% Load the Schaefer 2018 atlas NIfTI file (adjust the file name and path as needed)
nifti_file = 'D:/ibp/schaefer/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii.gz';
nifti_data = load_untouch_nii(nifti_file);

% Extract the atlas data
atlas_data = nifti_data.img;

% Initialize variables
roi_coords = zeros(num_rois, 3);
roi_labels = unique(atlas_data(:));
roi_labels(roi_labels == 0) = []; % Remove the background label

% Calculate the centroid of each ROI
for i = 1:num_rois
    [x, y, z] = ind2sub(size(atlas_data), find(atlas_data == roi_labels(i)));
    roi_coords(i, :) = [mean(x), mean(y), mean(z)];
end

% Convert voxel coordinates to MNI coordinates (if necessary)
% Adjust the transformation matrix based on your specific NIfTI file
voxel_size = [1, 1, 1]; % Assuming 1mm isotropic voxels
origin = [90, 126, 72]; % Example origin for FSL MNI152 template
roi_coords_mni = bsxfun(@times, roi_coords - 1, voxel_size) + origin;

% Save the ROI coordinates
save('schaefer_2018_400_coords.mat', 'roi_coords_mni', 'roi_labels');

% Load the Schaefer atlas ROI coordinates
load('schaefer_2018_400_coords.mat');  % Load ROI coordinates (roi_coords_mni) and labels (roi_labels)

% Define colors for each cluster, supporting up to num_clusters clusters
cluster_colors = distinguishable_colors(num_clusters);

% Plot each ROI colored by its cluster label
figure;
hold on;
for i = 1:num_clusters
    cluster_indices = find(most_balanced_community_labels == unique_labels(i));
    scatter3(roi_coords_mni(cluster_indices, 1), roi_coords_mni(cluster_indices, 2), roi_coords_mni(cluster_indices, 3), ...
        36, cluster_colors(i,:), 'filled');
end

% Enhance visualization
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Brain Clusters Visual');
view(3);
axis equal;
grid on;
hold off;

% Combine all clusters into a single index vector in the desired order
cluster_order = [];
for i = 1:num_clusters
    cluster_order = [cluster_order; cluster_rois{i}(:)];
end

% Reorder the correlation matrix according to the cluster order
ordered_correlation_matrix = correlation_matrix(cluster_order, cluster_order);

% Plot the reordered correlation matrix
figure;
imagesc(ordered_correlation_matrix);
colorbar;
title('Correlation Matrix Ordered by Clusters');
xlabel('ROIs');
ylabel('ROIs');

% Change the colormap to 'jet'
colormap(jet);

% Add lines to separate the clusters visually
hold on;
% Calculate cumulative sum of cluster sizes to determine separation lines
cluster_sizes = cellfun(@numel, cluster_rois);
separation_lines = cumsum(cluster_sizes);

% Draw lines to separate clusters
for i = 1:length(separation_lines)-1
    line_position = separation_lines(i) + 0.5;
    plot([line_position line_position], ylim, 'k-', 'LineWidth', 1);
    plot(xlim, [line_position line_position], 'k-', 'LineWidth', 1);
end

hold off;

% Save the reordered correlation matrix to Excel
filename = 'correlation_matrix.xlsx';
output_dir = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files/result_matlab_storage';
filepath = fullfile(output_dir, filename);

if isnumeric(correlation_matrix)
    writematrix(correlation_matrix, filepath);
    disp(['Correlation matrix saved successfully to ', filepath]);
else
    error('The correlation_matrix must be a numeric matrix to save as an Excel file.');
end

% Define the file path where you want to save second_level_centers.mat
save_path = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files/result_matlab_storage/cluster1.mat';

% Save the variable second_level_centers to the specified file
save(save_path, 'cluster1');

% Define the file path where you want to save second_level_centers.mat
save_path = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files/result_matlab_storage/cluster2.mat';

% Save the variable second_level_centers to the specified file
save(save_path, 'cluster2');

% Define the file path where you want to save second_level_centers.mat
save_path = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files/result_matlab_storage/cluster3.mat';

% Save the variable second_level_centers to the specified file
save(save_path, 'cluster3');

% Define the file path where you want to save second_level_centers.mat
save_path = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files/result_matlab_storage/cluster4.mat';

% Save the variable second_level_centers to the specified file
save(save_path, 'cluster4');

% Define the file path where you want to save second_level_centers.mat
save_path = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files/result_matlab_storage/DLI_storage.mat';

% Save the variable second_level_centers to the specified file
save(save_path, 'DLI_storage');

