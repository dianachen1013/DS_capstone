% Start timing the script
total_tic = tic;

% Define the directory containing the data files
data_dir = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files';
pattern = 'sub-*_ses-*_rsfmri_BP_space-fsnative_atlas-schaefer-400_desc-timeseries.txt';
files = dir(fullfile(data_dir, pattern));
out_dir = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files/result_matlab_storage';

% Ensure the directory exists, if not, create it
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Check if files are found
if isempty(files)
    error('No files found matching the pattern.');
end

num_samples = length(files);
num_rois = 400;
window_length = 22; % 22 seconds
step_size = 1; % 1 TR step
TR = 1; % Repetition Time in seconds
window_length_TRs = window_length / TR; % Window length in TRs

% Initialize storage for all participants
all_participant_states = [];
participant_indices = [];
DLI_storage = []; % Initialize storage for DLI values
state_labels_storage = []; % Initialize storage for state labels

% Loop through each file to process data
for file_idx = 1:num_samples
    participant_tic = tic;
    
    % Load the data from the current file
    file_path = fullfile(files(file_idx).folder, files(file_idx).name);
    data = load(file_path);

    % Ensure the file has 450 columns
    if size(data, 2) ~= 450
        error('Expected 450 columns for ROIs, but found %d', size(data, 2));
    end

    % Remove the first 50 columns
    tc = data(:, 51:end);

    % Apply band-pass filtering to the first ROI data
    %low_cutoff = 0.01; % Lower cutoff frequency
    %high_cutoff = 0.1; % Upper cutoff frequency
    %sampling_rate = 1; % Sampling rate in Hz (assuming TR = 1 second)
    %tc_after = bandpass(tc, [low_cutoff high_cutoff], sampling_rate);

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


    %for roi = 1:num_rois
    %    roi_signal = tc(:, roi);
       % Calculate DLI
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
   

    % Store the DLI matrix for this participant
    DLI_storage = [DLI_storage; DLI_matrix];

    % First-level clustering using k-means for this participant
    clustering_tic = tic;
    num_clusters_first_level = 11; % Adjust this number as needed
    [first_level_labels, first_level_centers] = kmeans(DLI_matrix, num_clusters_first_level, 'Distance', 'cosine');
    clustering_time = toc(clustering_tic);
    
    % Store the clustered states for this participant
    all_participant_states = [all_participant_states; first_level_centers];
    participant_indices = [participant_indices; repmat(file_idx, num_clusters_first_level, 1)];
    
    % Store state labels for this participant
    state_labels_storage = [state_labels_storage; first_level_labels];
    
    % Print timing information for this participant
    fprintf('Processed participant %d in %.2f seconds (clustering time: %.2f seconds)\n', file_idx, toc(participant_tic), clustering_time);
end

% Save the first-level clustering results
first_level_clusters_path = fullfile(out_dir, 'first_level_clusters.mat');
save(first_level_clusters_path, 'all_participant_states', 'participant_indices', 'DLI_storage', 'state_labels_storage');

% Second-level clustering using k-means for all participants
second_clustering_tic = tic;
num_clusters_second_level = 300; % Adjust this number as needed
[second_level_labels, second_level_centers] = kmeans(all_participant_states, num_clusters_second_level, 'Distance', 'cosine');
second_clustering_time = toc(second_clustering_tic);

% Store second-level clustering results
second_level_clusters_path = fullfile(out_dir, 'second_level_clusters.mat');
save(second_level_clusters_path, 'second_level_centers', 'second_level_labels', 'participant_indices');

% Count the number of states in each second-level cluster
state_counts = histcounts(second_level_labels, num_clusters_second_level);

% Display the number of states in each second-level cluster
total_states = 0;
disp('Number of states in each second-level cluster:');
for cluster_idx = 1:num_clusters_second_level
    fprintf('Cluster %d: %d states\n', cluster_idx, state_counts(cluster_idx));
    total_states = total_states + state_counts(cluster_idx);
end

% Display the total number of states
fprintf('Total number of states in all clusters: %d\n', total_states);

% Store each group's each ROI's DLI values
group_dli_tic = tic;
group_dli_values = cell(num_clusters_second_level, 1);
for cluster_idx = 1:num_clusters_second_level
    group_dli_values{cluster_idx} = DLI_storage(second_level_labels == cluster_idx, :);
end
group_dli_time = toc(group_dli_tic);

% Save the group DLI values
group_dli_values_path = fullfile(out_dir, 'group_dli_values.mat');
save(group_dli_values_path, 'group_dli_values');

% End timing the script
total_time = toc(total_tic);

% Print timing information for each step
fprintf('First-level clustering time: %.2f seconds\n', clustering_time);
fprintf('Second-level clustering time: %.2f seconds\n', second_clustering_time);
fprintf('Grouping DLI values time: %.2f seconds\n', group_dli_time);
fprintf('Total script time: %.2f seconds\n', total_time);

% Display the paths where files are saved
fprintf('Files saved in the following paths:\n');
fprintf('First-level clustering results: %s\n', first_level_clusters_path);
fprintf('Second-level clustering results: %s\n', second_level_clusters_path);
fprintf('Group DLI values: %s\n', group_dli_values_path);

% Define the file path where you want to save second_level_centers.mat
save_path = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files/second_level_centers.mat';

% Save the variable second_level_centers to the specified file
save(save_path, 'second_level_centers');