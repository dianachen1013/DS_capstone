% Define the directory containing the data files
data_dir = 'D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing\valid_files';
pattern = 'sub-*_ses-*_rsfmri_BP_space-fsnative_atlas-schaefer-400_desc-timeseries.txt';
files = dir(fullfile(data_dir, pattern));

% Check if files are found
if isempty(files)
    error('No files found matching the pattern.');
end

% Filter out files from sub-05
files = files(~contains({files.name}, 'sub-05'));

% Recalculate the number of participants and samples after filtering
num_samples = length(files);
num_rois = 400;
num_sessions = 11; 
num_participants = num_samples / num_sessions;  % This should be an integer if sub-05 is correctly excluded

% Initialize the matrix to store the DLI data as 3D matrix
dli_3d = zeros(num_participants, num_rois, num_sessions);

% Initialize DLI storage matrix (454*22*10)*400 = 99880*400
total_windows = 454; % Total number of windows per session
DLI_storage = zeros(total_windows * num_sessions * num_participants, num_rois);

% Start overall timer
overall_tic = tic;

% Process each file
for file_idx = 1:num_samples
    % Update participant and session indices
    session_idx = mod(file_idx - 1, num_sessions) + 1;
    participant_idx = ceil(file_idx / num_sessions);

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

    % Store the average DLI for the session in the 3D matrix
    dli_3d(participant_idx, :, session_idx) = mean(DLI_matrix, 1);

    % Optionally, store all DLI data in DLI_storage if needed
    % This will store all the DLI matrices from different windows for further analysis
    start_idx = (file_idx - 1) * total_windows + 1;
    end_idx = start_idx + total_windows - 1;
    DLI_storage(start_idx:end_idx, :) = DLI_matrix;

    % Display processing time for this file
    fprintf('Processed file %d of %d (Participant: %d, Session: %d) in %.2f seconds\n', ...
        file_idx, num_samples, participant_idx, session_idx, toc(file_tic));
end

% Display overall processing time
fprintf('All files processed in %.2f seconds\n', toc(overall_tic));
%% Step 1: Define the Clusters
% Define clusters (cluster1, cluster2, cluster3, cluster4)
% Example: cluster1 = [1, 2, 3, 4, 5]; % Add your actual cluster indices
clusters = {cluster1, cluster2, cluster3, cluster4};
nClusters = numel(clusters);

%% Step 2: Calculate Average DLI Per Cluster
% Initialize a matrix to store the average DLI per cluster [participants, clusters, sessions]
average_dli_per_cluster = zeros(size(dli_3d, 1), nClusters, size(dli_3d, 3));

% Calculate the mean DLI for each cluster
for participant = 1:size(dli_3d, 1)
    for session = 1:size(dli_3d, 3)
        for cluster = 1:nClusters
            cluster_indices = clusters{cluster};  % Get the indices for the current cluster
            % Calculate the mean DLI for the current cluster, participant, and session
            average_dli_per_cluster(participant, cluster, session) = mean(dli_3d(participant, cluster_indices, session), 2);
        end
    end
end

%% Step 3: Perform Paired t-Tests
% Initialize a matrix to store p-values for the t-tests [sessions-1, clusters]
p_values = zeros(size(dli_3d, 3) - 1, nClusters);

% Perform paired t-tests comparing each session with the baseline (first session)
for cluster = 1:nClusters
    for session = 2:size(dli_3d, 3)
        % Extract data for all participants for the baseline and current session
        baseline_data = squeeze(average_dli_per_cluster(:, cluster, 1));  % Baseline session (1st session)
        comparison_data = squeeze(average_dli_per_cluster(:, cluster, session));  % Current session

        % Perform a paired t-test between baseline and the current session
        [~, p_values(session - 1, cluster)] = ttest(baseline_data, comparison_data);
    end
end

% Display the p-values
disp('Paired t-test p-values for each session compared to baseline for each cluster:');
disp(p_values);