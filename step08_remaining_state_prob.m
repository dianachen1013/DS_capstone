% Define parameters
data_dir = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files';
newDataDir = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/processed-data';
outputDir = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/output-data';
window_length = 22;
step_size = 1;
num_clusters = 3;

% Load the centroids variable
centroidsPath = 'D:/ibp/DLI-github/sample-data/实验数据/final_data/ChenDanQing/valid_files/result_matlab_storage/centroids.mat';
load(centroidsPath, 'centroids');

% Create directories if they do not exist
if ~exist(newDataDir, 'dir')
    mkdir(newDataDir);
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Define the available session numbers
session_numbers = [01, 02, 03, 04, 05, 06, 07, 08, 10, 11, 12];

% Initialize global min and max for color scale
global_min = Inf;
global_max = -Inf;

% Initialize matrix to store probabilities of remaining in the same state
same_state_probs = zeros(length(session_numbers), num_clusters);

% First loop to find global min and max values
for session_idx = 1:length(session_numbers)
    session_num = session_numbers(session_idx);
    session_id = sprintf('ses-%02d', session_num);

    % List all files for the current session
    pattern = ['sub-*_' session_id '_rsfmri_BP_space-fsnative_atlas-schaefer-400_desc-timeseries.txt'];
    files = dir(fullfile(data_dir, pattern));

    % Check if files are found
    if isempty(files)
        warning('No files found matching the pattern for %s.', session_id);
        continue;  % Skip this session if no files are found
    end

    num_samples = length(files);
    num_rois = 400;

    % Initialize variables to accumulate transition counts across participants
    average_transition_counts = zeros(num_clusters, num_clusters);

    % Process each file
    for file_idx = 1:num_samples
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
        left_brain = tc(:, 1:200);
        right_brain = tc(:, 201:400);

        % Calculate the global signal for each hemisphere
        global_signal_left = mean(left_brain, 2);
        global_signal_right = mean(right_brain, 2);

        % Initialize variables to store DLI results
        num_windows = floor((length(global_signal_left) - window_length) / step_size) + 1;
        DLI_matrix = zeros(num_windows, num_rois); % 400 ROIs

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

        % Perform k-means clustering on the DLI matrix using the provided centroids
        [idx, ~] = kmeans(DLI_matrix, num_clusters, 'Start', centroids);

        % Calculate the transition counts for the current participant
        transitions = zeros(num_clusters, num_clusters);
        for i = 1:(length(idx) - 1)
            transitions(idx(i), idx(i + 1)) = transitions(idx(i), idx(i + 1)) + 1;
        end

        % Accumulate transition counts across participants
        average_transition_counts = average_transition_counts + transitions;
    end

    % Calculate the average transition counts by dividing by the number of participants
    average_transition_counts = average_transition_counts / num_samples;

    % Calculate the overall transition probabilities
    total_transitions = sum(average_transition_counts(:));
    overall_transition_probs = average_transition_counts / total_transitions;

    % Update global min and max
    global_min = min(global_min, min(overall_transition_probs(:)));
    global_max = max(global_max, max(overall_transition_probs(:)));

    % Extract and store the probabilities of remaining in the same state
    for state = 1:num_clusters
        same_state_probs(session_idx, state) = overall_transition_probs(state, state);
    end
end

% Display the global min and max values for color scaling
fprintf('Global Min: %.4f, Global Max: %.4f\n', global_min, global_max);

% Display the probabilities of remaining in the same state for each session
fprintf('Probabilities of remaining in the same state for each session:\n');
disp(same_state_probs);


% Plot the probabilities of remaining in the same state for each session
figure;
hold on;

% Define colors for each state line
state_colors = lines(num_clusters);

% Plot each state's probabilities over the sessions
for state = 1:num_clusters
    plot(1:length(session_numbers), same_state_probs(:, state), '-o', ...
        'LineWidth', 2, 'Color', state_colors(state, :));
    
    % Add annotation to each curve
    text(length(session_numbers) + 0.1, same_state_probs(end, state), ...
        sprintf('State %d', state), 'FontSize', 12, 'Color', state_colors(state, :));
end

% Add labels and title
xlabel('Session Number');
ylabel('Probability of Remaining in the Same State');
title('Probabilities of Remaining in the Same State Across Sessions');

% Set x-axis ticks to match session numbers
xticks(1:length(session_numbers));
xticklabels(arrayfun(@(x) sprintf('Session %02d', x), session_numbers, 'UniformOutput', false));

% Add legend for clarity
legend(arrayfun(@(x) sprintf('State %d', x), 1:num_clusters, 'UniformOutput', false), 'Location', 'best');

% Adjust the plot for better visualization
grid on;
hold off;

% Save the plot if desired
saveas(gcf, fullfile(outputDir, 'same_state_probabilities_plot.png'));