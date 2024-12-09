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

% First loop to find global min and max values
for session_num = session_numbers
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

    % Display the overall transition probabilities for the current session
    fprintf('Overall Transition Probabilities for %s:\n', session_id);
    disp(overall_transition_probs);
end

% Display the global min and max values for color scaling
fprintf('Global Min: %.4f, Global Max: %.4f\n', global_min, global_max);

% Second loop to plot the heatmaps with consistent color scaling
for session_num = session_numbers
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

    % Step 1: Define a helper function directly within the script to convert hex to RGB
    hex2rgb = @(hex) reshape(sscanf(hex(2:end), '%2x%2x%2x'), 1, 3) / 255;

    % Step 2: Define your hex color codes
    hex_colors = {'#afd1e9', '#b0bee2', '#b299d4'};

    % Step 3: Convert the hex color codes to an RGB colormap
    custom_colormap = zeros(length(hex_colors), 3);
    for i = 1:length(hex_colors)
        custom_colormap(i, :) = hex2rgb(hex_colors{i});
    end

    % Step 4: Apply the custom colormap and plot the heatmap
    figure;
    colormap(custom_colormap); % Apply the custom colormap

    % Plot the heatmap for overall transition probabilities
    imagesc(overall_transition_probs);

    % Add colorbar to indicate the scale of probabilities
    colorbar;

    % Set the same color axis range for all heatmaps
    clim([global_min global_max]);

    % Add labels and title
    xlabel('To State');
    ylabel('From State');
    title(['Overall State Transition Probabilities for ' session_id]);

    % Set axis tick marks and labels
    set(gca, 'XTick', 1:num_clusters, 'YTick', 1:num_clusters);
    set(gca, 'XTickLabel', {'State 1', 'State 2', 'State 3'}, 'YTickLabel', {'State 1', 'State 2', 'State 3'});

    % Add numeric labels on the heatmap cells for clarity
    for i = 1:num_clusters
        for j = 1:num_clusters
            % Choose text color based on background intensity
            if overall_transition_probs(i, j) > 0.5 * max(overall_transition_probs(:))
                text_color = 'black';
            else
                text_color = 'white';
            end
            text(j, i, sprintf('%.2f', overall_transition_probs(i, j)), ...
                'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle', 'Color', text_color, 'FontWeight', 'bold');
        end
    end

    % Add grey dashed lines to separate the states
    hold on; % Keep the heatmap so that the lines are plotted on top
    for k = 1.5:(num_clusters-0.5)
        % Vertical lines
        line([k k], [0.5 num_clusters+0.5], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5);
        % Horizontal lines
        line([0.5 num_clusters+0.5], [k k], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1.5);
    end
    hold off;

    % Adjust the plot for better visualization
    axis square;
    set(gca, 'FontSize', 12);

    % Save the plot
    saveas(gcf, fullfile(outputDir, ['transition_probs_heatmap_' session_id '.png']));

    % Display the transition probabilities for the session
    fprintf('Transition probabilities for %s:\n', session_id);
    disp(overall_transition_probs);
end
