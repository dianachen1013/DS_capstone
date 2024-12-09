% Initialize a matrix to store the transition probabilities
transit_state_probs = zeros(length(session_numbers), num_clusters, num_clusters);

% Loop over each session
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

    % Extract and store the probabilities of transitions between different states
    for state_from = 1:num_clusters
        for state_to = 1:num_clusters
            if state_from ~= state_to
                transit_state_probs(session_idx, state_from, state_to) = overall_transition_probs(state_from, state_to);
            end
        end
    end
end

% Display the transition probabilities between different states
disp('Transition Probabilities between Different States for each Session:');
for session_idx = 1:length(session_numbers)
    fprintf('Session %02d:\n', session_numbers(session_idx));
    disp(squeeze(transit_state_probs(session_idx, :, :)));
end

% Optionally, plot the transition probabilities for one state to another state
figure;
hold on;
for state_from = 1:num_clusters
    for state_to = 1:num_clusters
        if state_from ~= state_to
            plot(1:length(session_numbers), squeeze(transit_state_probs(:, state_from, state_to)), '-o', ...
                'LineWidth', 2, 'DisplayName', sprintf('State %d -> State %d', state_from, state_to));
        end
    end
end

% Add labels and title
xlabel('Session Number');
ylabel('Transition Probability');
title('Transition Probabilities between Different States Across Sessions');

% Set x-axis ticks to match session numbers
xticks(1:length(session_numbers));
xticklabels(arrayfun(@(x) sprintf('Session %02d', x), session_numbers, 'UniformOutput', false));

% Add a legend
legend('show');

% Adjust the plot for better visualization
grid on;
hold off;

% Save the plot if desired
saveas(gcf, fullfile(outputDir, 'transit_state_probabilities_plot.png'));
