

% Define the directory containing the data files
dataDir = 'D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing\valid_files';

% List all files in the directory to verify the contents
allFiles = dir(fullfile(dataDir, '*.txt'));
disp('All files in the directory:');
for i = 1:length(allFiles)
    disp(allFiles(i).name);
end

% List all files matching the pattern for session 1
pattern = 'sub-*_ses-01_rsfmri_BP_space-fsnative_atlas-schaefer-400_desc-timeseries.txt';
files = dir(fullfile(dataDir, pattern));
disp('Files matching the pattern for session 1:');
for i = 1:length(files)
    disp(files(i).name);
end

% Initialize variables to store combined results
combined_durations = [];
combined_up_counts = [];
combined_down_counts = [];

% Check if files are found
if isempty(files)
    error('No files found matching the pattern for session 1.');
end

for fileIdx = 1:length(files)
    % Extract subject and session information from the filename
    [~, fileName, ~] = fileparts(files(fileIdx).name);
    tokens = regexp(fileName, 'sub-(\d+)_ses-(\d+)', 'tokens');
    sampleID = str2double(tokens{1}{1});
    sessionID = str2double(tokens{1}{2});

    % Display the current file being processed
    disp(['Processing file: ', fileName]);

    % Load the data from the current file
    filePath = fullfile(files(fileIdx).folder, files(fileIdx).name);
    data = load(filePath);

    % Ensure the file has 450 columns
    if size(data, 2) ~= 450
        error('Expected 450 columns for ROIs, but found %d', size(data, 2));
    end

    % Remove the first 50 columns
    tc = data(:, 51:end);

    % Extract the first ROI's data (column 1 after removing the first 50 columns)
    first_roi_data = tc(:, 1);

    % Define the window length and step size
    window_length = 22; % Adjust this value based on the TR and analysis requirement
    step_size = 1; % Typically 1 TR step

    % Initialize variables to store DLI results
    num_windows = floor((length(first_roi_data) - window_length) / step_size) + 1;
    DLI = zeros(num_windows, 1);

    % Split the data into left and right hemispheres for DLI calculation
    left_brain = tc(:, 1:200);
    right_brain = tc(:, 201:400);

    for win = 1:num_windows
        % Define the window range
        win_start = (win - 1) * step_size + 1;
        win_end = win_start + window_length - 1;

        % Extract the window data
        window_left_brain = left_brain(win_start:win_end, :);
        window_right_brain = right_brain(win_start:win_end, :);
        window_roi_signal = first_roi_data(win_start:win_end);

        % Calculate the global signal for each hemisphere within the window
        global_signal_left = mean(window_left_brain, 2);
        global_signal_right = mean(window_right_brain, 2);

        % Calculate DLI for the first ROI within the window
        corr_left = corr(window_roi_signal, global_signal_left);
        corr_right = corr(window_roi_signal, global_signal_right);
        DLI(win) = corr_left - corr_right;
    end

    % Display the DLI values for the current file
    %disp('DLI values:');
    %disp(DLI);

    % Calculate basic statistics for the DLI values
    mean_DLI = mean(DLI);
    std_DLI = std(DLI); % This is the Laterality Fluctuation (LF)

    % Display the mean and standard deviation
    fprintf('Mean DLI: %.4f\n', mean_DLI);
    fprintf('Standard Deviation of DLI (LF): %.4f\n', std_DLI);

    % Calculate Mean Laterality Index (MLI)
    MLI = mean(DLI);

    % Calculate Laterality Fluctuation (LF) as standard deviation of DLI
    LF = std(DLI);

    % Find the indices where DLI is outside the standard deviation range
    outside_std_indices = find(DLI > MLI + LF | DLI < MLI - LF);

    % Initialize counters for upward and downward movements
    up_count = 0;
    down_count = 0;

    % Group continuous periods where DLI is outside the standard deviation range
    durations = [];
    if ~isempty(outside_std_indices)
        % Initialize variables
        start_idx = outside_std_indices(1);
        duration_start = start_idx;

        for i = 2:length(outside_std_indices)
            if outside_std_indices(i) ~= outside_std_indices(i-1) + 1
                % End of a continuous period
                end_idx = outside_std_indices(i-1);
                durations = [durations; sampleID, sessionID, duration_start, end_idx, end_idx - duration_start + 1, DLI(duration_start), DLI(end_idx)];
                duration_start = outside_std_indices(i);
                %Count upward and downward movements
                if DLI(duration_start) < MLI
                    up_count = up_count + 1;
                else
                    down_count = down_count + 1;
                end
                duration_start = outside_std_indices(i);
            end
        end
        % Add the last period
        end_idx = outside_std_indices(end);
        durations = [durations; sampleID, sessionID, duration_start, end_idx, end_idx - duration_start + 1, DLI(duration_start), DLI(end_idx)];
        % Count upward and downward movements
        if DLI(duration_start) < MLI
            up_count = up_count + 1;
        else
            down_count = down_count + 1;
        end
    end

    if ~isempty(durations)
        % Create a table with the results for the current file
        outliers_table = array2table(durations, 'VariableNames', {'SampleID', 'SessionID', 'Start_Time_Window', 'End_Time_Window', 'Duration', 'DLI_Start', 'DLI_End'});

        % Add upward and downward counts to the table
        outliers_table.Up_Count = repmat(up_count, height(outliers_table), 1);
        outliers_table.Down_Count = repmat(down_count, height(outliers_table), 1);

        % Display the table for the current file
        disp(['Periods where DLI is outside the standard deviation range for file: ', fileName]);
        disp(outliers_table);

        % Save the table to a CSV file
        writetable(outliers_table, [fileName, '_DLI_Outliers_Grouped.csv']);

        % Append the results to the combined durations
        combined_durations = [combined_durations; durations];
        combined_up_counts = [combined_up_counts; repmat(up_count, size(durations, 1), 1)];
        combined_down_counts = [combined_down_counts; repmat(down_count, size(durations, 1), 1)];
    end
end

% Ensure that combined_durations is not empty
if ~isempty(combined_durations)
    % Create a table with the combined results
    combined_outliers_table = array2table(combined_durations, 'VariableNames', {'SampleID', 'SessionID', 'Start_Time_Window', 'End_Time_Window', 'Duration', 'DLI_Start', 'DLI_End'});
    combined_outliers_table.Up_Count = combined_up_counts;
    combined_outliers_table.Down_Count = combined_down_counts;

    % Display the combined table
    disp('Combined periods where DLI is outside the standard deviation range:');
    disp(combined_outliers_table);

    % Save the combined table to a CSV file
    writetable(combined_outliers_table, 'Combined_DLI_Outliers_Grouped_Session1.csv');
else
    disp('No periods where DLI is outside the standard deviation range found for any file.');
end

% Assuming 'combined_outliers_table' is your table with the results

% Extract durations
durations = combined_outliers_table.Duration;

% Plot histogram
figure;
histogram(durations, 'BinWidth', 1);
xlabel('Duration (time windows)');
ylabel('Frequency');
title('Histogram of Durations where DLI is Outside the Standard Deviation');
grid on;


% Assuming 'combined_outliers_table' is your table with the results

% Create a box plot of durations by sample
figure;
boxplot(combined_outliers_table.Duration, combined_outliers_table.SampleID);
xlabel('Sample ID');
ylabel('Duration (time windows)');
title('Box Plot of Durations where DLI is Outside the Standard Deviation by Sample');
grid on;