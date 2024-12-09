% Define the directory containing the data files
data_dir = 'D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing';
pattern = 'sub-*_ses-*_rsfmri_BP_space-fsnative_atlas-schaefer-400_desc-timeseries.txt';
files = dir(fullfile(data_dir, pattern));

% Check if files are found
if isempty(files)
    error('No files found matching the pattern.');
end

% Initialize variables
required_sessions = [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12];
subjects_sessions = containers.Map('KeyType', 'int32', 'ValueType', 'any');

% Process each file to extract subject and session information
for file_idx = 1:length(files)
    file_name = files(file_idx).name;
    tokens = regexp(file_name, 'sub-(\d+)_ses-(\d+)', 'tokens');
    if ~isempty(tokens)
        subject_id = str2double(tokens{1}{1});
        session_id = str2double(tokens{1}{2});
        
        if isKey(subjects_sessions, subject_id)
            subjects_sessions(subject_id) = unique([subjects_sessions(subject_id), session_id]);
        else
            subjects_sessions(subject_id) = session_id;
        end
    end
end

% Find subjects with all required sessions
subject_ids = keys(subjects_sessions);
valid_subjects = [];
for i = 1:length(subject_ids)
    sessions = subjects_sessions(subject_ids{i});
    if all(ismember(required_sessions, sessions))
        valid_subjects = [valid_subjects, subject_ids{i}];
    end
end

% Define the directory for abandoned files
abandoned_dir = fullfile(data_dir, 'abandoned_files');
if ~exist(abandoned_dir, 'dir')
    mkdir(abandoned_dir);
end

% Define the directory for valid files
valid_dir = fullfile(data_dir, 'valid_files');
if ~exist(valid_dir, 'dir')
    mkdir(valid_dir);
end

% Initialize counters
abandoned_count = 0;
valid_count = 0;

% Copy files to the appropriate directory based on subject validity
for file_idx = 1:length(files)
    file_name = files(file_idx).name;
    tokens = regexp(file_name, 'sub-(\d+)_ses-(\d+)', 'tokens');
    if ~isempty(tokens)
        subject_id = str2double(tokens{1}{1});
        session_id = str2double(tokens{1}{2});
        if ismember(subject_id, valid_subjects) && ismember(session_id, required_sessions)
            copyfile(fullfile(files(file_idx).folder, file_name), fullfile(valid_dir, file_name));
            valid_count = valid_count + 1;
        else
            copyfile(fullfile(files(file_idx).folder, file_name), fullfile(abandoned_dir, file_name));
            abandoned_count = abandoned_count + 1;
        end
    end
end

% Display the number of files in each directory
fprintf('Number of files in valid_files: %d\n', valid_count);
fprintf('Number of files in abandoned_files: %d\n', abandoned_count);

% Save the list of valid subjects to a file
save(fullfile(data_dir, 'valid_subjects.mat'), 'valid_subjects');

% Verify that the session numbers are correct for each valid subject
disp('Valid session numbers for each subject:');
for i = 1:length(valid_subjects)
    subject_id = valid_subjects(i);
    sessions = subjects_sessions(subject_id);
    fprintf('Subject %d: ', subject_id);
    disp(sessions);
end
