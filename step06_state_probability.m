% Load the VV matrix from the file
load('D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing\valid_files\result_matlab_storage\VV.mat'); 

% Use only the 7th column of VV
VV_column = VV(:, 7);

% Calculate the number of states
num_states = max(VV_column); % Assuming states are numbered from 1 to max state number

% Initialize matrices to count transitions and remain counts
transition_counts = zeros(num_states, num_states);
remain_counts = zeros(1, num_states);

% Loop through each row to count transitions and remains
for i = 1:length(VV_column)-1
    current_state = VV_column(i);
    next_state = VV_column(i+1);
    if current_state == next_state
        remain_counts(current_state) = remain_counts(current_state) + 1;
    else
        transition_counts(current_state, next_state) = transition_counts(current_state, next_state) + 1;
    end
end

% Combine transition counts and remain counts into a single matrix
combined_matrix = transition_counts;
for i = 1:num_states
    combined_matrix(i, i) = remain_counts(i);
end

% Calculate the total count of all transitions and remains
total_transitions = sum(combined_matrix(:));

% Calculate the overall transition probabilities (sum should be 1)
transition_probabilities = combined_matrix / total_transitions;

% Display the combined matrix of counts
disp('Combined transition and remaining counts matrix:');
disp(combined_matrix);

% Display the transition probabilities
disp('Transition probabilities (total sum should be 1):');
disp(transition_probabilities);

% Verify that the sum of all probabilities is 1
total_probability = sum(transition_probabilities(:));
disp('Total probability (should be 1):');
disp(total_probability);
