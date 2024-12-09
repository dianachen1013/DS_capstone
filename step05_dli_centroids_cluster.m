% Load the centroids and DLI_storage matrices from the specified paths
load('D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing\valid_files\result_matlab_storage\centroids.mat');
load('D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing\valid_files\result_matlab_storage\DLI_storage.mat');

% Load the spatial clusters (assuming these are stored as separate variables)
load('D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing\valid_files\result_matlab_storage\cluster1.mat');
load('D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing\valid_files\result_matlab_storage\cluster2.mat');
load('D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing\valid_files\result_matlab_storage\cluster3.mat');
load('D:\ibp\DLI-github\sample-data\实验数据\final_data\ChenDanQing\valid_files\result_matlab_storage\cluster4.mat');

% Perform k-means clustering with custom initial centroids
opts = statset('MaxIter', 1000);  % Set options for k-means, if needed
[idx, C] = kmeans(DLI_storage, 3, 'Start', centroids, 'Options', opts);

% Initialize variables to store cosine distances and average DLIs
cosine_distances = zeros(size(DLI_storage, 1), 3);
average_dli = zeros(3, size(DLI_storage, 2));
num_rows_in_states = zeros(3, 1);

% Calculate cosine distances to each centroid
for i = 1:size(DLI_storage, 1)
    for j = 1:3
        cosine_distances(i, j) = 1 - dot(DLI_storage(i, :), centroids(j, :)) / ...
            (norm(DLI_storage(i, :)) * norm(centroids(j, :)));
    end
end

% Calculate the average DLI for each group (state)
for state = 1:3
    state_indices = (idx == state);
    num_rows_in_states(state) = sum(state_indices);
    if num_rows_in_states(state) > 0
        average_dli(state, :) = mean(DLI_storage(state_indices, :), 1);
    end
end

% Initialize matrices to store mean DLI of spatial clusters for each state
spatial_clusters_mean_dli = zeros(3, 4);

% Define the clusters
clusters = {cluster1, cluster2, cluster3, cluster4};

% Calculate the mean DLI for each spatial cluster in each state
for state = 1:3
    for cluster_num = 1:4
        cluster_rois = clusters{cluster_num};
        spatial_clusters_mean_dli(state, cluster_num) = mean(average_dli(state, cluster_rois));
    end
end

% Display results
disp('Number of rows in each state:');
disp(num_rows_in_states);

disp('Average DLI of each ROI in each state:');
disp(average_dli);

disp('Mean DLI of spatial clusters in each state:');
disp(spatial_clusters_mean_dli);


% Assuming spatial_clusters_mean_dli is the matrix containing your data
% Each row represents a state, and each column represents a cluster

% Define colors for the clusters
cluster_colors = [
    0, 0.4470, 0.7410;  % Color for Cluster 1 (blue)
    0.4660, 0.6740, 0.1880;  % Color for Cluster 2 (green)
    0.8500, 0.3250, 0.0980;  % Color for Cluster 3 (orange)
    0.6350, 0.0780, 0.1840;  % Color for Cluster 4 (red)
];

% Create a new figure with adjusted size
figure('Position', [100, 100, 800, 400]);

% Plot the first state
state = 1;
hold on;

% Plot each cluster in a different color with wider bars
for cluster = 1:size(spatial_clusters_mean_dli, 2)
    % Extract the data for the current cluster
    cluster_data = spatial_clusters_mean_dli(state, cluster);
    
    % Plot the data as horizontal bars with increased width
    barh(cluster, cluster_data, 'FaceColor', cluster_colors(cluster, :), 'EdgeColor', 'none', 'BarWidth', 0.8);
    
    % Add text labels next to the bars
    %if cluster_data > 0
    %    text(cluster_data + 0.01, cluster, ['Cluster ', num2str(cluster)], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', cluster_colors(cluster, :));
    %else
    %    text(cluster_data - 0.01, cluster, ['Cluster ', num2str(cluster)], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'Color', cluster_colors(cluster, :));
    %end
end

% Add a vertical line at 0
xline(0, 'k', 'LineWidth', 1.5);

% Set the y-axis labels to represent the clusters
yticks(1:size(spatial_clusters_mean_dli, 2));
yticklabels({'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'});

% Set the x-axis limits for better visualization
xlim([-0.1, 0.1]);  % Adjusted limits for better visibility

% Set the title for the state
title('State 1');

% Add labels R and L for right and left
text(-0.05, size(spatial_clusters_mean_dli, 2) + 0.5, 'R', 'FontSize', 12, 'HorizontalAlignment', 'center');
text(0.05, size(spatial_clusters_mean_dli, 2) + 0.5, 'L', 'FontSize', 12, 'HorizontalAlignment', 'center');

hold off;

% Add a legend
legend({'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4'}, 'Location', 'bestoutside');


