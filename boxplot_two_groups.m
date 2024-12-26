function [] = boxplot_two_groups(features)
    % boxplot_two_groups creates a boxplot for two groups from a matrix of 20 features.
    % Inputs:
    %   features - A 20xN matrix where each row is a signal (N is the number of samples per signal)
    
    % Validate that the input matrix has 20 rows
    if size(features, 1) ~= 20
        error('Input matrix must have exactly 20 rows.');
    end
    
    % Split the matrix into two groups of 10 rows each
    group1 = features(1:10, :);  % First 10 rows
    group2 = features(11:20, :); % Last 10 rows
    
    % Combine the data for boxplot
    combined_data = [group1; group2]; % 20 rows in total
    
    % Create grouping labels for the boxplot
    group_labels = [ones(10, 1); 2 * ones(10, 1)]; % 1 for group1, 2 for group2
    
    % Reshape the data and group labels for the boxplot
    reshaped_data = combined_data(:);
    reshaped_labels = repmat(group_labels, size(features, 2), 1);
    fprintf('Feature columns size : %d\n', size(features, 2));

    % Plot the boxplot
    boxplot(reshaped_data, reshaped_labels);
    
    % Add title and labels
    title('Boxplot for the 2 classes');
    xlabel('Group');
    ylabel('Feature Value');
end
